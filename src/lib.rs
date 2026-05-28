//! Count reads from multiple BAM files overlapping BED intervals.
//!
//! ## Algorithm
//!
//! BED regions are sorted by chrom+start. For each BAM file, we scan each
//! chromosome's reads in a single linear pass (one BAM index seek per chrom),
//! using an event-sorted sweep to maintain the set of active regions as we
//! advance through the BAM. This gives O(G × N/G × log(N/G)) total complexity
//! where G is the chromosome count and N is the region count — essentially
//! linear in the number of reads times the average region depth.
//!
//! This is far more efficient than one BAM random-access query per region
//! (which was the naive approach and cost 50k seeks for 50k regions).
//!
//! ## Read filtering (defaults match bedtools multicov)
//!
//! - Skip UNMAP (0x4), SECONDARY (0x100), QCFAIL (0x200), DUP (0x400).
//! - Skip reads with MAPQ below the minimum.
//! - `-s` / `-S` strand filters compare the read's strand to the BED strand
//!   (column 6). Reads where the BED record has no strand column, or the
//!   strand is `.`, are included regardless of the strand filter.
//!
//! ## Reference
//!
//! `BEDTools multicov` — Quinlan & Hall (2010). Bioinformatics 26(6): 841–842.
//! DOI: 10.1093/bioinformatics/btq033

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use noodles::bam;
use noodles::core::{Position, Region};
use rsomics_common::{Result, RsomicsError};

/// Read filtering and counting options.
// Each bool maps directly to one CLI flag; collapsing to an enum would be more
// confusing than it's worth for a plain option struct.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug)]
pub struct MulticovOpts {
    /// Minimum mapping quality for a read to be counted.
    pub min_mapq: u8,
    /// Include duplicate reads (FLAG 0x400). Default: false.
    pub include_dups: bool,
    /// Include failed-QC reads (FLAG 0x200). Default: false.
    pub include_failed_qc: bool,
    /// Only count properly-paired reads (FLAG 0x2 must be set). Default: false.
    pub proper_pairs_only: bool,
    /// Minimum overlap as a fraction of the BED region length.
    /// Default: 1e-9 (any 1-bp overlap qualifies).
    pub min_overlap_frac: f64,
    /// Require reciprocal overlap: the read must also overlap the BED region
    /// by at least `min_overlap_frac` of the *read* length. Default: false.
    pub reciprocal: bool,
    /// Strand filter: None = any; Some(true) = same strand; Some(false) = opposite strand.
    pub strand_filter: Option<bool>,
}

impl Default for MulticovOpts {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            include_dups: false,
            include_failed_qc: false,
            proper_pairs_only: false,
            min_overlap_frac: 1e-9,
            reciprocal: false,
            strand_filter: None,
        }
    }
}

struct BedRegion {
    chrom: String,
    start: u64,           // 0-based
    end: u64,             // exclusive
    strand: Option<bool>, // Some(true)=+, Some(false)=-, None=unspecified
    raw: Vec<u8>,
}

fn parse_bed_line(line: &[u8]) -> Option<BedRegion> {
    if line.is_empty() || line[0] == b'#' {
        return None;
    }
    let mut it = line.split(|&c| c == b'\t');
    let chrom = std::str::from_utf8(it.next()?).ok()?.to_string();
    let start = parse_u64(it.next()?)?;
    let end = parse_u64(it.next()?)?;
    if start >= end {
        return None;
    }
    // Skip name (col 4) and score (col 5) to reach strand (col 6).
    let _name = it.next();
    let _score = it.next();
    let strand = it.next().and_then(|s| match s {
        b"+" => Some(true),
        b"-" => Some(false),
        _ => None,
    });
    Some(BedRegion {
        chrom,
        start,
        end,
        strand,
        raw: line.to_vec(),
    })
}

fn parse_u64(b: &[u8]) -> Option<u64> {
    if b.is_empty() {
        return None;
    }
    let mut n: u64 = 0;
    for &c in b {
        let d = c.wrapping_sub(b'0');
        if d > 9 {
            return None;
        }
        n = n.checked_mul(10)?.checked_add(u64::from(d))?;
    }
    Some(n)
}

fn load_bed(bed_path: &Path) -> Result<Vec<BedRegion>> {
    let bytes = std::fs::read(bed_path).map_err(|e| RsomicsError::InvalidInput(e.to_string()))?;
    let mut regions = Vec::new();
    for raw in bytes.split(|&b| b == b'\n') {
        let line = match raw.last() {
            Some(b'\r') => &raw[..raw.len() - 1],
            _ => raw,
        };
        if let Some(r) = parse_bed_line(line) {
            regions.push(r);
        }
    }
    Ok(regions)
}

/// Count overlapping reads per BED region across all `bam_paths`.
///
/// Each BAM must have an accompanying `.bam.bai` index — run
/// `samtools index` or `rsomics-bam-index` if missing.
///
/// Returns the number of BED regions emitted.
pub fn multicov(
    bed_path: &Path,
    bam_paths: &[impl AsRef<Path>],
    opts: &MulticovOpts,
    output: &mut dyn Write,
) -> Result<u64> {
    let regions = load_bed(bed_path)?;
    let n_bams = bam_paths.len();
    let mut counts: Vec<Vec<u64>> = vec![vec![0u64; n_bams]; regions.len()];

    for (bam_idx, bam_path) in bam_paths.iter().enumerate() {
        count_bam_sweep(bam_path.as_ref(), &regions, opts, &mut counts, bam_idx)?;
    }

    emit(output, &regions, &counts)
}

/// Per-chrom region data for the sweep.
struct ChromRegion {
    ri: usize,  // index into the global regions Vec
    start: u64, // 0-based
    end: u64,   // exclusive
    reg_len: f64,
    strand: Option<bool>,
}

/// Count reads using a per-chromosome linear sweep instead of per-region seeks.
///
/// Regions are grouped by chromosome, sorted by start position. For each
/// chrom we do one BAM index seek to the first region's start, then scan
/// linearly: active regions (end > `read_start`) are kept; expired ones removed.
/// Each read is tested only against the active region set.
#[allow(clippy::too_many_lines)]
fn count_bam_sweep(
    bam_path: &Path,
    regions: &[BedRegion],
    opts: &MulticovOpts,
    counts: &mut [Vec<u64>],
    bam_idx: usize,
) -> Result<()> {
    let index_path = bam_path.with_extension("bam.bai");
    let index = bam::bai::fs::read(&index_path).map_err(|e| {
        RsomicsError::InvalidInput(format!(
            "cannot open BAM index {}: {e} — run `samtools index` first",
            index_path.display()
        ))
    })?;

    let file = File::open(bam_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bam_path.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    // Group regions by chrom, sorted by start within each chrom.
    let mut chrom_map: HashMap<&str, Vec<ChromRegion>> = HashMap::new();
    for (ri, reg) in regions.iter().enumerate() {
        chrom_map
            .entry(reg.chrom.as_str())
            .or_default()
            .push(ChromRegion {
                ri,
                start: reg.start,
                end: reg.end,
                reg_len: (reg.end - reg.start) as f64,
                strand: reg.strand,
            });
    }
    for v in chrom_map.values_mut() {
        v.sort_unstable_by_key(|r| r.start);
    }

    // Build skip-flags mask.
    let mut skip_flags: u16 = 0x0104; // UNMAP(4) | SECONDARY(256)
    if !opts.include_failed_qc {
        skip_flags |= 0x200; // QCFAIL
    }
    if !opts.include_dups {
        skip_flags |= 0x400; // DUP
    }

    let mut record = bam::Record::default();

    // Process each chromosome independently.
    for (chrom, mut chrom_regions) in chrom_map {
        // Sort by start so we can sweep front to back.
        chrom_regions.sort_unstable_by_key(|r| r.start);

        // The scan start is the leftmost region's start (1-based for noodles).
        let scan_start = chrom_regions[0].start + 1;
        let Ok(pos_start) = Position::try_from(scan_start as usize) else {
            continue;
        };
        // Scan to the rightmost region's end.
        let scan_end = chrom_regions.iter().map(|r| r.end).max().unwrap_or(0);
        let Ok(pos_end) = Position::try_from(scan_end as usize) else {
            continue;
        };

        let noodles_region = Region::new(chrom.as_bytes(), pos_start..=pos_end);
        let Ok(mut query) = reader.query(&header, &index, &noodles_region) else {
            continue; // chrom absent from BAM header
        };

        // next_region_idx advances as the sweep position moves forward.
        let mut next_region_idx = 0usize;
        // active: regions whose end > current read_start.
        let mut active: Vec<&ChromRegion> = Vec::new();

        loop {
            let n = query.read_record(&mut record).map_err(RsomicsError::Io)?;
            if n == 0 {
                break;
            }

            let flags = record.flags().bits();
            if (flags & skip_flags) != 0 {
                continue;
            }
            if opts.proper_pairs_only && (flags & 0x2) == 0 {
                continue;
            }
            let mq = record.mapping_quality().map_or(0, |q| q.get());
            if mq < opts.min_mapq {
                continue;
            }

            let Some(aln_start_pos) = record.alignment_start().transpose().ok().flatten() else {
                continue;
            };
            let read_start = aln_start_pos.get() as u64 - 1; // 0-based
            let read_span = ref_span(&record)?;
            if read_span == 0 {
                continue;
            }
            let read_end = read_start + read_span;

            // Advance next_region_idx: add regions whose start <= read_end.
            while next_region_idx < chrom_regions.len()
                && chrom_regions[next_region_idx].start < read_end
            {
                active.push(&chrom_regions[next_region_idx]);
                next_region_idx += 1;
            }

            // Remove expired regions (end <= read_start means no overlap).
            active.retain(|r| r.end > read_start);

            // Check each active region for overlap.
            for reg in &active {
                let lo = read_start.max(reg.start);
                let hi = read_end.min(reg.end);
                if hi <= lo {
                    continue;
                }
                let overlap = (hi - lo) as f64;

                if overlap / reg.reg_len < opts.min_overlap_frac {
                    continue;
                }
                if opts.reciprocal && overlap / (read_span as f64) < opts.min_overlap_frac {
                    continue;
                }

                // Strand filter.
                if let Some(require_same) = opts.strand_filter
                    && let Some(bed_plus) = reg.strand
                {
                    let read_forward = (flags & 0x10) == 0;
                    let strands_same = bed_plus == read_forward;
                    if require_same != strands_same {
                        continue;
                    }
                }

                counts[reg.ri][bam_idx] += 1;
            }
        }
    }

    Ok(())
}

fn ref_span(record: &bam::Record) -> Result<u64> {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let mut span: u64 = 0;
    for op in record.cigar().iter() {
        let op = op.map_err(RsomicsError::Io)?;
        match op.kind() {
            Kind::Match
            | Kind::Deletion
            | Kind::Skip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch => {
                span += op.len() as u64;
            }
            _ => {}
        }
    }
    Ok(span)
}

fn emit(output: &mut dyn Write, regions: &[BedRegion], counts: &[Vec<u64>]) -> Result<u64> {
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut ib = itoa::Buffer::new();
    let mut emitted: u64 = 0;

    for (ri, reg) in regions.iter().enumerate() {
        out.write_all(&reg.raw).map_err(RsomicsError::Io)?;
        for &c in &counts[ri] {
            out.write_all(b"\t").map_err(RsomicsError::Io)?;
            out.write_all(ib.format(c).as_bytes())
                .map_err(RsomicsError::Io)?;
        }
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        emitted += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(emitted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_bed_line_basic() {
        let line = b"chr1\t100\t200\tname\t0\t+";
        let r = parse_bed_line(line).unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, 200);
        assert_eq!(r.strand, Some(true));
    }

    #[test]
    fn parse_bed_line_minus_strand() {
        let line = b"chr2\t300\t400\tgene\t0\t-";
        let r = parse_bed_line(line).unwrap();
        assert_eq!(r.strand, Some(false));
    }

    #[test]
    fn parse_bed_line_no_strand() {
        let line = b"chr1\t0\t1000";
        let r = parse_bed_line(line).unwrap();
        assert_eq!(r.strand, None);
    }

    #[test]
    fn parse_bed_line_skips_header() {
        assert!(parse_bed_line(b"#track name=test").is_none());
    }

    #[test]
    fn parse_bed_line_degenerate() {
        // start == end → skip
        assert!(parse_bed_line(b"chr1\t100\t100").is_none());
    }
}
