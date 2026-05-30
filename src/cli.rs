use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_multicov::{MulticovOpts, multicov};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

// Each bool maps directly to one bedtools multicov CLI flag.
#[allow(clippy::struct_excessive_bools)]
#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-multicov", disable_help_flag = true)]
pub struct Cli {
    /// BAM files to count reads from (one or more).
    #[arg(long = "bams", required = true, num_args = 1..)]
    pub bams: Vec<PathBuf>,

    /// BED file defining the regions.
    #[arg(long = "bed", required = true)]
    pub bed: PathBuf,

    /// Minimum mapping quality. Default 0.
    #[arg(long = "min-mapq", default_value_t = 0)]
    pub min_mapq: u8,

    /// Include duplicate reads (FLAG 0x400).
    #[arg(short = 'D')]
    pub include_dups: bool,

    /// Include failed-QC reads (FLAG 0x200).
    #[arg(short = 'F')]
    pub include_failed_qc: bool,

    /// Only count properly-paired reads (FLAG 0x2 must be set).
    #[arg(short = 'p')]
    pub proper_pairs_only: bool,

    /// Minimum overlap as a fraction of the BED region. Default 1e-9.
    #[arg(short = 'f', default_value_t = 1e-9)]
    pub min_overlap_frac: f64,

    /// Require reciprocal overlap: read must also overlap >= `-f` of its own length.
    #[arg(short = 'r')]
    pub reciprocal: bool,

    /// Require same strandedness (BED column 6 vs read FLAG 0x10).
    #[arg(short = 's')]
    pub same_strand: bool,

    /// Require opposite strandedness (BED column 6 vs read FLAG 0x10).
    #[arg(short = 'S')]
    pub opposite_strand: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let strand_filter = if self.same_strand {
            Some(true)
        } else if self.opposite_strand {
            Some(false)
        } else {
            None
        };

        let opts = MulticovOpts {
            min_mapq: self.min_mapq,
            include_dups: self.include_dups,
            include_failed_qc: self.include_failed_qc,
            proper_pairs_only: self.proper_pairs_only,
            min_overlap_frac: self.min_overlap_frac,
            reciprocal: self.reciprocal,
            strand_filter,
        };

        let stdout = io::stdout();
        let mut out = stdout.lock();
        multicov(&self.bed, &self.bams, &opts, &mut out)?;
        Ok(())
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Count reads from multiple BAMs at BED intervals (bedtools multicov equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["--bams <bam>... --bed <bed> [OPTIONS]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "bams",
                aliases: &[],
                value: Some("<path>..."),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input BAM files (indexed, one or more)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "bed",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "BED file defining the regions to count",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "min-mapq",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("u8"),
                required: false,
                default: Some("0"),
                description: "Minimum mapping quality",
                why_default: None,
            },
            FlagSpec {
                short: Some('D'),
                long: "include-dups",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Include duplicate reads (FLAG 0x400)",
                why_default: None,
            },
            FlagSpec {
                short: Some('F'),
                long: "include-failed-qc",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Include failed-QC reads (FLAG 0x200)",
                why_default: None,
            },
            FlagSpec {
                short: Some('p'),
                long: "proper-pairs-only",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Only count properly-paired reads",
                why_default: None,
            },
            FlagSpec {
                short: Some('f'),
                long: "min-overlap-frac",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("1e-9"),
                description: "Minimum overlap as fraction of BED region length",
                why_default: None,
            },
            FlagSpec {
                short: Some('r'),
                long: "reciprocal",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Require reciprocal overlap (-f fraction of read length too)",
                why_default: None,
            },
            FlagSpec {
                short: Some('s'),
                long: "same-strand",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Require same strandedness (BED col 6 vs read FLAG 0x10)",
                why_default: None,
            },
            FlagSpec {
                short: Some('S'),
                long: "opposite-strand",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Require opposite strandedness",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Count reads from three BAMs at BED intervals",
            command: "rsomics-bed-multicov --bams s1.bam s2.bam s3.bam --bed peaks.bed",
        },
        Example {
            description: "Same strand, min MAPQ 20",
            command: "rsomics-bed-multicov --bams sample.bam --bed genes.bed -s --min-mapq 20",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
