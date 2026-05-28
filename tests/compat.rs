//! Compatibility tests: rsomics-bed-multicov vs bedtools multicov.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> PathBuf {
    env!("CARGO_BIN_EXE_rsomics-bed-multicov").into()
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn run_ours(bams: &[&str], bed: &str, extra: &[&str]) -> String {
    let dir = golden_dir();
    let bam_paths: Vec<PathBuf> = bams.iter().map(|b| dir.join(b)).collect();
    let bed_path = dir.join(bed);

    let mut cmd = Command::new(bin());
    for bam in &bam_paths {
        cmd.arg("--bams").arg(bam);
    }
    cmd.arg("--bed").arg(&bed_path);
    for arg in extra {
        cmd.arg(arg);
    }

    let out = cmd.output().expect("failed to run rsomics-bed-multicov");
    assert!(
        out.status.success(),
        "rsomics-bed-multicov failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn bedtools_version() -> Option<String> {
    Command::new("bedtools")
        .arg("--version")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
}

fn run_bedtools(bam: &str, bed: &str, extra: &[&str]) -> Option<String> {
    let dir = golden_dir();
    let bam_path = dir.join(bam);
    let bed_path = dir.join(bed);

    let mut cmd = Command::new("bedtools");
    cmd.args(["multicov", "-bams"])
        .arg(&bam_path)
        .arg("-bed")
        .arg(&bed_path);
    for arg in extra {
        cmd.arg(arg);
    }

    let out = cmd.output().ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8(out.stdout).unwrap())
}

#[test]
fn golden_output_matches() {
    let got = run_ours(&["test.bam"], "regions.bed", &[]);
    let expected = fs::read_to_string(golden_dir().join("expected.tsv")).unwrap();
    assert_eq!(
        got.trim(),
        expected.trim(),
        "golden mismatch:\ngot:\n{got}\nexpected:\n{expected}"
    );
}

#[test]
fn output_matches_bedtools() {
    let Some(version) = bedtools_version() else {
        eprintln!("bedtools not found — skipping compat test");
        return;
    };
    eprintln!("bedtools version: {version}");

    let Some(expected) = run_bedtools("test.bam", "regions.bed", &[]) else {
        eprintln!("bedtools multicov failed — skipping");
        return;
    };
    let got = run_ours(&["test.bam"], "regions.bed", &[]);
    assert_eq!(
        got.trim(),
        expected.trim(),
        "mismatch vs bedtools:\ngot:\n{got}\nexpected:\n{expected}"
    );
}

#[test]
fn same_strand_filter() {
    // -s: only same strand reads.
    // region1 (+): read1 (forward, same) + read2 (forward, same) = 2; read4 (reverse, diff) excluded
    // region2 (+): read3 (forward) = 1
    // region3 (-): read5 (forward, opposite to -) = 0
    let got = run_ours(&["test.bam"], "regions.bed", &["-s"]);
    let lines: Vec<&str> = got.trim().lines().collect();
    assert_eq!(lines.len(), 3);
    // region1: 2 reads (forward only)
    let c1: u64 = lines[0].split('\t').next_back().unwrap().parse().unwrap();
    assert_eq!(c1, 2, "region1 same-strand count: {}", lines[0]);
    // region2: 1 read
    let c2: u64 = lines[1].split('\t').next_back().unwrap().parse().unwrap();
    assert_eq!(c2, 1, "region2 same-strand count: {}", lines[1]);
    // region3 (-): read5 is forward → opposite strand to - → 0
    let c3: u64 = lines[2].split('\t').next_back().unwrap().parse().unwrap();
    assert_eq!(c3, 0, "region3 same-strand count: {}", lines[2]);
}

#[test]
fn min_mapq_filter() {
    // All test reads have MAPQ=60; --min-mapq 61 should yield all zeros.
    let got = run_ours(&["test.bam"], "regions.bed", &["--min-mapq", "61"]);
    for line in got.trim().lines() {
        let count: u64 = line.split('\t').next_back().unwrap().parse().unwrap();
        assert_eq!(count, 0, "expected 0 with high MAPQ threshold: {line}");
    }
}
