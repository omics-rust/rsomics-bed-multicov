use criterion::{Criterion, criterion_group, criterion_main};
use std::io;
use std::path::Path;

use rsomics_bed_multicov::{MulticovOpts, multicov};

fn bench_multicov(c: &mut Criterion) {
    // Tier-1 fixture; skip if absent.
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let bam = dir.join("test.bam");
    let bed = dir.join("regions.bed");
    if !bam.exists() || !bed.exists() {
        return;
    }

    let bam_paths = vec![bam];
    let opts = MulticovOpts::default();

    c.bench_function("multicov_golden", |b| {
        b.iter(|| {
            let mut out = io::sink();
            multicov(&bed, &bam_paths, &opts, &mut out).unwrap();
        });
    });
}

criterion_group!(benches, bench_multicov);
criterion_main!(benches);
