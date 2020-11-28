use criterion::{criterion_group, Criterion};
use std::time::Duration;

use jarnik::{
    get_arclength_coeffs
};

fn test_arclength_large_rad(c: &mut Criterion) {
    let mut group = c.benchmark_group("Arc Length Ratio");
    group.sample_size(10);
    group.measurement_time(Duration::new(100, 0));

    group.bench_function("Large R",  move |b| {
        b.iter(|| get_arclength_coeffs(1e8 as u64, (1e8 + 1e5) as u64));
    });

    group.finish();
}

criterion_group!(benches, test_arclength_large_rad);