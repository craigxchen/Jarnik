use criterion::criterion_main;

mod benchmarks;

criterion_main! {
    benchmarks::rad_from_1to1e5::benches,
    benchmarks::rad_1e8::benches,
}