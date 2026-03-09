use fast_ode::optimizer_tests::*;

#[bench]
fn bench_direct(b: &mut Bencher) {
    b.iter(direct_euler);
}

#[bench]
fn bench_vec_inplace(b: &mut Bencher) {
    b.iter(vec_inplace_euler);
}

#[bench]
fn bench_vec_collect(b: &mut Bencher) {
    b.iter(vec_collect_euler);
}

#[bench]
fn bench_tuple(b: &mut Bencher) {
    b.iter(tuple_euler);
}

#[bench]
fn bench_array(b: &mut Bencher) {
    b.iter(array_euler);
}

#[bench]
fn bench_struct(b: &mut Bencher) {
    b.iter(struct_euler);
}

#[bench]
fn bench_fastest(b: &mut Bencher) {
    b.iter(fastest);
}

#[bench]
fn bench_preparer(b: &mut Bencher) {
    b.iter(preparer);
}

#[bench]
fn bench_weird_fast(b: &mut Bencher) {
    b.iter(weird_fast);
}

#[bench]
fn bench_weird_slow(b: &mut Bencher) {
    b.iter(weird_slow);
}
