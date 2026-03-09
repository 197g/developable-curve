//! This module contains some experiments that show what the rustc and the llvm can or cannot optimize
use core::hint::black_box;

// #[inline]
// pub fn dgl_x(x: f64, v: f64) -> f64 {
//     x.cos() - v.sin()
// }
// #[inline]
// pub fn dgl_v(x: f64, v: f64) -> f64 {
//     x.tan() - v.exp()
// }
#[inline]
pub fn dgl_x(_x: f64, v: f64) -> f64 {
    v
}
#[inline]
pub fn dgl_v(x: f64, _v: f64) -> f64 {
    -x
}

pub fn direct_euler() -> (f64, f64) {
    let mut x = 1.;
    let mut v = 0.;
    let dt = 0.01;
    for _ in 0..10000 {
        let new_x = x + dgl_x(x, v) * dt;
        let new_v = v + dgl_v(x, v) * dt;
        x = new_x;
        v = new_v;
    }
    (x, v)
}

fn vec_dgl(state: &[f64]) -> Vec<f64> {
    let x = state[0];
    let v = state[1];
    vec![dgl_x(x, v), dgl_v(x, v)]
}
pub fn vec_inplace_euler() -> (f64, f64) {
    let mut state = vec![1., 0.];
    let dt = 0.01;
    for _ in 0..10000 {
        let diff = vec_dgl(&state);
        // diff.into_iter().zip(state.iter_mut()).for_each(|(diff, state)| *state = *state + diff * dt ); // same performance, [https://doc.rust-lang.org/std/iter/trait.Iterator.html#method.for_each] write: This is equivalent to using a for loop on the iterator
        for (diff, state) in diff.into_iter().zip(state.iter_mut()) {
            *state += diff * dt;
        }
    }
    (state[0], state[1])
}

pub fn vec_collect_euler() -> (f64, f64) {
    let mut state = vec![1., 0.];
    let dt = 0.01;
    for _ in 0..10000 {
        let diff = vec_dgl(&state);
        state = diff
            .into_iter()
            .zip(state.iter_mut())
            .map(|(diff, state)| *state + diff * dt)
            .collect::<Vec<f64>>();
    }
    (state[0], state[1])
}

fn tuple_dgl(state: (f64, f64)) -> (f64, f64) {
    let x = state.0;
    let v = state.1;
    (dgl_x(x, v), dgl_v(x, v))
}
pub fn tuple_euler() -> (f64, f64) {
    let mut state = (1., 0.);
    let dt = 0.01;
    for _ in 0..10000 {
        let diff = tuple_dgl(state);
        state.0 += diff.0 * dt;
        state.1 += diff.1 * dt;
    }
    state
}

fn array_dgl(state: [f64; 2]) -> [f64; 2] {
    let x = state[0];
    let v = state[1];
    [dgl_x(x, v), dgl_v(x, v)]
}
pub fn array_euler() -> (f64, f64) {
    let mut state = [1., 0.];
    let dt = 0.01;
    for _ in 0..10000 {
        let diff = array_dgl(state);
        for (diff, state) in diff.iter().zip(state.iter_mut()) {
            *state += diff * dt;
        }
    }
    (state[0], state[1])
}

#[derive(Copy, Clone)]
struct State {
    x: f64,
    v: f64,
}
fn struct_dgl(state: State) -> State {
    let x = state.x;
    let v = state.v;
    State {
        x: dgl_x(x, v),
        v: dgl_v(x, v),
    }
}
pub fn struct_euler() -> (f64, f64) {
    let mut state = State { x: 1., v: 0. };
    let dt = 0.01;
    for _ in 0..10000 {
        let diff = struct_dgl(state);
        state.x += diff.x * dt;
        state.v += diff.v * dt;
    }
    (state.x, state.v)
}

pub fn fastest() -> (f64, f64) {
    let mut state = vec![1., 0.];
    for _ in 0..10000 {
        let diff = vec![state[1], -state[0]];
        for (diff, state) in diff.into_iter().zip(state.iter_mut()) {
            *state += diff * 0.01;
        }
    }
    (state[0], state[1])
}

fn looper(mut state: Vec<f64>) -> Vec<f64> {
    //let mut state = input.clone();
    //let mut state = vec![input[0], input[1]];

    for _ in 0..10000 {
        let diff = vec_dgl(&state);
        //let diff = vec![state[1], -state[0]];
        // let x = state[0];
        // let v = state[1];
        // let diff = vec![v, -x];
        for (diff, state) in diff.into_iter().zip(state.iter_mut()) {
            *state += diff * 0.01;
        }
    }
    state
}
pub fn preparer() -> (f64, f64) {
    // let mut state = test::black_box(vec![1., 0.]);
    // let state = test::black_box(looper)(state);
    let state = vec![black_box(1.), black_box(0.)];
    let state = looper(state);
    (state[0], state[1])
}

pub fn weird_fast() -> Vec<f64> {
    let mut state = vec![black_box(1.), black_box(1.)];
    for _ in 0..10000 {
        let diff = vec![state[0], state[1]];
        for (diff, state) in diff.into_iter().zip(state.iter_mut()) {
            *state += diff * 0.01;
        }
    }
    state
}

pub fn weird_slow() -> Vec<f64> {
    let mut state = vec![black_box(1.), black_box(1.)];
    for _ in 0..10000 {
        // let a = state[0];
        // let b = state[1];
        // let diff = vec![-b, a];
        let diff = vec_dgl(&state);
        for (diff, state) in diff.into_iter().zip(state.iter_mut()) {
            *state += diff * 0.01;
        }
    }
    state
}

#[test]
fn compare() {
    assert_eq!(direct_euler(), vec_inplace_euler());
    assert_eq!(direct_euler(), vec_collect_euler());
    assert_eq!(direct_euler(), array_euler());
    assert_eq!(direct_euler(), tuple_euler());
    assert_eq!(direct_euler(), struct_euler());
    //assert_eq!(direct_euler(), lib_euler());
}
