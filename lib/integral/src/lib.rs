use gauss_quad::GaussLegendre;
use glam::Vec3;

pub fn integrate_curve(tangent: impl Fn(f32) -> Vec3, t: [f32; 2]) -> Vec3 {
    let [left, right] = t.map(f64::from);
    let scheme = GaussLegendre::new(5).unwrap();

    Vec3::new(
        scheme.integrate(left, right, |t| f64::from(tangent(t as f32).x)) as f32,
        scheme.integrate(left, right, |t| f64::from(tangent(t as f32).y)) as f32,
        scheme.integrate(left, right, |t| f64::from(tangent(t as f32).z)) as f32,
    )
}
