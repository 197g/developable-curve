use super::*;

pub struct Spiral {
    /// Radius of curvature.
    pub radius: f32,
    /// Z-offset per rotation.
    pub pitch: f32,
}

impl Curve for Spiral {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let (y, x) = (f64::from(t) * core::f64::consts::PI * 2.0f64).sin_cos();

        let f0 = f64::from(self.radius);
        let f1 = f64::from(self.radius) * core::f64::consts::PI * 2.0f64;
        let f2 = f64::from(self.radius) * (core::f64::consts::PI * 2.0f64).powi(2);
        let f3 = f64::from(self.radius) * (core::f64::consts::PI * 2.0f64).powi(3);

        DenormalTangentFrame {
            base: DVec3::new(x * f0, y * f0, f64::from(t) * f64::from(self.pitch)),
            tangent: DVec3::new(-y * f1, x * f1, f64::from(self.pitch)),
            derivative: DVec3::new(-x * f2, -y * f2, 0.0),
            third_derivative: DVec3::new(y * f3, -x * f3, 0.0),
        }
    }
}
