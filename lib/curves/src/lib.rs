use glam::Vec3;

pub use dc_theory::*;

pub trait Curve {
    fn at(&self, t: f32) -> DenormalTangentFrame;

    fn sample_at(&self, ts: &[f32]) -> Vec<DenormalTangentFrame> {
        ts.iter().map(|&t| self.at(t)).collect()
    }
}

impl<T: Curve + ?Sized> Curve for &'_ T {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        (**self).at(t)
    }
}

pub fn normal_ode(curve: impl Curve, parameter: impl Fn(f32) -> f32) -> impl Fn(Vec3, f32) -> Vec3 {
    move |normal: Vec3, t: f32| {
        let frame = curve.at(t);
        let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { normal });
        let lambda = parameter(t);
        dev.derivative_base + lambda * dev.derivative_free
    }
}

pub fn normal_and_flat_ode(
    curve: impl Curve,
    parameter: impl Fn(f32) -> f32,
) -> impl Fn(Vec3, f32) -> CurveDescription {
    move |normal: Vec3, t: f32| {
        let frame = curve.at(t);
        let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { normal });
        let lambda = parameter(t);
        let dt_normal = dev.derivative_base + lambda * dev.derivative_free;
        let speed = frame.tangent.length();

        CurveDescription {
            tangent: dev.frame.tangent,
            dt_normal,
            curvature: dev.surface_curvature,
            speed,
        }
    }
}

/// Steer the surface and horizontal direction by defining an angle between the tangent and
/// horizontal direction along the curve.
pub fn normal_and_angle_ode(
    curve: impl Curve,
    parameter: impl Fn(f32) -> f32,
) -> impl Fn(Vec3, f32) -> CurveDescription {
    move |normal: Vec3, t: f32| {
        let frame = curve.at(t);
        let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { normal });
        let target_angle = parameter(t);
        // angle(horizontal, frame.tangent) = atan2(<normal, frame.derivative>, lambda)
        //
        // See `dc-integral/src/lib.rs` for the derivation of this formula where lambda is the
        // parameter from the above formula. Now let's derive that lambda. Note how we
        // automatically get `lambda = 0` at the direction discontinuity.
        let lambda = target_angle.tan() * dev.normal.dot(frame.derivative);
        // ^ LLM anecdote: this was oneshot before the derivation. It badly fumbled the
        // derivation itself though, forgetting the square root in the identity or forgetting
        // that subtract `1` changes the numerator..
        let dt_normal = dev.derivative_base + lambda * dev.derivative_free;
        let speed = frame.tangent.length();

        CurveDescription {
            tangent: dev.frame.tangent,
            dt_normal,
            curvature: dev.surface_curvature,
            speed,
        }
    }
}

pub struct Circle {
    pub radius: f32,
}

impl Curve for Circle {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let base = Vec3::new(t.cos(), t.sin(), 0.0) * self.radius;
        let tangent = Vec3::new(-t.sin(), t.cos(), 0.0) * self.radius;
        let derivative = Vec3::new(-t.cos(), -t.sin(), 0.0) * self.radius;
        let binormal = Vec3::new(0.0, 0.0, 1.0);

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            binormal,
        }
    }
}

pub struct HermiteSpline<const N: usize> {
    /// The points the curve should pass through.
    pub points: [Vec3; N],
}

impl Curve for HermiteSpline<2> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];

        let base = p0.lerp(p1, t);
        let tangent = p1 - p0;
        let derivative = Vec3::ZERO;
        let binormal = Vec3::ZERO;

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            binormal,
        }
    }
}

impl Curve for HermiteSpline<3> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];
        let p2 = self.points[2];

        let base = p0 * (1.0 - t).powi(2) + p1 * 2.0 * (1.0 - t) * t + p2 * t.powi(2);
        let tangent = (p1 - p0) * 2.0 * (1.0 - t) + (p2 - p1) * 2.0 * t;
        let derivative = (p2 - 2.0 * p1 + p0) * 2.0;
        let binormal = Vec3::ZERO;

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            binormal,
        }
    }
}

// This is where it gets interesting, this curve may have a curl.
impl Curve for HermiteSpline<4> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];
        let p2 = self.points[2];
        let p3 = self.points[3];

        // TODO: this is a stupid way of evaluating this.
        let base = p0 * (1.0 - t).powi(3)
            + p1 * 3.0 * (1.0 - t).powi(2) * t
            + p2 * 3.0 * (1.0 - t) * t.powi(2)
            + p3 * t.powi(3);
        let tangent = (p1 - p0) * 3.0 * (1.0 - t).powi(2)
            + (p2 - p1) * 6.0 * (1.0 - t) * t
            + (p3 - p2) * 3.0 * t.powi(2);
        let derivative = (p2 - 2.0 * p1 + p0) * 6.0 * (1.0 - t) + (p3 - 2.0 * p2 + p1) * 6.0 * t;

        // At first gpt-41 said this was `0`, but that is only the case for `t=0` and `t=1`. The
        // binormal is the derivative of the derivative, which is a constant in this case.
        // ^ This is its own autocomplete when I removed the line. LMAO.
        let binormal = (p3 - 3.0 * p2 + 3.0 * p1 - p0) * 6.0;

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            binormal,
        }
    }
}
