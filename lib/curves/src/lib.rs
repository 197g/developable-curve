use dc_integral::CurveSegment;
use glam::{DVec3, Vec3};

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

pub fn normal_ode(curve: impl Curve, parameter: impl Fn(f32) -> f32) -> impl Fn(DVec3, f64) -> DVec3 {
    move |normal: DVec3, t: f64| {
        let frame = curve.at(t as f32);
        let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { axis: normal });
        let lambda = f64::from(parameter(t as f32));
        dev.derivative_base + lambda * dev.derivative_free
    }
}

pub fn normal_and_flat_ode(
    curve: impl Curve,
    parameter: impl Fn(f32) -> f32,
) -> impl Fn(DVec3, f64) -> CurveDescription {
    move |normal: DVec3, t: f64| {
        let dev =
            SurfaceDevelopment::from_frame_and_normal(curve.at(t as f32), SurfaceNormal { axis: normal });
        let lambda = f64::from(parameter(t as f32));
        let dt_normal = dev.derivative_base + lambda * dev.derivative_free;

        CurveDescription {
            tangent: dev.frame.tangent,
            dt_tangent: dev.frame.derivative,
            dt_normal,
            angle: None,
            dt_normal_normalized: None,
        }
    }
}

/// Steer the surface and horizontal direction by defining an angle between the tangent and
/// horizontal direction along the curve.
pub fn normal_and_tan_ode(
    curve: impl Curve,
    // Provides intended `angle` at a point `t` on the curve where `0` is orthogonal to the left
    // side of the tangent in the normal plane.
    parameter: impl Fn(f32) -> f32,
) -> impl Fn(DVec3, f64) -> CurveDescription {
    move |normal: DVec3, t: f64| {
        let frame = curve.at(t as f32);
        let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { axis: normal });

        // Set so that 0 refers to the same side as the 0 of the other parameterization.
        let target_angle = f64::from(parameter(t as f32)) + std::f64::consts::PI * 0.5;

        // angle(horizontal, frame.tangent) = atan2(<normal, frame.derivative>, lambda)
        //
        // See `dc-integral/src/lib.rs` for the derivation of this formula where lambda is the
        // parameter from the above formula. Now let's derive that lambda. Note how we
        // automatically get `lambda = 0` at the direction discontinuity. Plus see that we get the
        // same result for `-target_angle` and the signum there chooses the orientation.
        let tan_angle = target_angle.tan();

        let angle_or_zero = if tan_angle == 0.0 {
            1e-6
        } else if tan_angle.abs() < 1e-6 {
            1e-6 * tan_angle.signum()
        } else {
            tan_angle
        };

        let lambda = dev.normal.dot(frame.derivative).abs() / angle_or_zero;

        // Note this may get the sign of the parameter incorrect. Now look at the direction of the
        // surface development to determine the right sign. For using the existing sign, tangent,
        // dir, normal should form a right-handed system. The curve has a 'signum' already but that
        // is for the tangent normal system.
        let signum = dev
            .normal
            .cross(frame.tangent)
            .dot(dev.derivative_free)
            .signum();

        // ^ LLM anecdote: oneshot incorrectly. Previously believed to be correct though but it
        // multiplied instead of divided. Stupid machine, stupid me for trusting it too much and
        // getting myself confused in the process.
        //
        // It badly fumbled the derivation itself already, forgetting the square root in the
        // cos(x)-identity or forgetting that subtract `1` changes the numerator..
        let dt_normal = dev.derivative_base + lambda * signum * dev.derivative_free;
        let dt_normalized = frame.tangent.normalize();

        CurveDescription {
            tangent: dev.frame.tangent,
            dt_tangent: dev.frame.derivative,
            dt_normal,
            angle: Some(target_angle),
            dt_normal_normalized: Some(dt_normalized.normalize()),
        }
    }
}

pub fn stitch(end: CurveSegment, curve: impl Curve) -> (CurveSegment, f32) {
    // Determine `parameter` such that the horizontals match up (given that the normals will match
    // up as well).
    let parameter = {
        let start_frame = curve.at(0.0);
        let development = SurfaceDevelopment::from_frame_and_normal(start_frame, end.normal);

        let raw_angle = end.horizontal.angle_between(development.frame.tangent);
        let signum = development
            .frame
            .tangent
            .cross(end.horizontal)
            .dot(development.normal)
            .signum();

        signum * raw_angle - std::f64::consts::PI * 0.5
    };

    let ode = normal_and_tan_ode(curve, |_| parameter as f32);
    let basis = CurveSegment::initial(end.normal.axis, ode);

    let start = CurveSegment {
        flat_position: end.flat_position,
        flat_direction: end.flat_direction,
        flat_curvature: end.flat_curvature,
        ..basis
    };

    (start, parameter as f32)
}

pub struct Circle {
    pub radius: f32,
}

impl Curve for Circle {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let t = f64::from(t);

        let base = DVec3::new(t.cos(), t.sin(), 0.0) * f64::from(self.radius);
        let tangent = DVec3::new(-t.sin(), t.cos(), 0.0) * f64::from(self.radius);
        let derivative = DVec3::new(-t.cos(), -t.sin(), 0.0) * f64::from(self.radius);
        let binormal = DVec3::new(0.0, 0.0, 1.0);

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            binormal,
        }
    }
}

pub struct BezierSpline<const N: usize> {
    /// The points the curve should pass through.
    pub points: [Vec3; N],
}

impl Curve for BezierSpline<2> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];

        let base = p0.lerp(p1, t).as_dvec3();
        let tangent = (p1 - p0).as_dvec3();
        let derivative = DVec3::ZERO;
        let binormal = DVec3::ZERO;

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            binormal,
        }
    }
}

impl Curve for BezierSpline<3> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];
        let p2 = self.points[2];

        let base = p0 * (1.0 - t).powi(2) + p1 * 2.0 * (1.0 - t) * t + p2 * t.powi(2);
        let tangent = (p1 - p0) * 2.0 * (1.0 - t) + (p2 - p1) * 2.0 * t;
        let derivative = (p2 - 2.0 * p1 + p0) * 2.0;
        let binormal = DVec3::ZERO;

        DenormalTangentFrame {
            base: base.as_dvec3(),
            tangent: tangent.as_dvec3(),
            derivative: derivative.as_dvec3(),
            binormal,
        }
    }
}

// This is where it gets interesting, this curve may have a curl.
impl Curve for BezierSpline<4> {
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
            base: base.as_dvec3(),
            tangent: tangent.as_dvec3(),
            derivative: derivative.as_dvec3(),
            binormal: binormal.as_dvec3(),
        }
    }
}
