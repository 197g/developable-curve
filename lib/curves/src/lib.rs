use dc_integral::{CurveSegment, PipeFaceBase, TrianglePipeBase};
use glam::{DVec2, DVec3, Vec3};

mod affine;
mod spiral;
mod spline;

pub use affine::{Affine, Translate};
pub use dc_theory::*;
pub use spiral::Spiral;
pub use spline::BezierSpline;

pub trait Curve {
    fn at(&self, t: f32) -> DenormalTangentFrame;

    fn sample_at(&self, ts: &[f32]) -> Vec<DenormalTangentFrame> {
        ts.iter().map(|&t| self.at(t)).collect()
    }
}

impl dyn Curve {
    pub fn as_pipe(&self, base1: DVec3, opposing: SurfaceNormal, base2: DVec3) -> TrianglePipeBase {
        fn flat_around(dleft: DVec3, orthogonal_frame: [DVec3; 2], dright: DVec3) -> PipeFaceBase {
            // Defines the right and forward direction, *not* orthonormal.
            let [right, forward] = orthogonal_frame;

            // Establish the angles via atan2.
            let x = forward.normalize_or_zero();
            let negy = right.normalize_or_zero();

            let angle = |vec: DVec3| {
                let vec = vec.normalize_or_zero();
                let cos = vec.dot(x);
                let sin = vec.dot(-negy);
                sin.atan2(cos)
            };

            PipeFaceBase {
                base_left: DVec2::ZERO,
                base_right: DVec2::new(0.0, -right.length()),
                orientation_left: angle(dleft),
                orientation_right: angle(dright),
            }
        }

        fn flats_from_inner_normals(base: [DVec3; 3], opposing: [DVec3; 3]) -> [PipeFaceBase; 3] {
            let forwards = [
                opposing[2].cross(opposing[1]),
                opposing[0].cross(opposing[2]),
                opposing[1].cross(opposing[0]),
            ];

            let opposing_base = {
                let rights = [base[2] - base[1], base[0] - base[2], base[1] - base[0]];

                let forwards = [
                    opposing[0].cross(rights[0]),
                    opposing[1].cross(rights[1]),
                    opposing[2].cross(rights[2]),
                ];

                [
                    [rights[0], forwards[0]],
                    [rights[1], forwards[1]],
                    [rights[2], forwards[2]],
                ]
            };

            [
                flat_around(forwards[0], opposing_base[0], forwards[1]),
                flat_around(forwards[1], opposing_base[1], forwards[2]),
                flat_around(forwards[2], opposing_base[2], forwards[0]),
            ]
        }

        let frame = self.at(0.0);
        let normalb = frame.tangent.cross(base2 - frame.base);
        let normala = frame.tangent.cross(frame.base - base1);

        let [flat1, opposing_flat, flat2] = flats_from_inner_normals(
            [frame.base, base1, base2],
            [opposing.axis, normalb, normala],
        );

        TrianglePipeBase {
            base1,
            base2,
            opposing_normal: opposing.axis,
            flat1,
            flat2,
            opposing_flat,
        }
    }
}

impl<T: Curve + ?Sized> Curve for &'_ T {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        (**self).at(t)
    }
}

pub fn normal_ode(
    curve: impl Curve,
    parameter: impl Fn(f32) -> f32,
) -> impl Fn(DVec3, f64) -> DVec3 {
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
        let dev = SurfaceDevelopment::from_frame_and_normal(
            curve.at(t as f32),
            SurfaceNormal { axis: normal },
        );
        let lambda = f64::from(parameter(t as f32));
        let dt_normal = dev.derivative_base + lambda * dev.derivative_free;

        CurveDescription {
            tangent: dev.frame.tangent,
            dt_tangent: dev.frame.derivative,
            dt_normal,
            angle: None,
        }
    }
}

/// Steer the surface and ruling direction by defining an angle between the tangent and
/// ruling direction along the curve.
pub fn normal_and_angle_ode(
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
        // The direction of the ruling within the tangent plane gives us a unit vector based on
        // the coefficients hs, hc = target_angle.sin_cos() as
        //
        //     ruling = hc·unit_tangent + hs·derivative_free.
        //     dt_normal = a·unit_tangent + b·derivative_free
        //
        // Since the ruling is orthogonal to the tangent derivative this derives one linear
        // equation for the tangent derivative: <dt_normal, tangent> = -<normal, dt_tangent>
        //
        //      a = -<normal, frame.derivative> / ||tangent||
        //      b · hs = - a · hc
        //
        // b = -a / target_angle.tan()

        // angle(ruling, frame.tangent) = atan2(<normal, frame.derivative>, lambda)
        //
        // See `dc-integral/src/lib.rs` for the derivation of this formula where lambda is the
        // parameter from the above formula. Now let's derive that lambda. Note how we
        // automatically get `lambda = 0` at the direction discontinuity. Plus see that we get the
        // same result for `-target_angle` and the signum there chooses the orientation.
        let tan_angle = target_angle.tan();

        let angle_or_zero = if tan_angle == 0.0 {
            1e-9
        } else if tan_angle.abs() < 1e-9 {
            1e-9 * tan_angle.signum()
        } else {
            tan_angle
        };

        let a = -dev.normal.dot(frame.derivative) / frame.tangent.length();
        let b = -a / angle_or_zero;

        // ^ LLM anecdote: oneshot incorrectly. Previously believed to be correct though but it
        // multiplied instead of divided. Stupid machine, stupid me for trusting it too much and
        // getting myself confused in the process.
        //
        // It badly fumbled the derivation itself already, forgetting the square root in the
        // cos(x)-identity or forgetting that subtract `1` changes the numerator..
        let dt_normal = dev.derivative_base + b * dev.derivative_free;

        CurveDescription {
            tangent: dev.frame.tangent,
            dt_tangent: dev.frame.derivative,
            dt_normal,
            angle: Some(target_angle),
        }
    }
}

pub fn stitch(end: CurveSegment, curve: impl Curve) -> (CurveSegment, f32) {
    // Determine `parameter` such that the ruling match up (given that the normals will match
    // up as well).
    let parameter = {
        let start_frame = curve.at(0.0);
        let development = SurfaceDevelopment::from_frame_and_normal(start_frame, end.normal);

        let raw_angle = end.ruling.angle_between(development.frame.tangent);
        let signum = development
            .frame
            .tangent
            .cross(end.ruling)
            .dot(development.normal)
            .signum();

        signum * raw_angle - std::f64::consts::PI * 0.5
    };

    let ode = normal_and_angle_ode(curve, |_| parameter as f32);
    let basis = CurveSegment::initial(end.normal, ode);

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

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            third_derivative: DVec3::ZERO,
        }
    }
}
