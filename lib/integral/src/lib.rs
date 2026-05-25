use fast_ode::Coord;
use glam::{DVec2, DVec3};

use dc_theory::{CurveDescription, SurfaceNormal};

pub fn curve_ode(
    tangent: impl Fn(DVec3, f64) -> DVec3,
    base: DVec3,
    (start, end): (f64, f64),
) -> DVec3 {
    struct Ode<F: Fn(DVec3, f64) -> DVec3>(F);

    impl<F: Fn(DVec3, f64) -> DVec3> fast_ode::DifferentialEquation<3> for Ode<F> {
        fn ode_dot_y(&self, t: f64, y: &Coord<3>) -> (Coord<3>, bool) {
            let x = DVec3::from_array(y.0);
            let tangent = (self.0)(x, t);
            (Coord(tangent.to_array().map(f64::from)), true)
        }
    }

    let x0 = Coord::<3>(base.to_array().map(f64::from));

    let sol = fast_ode::solve_ivp(
        &Ode(tangent),
        (f64::from(start), f64::from(end)),
        x0,
        |_, _| true,
        1e-6,
        1e-6,
    );

    let x1 = match sol {
        fast_ode::IvpResult::FinalTimeReached(coord) => coord,
        // Deserves a warning, at least!
        fast_ode::IvpResult::StepTooSmall(_, coord) => coord,
        fast_ode::IvpResult::OdeRequestedExit(..)
        | fast_ode::IvpResult::CallbackRequestedExit(..) => {
            unreachable!("we do not request exit")
        }
    };

    DVec3::from_array(x1.0)
}

#[derive(Clone, Copy)]
pub struct CurveSegment {
    pub normal: SurfaceNormal,
    pub ruling: DVec3,
    pub flat_position: DVec2,
    pub flat_direction: f64,
    pub flat_curvature: f64,
    pub angle: f64,
}

fn warn_nonzero(actual: f64, what: &str) {
    if !(actual < 1e-8) {
        eprintln!("{what}: {actual}");
    }
}

impl CurveSegment {
    pub fn initial(normal: SurfaceNormal, ode: impl Fn(DVec3, f64) -> CurveDescription) -> Self {
        let descriptor = ode(normal.axis, 0.0);

        if let Some(angle) = descriptor.angle {
            Self::from_angle(normal.axis, descriptor, angle)
        } else {
            Self::from_parameter_with_unstable_angle_at_zero(normal.axis, descriptor)
        }
    }

    fn from_angle(normal: DVec3, frame: CurveDescription, target_angle: f64) -> Self {
        let forward = frame.tangent.normalize();
        let sideways = normal.cross(forward);
        let (s, c) = target_angle.sin_cos();

        let ruling = c * forward + s * sideways;
        let angle = target_angle;

        warn_nonzero(ruling.dot(frame.dt_normal), "ruling to dt normal");
        warn_nonzero(frame.tangent.dot(normal), "normal to tangent");
        warn_nonzero(ruling.dot(normal), "normal to ruling");

        assert!(
            (normal.dot(frame.dt_tangent) + frame.dt_normal.dot(frame.tangent)).abs() < 1e-8,
            "{:.8} / {:.8}",
            normal.dot(frame.dt_tangent),
            frame.dt_normal.dot(frame.tangent),
        );

        CurveSegment {
            normal: SurfaceNormal { axis: normal },
            ruling,
            flat_position: Default::default(),
            flat_direction: Default::default(),
            flat_curvature: 0.0,
            angle,
        }
    }

    fn from_parameter_with_unstable_angle_at_zero(
        normal: DVec3,
        end_descriptor: CurveDescription,
    ) -> Self {
        // The ruling is orthogonal to both.
        let pre_ruling = normal.cross(end_descriptor.dt_normal);

        // Solve (dt pre_ruling) × tangent = 0
        //   or (dt unit(pre_ruling)) × tangent = 0
        //
        // Observe dt (unit(pre_ruling)×tangent)
        //   = (dt unit(pre_ruling))×tangent + unit(pre_ruling)×(dt tangent)
        //
        // FIXME: something not right here. Choosing `v = 0` does indeed generate a cone for which
        // the derivative of the ruling (up to the tip) is `-tangent`. Indeed any skew cone can
        // be developed since the derivative of the ruling is always in the plane spanned by
        // the tangent and the ruling.
        //
        // We have (dt unit(pre_ruling))×tangent = 0 iff
        //   dt (unit(pre_ruling)×tangent) = unit(pre_ruling)×(dt tangent)
        //
        // On the left: dt (||tangent|| normal)
        //   = (dt ||tangent||) normal + ||tangent|| dt normal
        //
        // On the right: unit(pre_ruling)×(dt tangent)
        //   = unit(pre_ruling)×(dt ||tangent|| unit(frame.tangent) + ||tangent|| frame.normal)
        //   = (dt ||tangent||) normal + ||tangent|| unit(pre_ruling)×frame.normal
        //
        // So we have … iff ||tangent|| dt normal = ||tangent|| unit(pre_ruling)×frame.normal
        //   iff dt normal = unit(pre_ruling)×frame.normal
        //   iff dt normal = unit(|sign|·normal×dt normal)×frame.normal
        //   iff dt normal = |sign|(normal×dt normal)×frame.normal / ||dt normal||
        //   iff dt normal = |sign|((normal·frame.normal) dt normal - (dt normal·frame.normal) normal) / ||dt normal||
        //
        // Note that for dt normal = <normal, frame.derivative>·unit(frame.tangent)
        //   dt normal · frame.normal = frame.tangent·frame.normal = 0
        //   ||dt normal|| = |<normal, frame.derivative>| = |<normal, frame.normal>|
        //
        // FIXME: ugh, there is this sign on the lhs and an absolute on the RHS. wat. Does this
        // next step properly get rid of it?

        // There is probably a cheaper way to get this, do not pass the whole frame. Or do we?
        let signum = end_descriptor
            .tangent
            .cross(pre_ruling)
            .dot(normal)
            .signum();

        // Note: `<ruling, frame.tangent> = v · ||frame.tangent||`
        //
        // if you want to control this angle. Expanded:
        //
        // cos(ruling, frame.tangent) · ||ruling|| · ||frame.tangent||
        //     = <ruling, frame.tangent>
        //     = v · ||frame.tangent||
        //
        // v = cos(ruling, frame.tangent) · ||ruling||
        //     = cos(ruling, frame.tangent) · ||dt normal||
        //     = cos(ruling, frame.tangent) · sqrt(v² + <normal, frame.derivative>²)
        //
        // angle(ruling, frame.tangent) = acos(v / sqrt(v² + <normal, frame.derivative>²))
        //     = atan(|<normal, frame.derivative>| / v) ; by acos(x) = atan(sqrt(1-x²)/x)
        //     = atan2(|<normal, frame.derivative>|, v)
        //
        // So we have a discontinuity. If the normal is perpendicular to `frame.derivative` then the
        // `cos(angle) = +1/-1` so the two are parallel with no steering at all. Otherwise, we can
        // choose `v = 0` for a guaranteed tangent-perpendicular ruling line or any other
        // non-parallel angle with appropriate `v`.
        //
        // So now you're asking, can we control `v` so that the discontinuity never occurs? Not in
        // general if the frame.derivative is discontinuous. But also consider this an artifact of our
        // choice of ruling, the direction of which is discontinuous at the zero of `dt_normal`.
        // And indeed at the same point we get a length of `0|v=0`. So really we should maybe instead
        // be steering by the angle; and then calculating a corresponding `v` while having `v=0` and
        // using our angle regardless at the discontinuity?

        // I would prefer an acos2 with semantics
        //     (cos(a)·||A||·||B||, ||A||·||B||) -> arccos(a)
        // but this is good enough for now–we get to manually do the atan transform and it works out
        // cleaner.
        //
        // FIXME: if we were handed `derivative_free` we could avoid this angle calculation and
        // probably the signum itself, too. We could calculate `v` and the rest of this would fall out
        // from atan2. But also if we were handed the angle then we could avoid the ill-defined
        // calculation for that point entirely. Maybe having the angle as a parameter to the ODE is
        // better after all and calculate ruling as rotateAround(normal, angle).rotate(tangent)
        // which I assume should itself simplify (TBD).

        let ruling = pre_ruling * signum;
        // `angle_between` measures absolute angle and we want a signed one.
        let angle = ruling.angle_between(end_descriptor.tangent) * signum;

        CurveSegment {
            normal: SurfaceNormal { axis: normal },
            ruling,
            flat_position: Default::default(),
            flat_direction: Default::default(),
            flat_curvature: 0.0,
            angle,
        }
    }
}

pub fn curve_ode_with_curvature(
    tangent: impl Fn(DVec3, f64) -> CurveDescription,
    base: SurfaceNormal,
    flat_base: (DVec2, f64),
    (start, end): (f64, f64),
) -> CurveSegment {
    struct Ode<F: Fn(DVec3, f64) -> CurveDescription>(F);

    impl<F: Fn(DVec3, f64) -> CurveDescription> fast_ode::DifferentialEquation<6> for Ode<F> {
        fn ode_dot_y(&self, t: f64, ty: &Coord<6>) -> (Coord<6>, bool) {
            let [x, y, z, _, _, _] = ty.0;
            let normal = DVec3::new(x, y, z);
            let descriptor = (self.0)(normal, t);

            let [_, _, _, _cx, _cy, k] = ty.0;
            let [x, y, z] = descriptor.dt_normal.to_array().map(f64::from);

            let speed = f64::from(descriptor.tangent.length());
            let curvature = descriptor.curvature_to_normal(normal);
            // The unit speed curvature but `t` is not unit speed.
            let dkds = f64::from(curvature) * speed;

            // k describes the current heading.
            let (my, mx) = k.sin_cos();
            // Since the 2d curve does not have unit speed either, it must be adjusted itself, too.

            let dt = [x, y, z, mx * speed, my * speed, dkds];

            (Coord(dt), true)
        }
    }

    let x0 = Coord::<6>({
        let [x, y, z] = base.axis.to_array().map(f64::from);
        let [cx, cy] = flat_base.0.to_array().map(f64::from);
        [x, y, z, cx, cy, f64::from(flat_base.1)]
    });

    let ode = Ode(tangent);

    let sol = fast_ode::solve_ivp(
        &ode,
        (f64::from(start), f64::from(end)),
        x0,
        |_, _| true,
        1e-6,
        1e-6,
    );

    let x1 = match sol {
        fast_ode::IvpResult::FinalTimeReached(coord) => coord,
        // Deserves a warning, at least!
        fast_ode::IvpResult::StepTooSmall(_, coord) => coord,
        fast_ode::IvpResult::OdeRequestedExit(..)
        | fast_ode::IvpResult::CallbackRequestedExit(..) => {
            unreachable!("we do not request exit")
        }
    };

    let Coord([x, y, z, fx, fy, k]) = x1;
    let normal = DVec3::from_array([x, y, z]);

    // The ruling must be perpendicular to the plane normal and its derivative.
    // We are however free to choose a direction, let us pick a consistent one.
    let end_descriptor = (ode.0)(normal, end);

    let basis = if let Some(target_angle) = end_descriptor.angle {
        CurveSegment::from_angle(normal, end_descriptor, target_angle)
    } else {
        CurveSegment::from_parameter_with_unstable_angle_at_zero(normal, end_descriptor)
    };

    CurveSegment {
        flat_position: DVec2::from_array([fx, fy]),
        // We do not build a full frame..
        flat_direction: k,
        flat_curvature: end_descriptor.curvature_to_normal(normal),
        ..basis
    }
}
