use fast_ode::Coord;
use glam::{Vec2, Vec3};

pub fn curve_ode(
    tangent: impl Fn(Vec3, f32) -> Vec3,
    base: Vec3,
    (start, end): (f32, f32),
) -> Vec3 {
    struct Ode<F: Fn(Vec3, f32) -> Vec3>(F);

    impl<F: Fn(Vec3, f32) -> Vec3> fast_ode::DifferentialEquation<3> for Ode<F> {
        fn ode_dot_y(&self, t: f64, y: &Coord<3>) -> (Coord<3>, bool) {
            let x = Vec3::from_array(y.0.map(|v| v as f32));
            let tangent = (self.0)(x, t as f32);
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

    Vec3::from_array(x1.0.map(|v| v as f32))
}

#[derive(Clone, Copy)]
pub struct CurveDescription {
    pub tangent: Vec3,
    pub dt_normal: Vec3,
    pub curvature: f32,
    pub speed: f32,
}

#[derive(Clone, Copy)]
pub struct CurveSegment {
    pub normal: Vec3,
    pub horizontal: Vec3,
    pub flat_position: Vec2,
    pub flat_direction: f32,
    pub angle: f32,
}

pub fn curve_ode_with_curvature(
    tangent: impl Fn(Vec3, f32) -> CurveDescription,
    base: Vec3,
    flat_base: (Vec2, f32),
    (start, end): (f32, f32),
) -> CurveSegment {
    struct Ode<F: Fn(Vec3, f32) -> CurveDescription>(F);

    impl<F: Fn(Vec3, f32) -> CurveDescription> fast_ode::DifferentialEquation<6> for Ode<F> {
        fn ode_dot_y(&self, t: f64, ty: &Coord<6>) -> (Coord<6>, bool) {
            let [x, y, z, _, _, _] = ty.0;
            let x = Vec3::new(x as f32, y as f32, z as f32);
            let descriptor = (self.0)(x, t as f32);

            let [_, _, _, _cx, _cy, k] = ty.0;
            let [x, y, z] = descriptor.dt_normal.to_array().map(f64::from);

            let speed = f64::from(descriptor.speed);
            // The unit speed curvature but `t` is not unit speed.
            let dkds = f64::from(descriptor.curvature) * speed;

            // k describes the current heading.
            let (my, mx) = k.sin_cos();
            let dt = [x, y, z, mx * speed, my * speed, dkds];

            (Coord(dt), true)
        }
    }

    let x0 = Coord::<6>({
        let [x, y, z] = base.to_array().map(f64::from);
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
    let normal = Vec3::from_array([x, y, z].map(|v| v as f32));

    // The horizontal must be perpendicular to the plane normal and its derivative.
    // We are however free to choose a direction, let us pick a consistent one.
    let end_descriptor = (ode.0)(normal, end);
    let horizontal = end_descriptor.dt_normal.cross(normal);

    // There is probably a cheaper way to get this, do not pass the whole frame. Or do we?
    let signum = horizontal
        .cross(end_descriptor.tangent)
        .dot(normal)
        .signum();

    // Note: `<horizontal, frame.tangent> = v · ||frame.tangent||`
    //
    // if you want to control this angle. Expanded:
    //
    // cos(horizontal, frame.tangent) · ||horizontal|| · ||frame.tangent||
    //     = <horizontal, frame.tangent>
    //     = v · ||frame.tangent||
    //
    // v = cos(horizontal, frame.tangent) · ||horizontal||
    //     = cos(horizontal, frame.tangent) · ||dt normal||
    //
    // So `v = 0` for a guaranteed tangent-perpendicular horizontal line.

    // I would prefer an acos2 with semantics
    //     (cos(a)·||A||·||B||, ||A||·||B||) -> arccos(a)
    // but this is good enough for now.
    let angle = horizontal.angle_between(end_descriptor.tangent) * signum / std::f32::consts::PI;

    // Make this the right-hand coordinate system instead (tangent, horizontal, normal). This makes
    // it compatible with the curvature calculation.
    let horizontal = horizontal * signum;

    CurveSegment {
        normal,
        horizontal,
        flat_position: Vec2::from_array([fx, fy].map(|v| v as f32)),
        // We do not build a full frame..
        flat_direction: k as f32,
        angle,
    }
}
