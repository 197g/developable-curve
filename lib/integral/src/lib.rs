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
    pub dt_normal: Vec3,
    pub curvature: f32,
    pub speed: f32,
}

#[derive(Clone, Copy)]
pub struct CurveSegment {
    pub normal: Vec3,
    pub flat_position: Vec2,
    pub flat_direction: f32,
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
            let (mx, my) = k.sin_cos();
            let dt = [x, y, z, mx * speed, my * speed, dkds];

            (Coord(dt), true)
        }
    }

    let x0 = Coord::<6>({
        let [x, y, z] = base.to_array().map(f64::from);
        let [cx, cy] = flat_base.0.to_array().map(f64::from);
        [x, y, z, cx, cy, f64::from(flat_base.1)]
    });

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

    let Coord([x, y, z, fx, fy, k]) = x1;

    CurveSegment {
        normal: Vec3::from_array([x, y, z].map(|v| v as f32)),
        flat_position: Vec2::from_array([fx, fy].map(|v| v as f32)),
        // We do not build a full frame..
        flat_direction: k as f32,
    }
}
