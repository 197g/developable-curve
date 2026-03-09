use fast_ode::Coord;
use glam::Vec3;

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
