mod svg;

use dc_integral::{CurveDescription, CurveSegment, curve_ode_with_curvature};
use dc_theory::{Curve, DenormalTangentFrame, SurfaceDevelopment, SurfaceNormal};

fn run_surface_along(
    curve: &(dyn Curve + '_),
    ts: &[f32],
    var: &[f32],
    initial: SurfaceNormal,
) -> Vec<(DenormalTangentFrame, CurveSegment)> {
    let base = curve.at(0.0);

    // The initial surface normal must be orthogonal otherwise chaos.. Let's be loud.
    assert!(
        initial.normal.dot(base.tangent) < 1e-6,
        "Initial surface normal must be orthogonal to the tangent frame."
    );

    let points = curve.sample_at(ts);
    // now we scan and do integration steps along the way.
    assert_eq!(points.len(), var.len());

    let initial = CurveSegment {
        normal: initial.normal,
        flat_position: Default::default(),
        flat_direction: 0.0,
    };

    // We might like a choice here, not always calculate the flattening
    let normals: Vec<CurveSegment> = ts
        .iter()
        .zip(var)
        .scan((initial, 0.0), |state, (&ts, &v)| {
            let ode_base = SurfaceDevelopment::normal_and_flat_ode(curve, |_| v);

            let callback = move |pos, t| {
                let (dt_normal, curvature, speed) = ode_base(pos, t);

                CurveDescription {
                    dt_normal,
                    curvature,
                    speed,
                }
            };

            let normal_x0 = state.0.normal;
            let flat_x0 = (state.0.flat_position, state.0.flat_direction);
            let time_segment = (state.1, ts);
            let endpoint = curve_ode_with_curvature(callback, normal_x0, flat_x0, time_segment);

            state.0 = endpoint;
            state.1 = ts;

            Some(endpoint)
        })
        .collect();

    points.into_iter().zip(normals).collect()
}

fn main() -> Result<(), Box<dyn core::error::Error>> {
    let curve = dc_theory::Circle { radius: 1.0 };
    let ts = (0..=100)
        .map(|i| 2.0 * core::f32::consts::PI * i as f32 / 100.0)
        .collect::<Vec<_>>();

    {
        // A cylinder.
        let var = (0..=100).map(|_| 0.0).collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([1.0, 0.0, 0.0]);
        let surface = run_surface_along(&curve, &ts, &var, initial);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = crate::svg::to_svg(&surface)?;
        std::fs::write("/tmp/template-cylinder.svg", &svg)?;
    }

    {
        // A cone.
        println!("\n\n");
        let var = (0..=100).map(|_| 0.0).collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([0.96, 0.0, 0.2]);
        let surface = run_surface_along(&curve, &ts, &var, initial);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = crate::svg::to_svg(&surface)?;
        std::fs::write("/tmp/template-cone.svg", &svg)?;
    }

    {
        // Something more interesting.
        println!("\n\n");
        let var = (0..=100)
            .map(|i| (4.0 * core::f32::consts::PI * i as f32 / 100.0).sin() * 0.3)
            .collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([0.96, 0.0, 0.2]);
        let surface = run_surface_along(&curve, &ts, &var, initial);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = crate::svg::to_svg(&surface)?;
        std::fs::write("/tmp/template-neat.svg", &svg)?;
    }

    {
        // Something even more interesting.
        println!("\n\n");

        let ts = (0..=100).map(|i| i as f32 / 100.0).collect::<Vec<_>>();

        let var = (0..=100)
            .map(|i| (1.0 * core::f32::consts::PI * i as f32 / 100.0).sin() * 8.)
            .collect::<Vec<_>>();

        let c = 4.;
        let h = 2.;
        let curve = dc_theory::HermiteSpline::<4> {
            points: [
                glam::Vec3::from_array([c, c, -h]),
                glam::Vec3::from_array([-c, c, h]),
                glam::Vec3::from_array([-c, -c, -h]),
                glam::Vec3::from_array([c, -c, h]),
            ],
        };

        // The tangent here is [-2c, 0, 2h], it must be orthogonal to that.
        let initial = SurfaceNormal::from_array([0.5*h, 1.0, 0.5*c]);
        let surface = run_surface_along(&curve, &ts, &var, initial);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = crate::svg::to_svg(&surface)?;
        std::fs::write("/tmp/template-woah.svg", &svg)?;
    }

    Ok(())
}
