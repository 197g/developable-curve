use dc_curves::{
    Curve, DenormalTangentFrame, SurfaceNormal, normal_and_flat_ode, normal_and_tan_ode,
};

use dc_integral::{CurveSegment, curve_ode_with_curvature};

use dc_export::{svg, obj};

enum OdeParameterization {
    Derivative,
    Angle,
}

fn run_surface_along(
    curve: &(dyn Curve + '_),
    ts: &[f32],
    var: &[f32],
    initial: SurfaceNormal,
    ode: OdeParameterization,
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

    // FIXME: this should be done internally because `horizontal` is calculated by the derivative
    // of the normal of the segment, so this should involve the parameter choice `v`.
    let initial = CurveSegment {
        normal: initial.normal,
        horizontal: Default::default(),
        flat_position: Default::default(),
        flat_direction: 0.0,
        angle: 0.0,
    };

    // We might like a choice here, not always calculate the flattening
    let normals: Vec<CurveSegment> = ts
        .iter()
        .zip(var)
        .scan((initial, 0.0), |state, (&ts, &v)| {
            let (tmp0, tmp1);
            let ode_base = match ode {
                OdeParameterization::Derivative => {
                    tmp0 = normal_and_flat_ode(curve, |_| v);
                    &tmp0 as &dyn Fn(glam::Vec3, f32) -> _
                }
                OdeParameterization::Angle => {
                    tmp1 = normal_and_tan_ode(curve, |_| v);
                    &tmp1
                }
            };

            let normal_x0 = state.0.normal;
            let flat_x0 = (state.0.flat_position, state.0.flat_direction);
            let time_segment = (state.1, ts);

            let endpoint = curve_ode_with_curvature(ode_base, normal_x0, flat_x0, time_segment);

            state.0 = endpoint;
            state.1 = ts;

            Some(endpoint)
        })
        .collect();

    points.into_iter().zip(normals).collect()
}

fn main() -> Result<(), Box<dyn core::error::Error>> {
    let curve = dc_curves::Circle { radius: 1.0 };
    let ts = (0..=100)
        .map(|i| 2.0 * core::f32::consts::PI * i as f32 / 100.0)
        .collect::<Vec<_>>();

    {
        // A cylinder.
        let var = (0..=100).map(|_| 0.0).collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([1.0, 0.0, 0.0]);
        let surface =
            run_surface_along(&curve, &ts, &var, initial, OdeParameterization::Derivative);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = svg::to_svg(&surface)?;
        std::fs::write("/tmp/template-cylinder.svg", &svg)?;

        let obj = obj::to_obj(&surface)?;
        std::fs::write("/tmp/template-cylinder.obj", &obj)?;
    }

    {
        // A cone.
        println!("\n\n");
        let var = (0..=100).map(|_| 0.0).collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([0.96, 0.0, 0.2]);
        let surface =
            run_surface_along(&curve, &ts, &var, initial, OdeParameterization::Derivative);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = svg::to_svg(&surface)?;
        std::fs::write("/tmp/template-cone.svg", &svg)?;

        let obj = obj::to_obj(&surface)?;
        std::fs::write("/tmp/template-cone.obj", &obj)?;
    }

    {
        // Something more interesting.
        println!("\n\n");
        let var = (0..=100)
            .map(|i| (4.0 * core::f32::consts::PI * i as f32 / 100.0).sin() * 0.2)
            .collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([0.96, 0.0, 0.2]);
        let surface = run_surface_along(&curve, &ts, &var, initial, OdeParameterization::Angle);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = svg::to_svg(&surface)?;
        std::fs::write("/tmp/template-neat.svg", &svg)?;

        let obj = obj::to_obj(&surface)?;
        std::fs::write("/tmp/template-neat.obj", &obj)?;
    }

    for factor in [0.0, -0.3, -0.6] {
        // Something even more interesting.
        println!("\n\n");

        let ts = (0..=100).map(|i| i as f32 / 100.0).collect::<Vec<_>>();

        let var = (0..=100)
            .map(|i| {
                let rel = i as f32 / 100.0;
                (core::f32::consts::PI * rel).sin() * factor
            })
            .collect::<Vec<_>>();

        let c = 4.;
        let h = 2.;
        let curve = dc_curves::HermiteSpline::<4> {
            points: [
                glam::Vec3::from_array([c, c, -h]),
                glam::Vec3::from_array([-c, c, h]),
                glam::Vec3::from_array([-c, -c, -h]),
                glam::Vec3::from_array([c, -c, h]),
            ],
        };

        // The tangent here is [-2c, 0, 2h], it must be orthogonal to that.
        let initial = SurfaceNormal::from_array([0.5 * h, 4.0, 0.5 * c]);
        let surface =
            run_surface_along(&curve, &ts, &var, initial, OdeParameterization::Derivative);

        for (frame, segment) in &surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}, flat: {:.4?}",
                frame.base, frame.tangent, segment.normal, segment.flat_position,
            );
        }

        let svg = svg::to_svg(&surface)?;
        let name_particle = format!("_{factor}").replace('.', "p");
        std::fs::write(format!("/tmp/template-woah{name_particle}.svg"), &svg)?;

        let obj = obj::to_obj(&surface)?;
        std::fs::write(format!("/tmp/template-woah{name_particle}.obj"), &obj)?;
    }

    Ok(())
}
