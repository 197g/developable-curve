use dc_integral::curve_ode;
use dc_theory::{Curve, DenormalTangentFrame, SurfaceDevelopment, SurfaceNormal};

fn run_surface_along(
    curve: &(dyn Curve + '_),
    ts: &[f32],
    var: &[f32],
    initial: SurfaceNormal,
) -> Vec<(DenormalTangentFrame, SurfaceNormal)> {
    let base = curve.at(0.0);

    // The initial surface normal must be orthogonal otherwise chaos.. Let's be loud.
    assert!(
        initial.normal.dot(base.tangent) < 1e-6,
        "Initial surface normal must be orthogonal to the tangent frame."
    );

    let points = curve.sample_at(ts);
    // now we scan and do integration steps along the way.
    assert_eq!(points.len(), var.len());

    let normals: Vec<SurfaceNormal> = ts
        .iter()
        .zip(var)
        .scan((initial, 0.0), |state, (&ts, &v)| {
            let callback = SurfaceDevelopment::ode_integrator(curve, |_| v);
            let end = curve_ode(callback, state.0.normal, (state.1, ts));
            Some(SurfaceNormal { normal: end })
        })
        .collect();

    points.into_iter().zip(normals).collect()
}

fn main() {
    let curve = dc_theory::Circle { radius: 1.0 };
    let ts = (0..=100)
        .map(|i| 2.0 * core::f32::consts::PI * i as f32 / 100.0)
        .collect::<Vec<_>>();

    {
        // A cylinder.
        let var = (0..=100).map(|_| 0.0).collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([0.0, 0.0, 1.0]);
        let surface = run_surface_along(&curve, &ts, &var, initial);

        for (frame, normal) in surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}",
                frame.base, frame.tangent, normal.normal
            );
        }
    }

    {
        // A cone.
        println!("\n\n");
        let var = (0..=100).map(|_| 0.0).collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([-0.2, 0.0, 0.96]);
        let surface = run_surface_along(&curve, &ts, &var, initial);

        for (frame, normal) in surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}",
                frame.base, frame.tangent, normal.normal
            );
        }
    }

    {
        // Something more interesting.
        println!("\n\n");
        let var = (0..=100)
            .map(|i| (2.0 * core::f32::consts::PI * i as f32 / 100.0).sin() * 0.1)
            .collect::<Vec<_>>();

        let initial = SurfaceNormal::from_array([-0.2, 0.0, 0.96]);
        let surface = run_surface_along(&curve, &ts, &var, initial);

        for (frame, normal) in surface {
            println!(
                "base: {:.4?}, tangent: {:.4?}, normal: {:.4?}",
                frame.base, frame.tangent, normal.normal
            );
        }
    }
}
