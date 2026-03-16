use std::io::Write as _;

use dc_curves::{Curve as _, DenormalTangentFrame, normal_and_tan_ode, stitch};
use dc_integral::{CurveSegment, curve_ode_with_curvature};
use glam::Vec3;

#[derive(Debug, Clone)]
struct Parameterization {
    hermite: Vec<HermiteNode>,
    normal: Vec3,
}

#[derive(Debug, Clone, Copy)]
struct Parameter {
    r#where: f32,
    h: f32,
}

#[derive(Debug, Clone)]
struct HermiteNode {
    position: [f32; 3],
    tangent: [f32; 3],
    /// Hm, should we have on `parameter` on `Parameterization` where `loc` has an integer /
    /// fractional part instead?
    parameter: Vec<Parameter>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input: Parameterization = read_parameterization();

    let mut splines = Vec::new();
    for [a, b] in input.hermite.array_windows::<2>() {
        let spline = a.to_bezier(b);
        splines.push(spline);
    }

    let Some((first, tail)) = splines.split_first() else {
        panic!("At least two Hermite nodes are required to form a curve.")
    };

    let ode = normal_and_tan_ode(&first, |_| 0.0);
    let mut start = CurveSegment::initial(input.normal, ode);

    let mut segments: Vec<(DenormalTangentFrame, CurveSegment)> = vec![];

    for (idx, (curve, node)) in [first]
        .into_iter()
        .chain(tail)
        .zip(&input.hermite)
        .enumerate()
    {
        let first_parameter;
        if idx != 0 {
            let end_of_prev = start;
            (start, first_parameter) = stitch(end_of_prev, &curve);
            assert_eq!(start.normal, end_of_prev.normal);
            assert!(start.horizontal.angle_between(end_of_prev.horizontal) < 1.0e-6);
        } else {
            first_parameter = 0.0;
        }

        segments.push((curve.at(0.0), start));
        // We pad with a final (t; 1.0] evaluation in the end if that isn't present.
        let final_param_if_not_1 = node.parameter.last().and_then(|p| {
            if p.r#where < 1.0 {
                Some(Parameter {
                    r#where: 1.0,
                    h: p.h,
                })
            } else {
                None
            }
        });

        let first_parameter = Parameter {
            r#where: 0.0,
            h: first_parameter,
        };

        let final_param_if_empty = node.parameter.is_empty().then(|| Parameter {
            r#where: 1.0,
            h: first_parameter.h,
        });

        let mut iter = start;
        let mut ival_start = 0.0;
        let mut hval_start = first_parameter.h;

        let base_nodes: Vec<_> = core::iter::once(first_parameter)
            .chain(node.parameter.iter().copied())
            .chain(final_param_if_not_1)
            .chain(final_param_if_empty)
            .collect();

        let interpolate_to_p01 = base_nodes.array_windows::<2>().flat_map(|&[start, end]| {
            let n = ((end.r#where - start.r#where) / 0.01) as usize;
            let inner = (1..=n).map(move |i| {
                let t = start.r#where + i as f32 * (end.r#where - start.r#where) / n as f32;
                Parameter {
                    r#where: t,
                    h: linear_interpolate((start.r#where, end.r#where), (start.h, end.h))(t),
                }
            });

            std::iter::once(start)
                .chain(inner)
                .chain(std::iter::once(end))
        });

        for Parameter { r#where: loc, h } in interpolate_to_p01 {
            // Skip duplicated nodes.
            if !(loc > ival_start) {
                continue;
            }

            assert!(
                loc > 0.0 && loc <= 1.0,
                "Parameter location must be in (0, 1]"
            );

            let hinterpolate = linear_interpolate((ival_start, loc), (hval_start, h));
            let next = curve_ode_with_curvature(
                normal_and_tan_ode(&curve, hinterpolate),
                iter.normal,
                (iter.flat_position, iter.flat_direction),
                (ival_start, loc),
            );

            segments.push((curve.at(loc), next));

            iter = next;
            ival_start = loc;
            hval_start = h;
        }

        start = iter;
    }

    let obj = dc_export::obj::ObjConfig {
        tangent_scale: None,
        normalize_horizontal: true,
        ..Default::default()
    };

    let obj = obj.to_obj(&segments)?;
    std::io::stdout().write_all(&obj)?;

    Ok(())
}

fn read_parameterization() -> Parameterization {
    // For demonstration purposes, we will create a simple parameterization with two Hermite nodes.
    Parameterization {
        hermite: vec![
            HermiteNode {
                position: [0.0, 0.0, 0.0],
                tangent: [1.0, 0.0, 0.0],
                parameter: vec![Parameter {
                    r#where: 0.5,
                    h: 0.2,
                }],
            },
            HermiteNode {
                position: [1.0, 0.0, 0.2],
                tangent: [1.0, 0.1, 0.0],
                parameter: vec![Parameter {
                    r#where: 0.5,
                    h: 0.0,
                }],
            },
            HermiteNode {
                position: [2.0, 0.0, 0.0],
                tangent: [1.0, 0.1, 0.0],
                parameter: vec![Parameter {
                    r#where: 0.5,
                    h: -0.16,
                }],
            },
            HermiteNode {
                position: [3.0, 0.0, 0.2],
                tangent: [1.0, 0.2, 0.0],
                parameter: vec![Parameter {
                    r#where: 0.5,
                    h: 0.0,
                }],
            },
        ],
        normal: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn linear_interpolate(ival: (f32, f32), hval: (f32, f32)) -> impl Fn(f32) -> f32 {
    move |t| {
        let (x0, x1) = ival;
        let (h0, h1) = hval;

        if x1 == x0 {
            return h0; // Avoid division by zero; return h0 as a fallback.
        }

        h0 + (h1 - h0) * ((t - x0) / (x1 - x0))
    }
}

impl HermiteNode {
    fn to_bezier(&self, next: &HermiteNode) -> dc_curves::BezierSpline<4> {
        let p0 = self.position;
        let p1 = next.position;
        let t0 = self.tangent;
        let t1 = next.tangent;

        let points = [
            Vec3::from(p0),
            Vec3::from(p0) + Vec3::from(t0) / 3.0,
            Vec3::from(p1) - Vec3::from(t1) / 3.0,
            Vec3::from(p1),
        ];

        dc_curves::BezierSpline { points }
    }
}
