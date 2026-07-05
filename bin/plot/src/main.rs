use std::io::{self, Read as _, Write as _};

use clap::Parser;
use dc_curves::{DenormalTangentFrame, SurfaceNormal, normal_and_angle_ode, stitch};
use dc_export::svg;
use dc_integral::{CurveSegment, curve_ode_with_curvature};
use glam::Vec3;
use miniserde::{Deserialize, Serialize};

/// Develop a curve according to parameterization passed as JSON
#[derive(Parser)]
#[clap(about, version)]
pub struct Args {
    /// If not provided or explicitly '-', the program will read a JSON parameterization from stdin. For the JSON schema, see source code. Coming soon.
    #[clap(long = "file", short = 'f')]
    parameterization: Option<String>,
    /// If not provided write the result to stdout.
    #[clap(short, long)]
    output: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct Parameterization {
    hermite: Vec<HermiteNode>,
    nodes: Vec<Node>,
    normal: [f32; 3],
    parameter: Option<Vec<SurfaceParameter>>,
    pipe: Option<PipeParameterization>,
}

#[derive(Debug, Clone, Deserialize)]
struct PipeParameterization {
    base: PipeBase,
    develop: Vec<PipeParameter>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut io = IoResources::new(&args)?;

    let mut cfg_data = String::new();
    io.input.read_to_string(&mut cfg_data)?;
    let input: Parameterization = miniserde::json::from_str(&cfg_data)?;

    let mut splines: Vec<Box<dyn dc_curves::Curve>> = Vec::new();

    for [a, b] in input.hermite.array_windows::<2>() {
        let spline = a.to_bezier(b);
        splines.push(Box::new(spline));
    }

    for extra in input.nodes {
        match extra {
            Node { spiral, hermite }
                if [spiral.is_some(), hermite.is_some()]
                    .into_iter()
                    .filter(|&x| x)
                    .count()
                    != 1 =>
            {
                Err("Node must be one of `spiral`, `hermite`")?;
            }
            Node {
                spiral: Some(params),
                ..
            } => {
                if let Some(last) = splines.last() {
                    let spiral = dc_curves::Spiral {
                        radius: params.radius,
                        pitch: params.pitch,
                    };

                    let affine = dc_curves::Affine::with_aligned(last.as_ref(), 1.0, spiral, 0.0);
                    splines.push(Box::new(affine));
                } else {
                    splines.push(Box::new(dc_curves::Spiral {
                        radius: params.radius,
                        pitch: params.pitch,
                    }));
                }
            }
            Node {
                hermite: Some(params),
                ..
            } => {
                let Some([a, b]) = params.array_windows::<2>().nth(0) else {
                    // Basically an empty spline.
                    continue;
                };

                let spline = a.to_bezier(b);
                let isometry = splines.last().map(|last| {
                    dc_curves::Affine::with_aligned(last.as_ref(), 1.0, spline, 0.0).isometry
                });

                for [a, b] in params.array_windows::<2>() {
                    let spline = a.to_bezier(b);
                    if let Some(isometry) = isometry {
                        splines.push(Box::new(dc_curves::Affine {
                            inner: spline,
                            isometry,
                        }));
                    } else {
                        splines.push(Box::new(spline));
                    }
                }
            }
            _ => unreachable!(),
        }
    }

    let obj = dc_export::obj::ObjConfig {
        tangent_scale: None,
        normalize_ruling: true,
        comment: Some(cfg_data),
        ..Default::default()
    };

    let results = if let Some(params) = input.parameter {
        let normal = SurfaceNormal::from_array(input.normal);
        extrapolate_curve(splines, normal, params, obj)?
    } else if let Some(params) = input.pipe {
        extrapolate_pipe(splines, params, obj)?
    } else {
        todo!("Proper error");
    };

    let json = miniserde::json::to_string(&results);
    io.output.write_all(json.as_bytes())?;

    Ok(())
}

fn extrapolate_curve(
    splines: Vec<Box<dyn dc_curves::Curve>>,
    normal: SurfaceNormal,
    parameterization: Vec<SurfaceParameter>,
    obj: dc_export::obj::ObjConfig,
) -> Result<Results, Box<dyn std::error::Error>> {
    let Some((first, tail)) = splines.split_first() else {
        panic!("At least one segment is required to form a curve.")
    };

    let ode = normal_and_angle_ode(first.as_ref(), |_| 0.0);
    let mut start = CurveSegment::initial(normal, ode);

    let mut segments: Vec<(DenormalTangentFrame, CurveSegment)> = vec![];

    for (idx, curve) in [first].into_iter().chain(tail).enumerate() {
        let first_parameter;
        if idx != 0 {
            let end_of_prev = start;
            (start, first_parameter) = stitch(end_of_prev, curve.as_ref());
            assert_eq!(start.normal, end_of_prev.normal);
            assert!(start.ruling.angle_between(end_of_prev.ruling) < 1.0e-6);
        } else {
            first_parameter = 0.0;
        }

        let parameter = SurfaceParameter::extract(&parameterization, idx);

        if idx == 0 {
            segments.push((curve.at(0.0), start));
        }

        // We pad with a final (t; 1.0] evaluation in the end if that isn't present.
        let final_param_if_not_1 = parameter.last().and_then(|p| {
            if p.loc < 1.0 {
                Some(SurfaceParameter { loc: 1.0, h: p.h })
            } else {
                None
            }
        });

        let first_parameter = SurfaceParameter {
            loc: 0.0,
            h: first_parameter,
        };

        let final_param_if_empty = parameter.is_empty().then(|| SurfaceParameter {
            loc: 1.0,
            h: first_parameter.h,
        });

        let mut iter = start;
        let mut ival_start = 0.0;
        let mut hval_start = first_parameter.h;

        let base_nodes: Vec<_> = core::iter::once(first_parameter)
            .chain(parameter.into_iter())
            .chain(final_param_if_not_1)
            .chain(final_param_if_empty)
            .collect();

        let interpolate_to_p01 = base_nodes.array_windows::<2>().flat_map(|&[start, end]| {
            let n = ((end.loc - start.loc) / 0.010) as usize;
            let inner = (1..=n).map(move |i| {
                let t = start.loc + i as f32 * (end.loc - start.loc) / n as f32;
                SurfaceParameter {
                    loc: t,
                    h: linear_interpolate((start.loc, end.loc), (start.h, end.h))(t),
                }
            });

            inner.chain(std::iter::once(end))
        });

        for SurfaceParameter { loc, h } in interpolate_to_p01 {
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
                normal_and_angle_ode(curve.as_ref(), hinterpolate),
                iter.normal,
                (iter.flat_position, iter.flat_direction),
                (f64::from(ival_start), f64::from(loc)),
            );

            segments.push((curve.at(loc), next));

            iter = next;
            ival_start = loc;
            hval_start = h;
        }

        start = iter;
    }

    let obj = obj.surface(&segments)?;
    let svg = svg::to_svg(&segments, obj.scale)?;

    Ok(Results {
        obj: obj.contents,
        svg: svg.contents,
    })
}

fn extrapolate_pipe(
    splines: Vec<Box<dyn dc_curves::Curve>>,
    pipe: PipeParameterization,
    obj: dc_export::obj::ObjConfig,
) -> Result<Results, Box<dyn std::error::Error>> {
    const ZERO_FLAT: dc_integral::PipeFaceBase = dc_integral::PipeFaceBase {
        base_left: glam::DVec2::ZERO,
        base_right: glam::DVec2::ZERO,
        orientation_left: 0.0,
        orientation_right: 0.0,
    };

    let Some((first, tail)) = splines.split_first() else {
        panic!("At least one segment is required to form a curve.")
    };

    let mut start = dc_integral::TrianglePipeBase {
        base1: Vec3::from_array(pipe.base.east).as_dvec3(),
        base2: Vec3::from_array(pipe.base.west).as_dvec3(),
        opposing_normal: Vec3::from_array(pipe.base.north_normal).as_dvec3(),
        flat1: ZERO_FLAT,
        flat2: ZERO_FLAT,
        opposing_flat: ZERO_FLAT,
    };

    let mut segments: Vec<(DenormalTangentFrame, dc_integral::TrianglePipeBase)> = vec![];

    for (idx, curve) in [first].into_iter().chain(tail).enumerate() {
        let develop = PipeParameter::extract(&pipe.develop, idx);

        if idx == 0 {
            segments.push((curve.at(0.0), start));
        }

        let first_parameter = PipeParameter {
            loc: 0.0,
            relative_speed_a: 1.2,
            relative_speed_b: 1.15,
            yaw: 2.0 * 3.12,
        };

        // We pad with a final (t; 1.0] evaluation in the end if that isn't present.
        let final_param_if_not_1 = develop.last().and_then(|p| {
            if p.loc < 1.0 {
                Some(PipeParameter {
                    loc: 1.0,
                    ..first_parameter
                })
            } else {
                None
            }
        });

        let final_param_if_empty = pipe.develop.is_empty().then(|| PipeParameter {
            loc: 1.0,
            ..first_parameter
        });

        let mut iter = start;
        let mut ival_start = 0.0;
        let mut hval_start = first_parameter;

        let base_nodes: Vec<_> = core::iter::once(first_parameter)
            .chain(develop.into_iter())
            .chain(final_param_if_not_1)
            .chain(final_param_if_empty)
            .collect();

        let interpolate_to_p01 = base_nodes.array_windows::<2>().flat_map(|&[start, end]| {
            let n = ((end.loc - start.loc) / 0.010) as usize;
            let inner = (1..=n).map(move |i| {
                let t = start.loc + i as f32 * (end.loc - start.loc) / n as f32;
                PipeParameter {
                    loc: t,
                    relative_speed_a: linear_interpolate(
                        (start.loc, end.loc),
                        (start.relative_speed_a, end.relative_speed_a),
                    )(t),
                    relative_speed_b: linear_interpolate(
                        (start.loc, end.loc),
                        (start.relative_speed_b, end.relative_speed_b),
                    )(t),
                    yaw: linear_interpolate((start.loc, end.loc), (start.yaw, end.yaw))(t),
                }
            });

            inner.chain(std::iter::once(end))
        });

        for p @ PipeParameter { loc, .. } in interpolate_to_p01 {
            // Skip duplicated nodes.
            if !(loc > ival_start) {
                continue;
            }

            assert!(
                loc > 0.0 && loc <= 1.0,
                "Parameter location must be in (0, 1]"
            );

            let hinterpolate = |h: f32| -> PipeParameter {
                PipeParameter {
                    loc: h,
                    relative_speed_a: linear_interpolate(
                        (ival_start, loc),
                        (hval_start.relative_speed_a, p.relative_speed_a),
                    )(h),
                    relative_speed_b: linear_interpolate(
                        (ival_start, loc),
                        (hval_start.relative_speed_b, p.relative_speed_b),
                    )(h),
                    yaw: linear_interpolate((ival_start, loc), (hval_start.yaw, p.yaw))(h),
                }
            };

            let next = dc_integral::triangle_pipe_ode(
                |at: f64| {
                    let frame = curve.at(at as f32);
                    let params = hinterpolate(at as f32);

                    dc_integral::PipeDescription {
                        frame,
                        len_a: frame.tangent.length() * f64::from(params.relative_speed_a),
                        len_b: frame.tangent.length() * f64::from(params.relative_speed_b),
                        yaw: f64::from(params.yaw),
                    }
                },
                iter,
                (f64::from(ival_start), f64::from(loc)),
            );

            segments.push((curve.at(loc), next.as_next()));

            iter = next.as_next();
            ival_start = loc;
            hval_start = p;
        }

        start = iter;
    }

    let obj = obj.pipe(&segments)?;
    let svg = svg::pipe(&segments, obj.scale)?;

    Ok(Results {
        obj: obj.contents,
        svg: svg.contents,
    })
}

struct IoResources {
    input: Box<dyn io::Read>,
    output: Box<dyn io::Write>,
}

impl IoResources {
    fn new(args: &Args) -> Result<Self, std::io::Error> {
        let input: Box<dyn io::Read> = match &args.parameterization {
            Some(path) if path != "-" => Box::new(std::fs::File::open(path)?),
            _ => Box::new(std::io::stdin()),
        };

        let output: Box<dyn io::Write> = match &args.output {
            Some(path) if path != "-" => Box::new(std::fs::File::create(path)?),
            _ => Box::new(std::io::stdout()),
        };

        Ok(Self { input, output })
    }
}

#[derive(Debug, Serialize)]
struct Results {
    obj: String,
    svg: String,
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

#[derive(Debug, Clone, Copy, Deserialize)]
struct SurfaceParameter {
    loc: f32,
    h: f32,
}

impl SurfaceParameter {
    fn extract(this: &[Self], idx: usize) -> Vec<Self> {
        this.iter()
            .filter_map(|p| {
                if p.loc as usize == idx {
                    Some(SurfaceParameter {
                        loc: p.loc.fract(),
                        h: p.h,
                    })
                } else {
                    None
                }
            })
            .collect()
    }
}

#[derive(Debug, Clone, Copy, Deserialize)]
struct PipeBase {
    east: [f32; 3],
    west: [f32; 3],
    north_normal: [f32; 3],
}

#[derive(Debug, Clone, Copy, Deserialize)]
struct PipeParameter {
    loc: f32,
    relative_speed_a: f32,
    relative_speed_b: f32,
    yaw: f32,
}

impl PipeParameter {
    fn extract(this: &[Self], idx: usize) -> Vec<Self> {
        this.iter()
            .filter_map(|p| {
                if p.loc as usize == idx {
                    Some(PipeParameter {
                        loc: p.loc.fract(),
                        relative_speed_a: p.relative_speed_a,
                        relative_speed_b: p.relative_speed_b,
                        yaw: p.yaw,
                    })
                } else {
                    None
                }
            })
            .collect()
    }
}

#[derive(Debug, Clone, Deserialize)]
struct HermiteNode {
    position: [f32; 3],
    tangent: [f32; 3],
}

#[derive(Debug, Clone, Deserialize)]
struct SpiralNode {
    radius: f32,
    pitch: f32,
}

#[derive(Debug, Clone, Deserialize)]
struct Node {
    spiral: Option<SpiralNode>,
    hermite: Option<Vec<HermiteNode>>,
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
