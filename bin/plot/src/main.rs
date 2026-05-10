use std::io::{self, Read as _, Write as _};

use clap::Parser;
use dc_curves::{Curve as _, DenormalTangentFrame, normal_and_tan_ode, stitch};
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
    normal: [f32; 3],
    parameter: Vec<Parameter>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let mut io = IoResources::new(&args)?;

    let mut cfg_data = String::new();
    io.input.read_to_string(&mut cfg_data)?;
    let input: Parameterization = miniserde::json::from_str(&cfg_data)?;

    let mut splines = Vec::new();
    for [a, b] in input.hermite.array_windows::<2>() {
        let spline = a.to_bezier(b);
        splines.push(spline);
    }

    let Some((first, tail)) = splines.split_first() else {
        panic!("At least two Hermite nodes are required to form a curve.")
    };

    let ode = normal_and_tan_ode(&first, |_| 0.0);
    let mut start = CurveSegment::initial(Vec3::from_array(input.normal).as_dvec3(), ode);

    let mut segments: Vec<(DenormalTangentFrame, CurveSegment)> = vec![];

    for (idx, (curve, _node)) in [first]
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

        let parameter = Parameter::extract(&input.parameter, idx);

        if idx == 0 {
            segments.push((curve.at(0.0), start));
        }

        // We pad with a final (t; 1.0] evaluation in the end if that isn't present.
        let final_param_if_not_1 = parameter.last().and_then(|p| {
            if p.loc < 1.0 {
                Some(Parameter { loc: 1.0, h: p.h })
            } else {
                None
            }
        });

        let first_parameter = Parameter {
            loc: 0.0,
            h: first_parameter,
        };

        let final_param_if_empty = parameter.is_empty().then(|| Parameter {
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
                Parameter {
                    loc: t,
                    h: linear_interpolate((start.loc, end.loc), (start.h, end.h))(t),
                }
            });

            inner.chain(std::iter::once(end))
        });

        for Parameter { loc, h } in interpolate_to_p01 {
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
                (f64::from(ival_start), f64::from(loc)),
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
        comment: Some(cfg_data),
        ..Default::default()
    };

    let obj = obj.to_obj(&segments)?;
    let svg = svg::to_svg(&segments, obj.scale)?;

    let results = Results {
        obj: obj.contents,
        svg: svg.contents,
    };

    let json = miniserde::json::to_string(&results);
    io.output.write_all(json.as_bytes())?;

    Ok(())
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
struct Parameter {
    loc: f32,
    h: f32,
}

impl Parameter {
    fn extract(this: &[Self], idx: usize) -> Vec<Self> {
        this.iter()
            .filter_map(|p| {
                if p.loc as usize == idx {
                    Some(Parameter {
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

#[derive(Debug, Clone, Deserialize)]
struct HermiteNode {
    position: [f32; 3],
    tangent: [f32; 3],
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
