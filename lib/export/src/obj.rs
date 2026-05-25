use core::fmt::Write as _;

use dc_integral::CurveSegment;
use dc_theory::DenormalTangentFrame;

#[derive(Debug, Clone)]
pub struct ObjConfig {
    pub tangent_scale: Option<f32>,
    pub normal_scale: f32,
    pub normalize_ruling: bool,
    /// Rescales the model to a 3d printer build plate.
    pub buildplate_mm: f32,
    pub comment: Option<String>,
}

impl Default for ObjConfig {
    fn default() -> Self {
        Self {
            tangent_scale: Some(1.0),
            normal_scale: 1.0,
            normalize_ruling: false,
            buildplate_mm: 180.,
            comment: None,
        }
    }
}

pub fn to_obj(
    curve: &[(DenormalTangentFrame, CurveSegment)],
) -> Result<super::StrFileData, core::fmt::Error> {
    ObjConfig::default().to_obj(curve)
}

impl ObjConfig {
    pub fn to_obj(
        self,
        curve: &[(DenormalTangentFrame, CurveSegment)],
    ) -> Result<super::StrFileData, core::fmt::Error> {
        to_obj_with(curve, self)
    }
}

pub fn to_obj_with(
    curve: &[(DenormalTangentFrame, CurveSegment)],
    cfg: ObjConfig,
) -> Result<super::StrFileData, core::fmt::Error> {
    let mut string = String::new();

    if let Some(cmt) = cfg.comment {
        for line in cmt.lines() {
            writeln!(string, "# {line}")?;
        }
    }

    let (min, max) = curve.iter().fold(
        ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
        |(min, max), (frame, _)| {
            let [x, y, z] = frame.base.to_array();
            let min = [
                min[0].min(x - 1.0),
                min[1].min(y - 1.0),
                min[2].min(y - 1.0),
            ];
            let max = [
                max[0].max(x + 1.0),
                max[1].max(y + 1.0),
                max[2].max(z + 1.0),
            ];
            (min, max)
        },
    );

    let model_bounds = (max[0] - min[0]).max(max[1] - min[1]).max(max[2] - min[2]);

    let tangent_scale = f64::from(cfg.tangent_scale.unwrap_or(0.5));
    let normal_scale = f64::from(cfg.normal_scale);

    // Print optimized: 180mm build plate.
    let ruling_scale = 1.0f64;
    let model_scale = f64::from(cfg.buildplate_mm) / model_bounds;

    let write_frame = cfg.tangent_scale.is_some();
    let vertices_per_frame = if write_frame { 4 } else { 2 };

    for (frame, segment) in curve {
        let [tx, ty, tz] = frame.tangent.to_array().map(|x| x * tangent_scale);
        let [nx, ny, nz] = segment.normal.axis.to_array().map(|x| x * normal_scale);

        let radius_of_curvature = (1.1 * segment.flat_curvature).max(1.0).recip();

        let [hx, hy, hz] = if cfg.normalize_ruling {
            segment
                .ruling
                .normalize_or_zero()
                .to_array()
                .map(|x| x * ruling_scale.min(radius_of_curvature))
        } else {
            segment.ruling.to_array().map(|x| x * ruling_scale)
        };

        writeln!(
            &mut string,
            "v {:.4} {:.4} {:.4}",
            frame.base.x * model_scale,
            frame.base.y * model_scale,
            frame.base.z * model_scale
        )?;

        writeln!(
            &mut string,
            "v {:.4} {:.4} {:.4}",
            (frame.base.x + hx) * model_scale,
            (frame.base.y + hy) * model_scale,
            (frame.base.z + hz) * model_scale,
        )?;

        if write_frame {
            writeln!(
                &mut string,
                "v {:.4} {:.4} {:.4}",
                (frame.base.x + tx) * model_scale,
                (frame.base.y + ty) * model_scale,
                (frame.base.z + tz) * model_scale,
            )?;

            writeln!(
                &mut string,
                "v {:.4} {:.4} {:.4}",
                frame.base.x * model_scale + nx,
                frame.base.y * model_scale + ny,
                frame.base.z * model_scale + nz
            )?;
        }
    }

    // Assign some vertex normals, all vertices of a frame pointing in the plane normal direction.
    // So we only have one normal per frame.
    for (_, segment) in curve {
        let [nx, ny, nz] = segment.normal.axis.to_array();
        writeln!(&mut string, "vn {:.4} {:.4} {:.4}", nx, ny, nz)?;
    }

    // And connect all vertices of a frame to the next one with quads.
    for i in 0..curve.len() - 1 {
        let normal = i + 1; // 1-based indexing for .obj
        let next_normal = normal + 1;
        let base = i * vertices_per_frame + 1;

        // Let's make sure this is the right way up. Each ruling, tangent, normal triple is
        // right handed in that order.
        writeln!(
            &mut string,
            "f {}//{normal} {}//{normal} {}//{next_normal} {}//{next_normal}",
            base + 1,
            base,
            base + vertices_per_frame,
            base + 1 + vertices_per_frame,
        )?;
    }

    if write_frame {
        // And line elements for the tangents, normals and ruling.
        for i in 0..curve.len() {
            let base = i * vertices_per_frame + 1;
            writeln!(&mut string, "l {} {}", base, base + 2)?; // tangent
            writeln!(&mut string, "l {} {}", base, base + 3)?; // normal
        }
    }

    Ok(super::StrFileData {
        contents: string,
        scale: model_scale as f32,
    })
}
