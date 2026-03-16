use core::fmt::Write as _;

use dc_integral::CurveSegment;
use dc_theory::DenormalTangentFrame;

#[derive(Debug, Clone, Copy)]
pub struct ObjConfig {
    pub tangent_scale: Option<f32>,
    pub normal_scale: f32,
    pub normalize_horizontal: bool,
}

impl Default for ObjConfig {
    fn default() -> Self {
        Self {
            tangent_scale: Some(1.0),
            normal_scale: 1.0,
            normalize_horizontal: false,
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

    let tangent_scale = cfg.tangent_scale.unwrap_or(0.5);
    let normal_scale = cfg.normal_scale;
    let horizontal_scale = 1.0;

    let write_frame = cfg.tangent_scale.is_some();
    let vertices_per_frame = if write_frame { 4 } else { 2 };

    for (frame, segment) in curve {
        let [tx, ty, tz] = frame.tangent.to_array().map(|x| x * tangent_scale);
        let [nx, ny, nz] = segment.normal.to_array().map(|x| x * normal_scale);
        let [hx, hy, hz] = if cfg.normalize_horizontal {
            segment
                .horizontal
                .normalize_or_zero()
                .to_array()
                .map(|x| x * horizontal_scale)
        } else {
            segment.horizontal.to_array().map(|x| x * horizontal_scale)
        };

        writeln!(
            &mut string,
            "v {:.4} {:.4} {:.4}",
            frame.base.x, frame.base.y, frame.base.z
        )?;

        writeln!(
            &mut string,
            "v {:.4} {:.4} {:.4}",
            frame.base.x + hx,
            frame.base.y + hy,
            frame.base.z + hz
        )?;

        if write_frame {
            writeln!(
                &mut string,
                "v {:.4} {:.4} {:.4}",
                frame.base.x + tx,
                frame.base.y + ty,
                frame.base.z + tz
            )?;

            writeln!(
                &mut string,
                "v {:.4} {:.4} {:.4}",
                frame.base.x + nx,
                frame.base.y + ny,
                frame.base.z + nz
            )?;
        }
    }

    // Assign some vertex normals, all vertices of a frame pointing in the plane normal direction.
    // So we only have one normal per frame.
    for (_, segment) in curve {
        let [nx, ny, nz] = segment.normal.to_array();
        writeln!(&mut string, "vn {:.4} {:.4} {:.4}", nx, ny, nz)?;
    }

    // And connect all vertices of a frame to the next one with quads.
    for i in 0..curve.len() - 1 {
        let normal = i + 1; // 1-based indexing for .obj
        let next_normal = normal + 1;
        let base = i * vertices_per_frame + 1;

        // Let's make sure this is the right way up. Each horizontal, tangent, normal triple is
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
        // And line elements for the tangents, normals and horizontals.
        for i in 0..curve.len() {
            let base = i * vertices_per_frame + 1;
            writeln!(&mut string, "l {} {}", base, base + 2)?; // tangent
            writeln!(&mut string, "l {} {}", base, base + 3)?; // normal
        }
    }

    Ok(super::StrFileData { contents: string })
}
