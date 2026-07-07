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

trait WriteObj {
    type Err;
    fn append_vertex(&mut self, p: glam::DVec3) -> Result<(), Self::Err>;
    fn append_vn(&mut self, n: glam::DVec3) -> Result<(), Self::Err>;
    fn append_triangle(&mut self, idx: [usize; 3], normals: [usize; 3]) -> Result<(), Self::Err>;
    fn append_quad(&mut self, idx: [usize; 4], normals: [usize; 4]) -> Result<(), Self::Err>;
}

pub(crate) trait MkNever {
    type Output;

}

impl<R> MkNever for fn() -> R {
    type Output = R;
}

impl WriteObj for String {
    type Err = <fn() -> ! as MkNever>::Output;

    fn append_vertex(&mut self, p: glam::DVec3) -> Result<(), Self::Err> {
        let [nx, ny, nz] = p.to_array();
        writeln!(self, "v {:.4} {:.4} {:.4}", nx, ny, nz).unwrap();
        Ok(())
    }

    fn append_vn(&mut self, n: glam::DVec3) -> Result<(), Self::Err> {
        let [nx, ny, nz] = n.to_array();
        writeln!(self, "vn {:.4} {:.4} {:.4}", nx, ny, nz).unwrap();
        Ok(())
    }

    fn append_triangle(&mut self, idx: [usize; 3], normals: [usize; 3]) -> Result<(), Self::Err> {
        let [p0, p1, p2] = idx;
        let [n0, n1, n2] = normals;
        writeln!(self, "f {p0}//{n0} {p1}//{n1} {p2}//{n2}").unwrap();
        Ok(())
    }

    fn append_quad(&mut self, idx: [usize; 4], normals: [usize; 4]) -> Result<(), Self::Err> {
        let [p0, p1, p2, p3] = idx;
        let [n0, n1, n2, n3] = normals;
        writeln!(self, "f {p0}//{n0} {p1}//{n1} {p2}//{n2} {p3}//{n3}").unwrap();
        Ok(())
    }
}

pub fn to_obj(
    curve: &[(DenormalTangentFrame, CurveSegment)],
) -> Result<super::StrFileData, core::fmt::Error> {
    ObjConfig::default().surface(curve)
}

impl ObjConfig {
    pub fn surface(
        &self,
        curve: &[(DenormalTangentFrame, CurveSegment)],
    ) -> Result<super::StrFileData, core::fmt::Error> {
        let mut string = String::new();

        if let Some(cmt) = &self.comment {
            for line in cmt.lines() {
                writeln!(string, "# {line}")?;
            }
        }

        let model_bounds = Self::model_bounds(curve.iter().map(|(frame, _)| frame.base));
        let tangent_scale = f64::from(self.tangent_scale.unwrap_or(0.5));
        let normal_scale = f64::from(self.normal_scale);

        // Print optimized: 180mm build plate.
        let ruling_scale = 1.0f64;
        let model_scale = f64::from(self.buildplate_mm) / model_bounds;

        let write_frame = self.tangent_scale.is_some();
        let vertices_per_frame = if write_frame { 4 } else { 2 };

        for (frame, segment) in curve {
            let [tx, ty, tz] = frame.tangent.to_array().map(|x| x * tangent_scale);
            let [nx, ny, nz] = segment.normal.axis.to_array().map(|x| x * normal_scale);

            let radius_of_curvature = (1.1 * segment.flat_curvature).max(1.0).recip();

            let [hx, hy, hz] = if self.normalize_ruling {
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

    pub fn pipe(
        &self,
        segments: &[(DenormalTangentFrame, dc_integral::TrianglePipeBase)],
    ) -> Result<super::StrFileData, core::fmt::Error> {
        let mut string = String::new();

        if let Some(cmt) = &self.comment {
            for line in cmt.lines() {
                writeln!(string, "# {line}")?;
            }
        }

        let model_bounds = Self::model_bounds(
            segments
                .iter()
                .flat_map(|(frame, pipe)| [frame.base, pipe.base1, pipe.base2]),
        );

        let tangent_scale = f64::from(self.tangent_scale.unwrap_or(0.5));
        let normal_scale = f64::from(self.normal_scale);

        // Print optimized: 180mm build plate.
        let ruling_scale = 1.0f64;
        let model_scale = f64::from(self.buildplate_mm) / model_bounds;

        for (frame, segment) in segments {
            string.append_vertex(frame.base * model_scale);
            string.append_vertex(segment.base1 * model_scale);
            string.append_vertex(segment.base2 * model_scale);
        }

        for (frame, segment) in segments {
            let arm2 = segment.base2 - frame.base;
            let arm1 = segment.base1 - frame.base;

            let forward_normal = arm2.cross(arm1);
            let no = -segment.opposing_normal;
            let n1 = -arm1.cross(frame.tangent).normalize_or_zero();
            let n2 = -frame.tangent.cross(arm2).normalize_or_zero();

            string.append_vn(forward_normal);
            string.append_vn(no);
            string.append_vn(n1);
            string.append_vn(n2);
        }

        for i in 0..segments.len() {
            let base = i * 3 + 1;
            let normal = i * 4 + 1;

            string.append_triangle([base, base + 1, base + 2], [normal; 3]);
        }

        for i in 1..segments.len() {
            let prev = (i - 1) * 3 + 1;
            let prevn = (i - 1) * 4 + 1;
            let base = i * 3 + 1;
            let basen = i * 4 + 1;

            string.append_quad(
                [prev, base, base + 1, prev + 1],
                [prevn + 2, basen + 2, basen + 2, prevn + 2],
            );

            string.append_quad(
                [prev + 1, base + 1, base + 2, prev + 2],
                [prevn + 1, basen + 1, basen + 1, prevn + 1],
            );

            string.append_quad(
                [prev + 2, base + 2, base, prev],
                [prevn + 3, basen + 3, basen + 3, prevn + 3],
            );
        }

        Ok(super::StrFileData {
            contents: string,
            scale: model_scale as f32,
        })
    }

    fn model_bounds(iter: impl IntoIterator<Item = glam::DVec3>) -> f64 {
        let (min, max) = iter.into_iter().fold(
            ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]),
            |(min, max), coord| {
                let [x, y, z] = coord.to_array();
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

        (max[0] - min[0]).max(max[1] - min[1]).max(max[2] - min[2])
    }
}

pub fn to_obj_with(
    curve: &[(DenormalTangentFrame, CurveSegment)],
    cfg: ObjConfig,
) -> Result<super::StrFileData, core::fmt::Error> {
    cfg.surface(curve)
}
