//! Note: SVG's coordinate space is left-handed, so the y-axis is flipped compared to the usual
//! mathematical convention. We unflip it because we are not savages.
use core::fmt::Write as _;

use dc_integral::CurveSegment;
use dc_theory::DenormalTangentFrame;

pub struct SvgConfig {
    pub to_mm: f32,
    /// Disable lines between vertices of opposing borders.
    pub disable_cross_bars: bool,
    /// Include mount holes? Incompatible with and disables cross bars.
    pub mount_holes_mm: Option<f32>,
}

impl SvgConfig {
    pub fn pipe(
        self,
        curve: &[(DenormalTangentFrame, dc_integral::TrianglePipeBase)],
    ) -> Result<[super::StrFileData; 3], core::fmt::Error> {
        pipe(curve, self.apply_incompatibility())
    }

    fn apply_incompatibility(mut self) -> Self {
        self.disable_cross_bars |= self.mount_holes_mm.is_some();
        self
    }
}

trait SvgSink {
    type Err;

    fn path(&mut self, xy: impl IntoIterator<Item = [f64; 2]>) -> Result<(), Self::Err>;

    fn path_closed(&mut self, xy: impl IntoIterator<Item = [f64; 2]>) -> Result<(), Self::Err>;

    fn circle(&mut self, base: [f64; 2], sz: f64) -> Result<(), Self::Err>;
}

impl SvgSink for String {
    type Err = <fn() -> ! as super::obj::MkNever>::Output;

    fn path(&mut self, xy: impl IntoIterator<Item = [f64; 2]>) -> Result<(), Self::Err> {
        let mut xy = xy.into_iter();

        let Some([x, y]) = xy.next() else {
            return Ok(());
        };

        write!(self, r#"  <path d=""#).unwrap();
        write!(self, "M {x:.4} {y:.4} ").unwrap();

        for [x, y] in xy {
            write!(self, "L {x:.4} {y:.4} ").unwrap()
        }

        writeln!(self, r#"" stroke="black" fill-opacity="0.0" />"#).unwrap();

        Ok(())
    }

    fn path_closed(&mut self, xy: impl IntoIterator<Item = [f64; 2]>) -> Result<(), Self::Err> {
        let mut xy = xy.into_iter();

        let Some([x, y]) = xy.next() else {
            return Ok(());
        };

        write!(self, r#"  <path d=""#).unwrap();
        write!(self, "M {x:.4} {y:.4} ").unwrap();

        for [x, y] in xy {
            write!(self, "L {x:.4} {y:.4} ").unwrap()
        }

        writeln!(self, r#"Z" fill="black" />"#).unwrap();

        Ok(())
    }

    fn circle(&mut self, [x, y]: [f64; 2], sz: f64) -> Result<(), Self::Err> {
        write!(
            self,
            r#"  <circle cx="{x}" cy="{y}" r="{radius}" fill="black" />"#,
            radius = sz / 2.0
        )
        .unwrap();

        Ok(())
    }
}

/// Convert to an SVG.
///
/// The SVG coordinate space is left-handed since (0, 0) is top-left and increases right and down.
/// Since we do like a good right handed one instead put `svg_y = -y`, also we apply a scale so that
/// stroke strengths defaults work out in our favor.
pub fn to_svg(
    curve: &[(DenormalTangentFrame, CurveSegment)],
    to_mm: f32,
) -> Result<super::StrFileData, core::fmt::Error> {
    let viewport = ViewportEmbedding::new(
        to_mm,
        curve.iter().map(|(_, segment)| segment.flat_position),
    );

    let scale = f64::from(viewport.scale);
    let mut string = viewport.start_svg();

    let eof = string.split_off(string.find('>').unwrap() + 2);

    {
        writeln!(&mut string, r#"  <g>"#)?;

        string.path(curve.iter().map(|(_, segment)| {
            let [x, y] = segment.flat_position.to_array();
            [x * scale, -y * scale]
        }));

        for (_, segment) in curve.get(1..).into_iter().flatten() {
            let radius_of_curvature = (1.1 * segment.flat_curvature).max(1.0).recip();

            let [x, y] = segment.flat_position.to_array();
            let [x, y] = [x * scale, -y * scale];

            let (dir_y, dir_x) = segment.flat_direction.sin_cos();
            let angle_rotation = glam::DVec2::from_angle(segment.angle);
            let [dir_x, dir_y] = angle_rotation
                .rotate(glam::DVec2::new(dir_x, dir_y))
                .to_array();

            let [dx, dy] = [dir_x, dir_y].map(|x| x * radius_of_curvature * scale);
            let [dx, dy] = [dx, -dy];

            writeln!(
                &mut string,
                r#"  <path d="M {x} {y} l {dx} {dy}" stroke="black" fill="transparent" />"#
            )?;
        }

        writeln!(&mut string, r#"</g>"#)?;
    }

    string.extend(eof.chars());

    Ok(super::StrFileData {
        contents: string,
        scale: scale as f32,
    })
}

fn pipe(
    curve: &[(DenormalTangentFrame, dc_integral::TrianglePipeBase)],
    cfg: SvgConfig,
) -> Result<[super::StrFileData; 3], core::fmt::Error> {
    let access: [fn(&dc_integral::TrianglePipeBase) -> dc_integral::PipeFaceBase; 3] = [
        |segment| segment.opposing_flat,
        |segment| segment.flat1,
        |segment| segment.flat2,
    ];

    let mut output = [(); 3].map(|_| super::StrFileData {
        contents: String::new(),
        scale: cfg.to_mm,
    });

    for (access, output) in access.into_iter().zip(&mut output) {
        let viewport = ViewportEmbedding::new(
            cfg.to_mm,
            curve.iter().flat_map(|(_, segment)| {
                let face = access(segment);
                [face.base_right, face.base_left]
            }),
        );

        let scale = f64::from(viewport.scale);
        let mut string = viewport.start_svg();

        let eof = string.split_off(string.find('>').unwrap() + 2);

        let has_mask_path = cfg.mount_holes_mm.is_some();

        if has_mask_path {
            write!(
                &mut string,
                r#"  <mask id="maskHoles" mask-type="luminance">"#
            )
            .unwrap();

            write!(
                &mut string,
                r#"  <rect x="{x}" y="{y}" width="{w}" height="{h}" fill="white" />"#,
                x = viewport.min[0] * scale,
                y = -viewport.max[1] * scale,
                w = (viewport.max[0] - viewport.min[0]) * scale,
                h = (viewport.max[1] - viewport.min[1]) * scale,
            )
            .unwrap();
        }

        if let Some(diameter) = cfg.mount_holes_mm {
            let diameter = viewport.from_mm_to_inner() * f64::from(diameter);
            let offset_len = viewport.from_mm_to_inner() * f64::from(0.3) + 0.5 * diameter;

            let mut pts = curve.iter();
            let _ = pts.next();
            let _ = pts.next_back();

            for (_, segment) in pts {
                let face = access(segment);

                let offset = (face.base_right - face.base_left).normalize_or_zero();
                let offset = offset * offset_len;

                let [lx, ly] = (face.base_left + offset).to_array();
                let [rx, ry] = (face.base_right - offset).to_array();

                string.circle([lx * scale, -ly * scale], diameter * scale);
                string.circle([rx * scale, -ry * scale], diameter * scale);
            }
        }

        if has_mask_path {
            write!(&mut string, r#"  </mask>"#).unwrap();
            write!(&mut string, r#"  <g mask="url(#maskHoles)">"#).unwrap();
        }

        let left_forward = curve.iter().map(|(_, segment)| {
            let [x, y] = access(segment).base_left.to_array();
            [x * scale, -y * scale]
        });

        let right_forward = curve.iter().map(|(_, segment)| {
            let [x, y] = access(segment).base_right.to_array();
            [x * scale, -y * scale]
        });

        string.path_closed(left_forward.chain(right_forward.rev()));

        if has_mask_path {
            write!(&mut string, r#"  </g>"#).unwrap();
        }

        if !cfg.disable_cross_bars {
            let mut points = curve.iter();
            let _ = points.next();
            let _ = points.next_back();

            for (_, segment) in points {
                let face = access(segment);
                let [lx, ly] = face.base_left.to_array();
                let [rx, ry] = face.base_right.to_array();

                string.path([[lx * scale, -ly * scale], [rx * scale, -ry * scale]]);
            }
        }

        string.extend(eof.chars());
        output.contents = string;
    }

    Ok(output)
}

struct ViewportEmbedding {
    scale: f32,
    point_to_mm: f32,
    min: [f64; 2],
    max: [f64; 2],
}

impl ViewportEmbedding {
    fn new(point_to_mm: f32, pos: impl IntoIterator<Item = glam::DVec2>) -> Self {
        let (min, max) = pos.into_iter().fold(
            (
                [f64::INFINITY, f64::INFINITY],
                [f64::NEG_INFINITY, f64::NEG_INFINITY],
            ),
            |(min, max), pos| {
                let [x, y] = pos.to_array();
                let min = [min[0].min(x - 1.0), min[1].min(y - 1.0)];
                let max = [max[0].max(x + 1.0), max[1].max(y + 1.0)];
                (min, max)
            },
        );

        ViewportEmbedding {
            // This is chosen such that the width of lines works visually together with the other parts of
            // paining the figure. Whereas the real-world is then expressed by providing a physical size of
            // the viewBox itself.
            scale: 400.,
            point_to_mm,
            min,
            max,
        }
    }

    /// Convert a size, given in mm, to a point.
    fn from_mm_to_inner(&self) -> f64 {
        1.0 / f64::from(self.point_to_mm)
    }

    fn start_svg(&self) -> String {
        let scale = f64::from(self.scale);

        // On SVG: x is up, y is right.
        let width = self.max[0] - self.min[0];
        let height = self.max[1] - self.min[1];
        let [mx, my] = [self.min[0] * scale, -self.max[1] * scale];

        format!(
            r#"<svg version="1.1" viewBox="{mx:.4} {my:.4} {width:.4} {height:.4}" width="{width_mm}mm" height="{height_mm}mm" xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMidYMid">{}</svg>"#,
            "\n",
            width = width * scale,
            height = height * scale,
            width_mm = f64::from(self.point_to_mm) * width,
            height_mm = f64::from(self.point_to_mm) * height,
        )
    }
}
