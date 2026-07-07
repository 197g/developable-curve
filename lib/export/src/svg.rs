//! Note: SVG's coordinate space is left-handed, so the y-axis is flipped compared to the usual
//! mathematical convention. We unflip it because we are not savages.
use core::fmt::Write as _;

use dc_integral::CurveSegment;
use dc_theory::DenormalTangentFrame;

trait SvgSink {
    type Err;

    fn path(&mut self, xy: impl IntoIterator<Item = [f64; 2]>) -> Result<(), Self::Err>;

    fn delta_path(
        &mut self,
        base: [f64; 2],
        xy: impl IntoIterator<Item = [f64; 2]>,
    ) -> Result<(), Self::Err>;
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

    fn delta_path(
        &mut self,
        [x, y]: [f64; 2],
        xy: impl IntoIterator<Item = [f64; 2]>,
    ) -> Result<(), Self::Err> {
        write!(self, r#"  <path d="M {x} {y}"#).unwrap();

        for [dx, dy] in xy {
            write!(
                self,
                r#" l {dx} {dy}" stroke="black" fill="transparent" />"#
            )
            .unwrap();
        }

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

pub fn pipe(
    _: &[(DenormalTangentFrame, dc_integral::TrianglePipeBase)],
    scale: f32,
) -> Result<super::StrFileData, core::fmt::Error> {
    Ok(super::StrFileData {
        contents: "".to_string(),
        scale,
    })
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

    fn start_svg(&self) -> String {
        let scale = f64::from(self.scale);

        // On SVG: x is up, y is right.
        let width = (self.max[0] - self.min[0]).max(5.);
        let height = (self.max[1] - self.min[1]).max(5.);
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
