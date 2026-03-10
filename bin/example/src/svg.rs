//! Note: SVG's coordinate space is left-handed, so the y-axis is flipped compared to the usual
//! mathematical convention. We unflip it because we are not savages.
use core::fmt::Write as _;

pub fn to_svg(
    curve: &[(super::DenormalTangentFrame, super::CurveSegment)],
) -> Result<String, core::fmt::Error> {
    let scale = 40.0;

    let (min, max) = curve.iter().fold(
        (
            [f32::INFINITY, f32::INFINITY],
            [f32::NEG_INFINITY, f32::NEG_INFINITY],
        ),
        |(min, max), (_, segment)| {
            let [x, y] = segment.flat_position.to_array();
            let min = [min[0].min(x - 1.0), min[1].min(y - 1.0)];
            let max = [max[0].max(x + 1.0), max[1].max(y + 1.0)];
            (min, max)
        },
    );

    let [rx, ry] = [
        (max[0] - min[0]).max(5.) * scale,
        (max[1] - min[1]).max(5.) * scale,
    ];

    let [mx, my] = [min[0] * scale, min[1] * scale];

    let mut string = format!(
        r#"<svg version="1.1" viewBox="{mx:.4} {my:.4} {rx:.4} {ry:.4}" xmlns="http://www.w3.org/2000/svg" transform="scale(1,-1)" >\n</svg>"#
    );

    let eof = string.split_off(string.find('>').unwrap() + 2);

    {
        writeln!(&mut string, r#"  <g>"#)?;
        write!(&mut string, r#"  <path d=""#)?;

        if let Some((_, first)) = curve.first() {
            let [x, y] = first.flat_position.to_array();
            let [x, y] = [x, y].map(|x| x * scale);
            write!(&mut string, "M {x:.4} {y:.4} ")?;
        }

        for (_, segment) in curve.get(1..).into_iter().flatten() {
            let [x, y] = segment.flat_position.to_array();
            let [x, y] = [x, y].map(|x| x * scale);
            write!(&mut string, "L {x:.4} {y:.4} ")?;
        }

        writeln!(&mut string, r#"" stroke="black" fill="transparent" />"#)?;

        for (_, segment) in curve.get(1..).into_iter().flatten() {
            let [x, y] = segment.flat_position.to_array();
            let [x, y] = [x, y].map(|x| x * scale);
            let (dir_y, dir_x) = segment.flat_direction.sin_cos();
            let [dx, dy] = [dir_x, dir_y].map(|x| x * scale);
            let [dx, dy] = [-dy, dx];

            writeln!(
                &mut string,
                r#"  <path d="M {x} {y} l {dx} {dy}" stroke="black" fill="transparent" />"#
            )?;
        }

        writeln!(&mut string, r#"</g>"#)?;
    }

    string.extend(eof.chars());
    Ok(string)
}
