use core::fmt::Write as _;

pub fn to_obj(
    curve: &[(super::DenormalTangentFrame, super::CurveSegment)],
) -> Result<String, core::fmt::Error> {
    let mut string = String::new();

    let tangent_scale = 0.5;
    let normal_scale = 0.5;

    // The svg generated of flat is for the left-side of this though it does not matter.
    let horizontal_scale = -1.0;

    const VERTICES_PER_FRAME: usize = 4;

    for (frame, segment) in curve {
        let [tx, ty, tz] = frame.tangent.to_array().map(|x| x * tangent_scale);
        let [nx, ny, nz] = segment.normal.to_array().map(|x| x * normal_scale);
        let [hx, hy, hz] = segment.horizontal.to_array().map(|x| x * horizontal_scale);

        writeln!(
            &mut string,
            "v {:.4} {:.4} {:.4}",
            frame.base.x, frame.base.y, frame.base.z
        )?;

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

        writeln!(
            &mut string,
            "v {:.4} {:.4} {:.4}",
            frame.base.x + hx,
            frame.base.y + hy,
            frame.base.z + hz
        )?;
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
        let base = i * VERTICES_PER_FRAME + 1;

        // Let's make sure this is the right way up. Each horizontal, tangent, normal triple is
        // right handed in that order.
        writeln!(
            &mut string,
            "f {}//{normal} {}//{normal} {}//{next_normal} {}//{next_normal}",
            base + 3,
            base,
            base + VERTICES_PER_FRAME,
            base + 3 + VERTICES_PER_FRAME,
        )?;
    }

    // And line elements for the tangents, normals and horizontals.
    for i in 0..curve.len() {
        let base = i * VERTICES_PER_FRAME + 1;
        writeln!(&mut string, "l {} {}", base, base + 1)?; // tangent
        writeln!(&mut string, "l {} {}", base, base + 2)?; // normal
    }

    Ok(string)
}
