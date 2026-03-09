use glam::Vec3;

/// A *non*-normalized frame. Pretty much none of our definitions care about it.
pub struct DenormalTangentFrame {
    pub tangent: Vec3,
    /// Derivative of the tangent. Also called `normal` if `tangent` is constant unit length and
    /// curve is unit speed.
    pub derivative: Vec3,
    /// Describe the curl.
    pub binormal: Vec3,
}

pub struct SurfaceNormal {
    /// This should be *constant* length.
    pub normal: Vec3,
}

pub struct SurfaceDevelopment {
    pub normal: Vec3,
    /// A point fulfilling the constraints with the tangent frame to remain orthogonal to the
    /// tangent.
    pub derivative_base: Vec3,
    /// Derivative is free along a line normal to the normal. Add any multiple of this to the base
    /// and you get a valid derivative.
    pub derivative_free: Vec3,
}

impl SurfaceDevelopment {
    /// Assume frame and normal describe the line x×[0; inf) of the developable surface.
    ///
    /// Choose a derivative of the plane normal at that point such that the surface is developable.
    /// For this to hold the derivative is in the plane of the surface and constrained by a linear
    /// relation of the normal and the frame normal. Except for points where the curvature is `0`,
    /// the choice of derivative implies the direction of the surface: the plane normal's
    /// derivative is also inside it. If the curvature is `0` then we permit only the trivial
    /// solution; but the choice for the plane becomes arbitrary. You may want a continuous
    /// continuation at such points. Or; you don't, discontinuities are okay! Quite an exciting
    /// loophole from the theory.
    ///
    /// Also note I did not say anything about Darboux frames. They are a special case here.
    pub fn from_frame_and_normal(frame: &DenormalTangentFrame, SurfaceNormal { normal }: SurfaceNormal) -> Self {
        // The constraints on the derivative are as follows:
        // 1. The surface normal has constant length.
        // 2. Surface must remain orthogonal to the first surface direction u, the tangent.
        // 3. Surface must remain orthogonal to the second surface direction v.
        //
        // From(1), `<normal, normal> = c` follows `<derivative, normal> = 0`.
        // From(2), `<normal, frame.tangent> = 0` follows `dt<normal, frame.tangent> = 0`
        //   and `<derivative, frame.tangent> = -<normal, dt frame.tangent>`.
        //   The right hand side is a known quantity since the frame is given. So this is a linear
        //   constraint on the derivative; describing a plane orthogonal to the tangent.
        // From(3), `dt<normal, v> = 0` follows `<derivative, v> = -<normal, dt v> = 0` by the
        //   requirements on `v` (we will uphold that later when choosing `v`/`dt v`.
        //
        // To recap the derivative is
        // - in the plane spanned by the normal and `v` and
        // - orthogonal to `normal` and
        // - and in the plane defined by `frame.normal` and an offset
        //
        // The latter two planes are orthogonal so this defines a line. The open variable is the
        // coordinate along that line. Note neither required `t` to be from a unit length
        // parameterized curve.

        // `0.0` here is due to bad input, these two must be orthogonal (but not both unit length).
        // Note that `<dir, frame.tangent> = <dir, normal> = 0` implies adding it to any choice of
        // derivative is neutral to the constraints.
        let dir = normal.cross(frame.tangent).normalize_or_zero();

        // Orthogonality to the tangent implies a linear relation to the normal here, even if not
        // unit length. We have `frame.derivative = l' frame.tangent + l frame.normal` for some `l`
        // being the length of the tangent. (It must not vanish). The first term disappears under
        // the dot product.
        let offset = -normal.dot(frame.derivative);

        let base = offset / frame.tangent.length_squared() * frame.tangent;

        Self {
            normal,
            derivative_base: base,
            derivative_free: dir,
        }
    }
}
