use glam::DVec3;

/// A *non*-normalized frame. Pretty much none of our definitions care about it.
#[derive(Clone, Copy, Debug)]
pub struct DenormalTangentFrame {
    /// The point on the curve.
    pub base: DVec3,
    /// The tangent `f(t)'`. Note we do not require unit speed curves. This can have any (non-zero)
    /// length.
    pub tangent: DVec3,
    /// Derivative of the tangent. Also called `normal` if `tangent` is constant unit length and
    /// curve is unit speed.
    pub derivative: DVec3,
    /// Describe the curl.
    pub binormal: DVec3,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceNormal {
    /// This should be *constant* length.
    pub axis: DVec3,
}

impl SurfaceNormal {
    pub fn from_array(arr: [f32; 3]) -> Self {
        Self {
            axis: DVec3::from_array(arr.map(f64::from)).normalize(),
        }
    }
}

pub struct SurfaceDevelopment {
    pub frame: DenormalTangentFrame,
    pub normal: DVec3,
    /// A point fulfilling the constraints with the tangent frame to remain orthogonal to the
    /// tangent.
    pub derivative_base: DVec3,
    /// Derivative is free along a line normal to the normal. Add any multiple of this to the base
    /// and you get a valid derivative.
    pub derivative_free: DVec3,
    /// What direction is the basic frame oriented? We want to stay on a consistent size even if
    /// the direction of the curvature flips.
    pub signum: f32,
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
    pub fn from_frame_and_normal(
        frame: DenormalTangentFrame,
        SurfaceNormal { axis: normal }: SurfaceNormal,
    ) -> Self {
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

        let sign_of_curve = frame.tangent.cross(frame.derivative).dot(normal).signum();

        Self {
            frame,
            normal,
            derivative_base: base,
            derivative_free: dir,
            signum: sign_of_curve as f32,
        }
    }
}

#[derive(Clone, Copy)]
pub struct CurveDescription {
    pub tangent: DVec3,
    pub dt_tangent: DVec3,
    pub dt_normal: DVec3,
    /// Angle between the intended horizontal direction and the tangent. If set it must be
    /// consistent with the `dt_normal` direction and its sign determines the orientation. Please
    /// note that usually an angle of `+-90` is not possible.
    pub angle: Option<f64>,
    pub dt_normal_normalized: Option<DVec3>,
}

impl CurveDescription {
    pub fn curvature_to_normal(&self, normal: DVec3) -> f64 {
        // Projected curvature is derived from projected curve normal. We can recover the
        // derivative of the othonormal frame of the curve first:
        //
        //     with l(x) = |frame.tangent| = <frame.tangent, frame.tangent>^(1/2) and l' = dl/dt,
        //
        //     The derivative of the tangent splits into a component of unit direction, and the
        //     derivative of this unit direction vector in the direction of the normal. With `l`
        //     being the parameterization's curve speed the 'normal' is `l` times larger than its
        //     unit speed curve. So:
        //
        //     frame.derivative = (l · frame.utangent)' = l' frame.utangent + l · dt frame.utangent
        //     <frame.tangent, frame.derivative> = l' / l <frame.tangent, frame.tangent> = l' · l
        //
        //     frame.derivative = l · dt frame.utangent + (<frame.tangent, frame.derivative> / l) · frame.tangent / l
        //     frame.derivative = l · dt frame.utangent + (<frame.tangent, frame.derivative> / l²) · frame.tangent
        //     l · dt frame.utangent = frame.derivative - (<frame.tangent, frame.derivative> / l²) · frame.tangent
        //
        // We can measure the flattened curvature through the distance metric that is necessarily
        // preserved by our projection.
        //
        //   ĸ = || du r × ds (du r) || = || ds (du r) || since du r is a unit vector.
        //
        // For a non-unit speed curve we need to find their unit-speed equivalents. In particular:
        //
        //     ds normal = (dt / ds) · (ds / dt) ds normal = dt normal · (dt / ds) = dt normal / l
        //
        //     (ds / dt) = ||frame.tangent|| = l
        //
        // Thusly
        //
        //     ds (du r) = ds (frame.utangent × normal)
        //     = ds frame.utangent × normal + frame.utangent × ds normal
        //     = dt frame.utangent × normal / l + frame.utangent × dt normal / l
        //
        //     dt frame.utangent × normal = l · dt frame.utangent × normal / l
        //     = ((frame.derivative - (<frame.tangent, frame.derivative> / l²) · frame.tangent) × normal) / l
        //
        // NOTE: previously got caught in a GPT-4.1 rabbit hole. It one-shot:
        //
        // `kappa = (dt frame.tangent×normal) / <frame.tangent, frame.tangent>`.
        //
        // (without the above derivation) and used a dot instead of cross product. It did not
        // provide any explanation at all. This however caught me caught up in ignoring the term
        // involving the tangent at all which led to crazy curvatures where l' != 0.
        let v = normal.cross(self.tangent) / self.tangent.length();

        self.tangent.cross(normal).dot(self.dt_tangent)
            / self.tangent.length_squared()
            / self.tangent.length()
    }
}
