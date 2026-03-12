use glam::Vec3;

mod curves;

pub use curves::{Circle, Curve, HermiteSpline};

/// A *non*-normalized frame. Pretty much none of our definitions care about it.
#[derive(Clone, Copy)]
pub struct DenormalTangentFrame {
    /// The point on the curve.
    pub base: Vec3,
    /// The tangent `f(t)'`. Note we do not support unit speed curves. This can have any (non-zero)
    /// length.
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

impl SurfaceNormal {
    pub fn from_array(arr: [f32; 3]) -> Self {
        Self {
            normal: Vec3::from_array(arr).normalize(),
        }
    }
}

pub struct SurfaceDevelopment {
    pub frame: DenormalTangentFrame,
    pub normal: Vec3,
    /// A point fulfilling the constraints with the tangent frame to remain orthogonal to the
    /// tangent.
    pub derivative_base: Vec3,
    /// Derivative is free along a line normal to the normal. Add any multiple of this to the base
    /// and you get a valid derivative.
    pub derivative_free: Vec3,
    /// The curvature at the surface, not of the curve.
    pub surface_curvature: f32,
    /// What direction is the basic frame oriented? We want to stay on a consistent size even if
    /// the direction of the curvature flips.
    pub signum: f32,
}

impl SurfaceDevelopment {
    pub fn normal_ode(
        curve: impl Curve,
        parameter: impl Fn(f32) -> f32,
    ) -> impl Fn(Vec3, f32) -> Vec3 {
        move |normal: Vec3, t: f32| {
            let frame = curve.at(t);
            let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { normal });
            let lambda = parameter(t);
            dev.derivative_base + lambda * dev.derivative_free
        }
    }

    pub fn normal_and_flat_ode(
        curve: impl Curve,
        parameter: impl Fn(f32) -> f32,
    ) -> impl Fn(Vec3, f32) -> (DenormalTangentFrame, Vec3, f32, f32) {
        move |normal: Vec3, t: f32| {
            let frame = curve.at(t);
            let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { normal });
            let lambda = parameter(t);
            let dt_normal = dev.derivative_base + lambda * dev.derivative_free;
            let speed = frame.tangent.length();
            (dev.frame, dt_normal, dev.surface_curvature, speed)
        }
    }

    /// Steer the surface and horizontal direction by defining an angle between the tangent and
    /// horizontal direction along the curve.
    pub fn normal_and_angle_ode(
        curve: impl Curve,
        parameter: impl Fn(f32) -> f32,
    ) -> impl Fn(Vec3, f32) -> (DenormalTangentFrame, Vec3, f32, f32) {
        move |normal: Vec3, t: f32| {
            let frame = curve.at(t);
            let dev = SurfaceDevelopment::from_frame_and_normal(frame, SurfaceNormal { normal });
            let target_angle = parameter(t);
            // angle(horizontal, frame.tangent) = atan2(<normal, frame.derivative>, lambda)
            //
            // See `dc-integral/src/lib.rs` for the derivation of this formula where lambda is the
            // parameter from the above formula. Now let's derive that lambda. Note how we
            // automatically get `lambda = 0` at the direction discontinuity.
            let lambda = target_angle.tan() * dev.normal.dot(frame.derivative);
            // ^ LLM anecdote: this was oneshot before the derivation. It badly fumbled the
            // derivation itself though, forgetting the square root in the identity or forgetting
            // that subtract `1` changes the numerator..
            let dt_normal = dev.derivative_base + lambda * dev.derivative_free;
            let speed = frame.tangent.length();
            (dev.frame, dt_normal, dev.surface_curvature, speed)
        }
    }

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
        SurfaceNormal { normal }: SurfaceNormal,
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

        // Projected curvature is derived from projected curve normal. We can recover the
        // derivative of the othonormal frame of the curve first:
        //
        //     with l(x) = |frame.tangent| = <frame.tangent, frame.tangent>^(1/2) and l' = dl/dt,
        //
        //     frame.derivative = l · frame.normal + l' frame.tangent / l
        //     <frame.tangent, frame.derivative> = l' / l <frame.tangent, frame.tangent> = l' · l
        //
        //     frame.derivative = l · frame.normal + (<frame.tangent, frame.derivative> / l) · frame.tangent / l
        //     frame.derivative = l · frame.normal + (<frame.tangent, frame.derivative> / l²) · frame.tangent
        //     l · frame.normal = frame.derivative - (<frame.tangent, frame.derivative> / l²) · frame.tangent
        //
        // What we're interested in is the length of the curve normal vector projected onto the
        // plane (given by its normal). That projection is of the form `T + o · normal` where `o`
        // is `-<normal, T>`. We borrow a trick from curvature. A cross product preserves the
        // lengths for orthogonal vectors: ||A|| = ||A×B|| / ||B||. And we can compute the cross
        // product without doing the projection itself since B×B = 0. Since our plane normal has
        // unit length we may as well be computing
        //
        //     ||frame.normal × normal|| = || l · frame.normal × normal || / l
        //     = || frame.derivative × normal - (<frame.tangent, frame.derivative> / l²) · frame.tangent × normal || / l
        //     = || frame.derivative × normal || / l
        //
        // which is so absurdly clean I'm not even sure it is correct. Just one more division by
        // `l` to correct for the non-unit speed of the curve—which makes this even cleaner. WTF.
        //
        // `kappa = (dt frame.tangent×normal) / <frame.tangent, frame.tangent>`.
        //
        //     GPT-4.1 almost oneshot this (without the above derivation) but used a dot instead of
        //     cross product. It did not oneshot the explanation at all.
        let kappa = normal.cross(frame.derivative).length() / frame.tangent.length_squared();
        let sign_of_curve = frame.tangent.cross(frame.derivative).dot(normal).signum();

        Self {
            frame,
            normal,
            derivative_base: base,
            derivative_free: dir,
            surface_curvature: kappa * sign_of_curve,
            signum: sign_of_curve,
        }
    }
}
