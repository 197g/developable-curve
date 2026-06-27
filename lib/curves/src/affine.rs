pub use super::{Curve, DenormalTangentFrame};
use glam::{Affine3, DQuat, Vec3};

pub struct Affine<C> {
    pub inner: C,
    pub isometry: Affine3,
}

pub struct Translate<C> {
    pub inner: C,
    pub translation: Vec3,
}

impl<C: Curve> Affine<C> {
    /// Chain two curves such that their Frenet-Sorret frames will align.
    ///
    /// Note that it is not strictly necessary for normal continuity to align the upwards direction
    /// of these curves. As long the previous normal stays orthogonal to the tangent we may find
    /// some continuation. However! Note that a discontinuity when the frame binormal it parallel
    /// with the surface normal. In this case the new curve can not influence the normal direction
    /// at all and is instead on the plane between tangent and normal—indeed the tangent developable
    /// is an example of this.
    ///
    /// FIXME: Add constructor that aligns the surface normal to a target angle. We want to have
    /// curve segments with their interior controlling the parameterization by that angle instead of
    /// the other way around (in particular holding it constant). It seems reasonable that the use
    /// case is either you have computed the angle of the previous segment a-prior, in which case
    /// you can just use that number, or the angle is more important than the segment's orientation
    /// in 3D. (This does not reason about a separate utility for initialization of the normal at
    /// the first segment, that should be an extra method).
    pub fn with_aligned(before: &(impl Curve + ?Sized), end: f32, inner: C, start: f32) -> Self {
        fn view_for(frame: DenormalTangentFrame) -> DQuat {
            let up = frame
                .tangent
                .cross(frame.derivative)
                .normalize_or(glam::DVec3::Y);
            let forward = frame.tangent.normalize_or(-glam::DVec3::Z);
            DQuat::look_to_rh(forward, up)
        }

        let end_frame = before.at(end);
        let start_frame = inner.at(start);
        let offset = (end_frame.base - start_frame.base).as_vec3();

        let end_view = view_for(end_frame);
        let start_view = view_for(start_frame);

        let isometry = start_view.mul_quat(end_view.inverse());

        Affine {
            inner,
            isometry: Affine3::from_scale_rotation_translation(
                Vec3::ONE,
                isometry.as_quat(),
                offset,
            ),
        }
    }
}

impl<C: Curve> Curve for Affine<C> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let DenormalTangentFrame {
            base,
            tangent,
            derivative,
            third_derivative,
        } = self.inner.at(t);

        let tr = self.isometry.as_daffine3();

        DenormalTangentFrame {
            base: tr.transform_point3(base),
            tangent: tr.transform_vector3(tangent),
            derivative: tr.transform_vector3(derivative),
            third_derivative: tr.transform_vector3(third_derivative),
        }
    }
}

impl<C: Curve> Curve for Translate<C> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let DenormalTangentFrame {
            base,
            tangent,
            derivative,
            third_derivative,
        } = self.inner.at(t);

        DenormalTangentFrame {
            base: base + self.translation.as_dvec3(),
            tangent,
            derivative,
            third_derivative,
        }
    }
}
