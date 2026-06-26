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
