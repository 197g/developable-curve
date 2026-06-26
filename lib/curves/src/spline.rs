use super::*;

pub struct BezierSpline<const N: usize> {
    /// The points the curve should pass through.
    pub points: [Vec3; N],
}

impl Curve for BezierSpline<2> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];

        let base = p0.lerp(p1, t).as_dvec3();
        let tangent = (p1 - p0).as_dvec3();
        let derivative = DVec3::ZERO;
        let binormal = DVec3::ZERO;

        DenormalTangentFrame {
            base,
            tangent,
            derivative,
            third_derivative: binormal,
        }
    }
}

impl Curve for BezierSpline<3> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];
        let p2 = self.points[2];

        let base = p0 * (1.0 - t).powi(2) + p1 * 2.0 * (1.0 - t) * t + p2 * t.powi(2);
        let tangent = (p1 - p0) * 2.0 * (1.0 - t) + (p2 - p1) * 2.0 * t;
        let derivative = (p2 - 2.0 * p1 + p0) * 2.0;
        let binormal = DVec3::ZERO;

        DenormalTangentFrame {
            base: base.as_dvec3(),
            tangent: tangent.as_dvec3(),
            derivative: derivative.as_dvec3(),
            third_derivative: binormal,
        }
    }
}

// This is where it gets interesting, this curve may have a curl.
impl Curve for BezierSpline<4> {
    fn at(&self, t: f32) -> DenormalTangentFrame {
        let p0 = self.points[0];
        let p1 = self.points[1];
        let p2 = self.points[2];
        let p3 = self.points[3];

        // TODO: this is a stupid way of evaluating this.
        let base = p0 * (1.0 - t).powi(3)
            + p1 * 3.0 * (1.0 - t).powi(2) * t
            + p2 * 3.0 * (1.0 - t) * t.powi(2)
            + p3 * t.powi(3);
        let tangent = (p1 - p0) * 3.0 * (1.0 - t).powi(2)
            + (p2 - p1) * 6.0 * (1.0 - t) * t
            + (p3 - p2) * 3.0 * t.powi(2);
        let derivative = (p2 - 2.0 * p1 + p0) * 6.0 * (1.0 - t) + (p3 - 2.0 * p2 + p1) * 6.0 * t;

        // At first gpt-41 said this was `0`, but that is only the case for `t=0` and `t=1`. The
        // binormal is the derivative of the derivative, which is a constant in this case.
        // ^ This is its own autocomplete when I removed the line. LMAO.
        let binormal = (p3 - 3.0 * p2 + 3.0 * p1 - p0) * 6.0;

        DenormalTangentFrame {
            base: base.as_dvec3(),
            tangent: tangent.as_dvec3(),
            derivative: derivative.as_dvec3(),
            third_derivative: binormal.as_dvec3(),
        }
    }
}
