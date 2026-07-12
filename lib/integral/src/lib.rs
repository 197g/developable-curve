use fast_ode::Coord;
use glam::{DVec2, DVec3};

use dc_theory::{CurveDescription, DenormalTangentFrame, SurfaceNormal};

pub fn curve_ode(
    tangent: impl Fn(DVec3, f64) -> DVec3,
    base: DVec3,
    (start, end): (f64, f64),
) -> DVec3 {
    struct Ode<F: Fn(DVec3, f64) -> DVec3>(F);

    impl<F: Fn(DVec3, f64) -> DVec3> fast_ode::DifferentialEquation<3> for Ode<F> {
        fn ode_dot_y(&self, t: f64, y: &Coord<3>) -> (Coord<3>, bool) {
            let x = DVec3::from_array(y.0);
            let tangent = (self.0)(x, t);
            (Coord(tangent.to_array().map(f64::from)), true)
        }
    }

    let x0 = Coord::<3>(base.to_array().map(f64::from));

    let sol = fast_ode::solve_ivp(
        &Ode(tangent),
        (f64::from(start), f64::from(end)),
        x0,
        |_, _| true,
        1e-6,
        1e-6,
    );

    let x1 = match sol {
        fast_ode::IvpResult::FinalTimeReached(coord) => coord,
        // Deserves a warning, at least!
        fast_ode::IvpResult::StepTooSmall(_, coord) => coord,
        fast_ode::IvpResult::OdeRequestedExit(..)
        | fast_ode::IvpResult::CallbackRequestedExit(..) => {
            unreachable!("we do not request exit")
        }
    };

    DVec3::from_array(x1.0)
}

#[derive(Clone, Copy)]
pub struct CurveSegment {
    pub normal: SurfaceNormal,
    pub ruling: DVec3,
    pub flat_position: DVec2,
    pub flat_direction: f64,
    pub flat_curvature: f64,
    pub angle: f64,
}

fn warn_nonzero(actual: f64, what: &str) {
    if !(actual < 1e-8) {
        eprintln!("{what}: {actual}");
    }
}

impl CurveSegment {
    pub fn initial(normal: SurfaceNormal, ode: impl Fn(DVec3, f64) -> CurveDescription) -> Self {
        let descriptor = ode(normal.axis, 0.0);

        if let Some(angle) = descriptor.angle {
            Self::from_angle(normal.axis, descriptor, angle)
        } else {
            Self::from_parameter_with_unstable_angle_at_zero(normal.axis, descriptor)
        }
    }

    fn from_angle(normal: DVec3, frame: CurveDescription, target_angle: f64) -> Self {
        let forward = frame.tangent.normalize();
        let sideways = normal.cross(forward);
        let (s, c) = target_angle.sin_cos();

        let ruling = c * forward + s * sideways;
        let angle = target_angle;

        warn_nonzero(ruling.dot(frame.dt_normal), "ruling to dt normal");
        warn_nonzero(frame.tangent.dot(normal), "normal to tangent");
        warn_nonzero(ruling.dot(normal), "normal to ruling");

        assert!(
            (normal.dot(frame.dt_tangent) + frame.dt_normal.dot(frame.tangent)).abs() < 1e-8,
            "{:.8} / {:.8}",
            normal.dot(frame.dt_tangent),
            frame.dt_normal.dot(frame.tangent),
        );

        CurveSegment {
            normal: SurfaceNormal { axis: normal },
            ruling,
            flat_position: Default::default(),
            flat_direction: Default::default(),
            flat_curvature: 0.0,
            angle,
        }
    }

    fn from_parameter_with_unstable_angle_at_zero(
        normal: DVec3,
        end_descriptor: CurveDescription,
    ) -> Self {
        // The ruling is orthogonal to both.
        let ruling = normal.cross(end_descriptor.dt_normal);

        let forward = end_descriptor.tangent.normalize();
        let sideways = normal.cross(forward);

        let x = ruling.dot(forward);
        let y = ruling.dot(sideways);

        // Note: `<ruling, frame.tangent> = v · ||frame.tangent||`
        //
        // if you want to control this angle. Expanded:
        //
        // cos(ruling, frame.tangent) · ||ruling|| · ||frame.tangent||
        //     = <ruling, frame.tangent>
        //     = v · ||frame.tangent||
        //
        // v = cos(ruling, frame.tangent) · ||ruling||
        //     = cos(ruling, frame.tangent) · ||dt normal||
        //     = cos(ruling, frame.tangent) · sqrt(v² + <normal, frame.derivative>²)
        //
        // angle(ruling, frame.tangent) = acos(v / sqrt(v² + <normal, frame.derivative>²))
        //     = atan(|<normal, frame.derivative>| / v) ; by acos(x) = atan(sqrt(1-x²)/x)
        //     = atan2(|<normal, frame.derivative>|, v)
        //
        // So we have a discontinuity. If the normal is perpendicular to `frame.derivative` then the
        // `cos(angle) = +1/-1` so the two are parallel with no steering at all. Otherwise, we can
        // choose `v = 0` for a guaranteed tangent-perpendicular ruling line or any other
        // non-parallel angle with appropriate `v`.
        //
        // So now you're asking, can we control `v` so that the discontinuity never occurs? Not in
        // general if the frame.derivative is discontinuous. But also consider this an artifact of our
        // choice of ruling, the direction of which is discontinuous at the zero of `dt_normal`.
        // And indeed at the same point we get a length of `0|v=0`. So really we should maybe instead
        // be steering by the angle; and then calculating a corresponding `v` while having `v=0` and
        // using our angle regardless at the discontinuity?

        // The angle as measured by its 2D projections. (Suggested by Gemini).
        let angle = y.atan2(x);

        CurveSegment {
            normal: SurfaceNormal { axis: normal },
            ruling,
            flat_position: Default::default(),
            flat_direction: Default::default(),
            flat_curvature: 0.0,
            angle,
        }
    }
}

pub fn curve_ode_with_curvature(
    tangent: impl Fn(DVec3, f64) -> CurveDescription,
    base: SurfaceNormal,
    flat_base: (DVec2, f64),
    (start, end): (f64, f64),
) -> CurveSegment {
    struct Ode<F: Fn(DVec3, f64) -> CurveDescription>(F);

    impl<F: Fn(DVec3, f64) -> CurveDescription> fast_ode::DifferentialEquation<6> for Ode<F> {
        fn ode_dot_y(&self, t: f64, ty: &Coord<6>) -> (Coord<6>, bool) {
            let [x, y, z, _, _, _] = ty.0;
            let normal = DVec3::new(x, y, z);
            let descriptor = (self.0)(normal, t);

            let [_, _, _, _cx, _cy, k] = ty.0;
            let [x, y, z] = descriptor.dt_normal.to_array().map(f64::from);

            let speed = f64::from(descriptor.tangent.length());
            let curvature = descriptor.curvature_to_normal(normal);
            // The unit speed curvature but `t` is not unit speed.
            let dkds = f64::from(curvature) * speed;

            // k describes the current heading.
            let (my, mx) = k.sin_cos();
            // Since the 2d curve does not have unit speed either, it must be adjusted itself, too.

            let dt = [x, y, z, mx * speed, my * speed, dkds];

            (Coord(dt), true)
        }
    }

    let x0 = Coord::<6>({
        let [x, y, z] = base.axis.to_array().map(f64::from);
        let [cx, cy] = flat_base.0.to_array().map(f64::from);
        [x, y, z, cx, cy, f64::from(flat_base.1)]
    });

    let ode = Ode(tangent);

    let sol = fast_ode::solve_ivp(
        &ode,
        (f64::from(start), f64::from(end)),
        x0,
        |_, _| true,
        1e-6,
        1e-6,
    );

    let x1 = match sol {
        fast_ode::IvpResult::FinalTimeReached(coord) => coord,
        // Deserves a warning, at least!
        fast_ode::IvpResult::StepTooSmall(_, coord) => coord,
        fast_ode::IvpResult::OdeRequestedExit(..)
        | fast_ode::IvpResult::CallbackRequestedExit(..) => {
            unreachable!("we do not request exit")
        }
    };

    let Coord([x, y, z, fx, fy, k]) = x1;
    let normal = DVec3::from_array([x, y, z]);

    // The ruling must be perpendicular to the plane normal and its derivative.
    // We are however free to choose a direction, let us pick a consistent one.
    let end_descriptor = (ode.0)(normal, end);

    let basis = if let Some(target_angle) = end_descriptor.angle {
        CurveSegment::from_angle(normal, end_descriptor, target_angle)
    } else {
        CurveSegment::from_parameter_with_unstable_angle_at_zero(normal, end_descriptor)
    };

    CurveSegment {
        flat_position: DVec2::from_array([fx, fy]),
        // We do not build a full frame..
        flat_direction: k,
        flat_curvature: end_descriptor.curvature_to_normal(normal),
        ..basis
    }
}

/// Start of a triangular pipe development.
///
/// For a simplified constructor, see in `dc-curve`.
#[derive(Clone, Copy, Debug)]
pub struct TrianglePipeBase {
    pub base1: DVec3,
    pub base2: DVec3,
    pub opposing_normal: DVec3,
    pub flat1: PipeFaceBase,
    pub flat2: PipeFaceBase,
    pub opposing_flat: PipeFaceBase,
}

/// Start of a face of [`TrianglePipeBase`].
#[derive(Clone, Copy, Debug)]
pub struct PipeFaceBase {
    pub base_left: DVec2,
    pub base_right: DVec2,
    pub orientation_left: f64,
    pub orientation_right: f64,
}

/// Parameterization of a triangular pipe development.
///
/// FIXME: this is not a very good description yet. It's technically exhaustive but the parameter
/// _choice_ (insofar as it is one, I do not know) is really odd for physical design. Control by
/// lengths makes it difficult to avoid an inversion of the pipe. (We could measure this as the sign
/// of the dot-product of the tangent with the pipe cross section normal).
pub struct PipeDescription {
    pub frame: DenormalTangentFrame,
    /// The speed of the first (in counter-clockwise order) support curve.
    pub len_a: f64,
    /// The speed of the first (in counter-clockwise order) support curve.
    pub len_b: f64,
    /// The speed of change of normal.
    ///
    /// FIXME: is that really a good description? Both curvatures follow from it maybe there's an
    /// alternative parameterization with more direct control of an important property. The only
    /// direct equivalence that we have is that `0.0` is a flat opposing face and preserves _all_
    /// face's normals.. That's cumbersome for a lot of curves though.
    pub yaw: f64,
}

/// End point and path sketch of a pipe development.
pub struct TrianglePipeSegment {
    pub pipe: PipeSegment,
    pub flat1: PipeFaceSegment,
    pub flat2: PipeFaceSegment,
    pub opposing_flat: PipeFaceSegment,
}

/// The 3d description of the pipe frame itself.
pub struct PipeSegment {
    pub base1: DVec3,
    pub base2: DVec3,
    pub base: DVec3,
}

pub struct PipeFaceSegment {
    /// The 3d normal of this face.
    pub normal: SurfaceNormal,
    /// The left point of the flattened representation.
    pub base_left: DVec2,
    /// The right point of the flattened representation.
    pub base_right: DVec2,
    /// The direction where the face's left flank points.
    pub orientation_left: f64,
    /// The direction where the face's right flank points.
    pub orientation_right: f64,
    /// The flat curvature of the left flank.
    pub curvature_left: f64,
    /// The flat curvature of the right flank.
    pub curvature_right: f64,
}

impl TrianglePipeSegment {
    /// Use this endpoint as the start of another segment.
    pub fn as_next(&self) -> TrianglePipeBase {
        TrianglePipeBase {
            base1: self.pipe.base1,
            base2: self.pipe.base2,
            opposing_normal: self.opposing_flat.normal.axis,
            flat1: self.flat1.as_next(),
            flat2: self.flat2.as_next(),
            opposing_flat: self.opposing_flat.as_next(),
        }
    }
}

impl PipeFaceSegment {
    fn as_next(&self) -> PipeFaceBase {
        PipeFaceBase {
            base_left: self.base_left,
            base_right: self.base_right,
            orientation_left: self.orientation_left,
            orientation_right: self.orientation_right,
        }
    }
}

pub fn triangle_pipe_ode(
    tangent: impl Fn(f64) -> PipeDescription,
    tr: TrianglePipeBase,
    (start, end): (f64, f64),
) -> TrianglePipeSegment {
    struct Params {
        /// 9: The position of each curve (A, B)
        curves: [DVec3; 2],
        /// 3: The orientation of F (Fa and Fb follow from dt Y)
        opposing_normal: DVec3,
        /// 12: The locations of 2d sides
        flats: [[DVec2; 2]; 3],
        // 6: The orientation of 2d sides
        flat_orientation: [[f64; 2]; 3],
    }

    struct Ode<F: Fn(f64) -> PipeDescription>(F);

    impl<F: Fn(f64) -> PipeDescription> fast_ode::DifferentialEquation<27> for Ode<F> {
        /// See `docs/three-surface-problem.md` for derivation.
        fn ode_dot_y(&self, t: f64, y: &Coord<27>) -> (Coord<27>, bool) {
            let params = Params::read(y);
            let curve = (self.0)(t);

            let DenormalTangentFrame {
                base: y,
                tangent,
                derivative,
                third_derivative: _,
            } = curve.frame;

            let [a, b] = params.curves;

            // We usually expect this to be constant length `1` but it may not be over the ODE.
            let normalf = params.opposing_normal.normalize_or_zero();
            // The planes are spanned by the boundary curve tangent and the ruling direction.
            let normalb = tangent.cross(b - y).normalize_or_zero();
            let normala = tangent.cross(y - a).normalize_or_zero();

            // Dot product used as a coefficient and for calculating ||F×Fa||.
            let faf = normala.dot(normalf);
            let fbf = normalb.dot(normalf);

            // Turns out we need 1.0 / ||F×Fa|| twice and calculate this from the dot product. We
            // accept the pole as it will lead to disaster in another part of the formulate, too.
            // For now until we figure out how to remove it entirely.
            let normalize_faf = 1.0 / (1.0 - faf * faf).sqrt();
            let normalize_fbf = 1.0 / (1.0 - fbf * fbf).sqrt();

            // The directions are in the plane of both its surfaces. Normalization through
            // independent calculation of the lengths
            let dir_a = normalf.cross(normala) * normalize_faf;
            let dir_b = normalb.cross(normalf) * normalize_fbf;

            // From the basic theorem, the derivative of the surface normal is orthogonal to the
            // ruling direction within the surface plane. This does not inform us of its length,
            // which is instead constrained with the second derivative of each boundary curve.
            let dir_dtf = (b - a).cross(normalf).normalize_or_zero();

            // These two normalize on their own with the coefficient below.
            let dir_dtfa = (a - y).cross(normala);
            let dir_dtfb = (b - y).cross(normalb);

            let dtf = dir_dtf * curve.yaw;

            // Recall <dt F, dt Y> = -<dt² Y, F>
            //
            // FIXME: make this work even when the direction is orthogonal? This does not really
            // work anyways since it implies the above cross product defining normala/b to be zero
            // so it's moot. But we should not need it.. Maybe it cancels by substitution. Recall
            // that both of these components are essentially triple-products and so we can rotate
            // around the operations to find common components. Also dir_dtfa is defined by the
            // tangent in the dot product..
            //
            // -<derivative,normala> / triple(tangent, a - y, unit(tangent×(y-a)))
            // = -<derivative,normala> / <unit(tangent×(a-y)), tangent×(y-a)>
            // = <derivative,normala>/ ||tangent×(y-a)||
            //
            // Then this is multiplied onto dir_dtfa..
            //
            // And now if we replace normala with its normalize_or_zero variant the problem
            // disappears? But we do not need the value itself, only for a dot-product, so shorten
            // further before.

            // These coefficients express the relationship of `dir_dtfa` to `dt Fa`. We have no need
            // for explicit `dt Fa` so we keep the coefficient only for a dot product.
            let coeff_dtfa = -derivative.dot(normala) / dir_dtfa.dot(tangent);
            let coeff_dtfb = -derivative.dot(normalb) / dir_dtfb.dot(tangent);

            // What we are after is of course <dt Fa, dt A> and <dt Fb, dt B> only.
            //
            // That is:
            // <dt Fa, la· F×normala> = la·coeff_dtfa·<(a - y)×normala, F×normala>
            //
            // But:
            // <(a - y)×normala, F×normala>
            // = <a-y, F>·<normala,normala>-<a-y,normala>·<normala, F>
            // = <a-y, F>·<normala,normala>-<a-y,tangent×y-a>·<normala,F>
            // = <a-y, F>·<normala,normala>-0·<normala,F>
            // = <a-y, F>·<normala,normala>
            //
            // And hence:
            // la·coeff_dtfa·<(a - y)×normala, F×normala>
            // = la·<derivative,normala>/<normala,normala>·<(a - y)×normala, F×normala>
            // = la·<derivative,normala>/<normala,normala>·<a-y, F>·<normala,normala>
            // = la·<derivative,normala>·<a-y, F>
            // = la·<derivative,tangent×y-a>·<a-y, F>
            // = la·<y-a,derivative×tangent>·<a-y, F>
            // = la·<a-y,F - derivative×tangent>
            //
            // Is that cleaner?
            //
            // As it turns out we really want the dot product with `unit(dt A)` instead, i.e.
            // depending on how we defined it by the cross product.
            let dtfa_dira = dir_dtfa.dot(dir_a) * coeff_dtfa;
            let dtfb_dirb = dir_dtfb.dot(dir_b) * coeff_dtfb;

            // The curvature relative to surface normal is a triple product:
            //
            //     <F×dt A, dt² A>
            //
            // We can arrange that however we want. Most reasonable is dt A×dt² A since this is
            // simply to calculate, we'll see we do not need dt² A explicitly as dt A is itself
            // formed from a dot product and we have the right coefficients.
            //
            // We have, for some c0, that `dt A = c0 · F×Fa` and from `<F, dt A> = 0` we have:
            //
            // (F×Fa)×dt²A
            // = Fa·<F, dt²A> - F·<Fa, dt²A>
            // = F·<dt Fa, dtA> - Fa·<dt F, dtA>
            //
            // With c0
            //   = curve.len_a / ||F×Fa||
            // (from ||dt A|| = |curve.len_a|)
            // (note ||F×Fa|| = sqrt(1 - dot(F, Fa)²), both are unit-length normals)
            //
            // Anyways we want the dot-product of this with both normals for the triple product. For
            // any scalar f (but especially for `c0`):
            //
            // <Fa, f·(F×Fa)×dt²A> = f·(<Fa,F>·<dtFa,dtA>- <dtF,dtA>)
            // <F, f·(F×Fa)×dt²A> = f·(<dtFa,dtA> - <Fa,F>·<dtF,dtA>)
            //
            // Now, we only want the curvature or rather its dt/ds (where ds is the _implied_
            // unit-speed parameterization) adjusted term. These are related to the triple product
            // above by noting that `dt A` is `||dt A||·dA/ds` and `dt² A = ||dt A||²·dA/ds + …` so
            // we'd have to divide by `||dt A||³` here. Then we multiply again with one though. That
            // is, we really want `f* = c0/||dt A||²` here. But `|curve.len_a| = ||dt A||` by
            // definition and so:
            //
            // <F, f*·(F×Fa)×dt²A>
            //     = f*·(<dtFa,dtA> - <Fa,F>·<dtF,dtA>)
            //     = curve.len_a/||F×Fa||/||dt A||²·(<dtFa,dtA> - <Fa,F>·<dtF,dtA>)
            //     = curve.len_a/||F×Fa||/||dt A||²·curve.len_a·(<dtFa,dir a> - <Fa,F>·<dtF,dir a>)
            //     = curve.len_a²/||dt A||²/||F×Fa||·(<dtFa,dir a> - <Fa,F>·<dtF,dir a>)
            //     = (<dtFa,dir a> - <Fa,F>·<dtF,dir a>)/||F×Fa||
            //
            // This implies the curvature does not depend on curve.len_a (as intuition expects).
            // However, of course the derivative influences the directions of `Fa` so in an ODE
            // sense there still is an influence—just not on the linearization at this timestep.
            //
            // NOTE: interestingly `dir_a` already is `F×Fa/||F×Fa||` which contracts into the `f0a`
            // factor as `||F×Fa||² = 1.0 - <F, Fa>²`. And of course we again have another set of
            // rather symmetrical triple products `<dtFa,F×Fa>` and `<dtF,F×Fa>` where of course
            // each derivative of a normal is itself another rescaled cross product. However, I
            // don't readily see this simplifying our calculation. We need the by-product of `dt A`
            // itself and this adds more summation terms. Unless it cancels with `f0a` this seems
            // not worth it.

            assert!((normalf.length() - 1.0).abs() < 1e-6);
            assert!((normala.length() - 1.0).abs() < 1e-6);
            assert!((normalb.length() - 1.0).abs() < 1e-6);
            assert!((normalf.dot(dtf)).abs() < 1e-6);

            // Canceled from `curve.len_a / ||F×Fa|| / ||dt A||²` by adjusting the products.
            let f0a = normalize_faf;
            // Note: negative factor since dt B is from Fb×F instead. This is the simpler adjustment
            // to the formula above, keeping `raw_curve_at2…` symmetrical with the case of `dt A`.
            let f0b = -normalize_fbf;

            // Second component for the curvature calculation. This implies that the yaw _does_
            // relate to the curvatures but we can only steer them highly coupled.
            let dtf_dira = dtf.dot(dir_a);
            let dtf_dirb = dtf.dot(dir_b);

            let raw_curve_at1fa = f0a * (faf * dtfa_dira - dtf_dira);
            let raw_curve_at1f = f0a * (dtfa_dira - faf * dtf_dira);
            let raw_curve_at2f = f0b * (dtfb_dirb - fbf * dtf_dirb);
            let raw_curve_at2fb = f0b * (fbf * dtfb_dirb - dtf_dirb);

            // For the base curve we calculate the dot product from the explicit second derivative.
            let base_curve = tangent.cross(derivative);
            let raw_curve_at0fa = normala.dot(base_curve);
            let raw_curve_at0fb = normalb.dot(base_curve);

            let speedy = tangent.length();
            let speeda = curve.len_a;
            let speedb = curve.len_b;

            let speed_cory = 1.0 / tangent.length_squared();

            // 0.0 points towards +X.
            let orientation = |x: f64| {
                let (s, c) = x.sin_cos();
                DVec2::new(c, s)
            };

            // Fill in all the derivatives.
            let diff = Params {
                curves: {
                    // Clarify this explicit notation is not required anywhere else. This should
                    // simplify the development of a better parameterization by letting it transfer
                    // to this one under the hood just with better control.
                    let dta = dir_a * curve.len_a;
                    let dtb = dir_b * curve.len_b;
                    [dta, dtb]
                },
                opposing_normal: dtf,
                flats: {
                    let [of, ofa, ofb] = params.flat_orientation;

                    [
                        [orientation(of[0]) * speeda, orientation(of[1]) * speedb],
                        [orientation(ofa[0]) * speedy, orientation(ofa[1]) * speeda],
                        [orientation(ofb[0]) * speedb, orientation(ofb[1]) * speedy],
                    ]
                },
                flat_orientation: [
                    [raw_curve_at1f, raw_curve_at2f],
                    [raw_curve_at0fa * speed_cory, raw_curve_at1fa],
                    [raw_curve_at2fb, raw_curve_at0fb * speed_cory],
                ],
            };

            (diff.put(), true)
        }
    }

    impl Params {
        fn read(Coord(coeffs): &Coord<27>) -> Self {
            let (curves, coeffs) = coeffs.split_first_chunk::<6>().unwrap();
            let (normal, coeffs) = coeffs.split_first_chunk::<3>().unwrap();
            let (flats, coeffs) = coeffs.split_first_chunk::<12>().unwrap();
            let (orients, _) = coeffs.split_first_chunk::<6>().unwrap();

            Params {
                curves: {
                    let &[a, b] = curves.as_chunks::<3>().0.as_array().unwrap();
                    [DVec3::from_array(a), DVec3::from_array(b)]
                },
                opposing_normal: DVec3::from_array(*normal),
                flats: {
                    let &[y, a, b] = flats.as_chunks::<4>().0.as_array().unwrap();
                    let as_starts =
                        |[c0, c1, c2, c3]: [f64; 4]| [DVec2::new(c0, c1), DVec2::new(c2, c3)];
                    [as_starts(y), as_starts(a), as_starts(b)]
                },
                flat_orientation: *orients.as_chunks::<2>().0.as_array().unwrap(),
            }
        }

        fn put(&self) -> Coord<27> {
            let mut c = [0.0; 27];

            let coeffs = &mut c[..];
            let (curves, coeffs) = coeffs.split_first_chunk_mut::<6>().unwrap();
            let [a, b] = curves.as_chunks_mut::<3>().0.as_mut_array().unwrap();
            *a = self.curves[0].to_array();
            *b = self.curves[1].to_array();

            let (normal, coeffs) = coeffs.split_first_chunk_mut::<3>().unwrap();
            *normal = self.opposing_normal.to_array();

            let (flats, coeffs) = coeffs.split_first_chunk_mut::<12>().unwrap();
            let put_flats = |slice: &mut [f64], [vecl, vecr]: [DVec2; 2]| {
                *slice.as_chunks_mut::<2>().0.as_mut_array().unwrap() =
                    [vecl.to_array(), vecr.to_array()];
            };

            let [y, a, b] = flats.as_chunks_mut::<4>().0.as_mut_array().unwrap();
            put_flats(y, self.flats[0]);
            put_flats(a, self.flats[1]);
            put_flats(b, self.flats[2]);

            let (orients, _) = coeffs.split_first_chunk_mut::<6>().unwrap();
            *orients.as_chunks_mut::<2>().0.as_mut_array().unwrap() = self.flat_orientation;

            Coord(c)
        }
    }

    let x0 = Params {
        curves: [tr.base1, tr.base2],
        opposing_normal: tr.opposing_normal,
        flats: {
            let flat_to_pos = |fl: &PipeFaceBase| [fl.base_left, fl.base_right];
            [
                flat_to_pos(&tr.opposing_flat),
                flat_to_pos(&tr.flat1),
                flat_to_pos(&tr.flat2),
            ]
        },
        flat_orientation: {
            let flat_to_pos = |fl: &PipeFaceBase| [fl.orientation_left, fl.orientation_right];
            [
                flat_to_pos(&tr.opposing_flat),
                flat_to_pos(&tr.flat1),
                flat_to_pos(&tr.flat2),
            ]
        },
    }
    .put();

    let ode = Ode(tangent);

    let sol = fast_ode::solve_ivp(
        &ode,
        (f64::from(start), f64::from(end)),
        x0,
        |_, _| true,
        1e-6,
        1e-6,
    );

    let x1 = match sol {
        fast_ode::IvpResult::FinalTimeReached(coord) => coord,
        // Deserves a warning, at least!
        fast_ode::IvpResult::StepTooSmall(_, coord) => coord,
        fast_ode::IvpResult::OdeRequestedExit(..)
        | fast_ode::IvpResult::CallbackRequestedExit(..) => {
            unreachable!("we do not request exit")
        }
    };

    let params = Params::read(&x1);
    let end = (ode.0)(f64::from(end));

    eprintln!("{:?}", params.opposing_normal);

    let pipe = PipeSegment {
        base1: params.curves[0],
        base2: params.curves[1],
        base: end.frame.base,
    };

    const UNFINISHED: PipeFaceSegment = PipeFaceSegment {
        normal: SurfaceNormal { axis: DVec3::X },
        base_left: DVec2::ZERO,
        base_right: DVec2::ZERO,
        orientation_left: 0.0,
        orientation_right: 0.0,
        curvature_left: 0.0,
        curvature_right: 0.0,
    };

    TrianglePipeSegment {
        flat1: PipeFaceSegment {
            normal: SurfaceNormal {
                axis: end
                    .frame
                    .tangent
                    .cross(pipe.base - pipe.base1)
                    .normalize_or_zero(),
            },
            base_left: params.flats[1][0],
            base_right: params.flats[1][1],
            orientation_left: params.flat_orientation[1][0],
            orientation_right: params.flat_orientation[1][1],
            ..UNFINISHED
        },
        flat2: PipeFaceSegment {
            normal: SurfaceNormal {
                axis: end
                    .frame
                    .tangent
                    .cross(pipe.base2 - pipe.base)
                    .normalize_or_zero(),
            },
            base_left: params.flats[2][0],
            base_right: params.flats[2][1],
            orientation_left: params.flat_orientation[2][0],
            orientation_right: params.flat_orientation[2][1],
            ..UNFINISHED
        },
        opposing_flat: PipeFaceSegment {
            normal: SurfaceNormal {
                axis: params.opposing_normal,
            },
            base_left: params.flats[0][0],
            base_right: params.flats[0][1],
            orientation_left: params.flat_orientation[0][0],
            orientation_right: params.flat_orientation[0][1],
            ..UNFINISHED
        },
        pipe,
    }
}
