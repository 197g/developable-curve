Notation:
Using F for the normal of the developed surface
Using T for the tangent of a general, non-degenerate curve
Using ĸ for the krümmung of the curve
Using ŧ for the curl for typing convenience
Using dt X for `dX/dt` and proper fraction for other partial derivatives
Using <a, b> for the dot product

F = aN + bB
dt F = c(bN - aB) + …T

(1):
a = <F, N>
dt a = dt <F, N> = <dt F, N> + <F, dt N> = 
    cb + <F, ŧB + ĸT> =
    cb + <aN + bB, ŧB + ĸT> = 
    cb + bŧ =
    b(c + ŧ)

(2):
b = <F, B>
dt b = dt <F, B> = <dt F, B> + <F, dt B> =
    -ca + <F, -ŧN> =
    -ca + <aN + bB, -ŧN> =
    -ca - ŧa = 
    -a(c + ŧ)

The result of this looks wrong… but it's so simple? Anyways, this implies a
simple formula for the rotation. Note that the angular rate here is given
simply by `c + ŧ`. Where `c` is just the length of the derivative of the
surface normal without the tangent component (i.e. `derivative_free` in the
code). That angular rate is in relation to the given parameterization of the
curve, not angular rate with unit speed.

But also, this requires us to have an impossible angle of attack if we want to
cross the binormal. This is because the tangent component of `dt F` has indeed
length `-<dt T, F>/||T||` if we are to keep F and T orthogonal. `<dt T, F>` is
also directly related to `a` since the `T` component of `dt T` vanishes:

    <dt T, F> = <ĸN · <T, T>, F>

Measuring `dt T` under `F` removes its `T` component (from the derivative of
curvature speed), leaving us with `||T||²` times the unit-speed second
derivative, which defines `ĸ` through the relation with the unit normal `N`.
Anyways, we have

    -<dt T, F>/||T|| = -<ĸN, F>·||T|| = -ĸa · ||T||

When ĸ or a goes towards `0` but we keep the same orthogonal component to our
derivative vector, the angle goes steeper and steeper and the ruling direction,
orthogonal to the derivative, gets parallel with the tangent—leaving us without
a well-defined developed surface at the cusp. In all practical applications we
do not want this. Instead, a constant angle should be maintained. As the above
also tells us, the total signed change in angle ø is also readily computed from
the curve in this way. (What it also tells us: with `a = 0` or `ĸ = 0` we will
not have a change is angle in any case).

To do the inverse, with a ruling angle C, we get a differential equation. With
the parameterization of a = cos(ø) (from (1) due to: da/dø = -sin(ø) = -b)

    dt ø = ĸ · cos(ø) · ||T|| · C + ŧ

This would be simple to solve if ĸ, ||T||, ŧ were also constants, alas. For a
constant speed helix (ĸ, ŧ constant) we also have an interesting, far less
general, problem. We'll fold away the factors with u = ĸ · ||T|| · C.

    dt ø = u cos ø + ŧ

If we give up a constant angle, we can also derive a control equation if we
just define dt ø. (Remember C = tan angle if you ponder this further).

    C = (dt ø - ŧ) / (ĸ·cos ø·||T||)

This begs the question, can we define dt ø - ŧ such that it usefully reduces,
e.g. with cos ø? For ŧ=0 most certainly there are a lot of options.
