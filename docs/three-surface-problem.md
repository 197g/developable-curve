Given the curve Y(t) a most interesting question is to construct two additional
curves Y1, Y2 such that the surface between each pair of them is its own
developable surface. The idea here is that this will uniquely assemble from
three planar sheets into the three dimensional configuration (we'll demonstrate
later what condition is sufficient given three krümmungen ĸ,ĸ1,ĸ2).

A single plane is floppy, indeed given any particular curvature k of the 2d
boundary curve there are two dimensions of possible spatial curves allowable
for it. The second derivatives of Y are constrained into a line at each point
by `dt T·(F×T)=ĸ`. (The component of dt T in direction T is determined from
`||T|| · dt ||T||`). Then the parameterization of the ruling direction (or dt F
except for the degenerate case where <dt T, F> = 0) are again independent. If
we are given the ruling directions in 2d then we have exactly one dimension of
freedom at each point).

This directly demonstrates that two _independent_ such constraints uniquely
determine the spatial curvature and thus the whole curve. We calculate the meet
of the two linear constraints which uniquely exists as long as directions `F`
are not parallel. Then a set of three such surfaces give two such constraints
at each of their pairwise connections—so they form a unique triangular pipe.

Let us setup initial parameters for one constraint problem. It seems
interesting in an engineering sense if we could envelope some given surface
with a rigid scaffolding (e.g. a splint). Start from Y and A(0), B(0), and F(0)
where F is the normal of the surface opposite to `Y`. We'll trivially find

    Fa(0) = normalize dt Y×(A(0) - Y(0))
    Fb(0) = normalize dt Y×(B(0) - Y(0))

Then as a basic constraint we must have dt A || Fa×F, dt B || Fb×F, and finally
dt Y || Fa×Fb so that the respective surfaces are developable. From the same
argument it follows that dt F, dt Fa, dt Fb are orthogonal to their edges, that
is to: `A(0) - B(0)`, `A(0) - Y(0)` and `B(0) - Y(0)` respectively as well as
each normal itself. For dt Fa and dt Fb we also have the bounds from dt² Y:

    <dt² Y, Fa> = -<dt Y, dt Fa>

Except for the degeneracy where these two conditions, so the edge and dt Y, are
parallel this uniquely determines dt Fa and dt Fb (of course). This then
further links us to the other curves since, just the same, it holds that:

    <dt² A, Fa> = -<dt A, dt Fa>
    <dt² B, Fb> = -<dt B, dt Fb>
    <dt² A, F> = -<dt A, dt F>
    <dt² B, F> = -<dt A, dt F>

Note that this confines dt² A to a plane orthogonal to dt A and dt F (same for
B). As long as those are independent, the only available degree of freedom is
thus a scalar multiple of A. This is enough for our purposes of deriving a 2d
curve equivalent of all sides; for which we need only the flat curvature given
from the triple product of the forms `<dt A × dt² A, F>` (over the six relevant
combinations of A, B, Y and F, Fa, Fb). The first term here is independent of
the unknown component of dt² A, we'll just suppose zero, i.e. corresponding to
an unchanging length of curve speed, and redo it at each point. Hence, three
parameters is enough to ODE-integrate A, B, F, Fa, Fb and their 2d equivalents.
We need only the lengths of dt A, dt B, and dt F.
