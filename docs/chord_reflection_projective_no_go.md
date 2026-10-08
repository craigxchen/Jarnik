# Chord reflections and anchor swaps cannot supply a nonconformal endpoint map

This note records a no-go result for a specific proposed construction. It does
not prove the uniform endpoint bound. The statements concern **one fixed
reflection** or a fixed finite family of conformal operations. They do not
exclude a projective map built from genuinely new arithmetic data.

Let `q_1,...,q_m` be distinct rational points of the unit circle, with an
anchor `q_0=1`, and let `q(x)=(x+iL)/(x-iL)` be the cotangent chart. Swapping
the anchor with `q_a` and reversing order induces

```text
S_a(q)=q_a/q.
```

For two actual points `q_a,q_b`, reflection in their chord's perpendicular
bisector through the circle center exchanges them and induces

```text
S_ab(q)=q_a q_b/q.
```

Both maps are rational projective maps in the cotangent coordinate. For
example, if `s=(X_a X_b-L^2)/(X_a+X_b)`, then `q(s)=q_a q_b` and
`S_ab(q(x))=q((s x+L^2)/(x-s))`. The parameter `s` may be rational rather
than integral; clearing its denominator changes the coordinate display, not
the projective map.

Every composition of these operations has the form `q -> kappa q` or
`q -> kappa/q`, with `|kappa|=1`. It is therefore a circle rotation or
reflection. On a *whole* tuple, common rational Gaussian scaling clears all
denominators, and primitive Gaussian normalization identifies the result with
the original tuple or its conjugate up to a Gaussian unit. Thus the primitive
least squared radius and shortest angular span are exactly preserved. In
particular these compositions never meet the nonconformality hypothesis of
the conditional growth theorem.

There is also a rigidity result for selectively applying one reflection.
Let `U(q)=alpha q` and `V(q)=beta/q` be two fixed rational projective circle
isometries (the same proof works for any two fixed projective maps). Suppose
`m>=5` distinct source points satisfy

```text
T(q_i) in {U(q_i),V(q_i)}   for every i
```

for a single rational projective map `T`. Assign each point to one branch,
choosing either branch if the outputs coincide. At least three points have
the same assignment. A projective map agreeing with another projective map
at three distinct source points equals it identically. Hence `T=U` or
`T=V`; it is conformal. Any points assigned to the other branch must solve
`U(q)=V(q)`, or `q^2=beta/alpha`, so there are at most two such points.
The count five is sharp for the three-point uniqueness argument: with four
points a two-versus-two assignment need not force either map. For `L=1`,
use cotangents `1,-1,2,1/2`, take `U(x)=x`, `V(x)=-x`, and
`T(x)=(5x-4)/(5-4x)`. Then `T` fixes `1,-1`, sends `2` to `-2`, and
`1/2` to `-1/2`, while `T` is neither `U` nor `V`.

More generally, if each output is selected from a **fixed** family of `K`
projective conformal maps, then `m>=2K+1` forces `T` to equal one member by
pigeonhole and three-point uniqueness. This does not address a family of
reflection axes whose size grows with the tuple.

Adjoining physically reflected points to the original tuple is different
from replacing the whole tuple. Its least squared radius is `N D`, where
`N` is the original squared radius and `D` is the reduced Gaussian
denominator norm of the reflection coefficient. That exact cost and the
associated endpoint radial-count obstruction are already proved in
`mirror_union_least_radius_cost.md` and
`reflection_union_radial_count_obstruction.md`; no radius saving follows
from the selective-map observation above.

The four-point sharpness fixture is checked in
[the rational inversion checker](check_actual_point_double_inversion.py).
