# Independent audit of the seven-row quotient obstruction

The arguments in
[the seven-row note](seven_row_quadratic_quotient_obstruction.md) are
correct on the stated geometric reading.  A boundary point of a
four-label forgetful target has exactly
`2^(n-4)=2^(m-1)` stable boundary divisors in its pullback, so balance
gives

```text
[k(C):k(x_i)]=2^(m-1).
```

More generally, forgetting a set `S` pulls a chosen boundary divisor
back to `2^|S|` distinct degree-one boundary divisors.  Factoring the
map through the normalization of its image therefore gives

```text
r_S divides 2^|S|.
```

The deck-transformation argument also checks exactly.  Distinct
single-label involutions are separated by anchored coordinate
functions.  Two of them fix the double-forgetful field and generate a
group of order at most four, so they commute and generate a Klein four
group.  With three anchors outside any chosen set of at most `n-3`
labels, a nonempty product relation would make one involution fix its
missing coordinate as well as all the others.  Hence those
involutions are independent over `F_2`.  Four of them at `n>=7` would
contradict the rank-three bound for an elementary two-subgroup of the
automorphism group of a genus-one curve.

For seven labels and three such involutions, the field identities

```text
k(d)=k(C)^G,
k(a,d)=k(C)^(<tau_b,tau_c>)
```

have the asserted degrees.  Thus the image in
`P^1_a times P^1_d` has both projection degrees two: it is a genuine
correspondence of bidegree `(2,2)`, not a rational quadratic relation
`d=f(a)`.

The squareclass statement should be read geometrically, after passage
to an algebraic closure: independence over the original constant
field alone can be caused by nonsquare constants and can disappear
after base change.  With this convention, the genus computation and
the final fixed-curve application of Faltings are sound.  The note
correctly makes no claim of a family-uniform estimate or a global
uniform `sqrt(R)` lattice-arc bound.
