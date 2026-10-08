# Repeated central weights force three reflection pairs

Let six distinct complex numbers `q1,...,q6` have modulus one, and let
`a,b` be nonzero real numbers. Suppose

```text
u=(a,-a,a,-a,b,-b),
sum_i u_i q_i^j=0,       j=0,1,2.                    (1)
```

Then reflection in the diameter bisecting `q5,q6` pairs all six nodes.
More precisely, put `t=q5 q6`. Exactly one of the following holds:

```text
q1 q2=q3 q4=q5 q6=t,
q1 q4=q3 q2=q5 q6=t.                                (2)
```

No short-arc hypothesis or exclusion of `b/a=±1` is needed.

## Proof by real parts

Rotate all nodes by a common unit so that `q5,q6` become conjugates.
The rotation preserves (1). Write `r_i=Re(q_i)` in this rotated frame.
Taking real parts of the first and second moments cancels the last pair.
Since `Re(q_i^2)=2r_i^2-1`, it follows that

```text
r1+r3=r2+r4,       r1^2+r3^2=r2^2+r4^2.
```

Thus the two pairs have the same sum and product, and hence
`{r1,r3}={r2,r4}` as multisets. Distinct points on the unit circle with
the same real part are conjugates. Global distinctness rules out matching
a point with itself, so `{q2,q4}={conjugate(q1),conjugate(q3)}`.
Undoing the rotation gives (2). The matching is unique: if both matchings
held, cancellation would give `q2=q4`.

There is also a square-root-free algebraic version. Set
`phi(z)=z+t/z`. The unit-circle and conjugate moment identities give
the same sums and sums of squares of `phi(q1),phi(q3)` and
`phi(q2),phi(q4)`, because `phi(q5)=phi(q6)`. Therefore

```text
0=(phi(q1)-phi(q2))(phi(q1)-phi(q4))
 =(q1-q2)(q1-q4)(q1 q2-t)(q1 q4-t)/(q1^2 q2 q4).
```

Distinctness and nonvanishing leave precisely the two alternatives (2).
The coefficient ratio `b/a` cancels before this factorization, so the
exceptional ratios in a generic elimination argument cause no gap.

## Consequence for an integer-circle arc

Suppose the six nodes arise by common scaling from distinct Gaussian
integers `z_i` of radius `R`. Equation (2) supplies a nontrivial equality
between products of two points. If the points lie in an arc of angular
width `theta<pi`, the sharp
[multiplicative rectangle theorem](multiplicative_rectangle_separation.md)
therefore gives

```text
theta sqrt(R) >= 2 sqrt(2).                          (3)
```

Thus this repeated-weight pattern cannot occur in a smaller endpoint arc.
It explains the product collisions in the earlier
[fixed-weight elliptic family](six_point_fixed_moment_elliptic_obstruction.md).

A single opposite pair of coefficients does not suffice. The rational
cotangent parameters `(2,3,5,12,17)` have primitive central vector
`(126,-7,50,-169,841,-841)` and have no nontrivial equality between
unordered pair products. The checker below verifies this distinction
directly, including repeated factors in the product test.

Run `python3 docs/check_six_point_repeated_weight_reflection.py`. It checks
the exact moments and both reflection matchings on rational members of the
elliptic family, and checks the single-opposite-pair counterexample. The
real-part argument above is the proof.
