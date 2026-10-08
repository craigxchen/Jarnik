# Six Pell points attain the full-clipping divisor bound

The divisor theorem in the [reciprocal-stretch grid](mobius_reciprocal_stretch_grid.md)
gives `m<=2tau(4)=6` for a nonconformal determinant-one map with
`N'=N=Nclip`. The following infinite six-point family attains equality.
It has bounded normalized arc length tending to `4sqrt(2)`, but does not
meet the endpoint-constant-at-most-two hypothesis of the critical reduction.
It is a sharpness family for the conditional map theorem, not an unbounded
point-count construction.

Take a positive odd solution of

```text
k^2-2b^2=-1,       b=1,5,29,169,...,
```

and put `M=((1,b),(0,1))`, `delta=1`, and `N=b^2(b^2+4)`. In the
actual source anchor frame use the six rational Gaussian half-angle rows

```text
H_0=1,
H_1=-b+2i,
H_2=k-2b+i,       H_3=-k-2b+i,
H_4=b+k+i,        H_5=b-k+i.
```

Choose physical source reference `z_0=-2b+ib^2` and physical target
reference `w_0=2b+ib^2`. The source and target rows are respectively
`z_i=z_0 H_i/conjugate(H_i)` and
`w_i=w_0 (MH_i)/conjugate(MH_i)`. Exact substitution of the Pell equation
gives the following correspondence:

| row | source `z_i` | target `w_i` | `x_i=Norm(MH_i)/Norm(H_i)` | `n_i=2x_i` |
|---|---|---|---:|---:|
| 0 | `-2b+ib^2` | `2b+ib^2` | 1 | 2 |
| 1 | `2b+ib^2` | `-2b+ib^2` | 1 | 2 |
| 2 | `k+i(b^2+1)` | `-2k+i(b^2-2)` | 1/2 | 1 |
| 3 | `-k+i(b^2+1)` | `2k+i(b^2-2)` | 1/2 | 1 |
| 4 | `-2k+i(b^2-2)` | `k+i(b^2+1)` | 2 | 4 |
| 5 | `2k+i(b^2-2)` | `-k+i(b^2+1)` | 2 | 4 |

Thus the target set is the source set with a different labeling. All
six points are distinct and have squared norm `N`: for the middle pair,
`k^2+(b^2+1)^2=N`; for the outer pair,
`4k^2+(b^2-2)^2=N`. The source and target configurations are primitive.
Indeed `gcd(b,k)=1` by the Pell equation. At every prime of `N`, a pair
of rows occupies opposite Gaussian orientations, as verified below, so
no Gaussian prime divides all six rows. Their exact least squared norms
are therefore `N'=N`. Primitivity here is collective Gaussian
primitivity. Individual points can have ordinary coordinate content
greater than one; the two anchor points have content `b`.

The two anchor gauges are different. In a single source half-angle chart,
the induced self-permutation has representative
`R_(H_1) M=((-b,-b^2-2),(2,b))`, whose trace is zero and whose
projective order is two. The unipotent matrix `M` by itself is not a
self-permutation in that one fixed chart.

To check full clipping, write the shear coefficients as

```text
U=2-ib,       V=ib,       z_0=-bU.
```

The odd prime supports of `b` and `b^2+4` are disjoint. Every odd prime
on these supports is split: if `p|b`, then `k^2=-1 mod p`; if
`p|b^2+4`, then `(b/2)^2=-1 mod p`. If `p|b` and `t=v_p(b)`, `U` is
a unit at both orientations, `V` has valuation `t` at both, and the
coefficient interval is `[-t,t]`. The anchor `z_0` has valuation `t`
at each orientation. The rows `z_2,z_3`, which are negative conjugates,
have ordinary content prime to `p` and occupy opposite orientations;
the source interval is also `[-t,t]`. Its clipped width is `2t`.

If `p|b^2+4` and `e=v_p(b^2+4)`, `V` is a unit and `U` has valuation
`e` at exactly one orientation. Since `z_0=-bU`, the coefficient
interval and source interval are both `[-e,0]` or both `[0,e]`,
depending on the chosen orientation. The conjugate pair `z_0,z_1`
attains both source endpoints. Its clipped width is `e`. Because `b`
is odd, `N` has no factor two. Therefore

```text
Nclip=b^(2) (b^2+4)=N=N'.
```

The divisor theorem now has `H=4delta^2NN'/Nclip^2=4`. Its three
positive divisors `1,2,4` each occur twice in the table, so
`m=6=2tau(4)` exactly. For `b>=5`, all six points lie in the upper
half-plane. Their shortest containing arc has angular width
`2arctan(2k/(b^2-2))`, giving normalized length

```text
C_b=2 N^(1/4) arctan(2k/(b^2-2)) -> 4sqrt(2).
```

The target has the same arc. Thus the family is geometrically bounded
but has a larger constant than the at-most-two critical hypotheses.
The [exact checker](check_full_clipping_six_point_pell_sharpness.py)
verifies four Pell instances, every source/target row and stretch, the
Gaussian primitive radius, and the local clipped widths at all primes
of `N`.
