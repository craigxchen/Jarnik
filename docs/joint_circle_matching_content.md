# Matching height with the retained Gaussian common content

The general inequality `h_matching <= (3/2) log R + 3 log C`
discards a potentially large common Gaussian factor. This note retains
that factor exactly and distinguishes its generic calculation from
actual arithmetic. The resulting method proves
[standalone radius inflation](joint_involution_standalone_inflation.md)
for the joint map, using a smaller actual content bound in place of
the unjustified generic total.

## 1. The common factor of all matching products

Let `z_1,...,z_6` be distinct Gaussian integers of the same norm `R^2`,
in an arc of length at most `C sqrt(R)`. For a perfect matching `M`,
put `P_M=product_(ij in M)(z_i-z_j)`. All fifteen products have the
same complex phase up to sign. Their ratios are rational, so

```text
P_M=gamma m_M,       gcd_M m_M=1,       gamma in Z[i].       (1)
```

The integrality of `gamma` follows by an integer Bezout combination.
The same combination shows that the Gaussian gcd of all the `P_M`
is `gamma`, up to a unit. In particular,

```text
h_matching + log|gamma| <= (3/2) log R + 3 log C.           (2)
```

There is no assumption that `gamma` is real.

## 2. A local lower bound for that factor

Fix a split prime `p=pi*bar(pi)`. Write `alpha_i=v_pi(z_i)` and
`n=v_p(R^2)`. Equality of norms gives `v_bar(pi)(z_i)=n-alpha_i`
for every row. Put the six `alpha_i` in increasing order, with
repetitions allowed. The ultrametric inequality gives

```text
v_pi(P_M) >= sum_(ij in M) min(alpha_i,alpha_j).
```

The minimum over matchings is the sum of the smallest three entries:
these entries can be paired with the largest three, and no choice of
three entries has smaller sum. Applying the same argument to the
conjugate valuations proves

```text
v_pi(gamma) >= alpha_(1)+alpha_(2)+alpha_(3),
v_bar(pi)(gamma) >= 3n-alpha_(4)-alpha_(5)-alpha_(6).       (3)
```

Thus this rational prime contributes at least

```text
(1/2)[3n + sum_(i<=3) alpha_(i)
                 - sum_(i>=4) alpha_(i)] log p            (4)
```

to `log|gamma|`. Extra cancellation in a chord only increases the
bound. Contributions at other primes may be omitted, since they
are nonnegative.

If a seventh anchor is retained, replace `alpha_i` by the relative
valuations `d_i=alpha_i-alpha_0` in the difference of sums in (4).
For the least primitive circle at this prime,

```text
n=max(0,d_1,...,d_6)-min(0,d_1,...,d_6).                  (5)
```

Formula (5) is just simultaneous integrality of all seven phases;
it is the local least-radius construction, not an extra assumption
about the six output directions.

For a simple cut with `d_i` equal to `e` on a set `T` and zero off it,
and with `T` nonempty, (4)--(5) reduce to

```text
contribution to log|gamma| >= (| |T|-3 |/2) e log p.       (6)
```

Summing over the prime factors of a norm block `n_S` gives
`(| |T|-3 |/2) log n_S`. If `T` is empty, this layer is absent
from the least circle and contributes zero, not `3 log n_S/2`.

## 3. What the generic joint-map table would imply

In this paragraph the output circle retains the original marked anchor
as a seventh point. The exact polynomial minimum orders in
[joint_circle_lift_cut_valuations.json](joint_circle_lift_cut_valuations.json)
give fifty-six retained simple old cuts in the generic lift. The
numbers of output subsets of sizes `1,2,3,4,5,6` are respectively

```text
6, 12, 12, 12, 6, 8.
```

In particular the eight size-six cuts remain nontrivial relative to
the retained anchor. For a standalone six-point circle one can remove
their common allocation; those eight cuts then contribute zero, and
the corresponding generic total below is `24w`, not `36w`.

Consequently their contribution to (6), if these are the actual
primitive output allocations, is `36w+o(w)`. Combining this with
the actual matching-height lower bound `42w-o(w)` would give

```text
log R_out >= 52w-o(w)                                   (7)
```

for a bounded endpoint arc constant. The fixed original first three
directions have angular separations `exp(-16w+o(w))`, so such an
endpoint realization in these same directions would instead require
`log R_out<=32w+o(w)`. Thus the generic pattern cannot give the
desired endpoint realization.

**The premise about actual allocations is essential.** At a fixed old
prime, extra cancellation in the degree-seven covariants can alter
their primitive Gaussian phase valuations. The local table records
polynomial minimum orders; it does not rule out those increases.
Equations (1)--(6) are unconditional. Equation (7) is only the stated
conditional consequence. The completed local audit instead proves an
actual lower bound of `6.5w-o(w)` with the anchor retained, and
`8w-o(w)` even without it. These sufficient bounds use selected cuts
and exact complementary-cut symmetry; they do not assert that the
generic total survives. See
[the complete inflation proof](joint_involution_standalone_inflation.md).

## 4. An actual local loss of the generic contribution

The loss in the preceding warning really can occur. Take primitive
integer binary rows `P_i=X_i+i`, where

```text
(X_1,...,X_6)=(8,2,0,21,1,9),       pi=3+2i, N(pi)=13.
```

Only `P_1,P_4` are divisible by `pi`, each to order one, and every
`bar(P_i)` is a `pi`-unit. Among the fifteen old brackets, only
`Delta_14` is divisible by thirteen, also to order one. Thus the
old local data have precisely the cut `S={1,4}` without an extra
old pair residue at this prime.

The exact matching and linear-form values include

```text
(A,B,C,F,H,V)=(1008,-1080,208,1288,-872,1144).
```

In particular `V=13*88`, while `F,H` are thirteen-adic units.
The three degree-seven Gaussian lifts are

```text
Q_4=-27048-5240i,       Q_5=952+1120i,
Q_6=-8568-1232i.
```

Their `pi`-orders are `1,1,0`, and their `bar(pi)`-orders are all
zero. The output relative phase allocation is therefore
`(1,0,0,1,1,0)`, a balanced three-versus-three cut. Its allocation
contribution in (6) is zero, whereas the generic lift of `S={1,4}`
has output subset `{1,4}` and contribution `(1/2) log13`.

This example also rules out recovering the lost contribution solely
from mandatory output chord residues at that prime. The ordinary gcd
of its fifteen raw output matching determinants is a thirteen-adic
unit. Constructing the actual least integral circle from the six
Gaussian phases and the marked anchor, then taking the Gaussian gcd
of all matching chord products, gives

```text
v_pi(gamma)=v_bar(pi)(gamma)=0.
```

[check_joint_circle_matching_content.py](check_joint_circle_matching_content.py)
checks these claims with exact integer arithmetic, including the
actual equal-norm seven-point realization. This is a local example.
It does not have the required global nearuniform sixty-four-cut
profile and does not contradict a uniform endpoint bound. It shows
why the generic matching-content count cannot be used unchanged in
an argument for that bound.
