# All-anchor rotations and the exact least-radius charge

This note continues [central truncation](endpoint_central_truncation.md) and
the [primitive transition calculation](all_anchor_primitive_transition.md).
It treats a fixed number `m+1` of original directions `0,...,m`, with
`m>=2`. The result
is a limit on reanchoring as a height argument: the anchor changes supply
one coherent minimal integral circle, and the correction needed for it has
sub-block height. The exact pair-core ordering remains a real condition on
an arbitrary proposed one-anchor model, as the transition note shows; it is
automatic for the actual centrally truncated rows and is repairable by
negligible trimming. No endpoint uniformity claim follows here.

## 1. The pair transition is determined by the original Gram matrix

Let `P_0=1` and `P_i=h_0i=X_i+iY_i` for `1<=i<=m`. These are the actual
conjugate-primitive odd Gaussian chord numerators, with their literal unit
choices inherited from the common original unit class. For `a!=b`, put

```text
Z_ab=bar(P_a)P_b=S_ab+i Delta_ab,
g_ab=Norm(gcd_G(P_a,P_b))=gcd_Z(S_ab,Delta_ab),
h_ab=Z_ab/g_ab.                                          (1)
```

The second equality and the integrality and conjugate-primitivity of
`h_ab` follow prime by prime from conjugate-primitivity of `P_a,P_b`;
`g_0b=1` and `h_a0=bar(P_a)`. Equivalently,
`P_b bar(P_a)=g_ab h_ab`. This is the exact transition identity of the
primitive transition note. In particular, for an anchor `a` and two other
labels `b,c`, its row norm and bracket are

```text
N(h_ab)=N_a N_b/g_ab^2,
Im(bar(h_ab)h_ac)=N_a Delta_bc/(g_ab g_ac),             (2)
```

where `N_0=1` and `Delta_0j=Y_j`. Every denominator in (2) clears because
the rows on the left are Gaussian integers. Applying the norm/Pluecker/
one-Heron criterion to a reanchored array therefore gives identities in
the original integral lift, not new independent metric equations.

For the Boolean factors write

```text
A_i=product_(T contains i) H_T,  P_i=K_i A_i,
B=product_(nonempty T subset [m]) H_T,
D_ab=product_(T contains a,b) n_T  (a,b>0),
C_ab=product_(b in T,a notin T) H_T
     product_(a in T,b notin T) bar(H_T),               (3)
```

where membership of `0` in any `T` is false. Then
`bar(P_a)P_b=D_ab K_b bar(K_a) C_ab` for `a,b>0`; the analogous
formula for `a=0` has `D_0b=1`. The actual central truncation gives
`C_ab|h_ab` at every pair. With `K_ab=h_ab/C_ab` and
`G_ab=gcd_G(A_a,A_b)` (literal primewise generator), write
`gcd_G(P_a,P_b)=G_ab L_ab`. Then

```text
Norm(L_ab) K_ab=K_b bar(K_a),
g_ab=D_ab Norm(L_ab),
t_ab=Norm(L_ab) Im(h_ab)       (a,b>0).                (4)
```

Thus `Norm(L_ab)=gcd_Z(S_ab/D_ab,t_ab)` and divides `t_ab`.
The divisibility in the first line is the exact pair-core condition;
it is not implied by arbitrary small coefficients and primitive
one-anchor rows. For the actual rows it is the residual-allocation
ordering proved in the primitive transition note. Under the fixed-`m`
profile, `log|K_ab|`, `log Norm(L_ab)`, and `log|Im h_ab|` are all `o(w)`.

## 2. Every anchor has the same least integral radius

Let the original equal-norm Gaussian points be `z_0,...,z_m`, of squared
radius `R^2`, and let `G=gcd_G(z_0,...,z_m)`. At anchor `a`, a nonzero
integral starting point `u_a` produces all the phases
`u_a h_ab/bar(h_ab)` as Gaussian integers if and only if each
`bar(h_ab)` divides `u_a`. Hence the least possible squared radius is

```text
M_a=Norm(lcm_G(bar(h_ab):b!=a)).                         (5)
```

In fact the lcm in (5) is associated to `z_a/G`, so

```text
M_a=R^2/Norm(G)        for every a.                     (6)
```

Here is a valuation proof that includes arbitrary split-prime powers.
At a split prime `p=pi bar(pi)` write the `pi` exponent of `z_i` as
`a_i` and its `bar(pi)` exponent as `e-a_i`; the common norm forces
the same `e` for all `i`. The `pi` exponent of `bar(h_ab)` is
`(a_a-a_b)_+`, and its `bar(pi)` exponent is `(a_b-a_a)_+`.
Their maxima over `b` are respectively
`a_a-min_i a_i` and `max_i a_i-a_a`, exactly the exponents of
`z_a/G`. Inert and ramified prime powers have the same exponent in
all equal-norm points and disappear from both sides. Units only
alter the chosen generator. Thus reanchoring does not multiply or
raise the least radius: all anchors use the same canonical integral
circle `z_i/G`.

For the core alone define its canonical point at anchor `a` by

```text
B_a=product_(a in T) H_T product_(a notin T) bar(H_T),
B_0=bar(B).                                              (7)
```

Every nonempty cut separates `a` from at least one other label, so
`B_a=lcm_G(bar(C_ab):b!=a)` up to a unit, and
`Norm(B_a)=product_(T!=empty)n_T`. Moreover

```text
B_b/B_a=C_ab/bar(C_ab),
B_a/B_0=A_a/bar(A_a)       (a>0).                       (8)
```

These are rational phase identities; the displayed `B_a` themselves
are Gaussian integers. In particular the all-anchor core radius has
exact logarithmic height `sum_(T!=empty) log n_T`, not a sum of
different anchor-radius heights.

## 3. Exact correction lcm and the rotation cocycle

Let `L_a=lcm_G(bar(h_ab):b!=a)` with generators chosen so that
`L_a=z_a/G`, and suppose the actual pair-core divisibilities hold.
Then `B_a|L_a`, and write `L_a=B_a E_a`. At anchor zero set

```text
E=lcm_G(bar(K_1),...,bar(K_m)).                          (9)
```

There is an exact identity, up to a Gaussian unit,

```text
L_0=bar(B) E,          E_0=E,
E_a=E K_a/bar(K_a),
E_b/E_a=K_ab/bar(K_ab),
M_a=(product_(T!=empty)n_T) Norm(E).                    (10)
```

To prove the first identity, every `bar(P_i)=bar(K_i)bar(A_i)`
divides `bar(B) E`, so `L_0|bar(B)E`. Conversely `bar(B)|L_0`
because every nonempty block occurs in some row. Thus
`L_0=bar(B)E_0` with `E_0|E`. The all-pair core divisibilities imply
`B_a|L_a`, while (6) gives
`L_a/L_0=P_a/bar(P_a)`. Using (8) yields
`E_a=E_0 K_a/bar(K_a) in Z[i]`. Since each `K_a` is
conjugate-primitive, `bar(K_a)|E_0`. Consequently `E|E_0`, proving
`E_0=E`. The remaining formulas follow by taking ratios and norms.
This proof also shows exactly where arbitrary one-anchor mock-ups can
fail: without pair-core divisibility, `bar(B)` can overpay for an
outside row's residual prime and `E_0` can be strictly smaller than
the proposed lcm `E`.

The correction radius is therefore charged once, by the Gaussian lcm
of the original `K_i`, not once per anchor. In particular

```text
0 <= log Norm(E) <= sum_i log Norm(K_i)=sum_i log k_i=o(w),
log M_a=sum_(T!=empty)log n_T+log Norm(E)
       =(2^m-1)w+o(w)                                  (11)
```

for fixed `m`. If the empty core block has norm `n_empty` and the
original common-modulus correction has modulus `D`, central truncation
also gives the exact source-radius charge

```text
R^2=D^2 n_empty product_(T!=empty)n_T,
Norm(G)=D^2 n_empty/Norm(E).                            (12)
```

In the actual primewise truncation, `Norm(E)|D^2`. At a prime of
original exponent `e`, the retained core exponent is `e'=e-2r`.
If its cut is nonempty, the allocation range across the selected
rows is `e'+v` with `0<=v<=2r`; if its cut is empty, the range is
at most `2r`. This range is the exponent of the common least
circle in (6), while `e'` contributes to the nonempty core only
in the first case. Thus the prime exponent of `Norm(E)` is at most
`2r`, exactly the exponent of `D^2`. Equation (12) consequently
identifies the common source gcd as one empty block times the
unspent part of the common-modulus correction, in norm.

Thus the apparent loss of one block in the least circle is the
common empty block; all nonempty cut blocks and the correction lcm
are paid exactly once.

## 4. What common rotations can change

For clarity, the common Gaussian factor at anchor `a` also has an
exact core charge. Let `F_a=gcd_G(h_ab:b!=a)`, let `J_0=H_[m]`,
and let `J_a=bar(H_{ {a} })` for `a>0`. Since the core numerators
`C_ab` have gcd `J_a`, the actual pair-core divisibilities give

```text
J_a | F_a,       Norm(F_a)/Norm(J_a)
               | product_(b!=a) Norm(K_ab).            (13)
```

The common-rotation gcd identity and (2) also express its norm solely
through the original Gram entries:

```text
Norm(F_a)=gcd_Z(
  {N_a N_b/g_ab^2 : b!=a},
  {N_a Delta_bc/(g_ab g_ac) : b<c, b,c!=a}).
```

For any core prime other than the singleton cut at this anchor,
some row `b` omits it; the remaining valuation of the gcd must come
from that row's `K_ab`. At the singleton cut, any excess above its
core exponent must likewise come from the corrections. This proves
the second divisibility prime by prime. Therefore
`Norm(F_a)=n_{chi_a} exp(o(w))=exp(w+o(w))`, where
`chi_0=[m]` and `chi_a={a}` otherwise. The `m+1` full-set factors
use distinct core cuts, but they are factors of the one canonical
circle in (10); multiplying their radii would count that circle
repeatedly.

Fix the original oriented Gram matrix. By the common-rotation
classification, every integral lift is `P_i'=H'Q_i` with
`Norm(H')=d_0`, where `d_0=Norm(gcd_G(P_1,...,P_m))`
has logarithmic height `w+o(w)`. If `P_i'=zeta P_i`, then
`|zeta|=1`. Among lifts whose rows remain conjugate-primitive, the
primitive reanchored numerators obey

```text
h_ab'=h_ab                 (a,b>0),
h_a0'=bar(zeta) h_a0,      h_0b'=zeta h_0b.              (14)
```

The first identity follows from the fixed Gram entries in (1).
The other two show why the complete Gram matrix at a *different*
anchor is not fixed before the original common rotation is chosen:
it includes an edge to label `0`. Its common Gaussian factor is
therefore a derived quantity, not another free rotation parameter.
Once the original lift is chosen, (1), (2), and (10) force every
reanchored factor and rational phase. Their cocycles are algebraic
consequences of that lift.

For the actual near-real lift, the anchor-zero bounds give
`|Y_i|=exp(o(w))`; (4) and the pair-distance estimate give
`|Im h_ab|=exp(o(w))` at all other anchors. The row norms are
`exp(2^(m-1)w+o(w))`. The tiny-imaginary-coordinate uniqueness
criterion applies at every anchor, but these are uniqueness statements
for derived anchored data. They do not combine into a second
height inequality or exclude the one exceptional original rotation.

The [primewise checker](check_all_anchor_rotation_radius_cocycle.py)
verifies the pair-core exponent ordering, the exact common least
radius, the correction cocycle, and `Norm(E)|D^2` on all 6,573
endpoint-allocation fixtures with four labels and original prime
exponents two through eight. The identities above are proved for
arbitrary exponents and any fixed number of labels.
