# Elliptic Vieta heights and the original moving frame

An elementary binary-quadratic factorization already gives the same
leading upper bound for the normalized coordinate after one Vieta swap. Canonical
heights give a compatible fixed-curve, fixed-frame comparison. Neither
observation by itself improves the actual
`13w-log E_i` row-size allowance or the `5w-(log E_i)/2` source-anchor
multiplier allowance. The original rational frame has logarithmic height
`12w+o(w)`, and its height-comparison error is not controlled by the
coefficients of the invariant.

All elliptic statements below assume that the anticanonical curve is
smooth. Geometric integrality and exclusion of boundary nodes alone would not
exclude a cuspidal curve of geometric genus zero.

The persistent [exact checker](check_five_row_elliptic_vieta_frame_obstruction.py)
verifies the polynomial norm identity, 9,216 primitive factor pairs,
all labelled cut counts, exact frame contents, and the moving-frame
examples below. Its literal full-cut fixture tests arithmetic identities
and makes no endpoint or small-residue claim.

Every profile expression `k w+o(w)` in this note uses the regime
`eta -> 0` and `beta+sigma=o(w)` from the source hypotheses; the
normalized Vieta conclusion also uses `log C_Q=o(w)`. For fixed
`eta`, the profile error is instead `O(eta w+beta+sigma)`, with an
additional `log C_Q` in bounds where the coefficient height occurs.

## 0. The normalized size bound needs no elliptic theory

Fix three retained directions at `(1,0),(0,1),(1,1)` and denote the other
two normalized projective coordinates by `a,b`. Substituting these
three constant vectors into `Q` gives a bihomogeneous form `q(a,b)` of
bidegree `(2,2)` with coefficient sum at most `C_Q`.

For fixed primitive `b=(u,v)`, regard `q(a,b)` as a nonzero binary
quadratic in `a`. Each of its three coefficients has absolute value at
most `C_Q max(|u|,|v|)^2`. Removing their ordinary common content can
only decrease height. Its two rational projective roots `a,a'` give a
factorization of the primitive quadratic into primitive integer linear
forms, up to sign, by Gauss's lemma.

For two real linear coefficient vectors `l_1,l_2`, their product's
quadratic coefficient vector has Euclidean norm at least
`||l_1||_2 ||l_2||_2/sqrt(2)`. This follows by expanding the squared
norms: for `l_1=(a,b)` and `l_2=(c,d)`, exactly

```text
2||(ac,ad+bc,bd)||_2^2 - ||(a,b)||_2^2 ||(c,d)||_2^2
   = (ac+bd)^2+(ad+bc)^2.
```

Bounding a three-entry
coefficient norm by `sqrt(3)` times its maximum coefficient therefore
proves the explicit projective-height inequality

```text
h(a)+h(a') <= h(primitive quadratic)+(log 6)/2,
h(a') <= 2h(b)-h(a)+log C_Q+(log 6)/2.                (0)
```

The proof includes a root at infinity. It only requires that the row
quadratic is nonzero and has the two rational roots under consideration;
smoothness of the curve is unnecessary.

The actual full profile gives `h(a)=h(b)=4w+o(w)`. Hence
`log C_Q=o(w)` already implies `h(a')<=4w+o(w)`. This one-sided estimate
needs no coefficient-uniform Weierstrass or canonical-height comparison.
It is sufficient for the frame upper bound below, which nevertheless
remains weaker than the transferred-core row-size bound.

## 1. An exact conditional canonical-height estimate

Let `C` be a fixed smooth anticanonical section, choose a rational origin
`O`, and use the degree-one canonical height `q` and pairing from
[the compatibility note](five_row_elliptic_height_compatibility.md).
The restriction of each forgetful map `f_j:C -> P^1` has degree two.
Write its fiber divisor class as

```text
f_j^* O(1) ~ [O]+[T_j].
```

The corresponding deck involution is `iota_j(P)=T_j-P`. A canonical
height for this degree-two class is

```text
F_j(P)=q(P)+q(P-T_j).
```

For any two labels, exact quadratic-height algebra gives

```text
F_j(iota_i P)-F_j(P)
   = -4 <P-T_i/2, T_i-T_j>.                         (1)
```

Division by two here is in the rationalized Mordell-Weil group; no
rational half-point is required. For `j=i` the difference is exactly
zero, as it must be because that forgetful coordinate is unchanged.

Let `C_j` be a valid comparison constant for the chosen actual coordinate:

```text
|h(f_j(P))-F_j(P)| <= C_j.
```

Then

```text
|h(f_j(iota_i P))-h(f_j(P))|
 <= 4 sqrt(q(P-T_i/2) q(T_i-T_j)) + 2 C_j.           (2)
```

This is unconditional with these displayed constants. If the curve,
origin and maps are fixed, its right side is `O_C(sqrt(H(P))+1)`.
For moving `Q`, the assertion that it is `o(w)` requires, for example,
`q(T_i),q(T_j),C_j=o(w)` and `q(P)=O(w)`. The mere inequality
`log C_Q=o(w)` is not a substitute for those comparisons. A uniform
bound of these data by `O(log C_Q+1)` would suffice, but is not proved
or assumed here. The relevant general height properties are recalled in
[Silverman's survey, slides 7--10 and 15--16](https://legacy.slmath.org/attachments/workshops/301/HtSurveyMSRIJan06.pdf).

## 2. The exact frame term

Fix three retained rows with labels `a,b,c`, and use their projective
frame to read the varying row. Its normalized coordinate is one of the
other forgetful coordinates, so (2) applies. Let `M` be the primitive
integer matrix restoring the actual Cartesian row frame. For a primitive
integer vector `v`, put

```text
c_M(v)=gcd_Z((Mv)_1,(Mv)_2)>0,
s_M(v)=||Mv||_2/||v||_2.
```

The actual primitive row is `Mv/c_M(v)`, up to sign, and exactly

```text
log ||prim(Mv)||_2
 = log ||v||_2 + log s_M(v) - log c_M(v).              (3)
```

Consequently the old-to-new difference has the additional term

```text
log(s_M(v')/s_M(v)) - log(c_M(v')/c_M(v)).             (4)
```

Canonical height on `C` controls neither part of (4). In the height
machine, composing `f_j` with `M` leaves its degree-two divisor class
unchanged, but changes the chosen sections and their metrics. The
resulting bounded error depends on `M`, even for one fixed curve.

For the usual projective logarithmic height, the elementary inequalities
for `M` and its adjugate give

```text
|h(Mv)-h(v)| <= h(M)+log 2.
```

Thus a fully charged version of (2) for the actual rows adds
`2h(M)+2log 2`. In particular a fixed-curve, fixed-frame theorem is
valid; a coefficient-only error independent of the moving frame is not.

## 3. The source frame has height twelve

Write `Delta_uv=det(P_u,P_v)` and set

```text
d=gcd_Z(Delta_bc,Delta_ca)>0,
M=(Delta_bc P_a, Delta_ca P_b)/d.                    (5)
```

The matrix is primitive because each column row `P_a,P_b` is primitive.
It sends the three standard projective directions to `P_a,P_b,P_c`.
Its determinant has absolute value

```text
|det M|=|Delta_bc Delta_ca Delta_ab|/d^2.             (6)
```

In the full five-row profile, the common core of the two brackets in
`d` consists of exactly four cut norms (the cuts containing `a,b,c`).
The remaining core parts are coprime. Extra integer gcd divides the
product of the two bracket corrections. Therefore

```text
log d=4w+o(w),
h(M)=12w+o(w),
log|det M|=16w+o(w).                                 (7)
```

Its largest and smallest singular values are consequently
`exp(12w+o(w))` and `exp(4w+o(w))`.

For the old varying row `i`, an exact normalized vector is the primitive
reduction of

```text
v_raw=(Delta_ca Delta_ib, Delta_bc Delta_ai).
```

Its ordinary content has logarithm `12w+o(w)`, and hence
`log||v||=4w+o(w)`. If that content is `s`, then exactly

```text
Mv = [Delta_bc Delta_ca Delta_ab/(d s)] P_i.
```

For each cut, the minimum core exponent in these two matching products
equals the minimum in all three matchings of `a,b,c,i`. If that cut
contains `r` of these four labels, the common exponent is
`max(0,r-2)`. Summing over the 32 cuts gives twelve. After removing
this common core, the two residual core products are coprime; any
extra ordinary content divides the product of their four bracket
corrections. The checker verifies both the cutwise minima and this
correction divisibility on literal Gaussian rows.

Thus the old frame content is `c_M(v)=exp(8w+o(w))`, and its stretch
is `s_M(v)=exp(12w+o(w))`. The difference in (3) is precisely
`4w+o(w)` at the old row.

Applying the elementary upper bound (0), the frame upper bound gives only

```text
log||P_i'|| <= 16w+o(w)-log c_M(v').                 (8)
```

Without a new estimate for the finite content or the stretch, this is
weaker than the existing `13w-log E_i+o(w)` bound. Keeping the new row at
`8w+o(w)` would require control of the *new* expression
`log s_M(v')-log c_M(v')` at `4w+o(w)`. This is exactly the missing
frame term, not a consequence of (1).

## 4. A fixed-section obstruction to a frame-independent error

Take any fixed smooth section, one fixed rational interior point where
the chosen swap is unramified, and primitive old/new row vectors `v,v'`.
Their determinant is nonzero. Choose `U in SL_2(Z)` with `Uv=(1,0)`.
Write `Uv'=(a,b)`; then `b!=0`. The matrices

```text
S_N=[[1,N],[0,1]]
```

fix the old vector and send the new vector to `(a+Nb,b)`. Its primitive
height tends to infinity, since unimodularity preserves primitivity.
The curve, moduli points, canonical heights, ramification status and
invariant coefficients are all unchanged. Thus a row-height comparison
with an error depending only on `Q` is false even for fixed `Q`.

This example is not a full-profile endpoint configuration. It isolates
the necessity of the frame dependence; it does not disprove an estimate
that additionally uses new arithmetic information from that profile.

For a literal matrix check, take `v=(12,1), v'=(22,1)`, use
`U=[[0,1],[-1,12]]`, and postcompose with `V=[[1,0],[2,1]]`. Then

```text
V S_N U v =(1,2),
V S_N U v'=(1-10N,-8-20N),
det(V S_N U v,V S_N U v')=-10.
```

Both vectors are primitive and have odd norm. At `N=1,10,10000` the old
squared norm is always `5`, while the new squared norms are respectively
`865, 53065, 50003000065`. Exact integer checks also verified (5)--(6)
and the displayed old-row content identity on three five-row fixtures,
including a literal full-cut Gaussian tuple. These checks are saved in
[the executable certificate](check_five_row_elliptic_vieta_frame_obstruction.py).

No saving in the source-anchor multiplier or primitive pair-residue
exponent is established here. Canonical height compares global divisor
heights; it also does not ensure that a new Gaussian denominator is
supported on the old anchor's Gaussian primes.
