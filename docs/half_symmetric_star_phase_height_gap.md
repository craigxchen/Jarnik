# An isolated star phase excludes endpoint half-symmetric profiles

## Scope and result

This note adds the **actual endpoint arc** hypothesis to the exact Gaussian
factor setting. A half-symmetric cut profile has one phase direction that
can be isolated by multiplying all rows in each half. The resulting
height bounds contradict endpoint geometry for large radius. The theorem
allows arbitrary moving split-prime supports, unequal weights, and nested
prime powers; no near-critical good-pair selection is needed.

Let `M=2k` with `k>=5`, partition the labels into halves `A,B` of size
`k`, and let `z_i` be distinct Gaussian integers of common modulus `R`,
all with one literal common unit. Assume they lie on an arc of length at
most `C sqrt(R)`, with `0<C<=2`. Write `N=R^2` and `W=log N`. At every odd
split prime `p`, fix one Gaussian prime `pi_p` above `p`. After dividing
the complete common Gaussian gcd as described below, relabel the
residual allocation exponents as `0<=a_i(p)<=e_p`. Their threshold layers are
`S_(p,l)={i:a_i(p)>=l}`. Suppose every residual layer has one of the
forms

```text
regular: |S intersect A|=|S intersect B|=t,  1<=t<=k-1;
star:    S=A or S=B.                              (1)
```

Remove the complete common Gaussian gcd `d` from the allocation columns
first, while retaining it in the literal source factorization and in the
original radius. At a split prime with raw exponents `tilde a_i` and
total exponent `tilde e`, put `m=min_i tilde a_i`, `M_p=max_i tilde a_i`.
The removed factor is `pi_p^m bar(pi_p)^(tilde e-M_p)`; relabel the
residual column by `a_i=tilde a_i-m`, `e_p=M_p-m`. Each residual
split-prime column then has
`min_i a_i(p)=0` and `max_i a_i(p)=e_p`, so all its threshold layers are
nonconstant. Inert and ramified factors, along with the removed split
layers, belong to `d`. Let `W_0` be the total `log p` weight of regular
residual layers, `V` the total weight of star layers, and

```text
D=W-W_0-V>=0.                                   (2)
```

Here `D=log Norm(d)=2 log |d|` is exactly the logarithmic norm of the
complete common Gaussian multiplier. Constant split layers were absorbed
into `d`; inert and ramified factors are common at fixed circle norm and
were absorbed there as well.
In a primitive Gaussian tuple, `D=0`.

**Half-symmetric phase-height theorem.** For fixed `k>=5` and `C<=2`,
the radius `R` is bounded under (1). For `k=4`, no tuple satisfies the
endpoint arc condition at any radius when `C<=2`; the elementary low-row
cases are stated in Section 5. The constants for `k>=5` arise from the
fixed-target Roth theorem and may be ineffective. This is a restriction
on a specified profile, rather than a general point-count bound.

## 1. The star factor and exact phase isolation

The common-unit rows can be written literally as

```text
z_i=epsilon d product_p pi_p^a_i(p) bar(pi_p)^(e_p-a_i(p)),
epsilon in {1,-1,i,-i}.                           (3)
```

All units and the factor `d` cancel from the ratio of the two half-row
products. Primewise,

```text
sum_(i in A) a_i(p)-sum_(j in B) a_j(p)
 =sum_l (|S_(p,l) intersect A|-|S_(p,l) intersect B|). (4)
```

Every regular layer contributes zero, and every star layer
contributes `+k` or `-k`. Because threshold cuts at one prime are nested,
that prime cannot have both an `A` and a `B` star layer, or both a star
and a nonconstant equal-count regular layer. Define

```text
Gamma=product_(p with A-star layers) pi_p^(number of such layers)
      product_(p with B-star layers) bar(pi_p)^(number of such layers).
```

It is conjugate-primitive, `Norm(Gamma)=exp(V)`, and (4) gives the exact
identity

```text
(product_(i in A) z_i)/(product_(j in B) z_j)
              =(Gamma/bar(Gamma))^k.              (5)
```

This retains the actual source factors and the literal common unit;
there is no arbitrary unit correction in (5).

Lift the source arguments `theta_i` continuously along their containing
arc, and let `Delta` be its angular width. Then `Delta<=C/sqrt(R)` and

```text
|sum_(i in A) theta_i-sum_(j in B) theta_j|<=k Delta.
```

If `phi=arg(Gamma)`, equality (5) says that the left-hand difference is
`2k phi` modulo `2pi`. Hence for some integer `m`,

```text
|phi-m pi/k|<=Delta/2.                          (6)
```

The integer `m` records the branch exactly. No small-angle assumption
is needed to pass from the lifted row arguments to (6).

## 2. All endpoint pair norms cap the star height

For every pair, let `A_ij` be its actual conjugate-primitive quotient
numerator and `P_ij=Norm(A_ij)`. The literal common-unit chord formula
is `|z_i-z_j|=2 t_ij sqrt(N/P_ij)`, where `t_ij` is a positive integer.
The chord is strictly shorter than the positive containing arc gap, so

```text
log P_ij>W/2+lambda,      lambda=log(4/C^2)>=0. (7)
```

At a regular layer having `t` rows in each half, the fraction of cross
pairs separated is `2t(k-t)/k^2<=1/2`. The fraction of pairs within `A`
separated is

```text
t(k-t)/binom(k,2)<=alpha:=k/[2(k-1)].           (8)
```

Star layers separate every cross pair and no within-half pair; the common
factor separates neither. Averaging (7) across all `k^2` cross pairs,
then across all `binom(k,2)` pairs within `A`, yields respectively

```text
W_0/2+V>W/2+lambda,
alpha W_0>W/2+lambda.                         (9)
```

Since `W=W_0+V+D`, these give the strict two-sided star budget

```text
V>D+2lambda,        D+2lambda>=0,
V+D<W/k-[2(k-1)/k]lambda.                    (10)
```

The first inequality makes `Gamma` a nonunit even when `C=2`; the
second gives `V<W/k`. Both hold for arbitrary weights and nested layers
under (1). A common multiplier (`D>0`) strengthens both restrictions.

## 3. The support-uniform phase-height contradiction

The finite target set `{m pi/k}` in (6) consists of rational multiples
of `pi`. The fixed-target phase lemma proved in
[fixed_hadamard_roth_phase_obstruction.md](fixed_hadamard_roth_phase_obstruction.md)
applies to the conjugate-primitive nonunit `Gamma`. For any `epsilon>0`,
it supplies a constant `c_(k,epsilon)>0`, independent of the Gaussian
primes and their exponents, with

```text
dist(arg(Gamma),(pi/k)Z)
    >=c_(k,epsilon) |Gamma|^(-2-epsilon)
    =c_(k,epsilon) exp[-(1+epsilon/2)V].       (11)
```

Combining (6), `Delta<=C/sqrt(R)`, and (11) gives

```text
V>=W/(4+2epsilon)-K_(k,C,epsilon).             (12)
```

Choose `epsilon>0` with `1/(4+2epsilon)>1/k`, possible for every
`k>=5`. Equations (10) and (12) then bound `W`, hence `R`. This uses one
isolated algebraic target ray; no fixed support or independent-block
assumption enters the comparison. For `k=4`, all targets `m pi/4` have
rational slopes (`0`, `infinity`, `+1`, or `-1`). The stronger rational
part of the same phase lemma gives `V>=W/2-O_C(1)`, whereas (10) gives
`V<W/4`, again bounding `R`.

The constants in the Roth step depend on `k`. Consequently this result
does not itself give a bound uniform in an increasing point count, nor
does it show that an arbitrary endpoint tuple has a half-symmetric
partition. The exact profile hypothesis is the source of isolation.

## 4. An effective twelve-row gap

For `M=12`, `k=6`, every target in (6) is an axis or has slope
`+/-sqrt(3)` or `+/-1/sqrt(3)`. Let `Gamma=x+iy`, `P=x^2+y^2`.
At an axis target, the nonzero perpendicular integer coordinate gives
`|sin(phi-beta)|>=1/sqrt(P)`. At a nonaxis target, the perpendicular
coordinate has the form `(sqrt(3)y +/- x)/2` or
`(y +/- sqrt(3)x)/2`. Its conjugate factor has modulus at most
`2 sqrt(P)` by Cauchy--Schwarz, while the product is a nonzero integer
`3y^2-x^2` or `y^2-3x^2`. Therefore, for the nearest target `beta`,

```text
dist(phi,(pi/6)Z)>=|sin(phi-beta)|>=1/(4P).     (13)
```

Exact axis equality would make `Gamma` a Gaussian unit, and exact
quadratic equality has no nonzero integer solution. Combining (13) with
(6) gives `Delta>=1/(2P)` and hence a geometric lower bound on the
containing arc length `L=R Delta`:

```text
L>=R/(2P)>(1/2)R^(2/3)                      (14)
```

because (10) gives `P=exp(V)<R^(1/3)` when `C<=2`. Thus an endpoint arc
`L<=C sqrt(R)` requires the simple explicit bound

```text
R<(2C)^6.                                       (15)
```

Retaining the `lambda` term and the common-factor cost `D` in (10)
gives `P<exp(-D) R^(1/3)(C^2/4)^(5/3)`, strengthening (14) to

```text
L>(1/2) exp(D) R^(2/3)(4/C^2)^(5/3),
R<exp(-6D) C^26/2^14 <= C^26/2^14.             (16)
```

In particular, at `C<=1` the right-hand side is below one, so no
nonzero twelve-row integer-circle tuple satisfies (1) and the endpoint
arc condition. Equation (14) is the geometric content under the
theorem's endpoint-derived pair lower bounds: a twelve-row containing
arc in this profile must grow at least as `R^(2/3)`, above the endpoint
exponent `1/2`.

## 5. Effective low-row consequences

The same isolation and pair averages give sharper conclusions for small
halves. These statements still require the literal common unit and the
complete cut condition (1).

For `k=4` (`M=8`), the target grid is `m pi/4`: every ray is an axis or
a diagonal. A conjugate-primitive nonunit Gaussian integer has a
nonzero perpendicular integer coordinate on each such ray, so

```text
dist(arg(Gamma),(pi/4)Z)>=1/[sqrt(2) sqrt(P)],
Delta>=sqrt(2)/sqrt(P).                         (17)
```

Equation (10) says `P<exp(-D) R^(1/2)(C^2/4)^(3/2)`. Comparing (17)
with `Delta<=C/sqrt(R)` gives the effective bound

```text
R<exp(-2D) C^10/256<=4                 when C<=2. (18)
```

But `Gamma` contains an odd split prime, so `P>=5`, whereas the same
upper budget would require `P<sqrt(R)<2`. Thus no eight-row tuple of
this profile fits an endpoint arc with `C<=2`.

For `k=3` (`M=6`), the sharp regular-layer within-half separation
fraction is `2/3`, rather than the general upper bound `3/4` in (8).
Repeating the within-half average gives
`V+D<W/4-(3/2)lambda`. The target grid `m pi/3` consists of axes and
quadratic rays, and the integer-norm argument from Section 4 gives
`Delta>=1/(2P)`. Consequently

```text
L>4 exp(D) C^(-3) sqrt(R),
L<=C sqrt(R)  implies  C^4>4 exp(D).              (19)
```

No six-row tuple of this profile fits an endpoint arc with
`C<=sqrt(2)`.

For `k=2` (`M=4`), the only regular cut has one row in each half, so
the within-half separation fraction is one. Equation (10) specialized
to this exact fraction gives `V+D<W/2-lambda`. The targets `m pi/2`
are axes; their nonzero integer perpendicular coordinate yields
`Delta>=2/sqrt(P)`. Together these imply

```text
Delta>4 exp(D/2)/(C sqrt(R)),
Delta<=C/sqrt(R)  implies  C^2>4 exp(D/2).       (20)
```

Hence no four-row tuple of this exact profile fits `C<=2`. Actual
four-point Pell clusters do not contradict this: their allocation cuts
or unit classes need not satisfy (1).

## 6. Where a general-profile extension breaks

Changing prime weights, moving supports, or nesting layers while every
nonstar layer retains equal A/B counts leaves (5), (8)--(10), and the
proof intact. Unequal counts introduce an additional factor into (5):

```text
(product_A z_i)/(product_B z_i)
 =(Gamma/bar(Gamma))^k · (H/bar(H)),             (21)
```

where the signed exponents of `H` are the A/B count imbalances of the
exceptional layers. An `o(W)` exceptional norm weight does not imply
that `arg(H)` is `o(Delta)`: Gaussian phases can move while their norm
budget remains small. Thus (21) does not place `Gamma` near the fixed
root grid in (6). Applying a scalar phase bound only to the **combined**
quotient, if its reduced numerator is nonunit, charges a lower bound to
the raw height budget `kV+log Norm H` (or less after conjugate
cancellation). That bound does not beat `V<W/k`; a unit reduction
supplies no phase lower bound. A perturbative theorem would need a
support-uniform control of the exceptional factor's argument or a new
simultaneous height estimate after its cancellation.

The separate [norm-only bipartite construction](bipartite_pair_factor_compatibility_obstruction.md)
has regular cuts with exactly `k/2` rows in each half and one star cut.
It satisfies all Gaussian triangle and four-cycle identities and the
endpoint-derived pair-norm window, but its prime arguments were
uncontrolled. The phase-height theorem here shows why that exact
construction cannot occupy an endpoint arc at unbounded radius.

## Verification

Run `python3 docs/check_half_symmetric_star_phase_height.py`. The
standard-library checker verifies every equal-half cut incidence through
`k=6`, including the sharp `k=3` fraction. It also builds a literal
ten-row Gaussian tuple with varying regular cut sizes, a common unit and
common Gaussian multiplier, and checks the star product identity and all
forty-five pair factor/chord identities exactly. The source tuple in this
finite check is not claimed to satisfy an endpoint arc. The Roth and
integer-ray lower bounds remain the displayed proof inputs.
