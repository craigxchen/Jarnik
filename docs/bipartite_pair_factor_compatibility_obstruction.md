# A bipartite near-critical factor profile with exact cycle compatibility

## Scope

This note gives an actual primitive, common-unit Gaussian tuple, for every
`M=4r>=12`, whose **pair norms** obey the endpoint lower bound and whose
near-critical pair graph is exactly `K_(M/2,M/2)`. All Gaussian triangle and
four-cycle identities hold because the factors come from the tuple itself.
The source arguments and the pair imaginary coordinates are uncontrolled.
Thus this is an obstruction to extracting a good triangle, or a fixed power
from a good four-cycle, using the norm window and factor cocycles alone. It
is not an endpoint arc construction or a divisor-strip construction.

The later [half-symmetric phase-height theorem](half_symmetric_star_phase_height_gap.md)
excludes an unbounded endpoint realization of this exact profile,
even when its prime supports and weights vary. That proof uses the
actual short-arc phases omitted from the norm construction here.

Put `lambda=log(4/C^2)>=0`, with `0<C<=2`. In an actual endpoint cluster,
every pair has `log P_ab>W/2+lambda`, where `W=log N`, and the Plotkin-good
window ends at `W/2+W/M`. The construction below meets these **norm**
inequalities strictly, without claiming the angular inequalities that
produce them in an endpoint cluster.

## 1. The exact Gaussian tuple

Write `M=2k` with even `k>=6`, and partition the labels into halves `A,B`
of size `k`. Let `F` be all cuts

```text
S=U union V,       U subset A, V subset B,
|U|=|V|=k/2,       L=|F|=binom(k,k/2)^2.             (1)
```

Choose pairwise distinct odd split primes `p_S` for `S in F` and one more
such prime `p_*`. Fix one Gaussian prime `pi_p` above each `p`. For a real
parameter `T`, take `s_S(T)` to be the least **odd** positive integer with
`s_S log p_S>=T`. Then

```text
T<=s_S log p_S<T+2 log p_S.
```

Let `E=2 sum_(S in F) log p_S`. Choose a fixed odd exponent `s_*` so that

```text
V=s_* log p_*>E+2 lambda.                           (2)
```

Use the cut `A|B` at `p_*`. Define `x_i(S)=1_(i in S)` and
`x_i(*)=1_(i in A)`. The rows are the literal Gaussian integers

```text
z_i=product_(S in F) pi_(p_S)^(s_S x_i(S))
                    bar(pi_(p_S))^(s_S(1-x_i(S)))
     · pi_(p_*)^(s_* x_i(*)) bar(pi_(p_*))^(s_*(1-x_i(*))). (3)
```

All rows have unit `1` and common squared radius
`N=p_*^s_* product_S p_S^s_S`. Every cut is nonconstant, so the tuple has
Gaussian gcd one. For an ordered pair, its actual conjugate-primitive
quotient numerator is

```text
A_ij=product_(x_i>x_j) pi_p^s_p
     product_(x_i<x_j) bar(pi_p)^s_p,
z_i/z_j=A_ij/bar(A_ij),
P_ij=Norm(A_ij),       log P_ij=sum_(cuts separating i,j) s_p log p. (4)
```

The product runs over both `F` and `*`. Since all exponents are odd and at
least one cut separates each pair, every `P_ij` is nonsquare and every
`A_ij` is nonreal. Distinct pairs have distinct prime-support signatures
in `F` when `k>=6`, hence distinct ordinary norms even before the extra
cut. For cross pairs this follows from independent half-subset choices.
For pairs within one half, half-subsets of size at least three distinguish
any two different pair-separation patterns; independence distinguishes
within-half patterns from cross patterns and from the other half.

## 2. The good graph is exactly bipartite

For a uniformly chosen cut in `F`, a cross pair is separated with
probability `1/2`. A pair within one half is separated with probability

```text
alpha=k/[2(k-1)]=1/2+1/(M-2).                     (5)
```

Write `w_S=s_S log p_S=T+e_S`, with `0<=e_S<2 log p_S`, and
`E_T=sum e_S<E`. Thus `W_0=sum_S w_S=L T+E_T` and `W=W_0+V`. Every cross
pair has

```text
log P_cross=(L/2)T+e_cross+V,
log P_cross-W/2=V/2+e_cross-E_T/2
                 >(V-E)/2>lambda.                    (6)
```

Its distance from the near-critical upper edge is

```text
log P_cross-(W/2+W/M)
 =-L T/M + O_(F,p_S,p_*,C)(1)<0                 (7)
```

for all sufficiently large `T`. Every within-half pair has

```text
log P_within=L alpha T+e_within,
log P_within-(W/2+W/M)
 = [2L/(M(M-2))]T+O_(F,p_S,p_*,C)(1)>0.        (8)
```

It also lies above `W/2+lambda` for large `T`. Therefore all pairs
meet the endpoint-derived lower **norm** bound, while precisely the `k^2`
cross pairs meet the near-critical upper norm bound. Since
`k^2=M^2/4>=ceil(M(M-2)/4)`, the Plotkin quota is met by a triangle-free
good graph. Every cut, including `A|B`, is balanced, so the weighted
layer-cake Plotkin sum `sum_(i<j) log P_ij=M^2 W/4` is exactly maximal.
The dominant cuts have weights `T+O(1)`; the extra cut has fixed weight
`V`. This is a near-uniform, nested-prime-power profile, rather than a
freely assigned pair-distance matrix.

## 3. What a good four-cycle actually says

Take `a,c in A` and `b,d in B`. Let

```text
F_1=A_ab A_cd,       F_2=A_ad A_cb.
```

At each prime their signed Gaussian exponent is the same:

```text
(x_a-x_b)+(x_c-x_d)=(x_a-x_d)+(x_c-x_b).        (9)
```

Consequently there is one conjugate-primitive Gaussian integer `G` and
positive ordinary integers `q_1,q_2` with `F_1=q_1 G`, `F_2=q_2 G`.
For the cuts (1), the two `q` supports are disjoint. A cut contributes
`p_S^s_S` to one of them exactly when it separates both `a,c` and `b,d`;
the two possible alignments decide which side. The `*` cut contributes
to neither. Thus `gcd(q_1,q_2)=1` and

```text
q_2 F_1=q_1 F_2,
q_1^2/q_2^2=(P_ab P_cd)/(P_ad P_cb),
q_1 | Im(F_1),             q_2 | Im(F_2).       (10)
```

The last divisibilities are integer divisibilities, retaining the signed
imaginary coordinates and the literal common unit. If all four cross
factors *additionally* had `|Im A_ij|<=H`, then

```text
|Im F_j|<=2H N^(1/4+1/(2M)),
q_j<=2H N^(1/4+1/(2M))                          (11)
```

provided the corresponding imaginary product is nonzero. This is the
direct four-cycle congruence-height test; it has no growing consequence
for `H` on this profile.

Indeed, the fraction of `F` cuts contributing to each `q_j` is

```text
beta=alpha^2/2=M^2/[8(M-2)^2],
log q_j=beta L T+O_(F,p_S)(1)=beta W+O(1).      (12)
```

For `M>=12`, `beta<=0.18<1/4`; it tends to `1/8`. Even if `H` were
bounded, the upper height in (11) has leading exponent `W/4`, above the
coprime cancellation height `beta W`. Divisibility alone does not force
residue growth.

Nor does (9) force a power. Choose a cut with `a,c` both selected in `A`
and `b` selected, `d` unselected in `B`. Such a cut exists for `k>=6`.
The signed exponent of `G` there is `s_S`, which is odd. Hence `G` is not
a Gaussian square, or an `ell`-th power for any even `ell`. To rule out
any other **fixed** `ell>1`, replace the odd rounding above by the least
exponent `s_S=1 mod 2ell` with `s_S log p_S>=T`; its error is less than
`2ell log p_S`. Enlarge `E` and `V` by the same fixed factor in (2).
Then the displayed signed exponent is coprime to `ell`, while (6)--(8)
are unchanged in their leading terms. The four-cycle equality is a
cancellation cocycle, not a new uniform power-divisibility certificate.

## Limitation and verification

The construction realizes actual source factors, exact triangle and
cycle compatibility, nonsquare and distinct pair norms, and the complete
endpoint-derived **norm** window. It does not control Gaussian prime
arguments, source arc length, or the imaginary coordinates of the good
pair factors. Those missing angular conditions may rule out this profile
as an endpoint cluster and may strengthen four-cycle congruences.

Run `python3 docs/check_bipartite_pair_factor_compatibility.py`. Its exact
combinatorial and prime-exponent checks cover the eight- and twelve-row
signature distinction, all twelve-row triangles and good four-cycles,
the coprime `q` supports and an odd reduced exponent. High-precision
decimal logarithms check one finite twelve-row instance at `C=1/2`;
the inequalities for arbitrarily large `T` follow from (6)--(8).
