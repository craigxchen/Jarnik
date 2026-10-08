# Central truncation with a growing number of uniform Gaussian blocks

## Status

This is a quantitative strengthening of
[endpoint_allocation_rounding.md](endpoint_allocation_rounding.md).
An endpoint cluster with cardinality tending to infinity supplies
`k=floor(c log_2 M)` rows for any fixed `0<c<1/2`. Their prime powers
can be truncated centrally to exact endpoint allocations, with Gaussian
integer row factors of one common modulus. Every one of the exponentially
many cut blocks retains its full relative weight.

The primitive anchor chord numerators are divisible by the resulting
core row products. A simultaneous extraction controlling pair moments
makes both the correcting coefficients and the nonzero integer values
negligible in logarithmic height relative to **each individual block**.
This is an exact arithmetic reduction, not a proof of uniformity. In
particular, the core points themselves need not form a short arc.

All logarithms below are natural. Every asymptotic assertion concerns
`M -> infinity` and permits arbitrary initial prime exponents and prime
support.

## 1. Retaining many rows with small interior weight

Remove the original cluster's common Gaussian divisor and write

```text
z_i = epsilon_i product_p pi_p^(a_ip) bar(pi_p)^(e_p-a_ip),
N=R^2>1,  W=log N=sum_p e_p log p,
E_i=sum_p min(a_ip,e_p-a_ip) log p.
```

Here all rational primes are odd and split, `N(pi_p)=p`, and the
normalized points lie on an arc of length at most `C sqrt(R)`.
Fix `0<C<=sqrt(2)`; taking a smaller fixed positive `C` is sufficient
for the uniformity question by subdivision. The audited estimate in
the rounding note is

```text
sum_(i=1..M) E_i <= W sqrt(M)/2.                  (1)
```

Consequently at least `M/2` rows satisfy

```text
E_i <= tau W,       tau=M^(-1/2).                (2)
```

Take their largest Gaussian-unit class, of size `L>=M/8`. All later
rows are selected from this class, so their units are identical. No
further common divisor is removed before the construction below:
all weights continue to refer to this same original `W`.

The normalized original threshold signs

```text
f_i(p,t)=2 indicator(a_ip>=t)-1,   1<=t<=e_p,
weight(p,t)=log p/W
```

are unit vectors with pairwise nonpositive inner products. Their
distance identity is

```text
mu_ij=<f_i,f_j>=1-2d_ij/W,
d_ij=sum_p |a_ip-a_jp| log p.                    (3)
```

## 2. Simultaneously extracting patterns and pair distances

Fix `2<=k<=L`. For an ordered random tuple of `k` distinct rows and
a nonempty subset `T` of its positions, put

```text
mu_T=E_layers product_(j in T) f_j.
```

The obtuse Bessel argument from
[uniform_profile_extraction.md](uniform_profile_extraction.md) gives

```text
E_rows mu_T^2 <= 2/(L-|T|+1).                   (4)
```

There are two useful versions of the selection. Pattern control alone
gives a tuple with relative error

```text
eta_basic=(2^k-1)sqrt(2/(L-k+1)).                (5)
```

For the stronger simultaneous conclusion, set

```text
A=sum_(nonempty T) mu_T^2,
P=sum_(i<j) mu_ij^2,
a=2(2^k-1)/(L-k+1),
b=2 binom(k,2)/(L-1).
```

Equation (4) gives `E A<=a` and `E P<=b`. Thus some tuple has
`A/a+P/b<=2`, and therefore both `A<=2a` and `P<=2b`.
Fourier inversion and Cauchy--Schwarz then give

```text
eta=2(2^k-1)/sqrt(L-k+1),
|2^k P_layer(epsilon)-1|<=eta  for every sign pattern,
|mu_ij|<=zeta=2sqrt(binom(k,2)/(L-1))  for every pair.    (6)
```

In particular the selected original pair distances satisfy

```text
W/2 <= d_ij <= (1+zeta)W/2.                    (7)
```

We use this slightly larger `eta` in the rest of the note. The
central truncation and its block estimates also hold with (5), but
the small integer values in Section 6 use the separate pair bound.

## 3. Endpoint blocks and exact central truncation

Number the selected rows `0,...,k-1`, and put `B=2^(k-1)`.
Round each `a_ip` to the nearest endpoint `b_ip in {0,e_p}`,
resolving a tie arbitrarily. By (2), each inherited pattern changes
in absolute weight by at most `k tau W`. Consequently

```text
eta_end=eta+2^k k tau.                           (8)
```

At each prime orient `gamma_p` opposite to row zero's rounded
allocation. Let `S subset {1,...,k-1}` be the set of rows whose
rounded allocation differs from row zero, and form the Gaussian
integer block

```text
H_S=product_(primes with pattern S) gamma_p^e_p.
```

The empty block is included. The two complementary sign patterns
combine to give

```text
w_S=log N(H_S),
|w_S-W/B| <= eta_end W/B.                        (9)
```

The blocks are pairwise coprime, individually coprime to their
conjugates, and coprime to the conjugates of the other blocks.

For each prime now define

```text
u_ip=min(a_ip,e_p-a_ip),
r_p=max_(selected i) u_ip,
e'_p=e_p-2r_p >= 0.
```

Replace a rounded low endpoint by `b'_ip=0` and a rounded high
endpoint by `b'_ip=e'_p`. Put `c_ip=a_ip-b'_ip`. These are
nonnegative integers satisfying `c_ip<=2r_p`: on a low row
`c_ip=u_ip`, and on a high row `c_ip=2r_p-u_ip`. Thus the following
are actual Gaussian integers:

```text
corez_i=epsilon_i product_p pi_p^b'_ip bar(pi_p)^(e'_p-b'_ip),
B_i=product_p pi_p^c_ip bar(pi_p)^(2r_p-c_ip),
z_i=B_i corez_i.                                (10)
```

All correcting factors have the same modulus. More precisely,

```text
D=product_p p^r_p  is a positive integer,
|B_i|=D,
s=log D=sum_p r_p log p <= sum_(selected i) E_i <= k tau W,
R_core=R/D,
W_core=log(R_core^2)=W-2s.                       (11)
```

This is a common-radius factorization by Gaussian integers. It
requires neither rational row denominators nor the removal of whole
prime supports.

If a tie occurs, `u_ip=e_p/2` for that row, hence `r_p=e_p/2` and
`e'_p=0`. Both choices of rounded endpoint then give `b'_ip=0` and
`c_ip=a_ip`, so (10) is independent of the tie. More generally any
prime with `e'_p=0` disappears from every core point and belongs
entirely to the factors `B_i`. These cases cause no division by zero
or negative Gaussian exponents.

## 4. Every block survives with its full relative weight

Define the core blocks by retaining the original grouping and
orientation:

```text
H'_S=product_(primes with pattern S) gamma_p^e'_p,
H_S/H'_S=product_(primes with pattern S) gamma_p^(2r_p),
ell_S=log N(H_S/H'_S)=2 sum_(primes with pattern S) r_p log p.
```

The removed divisor is a Gaussian square. It may remove a whole
prime power, but its weight is charged only at the actual removed
exponent. In particular,

```text
ell_S>=0,     sum_S ell_S=2s,
w'_S=log N(H'_S)=w_S-ell_S.                     (12)
```

Relative to the original block scale `W/B`, (9)--(12) imply

```text
|w'_S-W/B| <= eta_core_old W/B,
eta_core_old=eta_end+2Bs/W
            <= eta_end+2^k k tau.              (13)
```

There is a sharper estimate at the new intrinsic weight before
empty-block removal, which can bypass endpoint rounding altogether.
The surviving layers at prime `p` are exactly
`r_p<t<=e_p-r_p`. On all these layers every selected row has its
rounded endpoint sign, because a low allocation is at most `r_p`
and a high allocation is at least `e_p-r_p`. Thus the core layers
are literally a subset of the original threshold layers.

Let `q_S` be the original merged cut weight, combining a pattern
and its complement, and let `v_S` be the weight of the removed
original layers with this cut. Then

```text
|q_S-W/B|<=eta W/B,
v_S>=0,    sum_S v_S=2s,
w'_S=q_S-v_S.
```

The removed-layer masses `v_S` need not equal the block-divisor
masses `ell_S`: the outer layers of one prime may carry cuts
different from that prime's rounded cut. If `2s<W`, the sharper
bound is

```text
|w'_S-W_core/B| <= eta_core W_core/B,
eta_core=[eta+2(B-1)s/W]/[1-2s/W].              (14)
```

Indeed

```text
w'_S-W_core/B=(q_S-W/B)-(v_S-2s/B),
|v_S-2s/B|<=2s(1-1/B),
```

where the last bound uses `B>=2` and `0<=v_S<=2s`.
The formula (14) is an upper bound for the actual error, not an
assertion that a pattern achieves it. The same argument also gives
the sharpened old-scale error `eta+2Bs/W`; estimate (13) remains
useful for comparing each individual rounded Gaussian block with
its own retained divisor.

Now take a fixed `0<c<1/2` and

```text
k=floor(c log_2 M),      M -> infinity.
```

For sufficiently large `M`, `2<=k<=L`, and (6), (8), (11), and
(14) give

```text
eta_end=O(k 2^k/sqrt(M))=o(1),
eta_core=O(k 2^k/sqrt(M))=o(1),
W_core/W=1-O(k/sqrt(M)),
max_S ell_S/(W/B) <= 2Bk/sqrt(M)=o(1).          (15)
```

Hence every one of the `B` core blocks is nonunit once these errors
are less than one, and every block retains a `1-o(1)` fraction of
its original logarithmic norm. This conclusion includes the empty
block. If `C<sqrt(2)`, the strict pair-distance bound also gives
`W>=4(M-1)log(sqrt(2)/C)`. Thus each block's logarithmic norm tends
to infinity in the uniformity reduction, where such a fixed `C`
can always be chosen.

## 5. Exact ratios and the inherited angular scale

For `1<=i<=k-1`, put

```text
A'_i=product_(S contains i) H'_S.
```

Because the selected rows have one common unit,

```text
corez_i/corez_0=A'_i/bar(A'_i).                 (16)
```

For completeness, the correction ratio in (10) can itself be
represented without denominators. Define

```text
C_i=product_p pi_p^((c_ip-r_p)_+)
                  bar(pi_p)^((r_p-c_ip)_+).
```

Then, exactly,

```text
B_i/D=C_i/bar(C_i),
log N(C_i)=sum_p |c_ip-r_p| log p=s-E_i<=s,
V'_i=C_i bar(C_0),     |V'_i|<=D,
z_i/z_0=(V'_i A'_i)/bar(V'_i A'_i).             (17)
```

The first equality uses `p=pi_p bar(pi_p)` with no Gaussian unit
ambiguity. In particular, it is `B_i/B_0` that must be represented
as a conjugate quotient; using `B_i bar(B_0)` directly as the
numerator would square the desired ratio.

The original points supply angular diameter at most

```text
Delta <= C R^(-1/2)=C exp(-W/4).                (18)
```

This is the scale used below. Small `s/W` does not make the phases
of `B_i` close. Equations (10) and (17) therefore do not imply that
the `corez_i` lie on a short arc, either at `R_core` or at `R`.

## 6. Primitive chord numerators give stronger small-value equations

The original anchor ratios have the Gaussian integer numerators

```text
h_i=product_p pi_p^((a_ip-a_0p)_+)
                  bar(pi_p)^((a_0p-a_ip)_+),
z_i/z_0=h_i/bar(h_i),
gcd(h_i,bar(h_i))=1,
log N(h_i)=d_0i.                               (19)
```

The core product `A'_i` divides `h_i` in the Gaussian integers.
To check this prime by prime, if rows `0,i` have the same rounded
endpoint then the prime does not occur in `A'_i`. If their rounded
endpoints differ, its exponent in `h_i` has the prescribed
orientation and equals

```text
|a_ip-a_0p|=e_p-u_ip-u_0p >= e_p-2r_p=e'_p.
```

This remains correct when a tie gives `e'_p=0`. Consequently

```text
h_i=K_i A'_i,          K_i in Z[i] \ {0}.        (20)
```

At a prime with different rounded endpoints, the exponent left in
`K_i` is `2r_p-u_ip-u_0p`; at a prime with the same rounded endpoint
it is `|u_ip-u_0p|`. Both are between `0` and `2r_p`. Therefore

```text
log|K_i| <= s <= kW/sqrt(M).                    (21)
```

Combining (7), (19), and (21) also gives the more direct row-size
estimate `W/4-s<=log|A'_i|<=(1+zeta)W/4`.

Moreover `K_i A'_i` is coprime to its conjugate. In particular
`K_i` is coprime to `bar(K_i)` and to `bar(A'_i)`, although `K_i`
and `A'_i` may share primes in the same orientation.

For any nonzero complex `h` with `h/bar(h)=exp(i theta)`,
`|Im h|=|h| |sin(theta/2)|`, independent of the choice of `theta`
modulo `2pi`. Applying (18)--(19) and using distinctness gives

```text
0 < |Im h_i|
  <= (C/2) exp[(d_0i-W/2)/2]
  <= (C/2) exp(zeta W/4).                       (22)
```

Thus the corrected row equations can be taken to be the actual
primitive anchor residue equations, rather than potentially large
integer multiples of them. Write

```text
K_i=p_i+i q_i,   A'_i=X_i+i Y_i,
q_i X_i+p_i Y_i=t_i,    t_i=Im h_i in Z \ {0}.  (23)
```

Their exact bounds, with our fixed `C<=sqrt(2)<2`, are

```text
log max(1,|p_i|,|q_i|) <= s,
log max(1,|t_i|) <= zeta W/4.                   (24)
```

Since `zeta=O(k/sqrt(M))`, (14)--(15) and (24) yield

```text
max_(i,S) [log max(1,|p_i|,|q_i|,|t_i|)]/log N(H'_S)
 = O(k 2^k/sqrt(M)) = o(1)                     (25)
```

for every fixed `c<1/2`. The elementary coefficient factors `B_i`
also have logarithmic modulus negligible relative to every block.

The separate pair-moment extraction matters here. Pattern control
alone yields an error of size `eta_end W` in the row-product size.
Dividing that error by `W/B` introduces another factor `B`; it does
not give (25) throughout `c<1/2`. Equations (7), (19), and (22)
avoid this extra factor.

## 7. What the reduction does and does not supply

The result supplies exponentially many independent, almost equal
Gaussian blocks, exact Gaussian integer row factors of one common
modulus, and primitive nonzero integer linear-form equations whose
coefficients and values are small relative to each block. All
statements hold with arbitrary growing prime powers. A prime can
still divide both a core block and a correcting coefficient; only
its removed valuation mass has been charged. There is no claim
that the correcting coefficients have disjoint prime support.

Deleting the empty block, if useful later, leaves (16), (19), and
(23) unchanged, since none of the `A'_i` contains it. The remaining
core log norm is `W_core-w'_empty`, while (18) still refers to the
original `W`. This optional deletion cannot justify replacing the
inherited angular estimate by an assertion about a new core arc.

An arithmetic incompatibility theorem for these simultaneous
equations is still missing. In particular, neither (25) nor
conjugate coprimality by itself rules out the system. This note
strengthens the reduction available for such a theorem; it does
not establish the requested uniform endpoint bound.

## Verification

An independent proof audit checked the two-objective selection,
retained-layer normalization, Gaussian factorizations, and primitive
small-value estimates. Exact primewise checks covered all 138,464
allocation tuples with two through five selected rows and exponents
one through eight. See Section 8 of
[the independent audit](uniformity_audit_parallel.md).
These finite checks supplement the general proof; this result has
not been formalized in Lean.
