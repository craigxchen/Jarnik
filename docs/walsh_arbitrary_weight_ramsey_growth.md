# Arbitrary-weight Walsh growth from four-row phases and Ramsey theory

Every actual repeated-Walsh endpoint configuration satisfies

```text
W0/(M log M) -> infinity,
M=o_C(log R/loglog R).                                    (1)
```

There is no affinity, prime-size, or slack-concentration hypothesis.
A quantitative version proved below is

```text
W0 >= c M log M (log_2^* M)^(1/3),
M=O_C(log R/[loglog R (log_2^* R)^(1/3)]).                 (2)
```

The parameter refinement at the end improves the quantitative saving
from `(log_2^* R)^(1/3)` to
`sqrt(log_2^* R)/log_2(log_2^* R)`.

These bounds hold for sufficiently large arguments depending on C;
their thresholds are enormous. They concern the repeated-Walsh model,
not arbitrary circle configurations, and do not prove a uniform count.

## Model and previously proved identities

Let `M=2^t`, with `b>=5` physical copies of every nonzero Walsh label,
one distinct physical-column flip per row, and an arbitrary assignment
`a:F_2^t -> F_2^t minus {0}` respecting copy capacity. Physical split
primes are distinct, and their logarithms have total W0. Write F for
the flipped-prime log mass, D for the common Gaussian-content log norm,
so `F<=W0`, `D>=0`, and `log R^2=W0+D`. Row units are allowed.
Assume an actual source arc of length at most `C sqrt(R)`.

Use the notation and normalization of the
[pair-slack proof](walsh_pair_slack_fixed_family_obstruction.md):

```text
kappa=4log C-log 4-D,              aC=max(0,4log C-log 4),
q_a=W_a-4F_a/M,                  Q=sum q_a=W0-4F/M,
u_d=kappa-qhat(d)>=0             (d!=0),
U=sum_(d!=0)u_d=Q+(M-1)kappa <=W0+(M-1)aC,
|kappa|<=W0/(M-1)+aC.                                    (3)
```

For any row subspace V of size h, nonzero restriction ell, and row
coset P, the exact necessary phase inequality is

```text
Psi >=2(2Fmatch-Fp)+(4/M)(F-hFcell)-4log(h/2),
Psi=sum_(d in V minus {0})(-1)^(ell(d)+1)u_d.              (4)
```

For an aligned coset it also implies

```text
(M-4)[2Fp-4log(h/2)] <=M Psi+4W0-4kappa-4F.              (5)
```

The [residual proof](walsh_full_support_residual_mass.md) checks all
row-coset averages with the original ambient M and radius retained.

## Four-row averaging bounds every slack value

For fixed V,ell, average (4) over all `M/h` row cosets. Since
`E Fp=hF/M` and `E Fmatch=hFcell/M`, the Fcell terms cancel:

```text
Psi >=-(2h-4)F/M-4log(h/2).                               (6)
```

Take `V=span(d,e)` for distinct nonzero d,e and choose
`ell(d)=0`, `ell(e)=1`. With h=4, (6) gives

```text
u_d <=u_e+u_(d+e)+E,       E=4F/M+4log 2.                (7)
```

Sum over all `e!=0,d`. Each positive sum is `U-u_d`, so

```text
M u_d <=2U+(M-2)E,
u_d <=6W0/M+2aC+4log 2=:B.                               (8)
```

Unlike a bounded average alone, (8) controls every direction.

## Fixed-color Ramsey gives the qualitative improvement

Suppose `W0<=K Mlog M` along arbitrarily large M. We may enlarge K
to a fixed integer `K>=1`. Distinctness of the `r=b(M-1)` physical primes gives
`W0>=log((r+1)!) >=(bM/4)log M` for M sufficiently large, hence
`b<=4K`. From (8), for large M,

```text
0<=u_d<=(6K+1)log M.
```

Color every nonzero direction by `floor(8u_d/log M)`, using at most
`48K+9` colors.
The finite vector-space Ramsey theorem gives, for every fixed s, an
s-dimensional linear subspace H whose nonzero vectors all have the
same color when t is sufficiently large. The theorem concerns
one-dimensional linear subspaces; over `F_2` these are precisely
individual nonzero vectors. See Frederickson--Yepremyan,
[Theorem 1.1](https://arxiv.org/html/2308.13489v1).

Choose s first, depending only on K, with `s>=4096`, and put

```text
N=2^s,
h=2^floor(log_2[s/(8(log_2(s+2))^2)]).
```

Choose s large enough that `h>=32K` and
`N^(3/4)>(4K+1)h`. Both requirements hold because h tends to infinity
and `h<=s`. Once H exists, there is `c>=0` such that

```text
c<=u_d<c+(log M)/8       for all d in H minus {0},
c<=(6K+1)log M.                                         (9)
```

There are at most `sqrt(M)` rows whose flipped prime log is below
`(log M)/2`. At most `b(M/N-1)` rows have assigned labels vanishing
on H. Averaging their union over the `M/N` row cosets of H, and
taking M large enough that `N<=sqrt(M)`, produces a coset T with
at most `b+1` bad rows. Replace their restricted labels by arbitrary
nonzero functionals on H. The
[nonlinear aligned-flat theorem](walsh_arbitrary_nonlinear_aligned_flat_extraction.md)
applies to this N-row field without a capacity assumption and retains
at least `N^(3/4)` rows in disjoint aligned h-point flats. Discard
every flat meeting a patched row. Since at most `(b+1)h` rows are
lost, some actual unpatched flat survives. On this flat,

```text
Fp>=(h/2)log M.
```

Its nonzero character signs sum to one. Thus (9) gives
`Psi<=c+(h/8)log M`. Substitute into (5), using (3) and
`W0<=K Mlog M`, and divide by `Mlog M`. With s and h fixed,
letting M grow yields

```text
(7/8)h <=10K+1+o(1),
```

contrary to `h>=32K`. This proves the first assertion in (1).
Notice the order of quantifiers: K, the number of colors, and s
are fixed before M grows. No fixed-rate dimension extraction from
an unbounded number of colors has been assumed.

## Quantitative version with iterated logarithms

Define `log_2^* x` as the number of successive base-two logarithms
needed to reach a value at most one. The same primary source's
[Theorem 1.3](https://arxiv.org/html/2308.13489v1) bounds the required
Ramsey dimension `R_2(s;r)` by a base-two tower of height
`(r-1)(s-1)+1` with top `3s`.

For integer `K>=2^24`, choose `s=K^2`. The displayed formula for h
gives `32K<=h<=s`: it suffices to use
`log_2(K^2+2)<=2log_2 K+1` and
`K>=512(2log_2 K+1)^2`, which holds at `2^24` and continues to hold
by monotonicity. There are `r<=48K+9<=49K` colors. The tower estimate
therefore implies

```text
log_2^* R_2(s;r) <=50K^3                                (10)
```

for these K. Indeed, an upper bound is
`49K^3+2+log_2^*(3K^2)<=49K^3+K+3<=50K^3`;
use `3K^2<=2^K` and the large K cutoff. Put `Lstar=log_2^* M` and choose
`K=floor((Lstar/200)^(1/3))`. For sufficiently large M, K meets
the cutoff, whereas the source dimension `t=log_2 M` has
`log_2^* t=Lstar-1>50K^3`. Thus t exceeds the Ramsey threshold.
It also exceeds `2s`, so `M>=N^2` and the bad-row argument applies.
The remaining fixed-C requirement is `log M>=2aC+4log 2`.

To check uniformity when K now varies, put `l=log M`. Equations
(3), (5), and (9) imply the exact upper bound

```text
h(7/8-4/M)
 <=10K+1+8K/M+4log h/l+4aC/(Ml).                         (11)
```

Here `|kappa|<=2Kl+aC`. For `M>=32` and
`l>=max(aC,4log h)`, the right side is at most `13K`, while the
left side is at least `(3/4)h>=24K`. All these conditions follow
eventually from `M>=2^(2K^2)` and the fixed-C cutoff. Hence
`W0>K Mlog M`, proving the first inequality in (2).

For the radius form let `L=log R^2>=W0`. If `M<sqrt(L)`, the asserted
bound is automatic for large L. Otherwise `log M>=(log L)/2` and
`log_2^* M>=log_2^* L-1` for large L. Substituting in (2) gives
`M=O_C(L/[log L (log_2^* L)^(1/3)])`. Replacing L by `log R^2`
changes only fixed factors and a bounded shift of the iterated log.
The same two-case argument proves the little-o statement in (1).

## A stronger class-weight regularity corollary

Let `delta q_a=q_a-Q/(M-1)` for nonzero a and `delta q_0=0`.
Its nonzero Fourier coefficients equal `U/(M-1)-u_d`; the zero
coefficient vanishes. Since `0<=u_d<=B`, their variance is at most
`B^2/4`: use `u_d^2<=B u_d` and complete the square in the mean.
Parseval followed by Cauchy--Schwarz gives

```text
sum_(a!=0)|delta q_a| <= B(M-1)/(2sqrt M).
```

Adding back `4F_a/M` and recentering at the physical class mean yields

```text
sum_(a!=0)|W_a-W0/(M-1)|
 <=B(M-1)/(2sqrt M)+8F/M = O_C(W0/sqrt M).                (12)
```

The last estimate uses the prime distinctness lower bound. This
strengthens the earlier [density-based regularity](walsh_phase_class_mass_regularity.md).
The common-content bound already follows directly from `U>=0`:
`D<=4log C-log 4+W0/(M-1)`; it is not a new consequence here.

The [bounded checker](check_walsh_arbitrary_weight_ramsey_growth.py)
checks the four-row averages, slack-max algebra, and Parseval identities.
The arbitrary-dimensional Ramsey step uses the cited theorem. Astra
found the main argument; Root, Sol and a separate Luna instance audited
the qualitative proof, and Root and Sol checked the quantitative step.

## Tighter quantitative consequence

The cubic exponent in (10) came from the convenient choice `s=K^2`.
The aligned-flat theorem actually needs only `s=O(K(log K)^2)`.
Keeping that dependence gives the stronger bounds, with
`T_M=log_2^* M` and `T_R=log_2^* R`,

```text
W0 >= c M log M sqrt(T_M)/log_2(T_M),
M = O_C(log R log_2(T_R)/[loglog R sqrt(T_R)]).           (13)
```

These assertions concern sufficiently large arguments, so both
iterated-log parameters exceed one. This is a refinement of the same
family theorem, not an extension to arbitrary circle profiles.

Here are explicit choices checking every varying parameter. Suppose
`W0<=K Mlog M`, where K is an integer at least `2^24`, and set

```text
m=ceil(log_2 K),        s=2^14 K m^2,
h=2^floor(log_2[s/(8(log_2(s+2))^2)]).
```

For `m>=24`, `15+2log_2 m<=2m`, hence
`log_2(s+2)<=15+log_2 K+2log_2 m<=3m`. Therefore

```text
h>=s/[16(log_2(s+2))^2]>= (1024/9)K >32K.               (14)
```

Also `s<=K^2`. At `K=2^24` this is direct. If `m>=25`, it follows
from `2^14 m^2<=2^(m-1)<K`; the first inequality holds at m=25
and persists because `2^(m-1)/m^2` increases for integer m>=3.
Thus the extraction requirements `s>=4096`, `h<=s`, and
`2^(3s/4)>(4K+1)h` all hold. For the last one, use `K<=s`,
`h<=s`, and `2^(3s/4)>5s^2` for `s>=32`.

With at most `r<=49K` colors, the same cited tower theorem gives

```text
log_2^* R_2(s;r)
 <=49Ks+2+log_2^*(3s)
 <=50Ks <2^20 K^2 m^2.                                  (15)
```

Here `log_2^*(3s)<=s` for s>=4, and `s+2<=Ks`.
Now put `T=log_2^* M` and take

```text
K=floor(sqrt(T)/[2^12 log_2 T]).                         (16)
```

For sufficiently large M, K meets the cutoff. Since
`m<=log_2 T`, the last quantity in (15) is at most `T/16`.
But `log_2^*(log_2 M)=T-1>T/16`, so the ambient dimension
`t=log_2 M` exceeds the Ramsey threshold. Also `s<=K^2<=T`,
so t>2s for large T, and hence `M>=2^(2s)`. All the bad-row,
extraction, and fixed-C cutoffs in the earlier proof apply uniformly.
The exact inequality (11) again gives the contradiction
`24K<=13K`. Consequently `W0>K Mlog M`.

The floor in (16) loses at most a factor two once its argument is at
least two. This proves the first bound in (13). The same two-case
radius inversion used after (11) proves the second: if
`M>=sqrt(log R^2)`, then `log M` is comparable to `loglog R`
and `log_2^* M` is at least `log_2^* R-O(1)`; otherwise the
claim is immediate. Bounded shifts change `sqrt(T)/log_2 T`
by at most a fixed factor for large T.
The upper comparison `log M<=loglog R^2` follows from
`M<=W0<=log R^2`, itself immediate for large M from the distinct-prime
lower bound. Thus no upper comparison is being inferred from the
case assumption alone.

The companion checker verifies the new elementary integer cutoffs and
the floor-dependent choices. It does not replace the dimensional
Ramsey theorem or the original phase and extraction proofs.
Root derived this refinement; Sol independently checked the tower
statement, integer cutoffs, varying-K restrictions, and radius inversion.
A fresh Astra instance independently confirmed the same points.
