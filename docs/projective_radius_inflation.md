# Small nonisometric projective changes inflate the least radius

## Status

A small projective change need not preserve the frozen factors in the
sparse-reconstruction theorem. This note treats that more general
possibility. In the extracted regime, every rational projective change
of negligible coefficient height either acts as a rotation or reflection
on the relative directions, or makes the transformed primitive anchor
numerators almost pairwise coprime. In the latter case their least
circle radius increases from the parent scale `exp(W/2)` to
`exp(mW/4+o(W))`, where there are `m` nonanchor rows.

For `m>=3`, the transformed normalized arc length grows exponentially.
Thus these low-height projective transformations cannot furnish the
desired height descent, even if they are allowed to change every frozen
block. Higher-height transformations and nonprojective constructions
are outside this result. Endpoint uniformity remains unproved.

## 1. Anchored data and the reanchored projective map

Use primitive Gaussian anchor numerators

```text
h_i=x_i+i t_i,       1<=i<=m,
gcd_G(h_i,bar(h_i))=1,
x_i>0,       t_i!=0,
h_i/bar(h_i)=z_i/z_0.
```

For `i!=j`, let `G_ij=gcd_G(h_i,h_j)` and write the exact primitive
transition from [all_anchor_primitive_transition.md](all_anchor_primitive_transition.md)
as

```text
h_j bar(h_i)/N(G_ij)=x_ij+i t_ij.
```

In particular,

```text
x_i t_j-x_j t_i=N(G_ij)t_ij,       t_ij!=0.       (1)
```

Take any invertible rational projective matrix, clear its denominators,
and divide the integer entries by their common gcd:

```text
A=[[a,b],[c,d]],       H=max(|a|,|b|,|c|,|d|)>=1.
```

It acts on each real coordinate direction `(x_i,t_i)`. Reanchor at
the image of `(1,0)` by multiplying the complex outputs by `a-i c`.
Define

```text
q=a^2+c^2>0,       r=ab+cd,       s=ad-bc!=0,
v_i=q x_i+(r+i s)t_i.
```

Thus the transformed relative phase is `v_i/bar(v_i)`, with anchor
phase one. Introduce the Gaussian integers

```text
u=q+s-i r,       v=q-s+i r,
2v_i=u h_i+v bar(h_i).                           (2)
```

The symbol `v` without a subscript denotes this fixed coefficient;
`v_i` denotes a transformed numerator.

## 2. Exact destruction of common Gaussian divisors

Let `G'_ij=gcd_G(v_i,v_j)`. Equation (1) gives

```text
t_j v_i-t_i v_j=q N(G_ij)t_ij,
G'_ij divides q N(G_ij)t_ij.                     (3)
```

On the other hand, conjugate primitivity of `h_i` and (2) imply

```text
gcd_G(G'_ij,G_ij) divides v,
gcd_G(G'_ij,bar(G_ij)) divides u.                (4)
```

For example, a divisor common to `G'_ij` and `G_ij` divides
`2v_i-u h_i=v bar(h_i)`, and it is coprime to `bar(h_i)`.
This proves the first assertion; the conjugate version proves the
second. Since `G_ij` and `bar(G_ij)` are coprime, (4) yields

```text
gcd_G(G'_ij,N(G_ij)) divides u v.
```

Combining this with (3), prime valuation by prime valuation, proves
the useful exact divisibility

```text
G'_ij divides q u v t_ij.                       (5)
```

Indeed, after removing `gcd_G(G'_ij,q t_ij)`, the remaining
quotient divides both `G'_ij` and `N(G_ij)`, and hence divides
`uv`. All prime powers are retained in this argument.

If `uv!=0`, the elementary coefficient bounds give

```text
q<=2H^2,       |u|,|v|<=sqrt(20) H^2,
|G'_ij|<=40 H^6 |t_ij|.                         (6)
```

Thus the large old common divisor has disappeared from the bound.
Only the projective matrix and the old primitive residue remain.

For a concrete fixed-height example, the shear
`A=[[1,1],[0,1]]` gives `v_i=(x_i+t_i)+i t_i` and

```text
gcd_G(v_i,v_j) divides (1+2i)t_ij.               (7)
```

Its common factors are controlled by the old residue and the single
Gaussian prime `1+2i`, regardless of the old shared conductor.

## 3. The two exceptional projective symmetries

Since `q>0` and `r,s` are real integers, `uv=0` means exactly

```text
r=0,       s=q or s=-q.
```

In the first case `v_i=q h_i`; in the second `v_i=q bar(h_i)`.
The relative directions are respectively unchanged or conjugated.
These are the rotation and reflection cases after reanchoring.

Equivalently, the two columns of `A` are orthogonal and have the
same squared length. To check this, use
`q(b^2+d^2)=r^2+s^2`. Such transformations preserve the least
primitive radius of the relative configuration. A reflection must
be included here even though it reverses the oriented core factors.

Every other rational projective matrix has `uv!=0`, so (5)--(6)
apply without any frozen-factor preservation hypothesis.

## 4. Primitive reduction and the exact least-radius construction

The raw `v_i` need not have primitive coordinates. Put

```text
c_i=gcd(|Re v_i|,|Im v_i|),       w_i=v_i/c_i.
```

Because `(x_i,t_i)` is primitive and the integer matrix taking it
to `v_i` is `[[q,r],[0,s]]`,

```text
c_i divides |q s|<=4H^4.                        (8)
```

For a direct proof, `c_i` divides both `qs x_i` and `qs t_i`;
an integer Bezout identity for `x_i,t_i` gives (8).

If the primitive coordinates of `w_i` are both odd, divide it by
`1+i`; otherwise leave it unchanged. Denote the result by `d_i`.
Then `gcd_G(d_i,bar(d_i))=1` and

```text
v_i/bar(v_i)=epsilon_i d_i/bar(d_i),
epsilon_i in {1,i},
|v_i|/(4sqrt(2)H^4) <= |d_i| <= |v_i|.          (9)
```

The units in (9) affect no denominator. If `L` is a Gaussian lcm
of all `d_i`, the least possible radius for these relative phases
including the anchor is exactly

```text
R_new=|L|.                                      (10)
```

This follows from the same integrality argument as the reconstruction
theorem: an integral anchor must be divisible by every `bar(d_i)`;
the anchor `bar(L)` works, and the resulting configuration is
primitive. The inserted units do not affect either implication.

Since `d_i|v_i`, equation (6) also bounds all `gcd_G(d_i,d_j)`.
For any list of nonzero Gaussian integers one has

```text
sum_i log|d_i|-sum_(i<j)log|gcd_G(d_i,d_j)|
 <= log|lcm_G(d_1,...,d_m)| <= sum_i log|d_i|.   (11)
```

This elementary inequality follows by sorting the valuations at
each Gaussian prime; it does not require pairwise coprimality.

## 5. Quantitative height and arc consequences

Suppose a common error parameter `E>=0` satisfies

```text
|log|h_i|-W/4|<=E,
log|t_i|<=E,       log|t_ij|<=E,
J=log H.
```

Assume `E+J=o(W)`. The positive real part of `v_i` is dominated
by `q x_i`: its error `r t_i` is exponentially smaller. In
particular, eventually

```text
|h_i|/4 <= |v_i| <= 6H^2 |h_i|.
```

Equations (6), (9), and (11) then give

```text
log R_new=mW/4+O(m^2(E+J+1)).                    (12)
```

All implied constants are absolute. This can be used for fixed or
growing `m`, provided the displayed error is controlled.

Let `Delta_new` be the length of the shortest angular interval
containing the transformed relative phases and the anchor. Its
upper bound follows from
`|Im v_i|=|s t_i|<=2H^2 exp(E)` and the lower modulus bound.
For its lower bound, use any `i`, the nonzero integers `s,t_i`,
and `|v_i|<=6H^2 exp(W/4+E)`. All phases approach the anchor,
so these sine estimates give

```text
log Delta_new=-W/4+O(E+J+1).                    (13)
```

The normalized length on the least circle is exactly
`C_new=Delta_new sqrt(R_new)`. Therefore

```text
log C_new=(m-2)W/8+O(m^2(E+J+1)).               (14)
```

In the centrally extracted families, one can take
`E=O(kW/sqrt(M))`, with `m=k-1`, `k=O(log M)`, and `W>>M`.
If `J=o(w)` for the individual block scale `w~W/2^(k-1)`, then

```text
m^2(E+J+1)=o(W).
```

Thus (12)--(14) become

```text
log R_new=mW/4+o(W),
log C_new=(m-2)W/8+o(W).                        (15)
```

For every fixed `m>=3` the normalized arc length consequently
tends to infinity. For growing `m` it does so as well. The new
least radius is larger than the original parent-circle scale
`exp(W/2)` by a positive exponential factor. Intrinsically
normalizing the original selected tuple can only lower its radius,
so it does not weaken this exclusion of a radius descent.

## 6. A necessary coefficient height for any endpoint-preserving map

The finite constants also exclude much larger matrices than those of
height negligible relative to a block. Assume `m>=3` and `uv!=0`.
Put `J=log H` and define

```text
A_m=(m^2+m+4)/4,
B_m=(3m^2+m+4)/2,
D_m=(m/2)log(16sqrt(2))
      +(binom(m,2)/2)log 40+log 3.
```

If

```text
2J+2E+log 8 <= W/4,                             (16)
```

then `|t_i|/|h_i|<=1/8` and
`|r t_i|/|h_i|<=1/4`. Hence the modulus lower bound in Section 5
holds with its stated constant. Equations (6), (9), and (11) give
the explicit inequality

```text
log R_new >= mW/4 - [m(m+1)/2]E - (3m^2+m)J
             -m log(16sqrt(2))-binom(m,2)log40.
```

The lower angular bound, which needs no dominance assumption, is
`Delta_new>=1/(3H^2 exp(W/4+E))`. Therefore

```text
log C_new >= (m-2)W/8-A_m E-B_m J-D_m.           (17)
```

If the transformed configuration is required to remain in an arc
with a fixed normalized constant `C`, equation (17) forces

```text
J >= [(m-2)W/8-A_m E-D_m-log C]/B_m.             (18)
```

If (16) fails instead, one already has
`J>W/8-E-(log 8)/2`. The unconditional finite lower bound is the
minimum of these two branch bounds. For each fixed `m>=3`,
`E=o(W)`, and fixed `C`, every nonisometric projective map whose
output remains at the endpoint scale must therefore satisfy

```text
log H >= [(m-2)/(4(3m^2+m+4))]W-o(W).           (19)
```

In the growing central reduction, `m^2E=o(W)` and `m^2=o(W)`.
Equations (18) and the failed-dominance alternative give the
stronger scale statement

```text
log H >= (1-o(1)) W/(12m).                      (20)
```

Here `m -> infinity`; for fixed `m` use the precise coefficient
in (19). Since `w~W/2^m`, the necessary height in (20) is much
larger than one block's logarithmic size. This lower bound concerns
matrices whose outputs actually remain endpoint clusters. It is not
an assertion that every matrix above that height succeeds.

## 7. Scope

The theorem permits the projective matrix to vary and to change
all core factors. Its coefficient height, the old residue heights,
and the content losses are all charged. The only low-height
exceptions are the actual rotational and reflectional symmetries,
which preserve the primitive radius.

This closes one possible use of overlapping frozen designs: a
nonisometric projective change of negligible height cannot turn
their reconstruction into a smaller endpoint cluster. It does not
exclude larger coefficient matrices with a separately proved
radius gain, or a transformation that is not rational projective
on the half-angle directions. No uniform count follows from the
present calculation.

## 8. Verification

An independent audit checked the exact divisor cancellation in (5),
the rational coordinate content and ramified factor in (8)--(9),
the least-radius construction, and both sides of the angular estimate.
The finite constants in (17)--(18), the failed-dominance branch,
and the growing-row threshold (20) were also independently audited.
The comparison deliberately uses the parent radius `exp(W/2)`;
the intrinsic selected radius may be smaller because of its empty
cut, which only strengthens the no-descent conclusion.

Exact Gaussian-integer checks covered 14,385 nonisometric matrix/pair
instances, with primitive input coordinates `1<=x<=60`,
`0<|t|<=12`, matrix entries between minus seven and seven, and
random seed 17831. Every case satisfied (5), `c_i|qs`, and the
conjugate-primitivity assertion after removing the possible ramified
factor. These checks supplement the valuation proof; no new Lean
formalization is claimed.
