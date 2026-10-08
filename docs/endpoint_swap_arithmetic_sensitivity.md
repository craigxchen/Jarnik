# A two-unit change can restore endpoint scale in an endpoint-swapping map

For the established four-point Pell family, two nearby parameters in the
same nonconformal endpoint-swapping involution have different radius
growth. The choice `tau=x^2` makes normalized image width tend to infinity.
The choice `tau=x^2+2` has an exact primitive radius and normalized width
bounded by fifteen. The input width tends to zero in both cases.

Thus the second map preserves endpoint scale with a fixed constant. Its
constant being greater than two is not a failure of endpoint scale, and
its increasing radius is not by itself a failure of the conditional
projective route when radii are allowed in either order. This is an exact
four-point family calculation, not an existence theorem for arbitrary
large configurations and not a new general counting bound.

## 1. The source and the involutions

Use the [Pell cotangent family](integer_cotangent_lcm_height_target.md):

```text
U+V sqrt(5)=(9+4 sqrt(5))^n,       n=1 mod 10, n>=1,
U^2-5V^2=1,
a=1+10V^2+2UV,   b=1+10V^2-2UV,
c=13+130V^2+38UV,
x=4UV=a-b,                    ab=x^2+1.
```

The four rational circle phases have half-angle representatives

```text
H_0=(1,0),   H_1=(x,1),
H_2=(u,8),   u=30UV+10V^2+1,
H_3=(w,3),   w=20V^2+16UV+2.
```

Here a nonzero pair `(r,s)` means the phase `(r+is)/(r-is)`.
All the half-angle ratios lie in `[0,1/x]`, with both endpoints present.
The exact primitive squared radius and angular span are

```text
N=5abc,                    Delta=2 arctan(1/x).
```

Since `N` has order `V^6` and `x` has order `V^2`, the normalized source
width `N^(1/4) Delta` has order `V^(-1/2)` and tends to zero.

For `tau>0`, put

```text
M_tau = [[x,tau],[1,-x]].
```

It satisfies `M_tau^2=(x^2+tau)I` and swaps `H_0,H_1` projectively.
On the half-angle interval its action is

```text
r -> (1-xr)/(x+tau r).
```

The denominator is positive, its derivative is strictly negative, and
it swaps `0` and `1/x`. Therefore every such map preserves the exact
angular span `Delta`. The matrices below are nonconformal: their two
column norms differ when `tau!=1`.

## 2. The choice `tau=x^2` destroys endpoint scale

The images of `H_0,H_1,H_2` are represented by

```text
B_0=x+i,       M H_1=(2x^2,0),
B=x(u+8x)+i(u-8x)=B_0 b+16x^2.
```

The Pell equation gives

```text
b=u-8x=2U(U-V)-1,
b=1 (mod V),       b=-1 (mod U).
```

Both `b` and `U` are odd, and `V` is even, so `gcd(b,2x)=1`.
It follows that the real and imaginary coordinates of `B` are coprime:
any common divisor would divide `b` and `16x^2`. Its real coordinate is
even and its imaginary coordinate odd. Hence `B` is coprime to its
Gaussian conjugate. The same is true of `B_0`.

Furthermore, `gcd_G(B_0,x)=1` and `B-B_0 b=16x^2`, so
`gcd_G(B_0,B)` divides sixteen. The norm of `B_0` is odd; this Gaussian
gcd is a unit. Anchor at the real image of `H_1`. Both reduced fractions
`B_0/bar(B_0)` and `B/bar(B)` now occur in the target, and their
Gaussian denominators are coprime. Every integral realization must
clear both. Its least squared radius therefore satisfies

```text
N_bad >= Norm(B_0) Norm(B)
      >= x^4 (u+8x)^2.
```

This uses Gaussian denominator coprimality, not a coprimality assertion
about their ordinary integer norms. As `2 arctan(1/x)>=1/x` for `x>=1`,
the normalized image width obeys

```text
C_bad=N_bad^(1/4) Delta >= sqrt(u+8x) > sqrt(10) V.
```

Thus this parameter fails to preserve endpoint scale on actual sources
whose normalized widths tend to zero. It has no consequence against a
different parameter or a different projective map.

## 3. The choice `tau=x^2+2` restores endpoint scale

Now `x^2+tau=2ab`. The two interior images factor exactly:

```text
M H_2 = b (B_0+16a),
M H_3 = 2a (B_0+3b).
```

After removing these ordinary coordinate contents, the target phases
are represented by

```text
1,       B_0=x+i,
B_1=(17a-b)+i,       B_2=(a+2b)+i.
```

Here `a,b` are odd and coprime. Indeed any common divisor of `a,b`
would divide `x=a-b` and `ab=x^2+1`. Thus `B_1` is already
conjugate-primitive, while `B_2` has exactly one ramified factor:
put `C_2=B_2/(1+i)`. Its primitive phase has a harmless Gaussian unit
times `C_2/bar(C_2)`; that unit does not alter its reduced denominator.
Set

```text
d=288a-31b,       f=(7a+3b)/2.
```

Direct expansion using `a^2-3ab+b^2=-1` gives

```text
Norm(B_1)=a d,       Norm(C_2)=b f.
```

The Gaussian denominator of `B_0` is already covered by those of
`B_1,C_2`. To see this with full multiplicities, split `B_0` into its
coprime Gaussian factors of norms `a` and `b`; this is possible because
`Norm(B_0)=ab` and `gcd(a,b)=1`. The first divides `a` and hence
`B_1=B_0+16a`; the second divides `b` and hence `B_2=B_0+3b`.
The latter remains a divisor after division by `1+i`, since its norm
is odd. Consequently `B_0` divides `B_1 C_2` up to a Gaussian unit.

There is no Gaussian gcd left between `B_1` and `C_2`. Any common
Gaussian prime would be odd and would divide

```text
B_1-B_2=16a-3b=c,       Norm(C_2)=b f.
```

Let `p` be its underlying rational prime. Since
`gcd(c,b)=gcd(16a,b)=1`, it must divide `f`. The identity
`c+2f=23a` then forces `p=23` or `p|a`; in the second case
`c=16a-3b` and `gcd(a,b)=1` force `p=3`. Both three and twenty-three
are inert in the Gaussian integers, so neither can divide `B_1`, whose
imaginary coordinate is one. This proves coprimality, including all
prime powers. The ramified prime is absent because `Norm(B_1)` is odd.

Anchoring again at the real image of `H_1`, the exact Gaussian
denominator lcm is therefore `bar(B_1) bar(C_2)` up to a unit. In
particular the **least primitive** squared radius is exactly

```text
N_good=a b d f.                                      (1)
```

This has order `V^8`; the source squared radius has order `V^6`.
Thus the repair does increase the least radius, but by an amount
compatible with the source's smaller angular width.

## 4. The repaired constant is bounded, with an exact limit

Put

```text
h=1+10V^2,       s=2UV,       h^2-5s^2=1,
x=2s,    a=h+s,    b=h-s,
d=257h+319s,      f=5h+2s.
```

Using `Delta<2/x=1/s` and (1), with `r=h/s`, gives

```text
C_good^4 < (r^2-1)(257r+319)(5r+2).
```

Since `s>=72`, `r^2=5+1/s^2<81/16`. All displayed factors are
positive and increasing in `r`, so

```text
C_good^4 < 65*3589*53/256 < 15^4.
```

Moreover `r -> sqrt(5)` and `s Delta -> 1`, giving the exact limit

```text
C_good -> [4(257 sqrt(5)+319)(5 sqrt(5)+2)]^(1/4)
        = 14.732989546... .
```

The parameter shift therefore changes an unbounded normalized width
into a bounded one. The strict test `N_target<=N` is not met, and
the image constant is not at most two. Neither restriction is part
of the general statement that the image remains at exponent one half
with some fixed constant.

## 5. Every fixed positive rational alignment preserves scale on this family

The repair is not limited to the single parameter above. Fix positive
integers `p,q`, put `T=x^2+1=ab`, and use

```text
M_pq = [[qx,pT-qx^2],[q,-qx]].
```

This is the rational parameter `tau=(p/q)T-x^2` with denominators
cleared. Its determinant is `-pqT`, and it swaps the two endpoints.
Even when its upper-right entry is negative, the denominator of its
fractional-linear action is positive throughout `[0,1/x]`: its values
at the two endpoints are `qx` and `pT/x`. Its derivative is negative.
Thus it preserves the exact angular span for every `p,q>0`. Its
column inner product is `qx(p-q)T`, so it is nonconformal when `p!=q`.

The four images have phase representatives

```text
1, B_0,
F_1=q B_0+8pa,       F_2=2q B_0+3pb,
```

because `M_pq H_2=b F_1` and `M_pq H_3=a F_2` exactly. Moreover

```text
F_1 F_2 = 24p^2 ab  (mod B_0).
```

Since `B_0` divides `ab`, it divides `F_1 F_2`. Consequently the
single Gaussian multiplier `bar(F_1)bar(F_2)` clears all four phases,
and the least target squared radius obeys

```text
N_pq <= Norm(F_1) Norm(F_2).                         (2)
```

This upper bound uses unreduced integral denominators; Gaussian common
content and ramification can only decrease the resulting primitive
radius. It requires no coprimality classification for varying `p,q`.

For an explicit fixed bound, set `H=max(p,q)`. The earlier inequality
`h/s<9/4` gives `a/x<13/8` and `b/x<5/8`. Therefore

```text
Norm(F_1) < 197 H^2 x^2,
Norm(F_2) < (1217/64) H^2 x^2.
```

Combining (2) with `Delta<2/x` proves

```text
C_pq^4 < (239749/4) H^4 < 16^4 H^4,
so C_pq < 16 max(p,q).
```

This proves endpoint-scale preservation for every fixed positive
rational alignment `p/q`, with a nonconformal map unless `p=q`.
The exact radius formula and sharper constant fifteen remain specific
to `p/q=2`. The distinction from the bad parameter is arithmetic:
`tau=x^2` corresponds to the moving ratio `2x^2/(x^2+1)`, whose
reduced numerator and denominator are unbounded. A bound depending
on that denominator would not be uniform along the Pell family.

## Verification and remaining question

The [exact checker](check_endpoint_swap_arithmetic_sensitivity.py)
constructs all four primitive source and target points at four Pell
indices. It compares full Gaussian denominator clearing and final gcd
removal with an independent ordinary all-edge denominator lcm. It
verifies all stated gcds, the two exact image factorizations, the bad
radius divisor, the good radius equality, interval preservation and
the rational bounds implying endpoint failure or preservation. It also
tests eight rational alignments at each index, including ratios below
one and the conformal control, checking their explicit common Gaussian
multiplier and the uniform bound for each fixed alignment.

Root derived the two parameter calculations. Sol independently checked
both proofs and the exact good radius; Astra independently constructed
the primitive bad targets and confirmed their full denominator cost.
Root and Luna independently derived the fixed-rational extension and
audited its common-multiplier proof.
Numerical arc constants are diagnostics only; the asymptotic conclusions
and uniform upper bound follow from the displayed exact inequalities.

The unresolved step for the main problem is to find a suitable parameter
for arbitrary large endpoint configurations, with the least target
radius controlled after **all** rows are retained. The factorizations
above use the specific Pell identities and four-row source data. They
do not establish that existence statement or a uniform lattice-point
count. The [general content and overlap formulas](endpoint_swap_content_and_overlap.md)
retain the exact missing denominator cost and explain why discarding
target overlap cannot certify endpoint scale from five points onward.
