# Weaker intrinsic sufficient tests

The uniform endpoint theorem remains unproved. This note weakens the
arithmetic inequalities that would suffice for it; it does not prove any
of those inequalities for arbitrary endpoint tuples. In particular the
eight-point balanced bonus need not have coefficient one. Any fixed
coefficient greater than `2/5` suffices after the correct intrinsic cap.

The normalization and sampling argument are those of
[the finite local test](finite_local_test_attack.md) and Appendix D of
[the multipoint continuation](multipoint_continuation.md). All conductor
weights are log norms, so `W=log(R^2)`. Fix the arc constant `C0=1/2`.
Common Gaussian units are sufficient throughout.

## 1. General even tuple size and the exact threshold

Fix an even integer `k>=6`. After dividing a selected tuple by its own
Gaussian gcd, let `W'` be its conductor weight, let

```
D_k'=sum_layers w (r-k/2)^2,
Q_k'=sum_(layers with r=k/2) w.
```

Every retained layer has `1<=r<=k-1`. Suppose that some fixed
`delta>0`, `B>=0`, and `W0>=0` satisfy the following presently unproved
arithmetic assertion on every such primitive endpoint tuple:

```
D_k' + delta Q_k' <= (k/4) W' + B,       W'>=W0.       (1)
```

Define the symmetric capped test

```
f(0)=f(k)=k/4,
f(r)=(r-k/2)^2 + delta 1[r=k/2],       1<=r<=k-1.
```

Its binomial mean is exactly

```
mu_f = k/4 + g,
g = [delta binom(k,k/2)-k(k-1)/2]/2^k.                 (2)
```

Indeed the uncapped quadratic has mean `k/4`, while the two constant
cuts lose a combined `2[(k/2)^2-k/4]/2^k`. Consequently (1) is a
sufficient local arithmetic target whenever

```
delta > delta_crit(k) = k(k-1)/(2 binom(k,k/2)).         (3)
```

Some exact thresholds are

| k | delta_crit(k) |
|---:|---:|
| 6 | 3/4 |
| 8 | 2/5 |
| 10 | 5/28 |
| 12 | 1/14 |
| 14 | 7/264 |
| 16 | 4/429 |

Thus for any fixed positive proposed balanced-bonus coefficient, some
fixed even tuple size has a positive binomial gap. This observation does
not transfer an arithmetic theorem from one tuple size to another.
The threshold is exact for this capped-test family: equality gives zero
binomial gap, not a uniform bound by this criterion.

## 2. Lifting, sampling, and the height threshold

Let `A` denote the weight of inherited layers constant on the selected
`k` rows, and let `W` be the ambient conductor. Then `W'=W-A`, and
the cap gives exactly

```
sum_layers w f(r) = D_k' + delta Q_k' + (k/4)A
                 <= (k/4)W+B.                         (4)
```

Gaussian gcd division multiplies the normalized endpoint constant by
`exp(-A/4)`, so the normalized tuple remains in the permitted arc class.
The inherited determinant inequality gives

```
(k^2/4)A <= D_k'+(k^2/4)A
           <= (k/4)W-k(k-1)log2.
```

It follows that

```
W' >= ((k-1)/k)W + 4(k-1)log2/k.                       (5)
```

Thus `W>=k W0/(k-1)` makes (1) applicable to every selected tuple.

For completeness, the general-size version of the existing sampling
lemma has exactly the same proof. Let `H=osc(f)`, take a uniform
`k`-subset of an even ambient `M`-point cluster, and put
`mu_f=E f(Bin(k,1/2))`. The ambient determinant supplies

```
sum_layers w (p-1/2)^2 <= W/(4M),
W >= 4(M-1)log2.
```

Sampling with and without replacement differs by at most
`binom(k,2)H/M` per layer. Symmetry gives `b_f'(1/2)=0`, and the
Bernstein second-derivative formula gives
`|b_f''|<=2k(k-1)H`. Therefore

```
|E sum_layers w f(r)-mu_f W|
    <= [3k(k-1)H/(4M)] W.                              (6)
```

Combining (2), (4), and (6) proves the conditional explicit bound

```
M <= [3k(k-1)H/4+B/(2log2)]/g,                         (7)
```

when `M>=k` and the ambient height threshold holds. Below that threshold,

```
M <= 1 + k W0/[4(k-1)log2].                            (8)
```

The cases `M<k` are already bounded. Deleting one point handles odd
cardinality; a largest common-unit class costs at most four; partitioning
an arc of constant `C` costs at most `ceil(2C)`. These are the same
normalization operations as in the eight-point criterion.
For `0<delta<=1`, `H=(k/2-1)^2-delta`.

## 3. Two weaker eight-point targets

Use `V_j` for minority cut sizes `j=1,2,3,4`, so
`W=V_1+V_2+V_3+V_4` and `D=9V_1+4V_2+V_3` on a primitive tuple.

First, setting `delta=1/2` gives the sufficient, unproved inequality

```
D+(1/2)V_4 <= 2W+B,
equivalently 7V_1+2V_2 <= V_3+(3/2)V_4+B.               (9)
```

Its capped test is `(2,9,4,1,1/2,1,4,9,2)`, with
`g=7/256` and `H=17/2`. Formula (7) becomes

```
M <= 13056 + 128B/(7log2).                             (10)
```

A different useful weakening removes two copies of singleton mass:

```
5V_1+2V_2 <= V_3+V_4+B.                               (11)
```

Its capped test is `(2,7,4,1,1,1,4,7,2)`, with
`g=5/128` and `H=6`. Hence (11), if proved uniformly above a fixed
height, yields

```
M <= 32256/5 + 64B/(5log2).                            (12)
```

Both (9) and (11) are pointwise weaker than the existing intrinsic
target `7V_1+2V_2<=V_3+V_4+B`, since all the `V_j` are nonnegative.
Neither new target is asserted to hold, and they are not ordered relative
to one another by this observation.

## 4. Literal source-gcd formulations

These targets can be written using actual Gaussian source gcds, with no
choice of independent block model. For a primitive eight-tuple define

```
H_i=gcd_G(z_l : l!=i),
A=Norm(product_i H_i),
L_ij=gcd_G(z_l : l not in {i,j})/(H_i H_j),
F=Norm(product_(i<j) L_ij),
G_j=Norm(product_(|J|=j) gcd_G(z_l : l not in J)), j=2,3,
J_3=exp(V_3),     Q_4=exp(V_4),     N=R^2.
```

The quotients `L_ij` are Gaussian integers by primitivity. The existing
source-gcd dictionary gives `A=exp(V_1)`, `F=exp(V_2)`, and
`G_2=A^7 F`. Threshold counting also gives

```
G_3=A^21 F^6 J_3,       N=A F J_3 Q_4.                 (13)
```

To check the extra identity, a layer with zero side of size `r<=3`
occurs in exactly `binom(8-r,3-r)` deleted-triple gcds. These counts
are `21,6,1`. Both Gaussian orientations are included, so all minority
sizes and nested prime-power multiplicities are accounted for.

The half-bonus target (9) is therefore equivalent to each of

```
A^16 F^6 <= exp(2B) N^2 Q_4,
G_2 G_3 <= exp(2B) N^3 A^11.                           (14)
```

The reduced-singleton target (11) is equivalent to

```
A^6 F^3 <= exp(B) N,
G_2^3 <= exp(B) N A^15.                                (15)
```

The prior target has `A^8 F^3` in place of `A^6 F^3`, or `A^13`
in place of `A^15` on the right in the second formulation. Thus the
new source-gcd estimate allows a larger deletion content precisely at
the previously uncontrolled singleton layers.

More generally the primitive inequality

```
A^a F^b <= exp(B) N                                    (16)
```

has the inherited capped test
`(2,a+1,b+1,1,1,1,b+1,a+1,2)`. Its binomial gap above two is

```
g=(8a+28b-127)/128.                                    (17)
```

Consequently every fixed pair `(a,b)` with `8a+28b>127` supplies a
sufficient arithmetic target. In particular `(6,3)` is sufficient.
As above, a uniform additive constant and a height threshold are part
of the unproved assertion; no growing source factor may be hidden in
either one. Equations (14)--(17) are exact reformulations and sufficient
criteria, not a new unconditional estimate on actual short arcs.

## 5. An actual mixed-layer class satisfying a sufficient test

There is a concrete application of (17) beyond profiles containing only
even layers. The point `(a,b)=(0,5)` has positive gap `13/128`; its
arithmetic target is simply `F^5<=exp(B)N`, and it does not charge
singleton mass at all.

Assume that the actual primitive eight-point tuple satisfies the
prime-by-prime aggregate parity hypothesis

```
sum_i a_i(p) is even for every split prime p.            (18)
```

The proved [aggregate-parity certificate](aggregate_parity_balanced_character_height.md)
retains an explicit nonnegative nested-layer cancellation `K` and gives,
at `C=1/2`,

```
V_4+10V_3 >= 5V_2+2K+140log2.                          (19)
```

Thus one has the unconditional estimate within this actual arithmetic
class

```
5log F-log N <= 9V_3-V_1-V_2-2K-140log2.               (20)
```

In particular the sufficient local test `F^5<=N` is proved whenever

```
9V_3 <= V_1+V_2+2K.                                   (21)
```

No restriction on individual prime sizes, their exponents, or equality
of block weights is used. A simple subclass consists of all tuples
satisfying (18) with no triple layers: singleton, pair, and balanced
layers are allowed, with arbitrary singleton mass. In this subclass
the stronger estimate follows directly from (19):

```
A F^6 <= 2^(-140) exp(-2K) N.                          (22)
```

This is a proved application of an existing arithmetic certificate to
a new sufficient test, not a new independent character inequality.
The earlier target `A^8 F^3<=exp(B)N` does not follow from (19) alone
when singleton mass is unrestricted. The reduced test identifies
exactly which singleton cost can instead be omitted on this subclass.

Condition (18), and likewise (21), is not established for an arbitrary
eight-row restriction of a larger endpoint cluster. The argument therefore
does not supply the required universally valid local test or a new
general point-count bound. In particular one must not perform uniform
row sampling while silently retaining a primewise parity condition
which sampling can destroy.

## Verification

The [exact checker](check_weaker_intrinsic_finite_tests.py) checked all
displayed binomial means and the
critical coefficients for `k=6,8,...,16`. An exhaustive prime-power
audit of the source-gcd identity (13) checked all 329 nonconstant sorted
allocation vectors of length eight with entries between zero and four
and minimum zero. It computed the actual gcd valuation sums over every
deleted singleton, pair, and triple, then compared them with the
threshold-layer formulas. These finite checks support the elementary
proofs and make no assertion that the fixtures are endpoint tuples.
