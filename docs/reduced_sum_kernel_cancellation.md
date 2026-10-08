# A sum-kernel pole can cancel completely at an actual circle prime

For every integer `e>=1` there are eight distinct, primitive,
common-unit Gaussian circle points for which

```text
v_17(D_sum)=4e,
v_17(T)=0,
v_17(denominator(per(S)))=0.                              (1)
```

Here `D_sum` is the exact common denominator of the normalized secant
permutation products, and `T` is the product of the old primitive chord
residues. The same complete cancellation occurs for every positive
Taylor remainder `R_d` from
[sum_kernel_denominator_budget.md](sum_kernel_denominator_budget.md).
All eight points lie in a fixed arc of angular width less than one radian.

The construction is not an endpoint family: two retained directions have
a fixed positive separation while the primitive radius tends to infinity.
It disproves a proposed **primewise** bound on cancellation in terms of
`v_17(T)`. It does not disprove a global inequality involving `log T`,
since residues at other primes can grow in this family.

## 1. Rational half-angle coordinates and the leading pole coefficient

For positive integers `a_i` divisible by four, use `h_i=a_i+i` and the directions
`h_i/bar(h_i)`. Put

```text
n_i=a_i^2+1,
B_ij=1/(1+a_i a_j),
S_ij=sqrt(n_i n_j)/(1+a_i a_j).
```

All permutation products of `S` are rational, and

```text
per(S)=(product_i n_i) per(B).                           (2)
```

Suppose a prime `p` divides only the off-diagonal denominator
`delta=1+a_1 a_2`, while all `n_i` and all other `1+a_i a_j` are
units. Separate the permutations according to whether they use zero,
one, or both directed edges `12,21`. Exactly,

```text
per(B)=C_0+delta^(-1) C_1+delta^(-2) C_2,
C_2=per(B restricted to rows and columns3,...,8).        (3)
```

The functions `C_j` have only unit denominators near the specified local
configuration. Thus the leading pole coefficient is the complementary
six-row permanent. Positivity over the real numbers does not prevent that
coefficient from vanishing modulo `p`.

In normalized pair coordinates, this gives a useful exact conditional
estimate. If the isolated pole has `v_p(x_12)=b` and
`c=v_p(per(S restricted to3,...,8))<b`, then the leading term has order
`c-2b`, strictly below the other terms, and

```text
v_p(per(S))=c-2b.                                      (3a)
```

Thus the denominator loss is exactly the valuation of the complementary
permanent until that valuation reaches the next pole scale. The old
residues alone do not bound this new coefficient, as the construction
below shows.

There is also a limited unconditional no-loss case. If the only positive
real-part valuations at a prime occur on four disjoint edges covering
all eight rows, the permutation transposing all four pairs is the unique
term with the largest pole. Its norm factors are units at that prime,
since `gcd(x_ij,q_ij)=1`. The reduced permanent denominator then retains
the full common-denominator valuation `2 sum_edges v_p(x_ij)`. With fewer
pole pairs, the complementary permanent is again a possible source of
cancellation. No statement about the frequency of these cases in an
endpoint cluster is proved here.

## 2. A simple root modulo seventeen

Use the residue classes

```text
(a_1,...,a_8)=(1,16,12,2,6,11,15,5) mod17.              (4)
```

All eight classes are distinct. Every `a_i^2+1` is a unit, and the only
off-diagonal zero of `1+a_i a_j` is the pair `12`. For the complement
`(12,2,6,11,15,t)`, exact arithmetic in `F_17` gives

```text
per(B_6)(5)=0,
d/dt per(B_6)(5)=2.                                    (5)
```

These finite identities are checked by the exact dynamic-programming
certificate linked below. The derivative includes both the final row
and final column, including the derivative `-2t/(1+t^2)^2` of the
diagonal entry.

Fix `a_1=52` and
`(a_3,a_4,a_5,a_6,a_7)=(12,36,40,28,32)`. These values are divisible
by four and have the residue classes (4). Choose positive `a_2`
with

```text
52a_2+1 ==17^(2e) mod17^(2e+1),       a_2==0 mod4.       (6)
```

The two congruences are compatible because `52` is a unit modulo
seventeen and the moduli are coprime. Equation (6) implies
`v_17(delta)=2e` exactly. Define the locally analytic rational function

```text
F_e(t)=delta^2 C_0(t)+delta C_1(t)+C_2(t).               (7)
```

At `t==5 mod17`, every denominator in (7) is a unit. By (5),
`F_e(5)==0 mod17` and `F_e'(5)==2 mod17`.

For completeness, the lifting step is elementary. If an integer `t_k`
satisfies `F_e(t_k)==0 mod17^k`, then

```text
t_(k+1)=t_k+17^k b,
b=-(F_e(t_k)/17^k)/F_e'(t_k) mod17
```

gives a root modulo `17^(k+1)`. All divisions in the rational function
are by local units. Iterating gives an integer `a_8==5 mod17` such that

```text
F_e(a_8)==0 mod17^(4e).                                (8)
```

The congruence `a_8==0 mod4` may be imposed after lifting. Equations
(2), (3), and (8) show that `per(S)` is seventeen-adically
integral, since the factor removed in (7) has valuation exactly `4e`.

## 3. An actual common-unit Gaussian realization

Take positive representatives in (6), and after the lift (8) impose
`a_8==0 mod4` by the Chinese remainder theorem. All eight `h_i=a_i+i`
have odd norm and are coprime to their conjugates. Moreover

```text
-i h_i=1-a_i i ==1 mod(1+i)^3.                          (9)
```

Fix the primary Gaussian-prime convention once and for all: an odd
Gaussian integer is primary when it is congruent to one modulo
`(1+i)^3`. Products and conjugates of primary integers are primary,
and each odd Gaussian associate class has a unique primary representative.
Equation (9) therefore makes the factorization unit of every `h_i`
equal to the same unit `i`. This convention is fixed throughout the
infinite family; no prime representative has to change with `e`.

Let `L` be the Gaussian least common multiple of the `h_i`, chosen
primary, and put

```text
z_i=bar(L) h_i/bar(h_i).                                (11)
```

Every `z_i` is a Gaussian integer, all have norm `N(L)`, and their
prime-factor units all equal `-1`. Divide their primary common Gaussian
gcd to obtain
a primitive circle. Distinctness is preserved because the `a_i` have
distinct residues modulo seventeen. The rational directions and the
normalized kernel `S` are unchanged by this common division.

## 4. Exact pole size and all old residues

For `g_ij=gcd_G(h_i,h_j)`, the primitive relative numerator is, up to a
sign,

```text
h_i bar(h_j)/N(g_ij)
 =[a_i a_j+1+i(a_j-a_i)]/N(g_ij).                        (12)
```

All `n_i` are seventeen-adic units, so every `N(g_ij)` is a unit there.
All differences `a_i-a_j` are units by (4). Consequently

```text
v_17(t_ij)=0 for every pair,
v_17(x_12)=2e,
v_17(x_ij)=0 for every other pair.                       (13)
```

The exact denominator theorem `D_sum=L_X^2` now gives
`v_17(D_sum)=4e`. Equations (8) and (13) prove (1): the entire local
pole of the common denominator disappears from the reduced permanent.

For every Taylor degree `d`, all coefficients of `P_d` are dyadic
rationals. This follows from
`(1-u)^(-1/2)=sum_k binom(2k,k)u^k/4^k`; the exponent-one factors
are geometric series. Since each `q_ij` is a seventeen-adic unit,
every `P_d(u)` is seventeen-adically integral. Thus every

```text
R_d=per(S)-P_d(u)
```

is integral at seventeen, while its natural termwise denominator still
has valuation `4e`. The real remainders are positive because the inputs
are in a proper arc and every `u_ij>0`.

## 5. Scope of the obstruction

The points have arguments `2 arctan(1/a_i)` before a common rotation.
All `a_i>=4`, so their angular span is less than one. The retained rows
`a_3=12` and `a_4=36` have the fixed positive separation
`2(arctan(1/12)-arctan(1/36))`.

Meanwhile `a_2` tends to infinity with `e`. The pair `h_1,h_2` has
relative primitive norm at least `(a_2^2+1)/2705`, because their Gaussian
gcd divides the fixed `52+i`. Every primitive circle realizing this pair
has radius tending to infinity. Thus the normalized arc constant grows;
this family is not at endpoint scale.

In particular there are no universal constants `c,B` giving the
primewise estimate

```text
v_17(D_sum/denominator(per(S))) <= c v_17(T)+B,
```

and the analogous statement fails for every `R_d`. The left side is
`4e` and the residue valuation on the right is zero. This says nothing
against an estimate using the **global** size `log T`: other residue
primes and other norm factors can grow. A successful global denominator
argument must use that information or the shrinking endpoint span.

## Exact verification

[check_reduced_sum_kernel_cancellation.py](check_reduced_sum_kernel_cancellation.py)
verifies (5), the Hensel lifts for `e=1,2,3,4`, the primary congruences,
and actual primitive odd-norm Gaussian circle realizations.
It checks every primitive pair coordinate and its seventeen-adic valuation.
For `e=1`, it additionally computes and reduces the full rational permanent
and the first Taylor remainder, independently of the modular calculation.
All checks pass. The lifting proof above, rather than the finite range of
the check, establishes the arbitrarily large cancellation order.
