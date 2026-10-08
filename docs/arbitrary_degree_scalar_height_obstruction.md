# No scalar-covolume height gain in any positive multidegree

For **every** positive integer degree vector `d=(d_1,...,d_m)` with
even total `D=2E` and `max d_i<=E`, define

```text
r(d) = [x^E] product_i(1+x+...+x^d_i)
       -[x^(E-1)] product_i(1+x+...+x^d_i),
s(d) = (1/2) #{S subset [m]: d(S)=E},
A(d) = D 2^(m-3)-sum_(S subset [m]) max(0,d(S)-E).
```

Then, for `m>=4`,

```text
A(d) >= 2 min(s(d),r(d)-1).                            (1)
```

If `m>=5` and `max d_i<E`, the inequality is strict. If `max d_i=E`,
then `r=1` and `A=0`, so there is no numerical-zero relation to exploit.
On at most three supported rows every nonzero invariant in a fixed
admissible multidegree is a bracket monomial and has nonzero value at
distinct directions.

This closes the arbitrary-multidegree version of the scalar/covolume
route left open in
[invariant_degree_variation.md](invariant_degree_variation.md). It
supersedes the degree-one/two restriction in
[all_mixed_degree_scalar_height_obstruction.md](all_mixed_degree_scalar_height_obstruction.md)
for this precise method. It does **not** prove the fair-profile exclusion:
there is no new lower bound here for an individual irreducible numerical
zero whose scalar restrictions all vanish.

The proof reduces arbitrary degrees to a cone with explicitly described
rays. Its only finite part consists of **32** possible residual degree
vectors. The certificate is
[check_arbitrary_degree_scalar_height_obstruction.py](check_arbitrary_degree_scalar_height_obstruction.py).

## 1. What inequality (1) says about the actual arithmetic method

Use the actual full-cut Gaussian model with fixed `m` and fixed multidegree,
inherited block
logarithms `w+o(w)`, and negligible correction and primitive pair-residue
logarithms. The previously proved unequal-degree graph-content formula
gives primitive numerical evaluation height `A(d)w+o(w)`. Complementary
weighted-balanced scalar congruences force an outside-kernel numerical
relation to have coefficient logarithm at least `2w-o(w)`.

Let `rho` be the actual scalar-map rank, and `K` its kernel. The rank need
not equal `s(d)`. The relation lattice has rank `r-1`, since its numerical
evaluation functional is nonzero at distinct directions. If evaluation
on `K` is nonzero, its first guaranteed outside relation has successive-
minima denominator `rho`; if evaluation on `K` vanishes, that denominator
is `rho-1`, when an outside relation exists at all. Therefore every
applicable denominator `b` obeys

```text
1<=b<=min(s,r-1).
```

The full-covolume guarantee has leading coefficient `A/b`. Inequality
(1) makes this at least two, so it cannot yield a strict comparison with
the complementary scalar threshold. This statement is independent of
unproved rank assertions about unequal-degree scalar maps.

The inherited weight amplification under support projection remains
unchanged. An independently obtained small irreducible factor might
still be excluded by a stronger arithmetic theorem. Neither a new
individual-factor lower bound nor an all-cut restriction-image theorem
is supplied here. Also, near-uniform extraction alone does not assert
the negligible-remainder hypotheses of the actual model in this section.

## 2. Forced bracket removal preserves height and dimension

Suppose `max d_i<E` and a pair has `d_i+d_j>E`. Put

```text
k=d_i+d_j-E>0,
d_i'=d_i-k=E-d_j,
d_j'=d_j-k=E-d_i,
d_l'=d_l   for l!=i,j,
E'=E-k.
```

All residual degrees remain positive and satisfy `max d_l'<E'`. In
fact the sum of all degrees outside the pair equals `E'`; there are at
least two such positive degrees when `m>=4`.

Every bracket graph of degree vector `d` has at least `k` copies of the
edge `ij`: the outside degrees can absorb at most `2E-d_i-d_j` of the
`d_i+d_j` stubs at this pair. Integral bracket straightening consequently
shows that every invariant is divisible by `Delta_ij^k`. Dividing gives
an isomorphism of invariant spaces, with inverse multiplication by this
bracket power. Thus

```text
r(d')=r(d).                                             (2)
```

The scalar cuts have an especially simple change. Every old balanced
cut separates `i,j`, and its weight decreases by `k` under the removal.
It is therefore still balanced for `d'`. Among new balanced cuts that
do not separate the pair, the only two oriented possibilities are the
pair itself and its full complement: `d_i'+d_j'=E'`, and every outside
degree is positive. Hence

```text
s(d')=s(d)+1.                                          (3)
```

To see height preservation without any actual-point existence
assumption, extend `A` to nonnegative real vectors. For independent
uniform signs,

```text
A(d)=2^(m-3)[D-2 E_epsilon |sum_i d_i epsilon_i|].       (4)
```

Write `x=d_i`, `y=d_j`, `h=x+y`, `t=2E-x-y=E'`, and let `Y` be the
signed sum over the other rows. Then `|Y|<=t<h`. Averaging only over
the two pair signs gives

```text
E_pair |Y+x epsilon_i+y epsilon_j|
  = (h+max(|Y|,|x-y|))/2.
```

The residual pair has sum `t`, the same difference `x-y`, and the same
outside rows. Its conditional mean is obtained by replacing `h` by `t`.
The mean absolute value therefore decreases by `(h-t)/2=k`, while `D`
decreases by `2k`. Formula (4) proves

```text
A(d')=A(d).                                             (5)
```

Repeat the operation and reorder the degrees. The positive integer `E`
decreases each time, so it stops with

```text
0<d_1<=...<=d_m,
d_(m-1)+d_m<=D/2.                                      (6)
```

Throughout this process `A,r` stay fixed and `s` increases. Proving (1)
for (6) therefore proves it for the original vector. No normalization
by a degree gcd is needed; in particular, dividing by a gcd is not used
to overlook a primitive vector of odd total degree.

If initially `max d_i=E`, its invariant space is generated by the star
bracket monomial with all other degrees joined to the largest row.
Thus `r=1`. The largest signed summand dominates the sum of all others,
so its mean absolute value is `E`; (4) gives `A=0`.

## 3. Explicit rays of the residual degree cone

Let `t_j=d_j-d_(j-1)>=0`, with `d_0=0`. The final inequality in (6) is

```text
sum_(j=1..m-4) (m-3-j)t_j >= t_(m-2)+2t_(m-1)+t_m.      (7)
```

The variable `t_(m-3)` has coefficient zero. The cone in (7) is generated
by its nonnegative singleton directions on the left and at that zero
coefficient, together with a direction pairing one positive left
coefficient with one of the three right coefficients. To prove this,
distribute each unit of required right-hand mass among the available
left-hand masses. Any unused left mass is a singleton direction; the
zero variable is independent. This constructs a nonnegative ray
decomposition explicitly.

In degree coordinates these rays have a trailing support of some size
`q>=4`, preceded by zeros. On that support the four possible forms are

```text
U_q = (1,...,1),
V_q = (1^(q-3),(q-3),(q-3),(q-3)),
P_q = (1^(q-2),(q-2)/2,(q-2)/2),
Z_q = (1^(q-1),q-3).                                  (8)
```

For `q=4` these coincide. The last three types with `q>=5` are exactly
the positive/negative pairings in (7). Rational half-coordinates in
`P_q` are allowed because this is a cone decomposition, not a new
integer invariant degree vector.

Formula (4) makes `A` positively homogeneous and superadditive:

```text
A(u+v)>=A(u)+A(v).                                     (9)
```

This is just the triangle inequality for each signed sum. Zero padding
also has the exact effect

```text
A_m(0^(m-q),u)=2^(m-q) A_q(u).                         (10)
```

All rays in (8) have positive `A`. For a full-support ray put `S_n` for
the sum of `n` independent signs. Direct averaging gives

```text
A_q(U_q) = 2^(q-3)(q-2 E|S_q|),
A_q(V_q) = (q-3)2^(q-3),
A_q(P_q) = 2^(q-3)(q-2-E|S_(q-2)|),
A_q(Z_q) = 2^(q-2)-2.                                 (11)
```

For `V_q`, the three large signs have absolute sums `q-3` or
`3(q-3)`, both dominating the small signed sum. For `P_q`, the two
large signs either cancel or have absolute sum `q-2`. For `Z_q`,
the small sum exceeds the large coordinate `q-3` only when all `q-1`
small signs agree, giving the last formula. These observations derive
all the expressions in (11), rather than positing ray heights.

For `q=4`, all four values are two. The `U_q` values for `q=5,6,7,8`
are respectively `5,18,42,116`. For `q>=9`, Cauchy--Schwarz gives
`E|S_q|<=sqrt(q)` and `q-2sqrt(q)>=2`, hence
`A_q(U_q)>=2^(q-2)`. Similarly, when `q>=6`,
`q-2-sqrt(q-2)>=2`, so `A_q(P_q)>=2^(q-2)`.
The `V_q` bound is immediate from (11); the remaining `q=5` values
are positive as well.

Consequently, for `m>=6`, the smallest `A` among the four full-support
rays is `2^(m-2)-2`, attained by `Z_m`. In any ray decomposition of a
positive vector, the coefficients of the full-support rays sum to
`d_1`, since their first coordinate is one and all other rays have first
coordinate zero. Equations (9)--(11) prove

```text
A(d)>=(2^(m-2)-2)d_1>=2^(m-2)-2,       m>=6.            (12)
```

The five-row all-ones ray has height five and must be treated separately.
For `m=5`, condition (6) makes all pair weights at most `E` and all
triple weights at least `E`. Summing the positive deficits for subset
sizes three, four and five gives `D+(3/2)D+(1/2)D=3D`. Thus on this
entire cone

```text
A(d)=D,                        m=5.                    (13)
```

For positive integer degrees of even total, `D>=6`, so the integer
analogue of (12), `A>=6=2^(5-2)-2`, still holds. No invalid height-six
bound is applied to the odd-total all-ones ray.

## 4. All row counts at least ten

The balanced subsets form an antichain because every degree is positive.
Count permutations whose initial chain contains a balanced subset.
Each permutation contains at most one, while a subset of size `j`
occurs as its initial set in `j!(m-j)!` permutations. Therefore

```text
sum_(balanced S) 1/binomial(m,|S|)<=1,
2s <= binomial(m,floor(m/2)).                           (14)
```

The central binomial coefficient divided by `2^m` is nonincreasing in
`m`: it decreases from an even index to the next odd index and is
unchanged from an odd index to the next even index. At `m=10` it equals
`252/1024=63/256`. Hence, for all `m>=10`,

```text
A(d)>=2^(m-2)-2 > (63/256)2^m
                     >=binomial(m,floor(m/2))>=2s.     (15)
```

The strict middle inequality follows from `2^m/256-2>=2`. This proves
more than (1) in these row counts, with no restriction on the degree
sizes and no invariant-dimension calculation.

## 5. A finite, derived boundary for five through nine rows

For each `5<=m<=9`, minimize `A/D` on the finite ray list (8), with
zero-padding rule (10). Superadditivity then gives `A(d)>=c_m D` on
the entire residual cone. The exact minima are

| m | `c_m` | Necessary bound on `D` if `A<=2s` | Residual vectors to check |
|---:|---:|---:|---:|
| 5 | 1 | `D<=10` | 4 |
| 6 | 7/4 | `D<=80/7` | 5 |
| 7 | 3 | `D<=35/3` | 4 |
| 8 | 31/6 | `D<=420/31` | 8 |
| 9 | 9 | `D<=14` | 11 |

The bounds use (14). Thus they are consequences of the proposed failure,
not an arbitrary search cutoff. In each row, enumerate the nondecreasing
positive integer vectors of even total within the displayed bound and
retain only those satisfying (6). There are exactly 32 vectors in all.

Applying the defining formulas for `r,s,A` gives `A>2 min(s,r-1)` for
every one of them. As a compact check, the minimum ratios
`A/min(s,r-1)` among those with nonzero denominator are

```text
m:                  5       6       7        8        9
minimum ratio:     8/3      4      29/6    132/19      8.
```

The only zero-denominator vector in this finite check is the five-row
all-twos vector; it has `s=0`, `r=6`, and `A=10>0`. The checker verifies
the ray minima and every vector using exact rational/integer arithmetic.
Vectors outside these degree bounds have `A>2s` already. This completes
the proof for every residual degree vector in these row counts.

## 6. Four rows and conclusion of the combinatorial proof

For four sorted positive degrees, (6) says that the sum of the largest
two is at most the sum of the smallest two. Ordering forces all four
to be equal, say to `c`. Then

```text
r=c+1,                    A=2c.
```

The dimension is the dimension of homogeneous degree-`c` polynomials
in two independent matching invariants. These formulas can also be read
directly from the generating function and cut sum. Since forced bracket
removal preserved `A,r`, every original nonsaturated four-row vector
satisfies

```text
A=2(r-1)>=2 min(s,r-1).                                (16)
```

Together with (2)--(5), (15), and the 32 finite cases, this proves (1)
for all positive integral multidegrees. It excludes a strict scalar
height gain throughout this family of methods; it imposes no new
restriction on the actual primitive residues beyond the hypotheses
already used to identify the height and complementary threshold.
