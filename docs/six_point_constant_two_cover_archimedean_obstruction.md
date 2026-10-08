# An archimedean obstruction for constant two-torsion covers

Assume the branch Galois group of the genus-five curve
`C: y^2=Delta(s,t)` contains `A_12`, and assume the triangle pole divisor
consists of ten rational closed points, each with a biquadratic residue
field of degree four. Then an unramified geometrically connected cover
obtained from any positive-rank **pointwise rational** subspace of J[2]
cannot satisfy the ordinary Runge pole-orbit surplus for the pulled-back
triangle pole function. Archimedean places alone already exhaust its pole
orbits. This remains true after arbitrary finite extensions of the field
of definition and after constant quadratic twists of its equations.

This is a conditional obstruction for constant elementary abelian
2-covers. It does not address the Q-defined degree-1024 multiplication-by-two
pullback with nonconstant group scheme, other pole functions, or weights
whose branch Galois group is smaller. The group hypothesis is established
for the first exact witness below, but not for every fair weight.

## 1. First witness: the branch Galois group is S_12

For weights `(1,-20,250,-1000,1445,-676)`, the exact modular certificates in
[six_point_branch_two_torsion_audit.md](six_point_branch_two_torsion_audit.md)
prove irreducibility and show cycle patterns `(1,11)` at 59 and `(5,7)`
at 73. The group G is consequently transitive and primitive. Indeed any
nontrivial block system has between two and six blocks, each of size less
than eleven. The 11-cycle must fix every block (its induced permutation
has order one or eleven), and it cannot act on a block of that size.
This contradicts its existence.

The fifth power of the `(5,7)` permutation is a 7-cycle with five fixed
points. Jordan's prime-cycle theorem therefore gives `G>=A_12`.
The precise primary-source statement checked here is Theorem 1.1 in
[Gareth A. Jones, Primitive permutation groups containing a cycle](https://arxiv.org/pdf/1209.5169):
a primitive group containing a prime-length cycle with at least three
fixed points contains the alternating group. Its hypotheses hold with
`n=12,p=7`. Finally the existing exact real-root count is two. Complex
conjugation is therefore five transpositions and is odd, giving `G=S_12`.
The conservative estimates below require only `G>=A_12`.

## 2. Field degree forced by a rational rank-m subspace

Let `L/Q` be the branch splitting field and

```
V=J[2]= {w in F_2^12 : sum_i w_i=0}/<1>.
```

Suppose a subspace `W<=V` of dimension `1<=m<=10` is pointwise rational
over a number field K. Let H be its pointwise stabilizer in G and put
`K0=L^H`. Then `K0` is contained in K. In particular,

```
[K0:Q]=[G:H],                  r_infinity(K)>=[K:Q]/2.       (1)
```

Choose even representatives for a basis of W as the columns of a
12-by-m binary matrix, and denote its rows by `v_i in F_2^m`. Then

```
the v_i affinely span F_2^m,            sum_i v_i=0.         (2)
```

Affine spanning follows because the all-one column together with the m
chosen columns is linearly independent; a relation would contradict the
independence of W in the quotient. The sum condition is their even parity.
A permutation fixes each column class precisely when there is one common
translation t such that `v_(permutation(i))=v_i+t` for every i.
Write `n_x=#{i:v_i=x}` and

```
T={t in F_2^m : n_(x+t)=n_x for all x}.
```

The pointwise stabilizer in S_12 has exactly

```
|H_S|=|T| product_x n_x!.                                  (3)
```

Each possible translation has that many bijections between matching row
fibers. Since T acts freely on row patterns with constant multiplicity,
`|T|` divides 12. It is a power of two, so `|T|<=4`.

For `m>=2`, we have

```
product_x n_x! <= (11-m)!.                                 (4)
```

If there are k occupied patterns and `k>=m+2`, the factorial product is
maximized by the multiplicities `(13-k,1,...,1)`, giving (4). If
`k=m+1`, the patterns are an affine basis. Combining both relations in
(2) shows that every multiplicity is even. Hence `m<=5` and the product
is at most `(12-2m)! 2^m`, which is at most `(11-m)!` for `m=2,3,4,5`.
These cases exhaust the possibilities by affine spanning.

Consequently, since `|G|>=12!/2`,

```
[K0:Q] >= 12!/[8(11-m)!]                         (m>=2).    (5)
```

For m=2 one can sharpen (3) to `|H_S|<=9!`. If `|T|=1`, this is (4).
If `|T|=2`, affine spanning requires both translation cosets in F_2^2,
so the four multiplicities are `(a,a,b,b)`, where `a,b>=1,a+b=6`.
Thus `|H_S|<=2(5!)^2<9!`. If `|T|=4`, all four multiplicities are three,
and `|H_S|=4(3!)^4<9!`. Therefore

```
[K0:Q]>=660                                      (m=2).    (6)
```

## 3. Ranks at least two already lose to geometric poles

The base pole divisor has 40 geometric points. An unramified degree-`2^m`
cover has `40*2^m` geometric poles of the pulled-back function; its number
of K-orbits is at most this. Equations (1) and (6) give, when m=2,

```
r_infinity(K)>=330>160.
```

For every `3<=m<=10`, equation (5) similarly gives

```
r_infinity(K)>=12!/[16(11-m)!] > 40*2^m.                    (7)
```

The weakest comparison is m=3, namely `742.5>320`; the ratio of successive
left/right bounds is `(11-m)/2>=1` over the remaining range. An extension
of K only increases the lower bound for its archimedean places, while
the geometric pole bound is unchanged. No finite-place hypothesis is used.

## 4. Rank one: preserve the ten biquadratic pole orbits

For a nonzero element of V its unordered even partition has type
`2+10`, `4+8`, or `6+6`. Its orbit under either A_12 or S_12 has size

```
66, 495, or 462, respectively.                             (8)
```

The A_12 orbits are as large as the S_12 orbits because each S_12
stabilizer contains an odd permutation, for example a transposition
inside one part. Hence `[K0:Q]>=66`.

Moreover K0 has no quadratic subfield. If `G=A_12`, any homomorphism
from G to a group of order two is trivial, since 3-cycles generate G.
If `G=S_12`, its only index-two subgroup is A_12 (transpositions are
conjugate and generate S_12); H contains an odd permutation, so H is
not contained in that subgroup. Galois correspondence proves the claim.

Let E be the biquadratic degree-four residue field of any one base pole.
Every nontrivial subfield of E contains a quadratic field. Thus
`E intersect K0=Q`. Since E/Q is Galois, E and K0 are linearly disjoint,
so the pole remains one closed point of degree four over K0. All ten
base pole orbits consequently remain ten over K0.

Write `e=[K:K0]`. After this extension each base closed point can split
into at most e closed points. A degree-two cover has at most two closed
points above any one base closed point. Therefore

```
number of pole K-orbits upstairs <=20e,
r_infinity(K)>=[K:Q]/2>=33e.                               (9)
```

This again prevents the strict Runge surplus, independent of finite
places, lift fields, and twists of the unramified double-cover equation.
Biquadratic residue fields are needed in this step; merely knowing that
the rational pole degree is four would not prove linear disjointness.

## 5. Applicability and checks

The standard Runge condition requires more pole orbits over K than the
places in S, and S includes all archimedean places. Equations (7)--(9)
make this impossible for the specified covers and pole function. The
statement includes any further field needed to define a point on the cover,
since K was arbitrary subject to W being pointwise rational. It does not
assume that every geometric cover has a pointwise rational deck group.

The required biquadratic-pole hypothesis is precisely the nonexceptional
case of [six_point_triangle_splitting_fields.md](six_point_triangle_splitting_fields.md).
The existence of forty geometric poles and their being disjoint from the
branch locus are the nonresonant triangle facts in the radius-content notes.
No claim is made here that all fair weights have `G>=A_12`, or that the
first witness verifies all ten pole field hypotheses merely because its
branch Galois group is S_12.

`check_six_point_constant_two_cover_budget.py` checks (4) by exhaustive
integer partitions of 12, checks the exact rank-two stabilizer maximum
by enumerating all four row multiplicities, and checks every numerical
comparison and the rank-one orbit sizes. It passes. The previously added
`check_six_point_branch_two_torsion.py` supplies the exact modular inputs;
`check_six_point_branch_arithmetic.py` supplies the exact real-root count.
The first witness and the conditional archimedean obstruction are distinct
claims, and neither closes the nonconstant-group covering route.
