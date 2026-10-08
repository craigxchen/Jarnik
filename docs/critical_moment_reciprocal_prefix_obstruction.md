# Reciprocal identities exclude every prescribed barycentric baseline

The [reciprocal Gram condition](critical_moment_reciprocal_gram.md)
excludes every prescribed barycentric baseline on a valid critical
root-order continuation, including every order-preserving nonnegative
shift after column eleven. In particular it excludes the construction in
[the exact cross-product relaxation](critical_moment_exact_crossproduct_relaxation.md).
The contradiction is uniform over its barycentric baseline choices,
positive block gaps, and fixed points; no minimum shift is needed.
A cross rectangle first excludes large shifts with arbitrary tail signs.
Retaining the critical pair root order gives a dimension-free bound
`a_12<46`, and the baseline spacing then gives a contradiction for every
shift size. None of these arguments requires the separately proposed order of
the eight reciprocal sums.

This is a strict separation of two template conditions. It does not
exclude every punctured orthogonal moment template or improve the general
lattice-circle point bound.

## The eleven-column prefix

Let the initial masks be

```
128,64,206,205,35,19,240,8,4,60,195,
```

with sign `+1` at a set bit and `-1` otherwise. Number columns from one,
and put `q_h=1/a_h` for strictly increasing positive magnitudes. Define

```
d_Ah=(S_0h-S_1h)/2,       d_Bh=(S_4h-S_5h)/2.
```

On this prefix, `d_A` is supported exactly at columns three and four,
with values `-1,+1`; `d_B` is supported at five and six, again with
values `-1,+1`. Thus the two supports are disjoint, and

```
sum_prefix d_Ah q_h=-p,       p=q_3-q_4>0,
sum_prefix d_Bh q_h=-r,       r=q_5-q_6>0,
sum_prefix d_Ah d_Bh q_h^2=0.                           (1)
```

The cross Gram equalities `C_ik=gamma-(beta_i-beta_k)^2/2`, for
`(i,k)=(0,4),(0,5),(1,4),(1,5)`, imply

```
C_04-C_05-C_14+C_15=(beta_0-beta_1)(beta_4-beta_5).
```

Dividing by four gives the necessary rectangle identity

```
sum_h d_Ah d_Bh q_h^2
    =(sum_h d_Ah q_h)(sum_h d_Bh q_h).                  (2)
```

For the prefix alone its left side is zero and its right side is `pr>0`.
Therefore no strictly decreasing positive prefix reciprocals satisfy all
the cross Gram conditions when the remaining reciprocals vanish.

## A necessary tail budget

For any finite extension of the prefix, define

```
T_1=sum_(h>11) q_h,       T_2=sum_(h>11) q_h^2.
```

The tail contributions `u,v,w` to the two linear sums and the quadratic
sum in (2) satisfy `|u|,|v|<=T_1` and `|w|<=T_2`. Hence (2) becomes
`w=(-p+u)(-r+v)`. In particular, whenever `T_1<min(p,r)`,

```
(p-T_1)(r-T_1) <= T_2.                                (3)
```

This inequality is valid without any assumptions about how the tail
masks are ordered. Let `N=n-11>0` be the number of tail columns and
`H=a_12` their smallest magnitude. Then `T_1<=N/H` and `T_2<=N/H^2`.
If `H>N/min(p,r)`, (3) consequently forces

```
(pH-N)(rH-N) <= N.                                    (4)
```

Thus the required reciprocal identity prevents all tail magnitudes from
escaping to infinity with the prefix fixed. This retains an explicit
size cost, rather than relying only on convergence to the inconsistent
prefix limit.

## Every prescribed barycentric baseline has fixed early magnitudes

The [barycentric first-moment construction](critical_moment_first_moment_barycentric.md)
starts each prefix increment at two and adds nonnegative corrections
only at heights `delta+e_i` or `delta-e_i`, where

```
delta=(0,0,0,0,-1,-1,1,1).
```

For the first five proper prefixes `t=0,1,2,3,4`, the distances from
`delta` in the seven nonanchor coordinates have `l1` norms

```
4,3,2,7,4.
```

None is an axial correction height. Therefore every baseline permitted
by that construction, irrespective of which later axial visits receive
the corrections, satisfies

```
a_1=2, a_2=4, a_3=6, a_4=8, a_5=10, a_6>=12.
```

The large-shift construction leaves these initial magnitudes fixed.
Consequently

```
p=1/24,       r>=1/60.                                (5)
```

Substituting (5) into (4) gives the particularly simple necessary bound

```
a_12 < 60N+40.                                       (6)
```

Indeed, if `H>=60N+40`, then `N/H<1/60<=min(p,r)`, and

```
(pH-N)(rH-N)
 >= ((60N+40)/24-N)((60N+40)/60-N)
 = N+10/9 > N,
```

contradicting (4).

In the 31-block construction every tail magnitude is its positive
baseline plus a shift `H_l=R exp(t_l)`, with `t_l>=0`. Thus `a_12>R`.
Its exact fixed-point argument imposes `R>=8nA=128n^2`, where `A=16n`.
This exceeds `60(n-11)+40` for every `n>=12`, contradicting (6).

This arbitrary-tail exclusion holds for all the positive gaps and all
the fixed points allowed by the large-shift construction. The stronger
root-order argument below also handles moderate shifts of the prescribed
baseline. Neither argument substitutes
the auxiliary cross-product coefficients for the actual middle
coefficients of a full moment solution.

## A dimension-free bound for valid root-order continuations

The preceding tail budget allowed arbitrary tail signs. The
[critical root-order theorem](critical_odd_moment_root_order.md) supplies
additional information for a full moment solution, or a candidate word
retaining that necessary order. For the cross pair `(1,6)`, the prefix
differing columns and row-one signs are

```
(2,-1), (4,-1), (5,+1), (6,+1), (7,-1).
```

These are the first five entries of the critical pattern
`-,-,+,+,-,-,+,+,...`. Consequently its remaining support signs begin
`-,+,+,-,-,+,+,-,...`; every partial sum of this tail lies in `[-1,1]`.

For decreasing positive reciprocal magnitudes `u_1>...>u_L` on this
tail, Abel summation gives

```
|sum_(j=1..L) sigma_j u_j| <= u_1 <= q_12.              (7)
```

Indeed, with `A_j=sum_(l<=j) sigma_l`, the sum is
`A_L u_L+sum_(j<L) A_j(u_j-u_(j+1))`; taking absolute values and
using `|A_j|<=1` proves (7). The tail length is even and its total sign
sum is zero, so the bound is strict when the tail is nonempty. An empty
tail contributes zero and causes no exception to the argument.

The prescribed early magnitudes give the prefix reciprocal moments

```
p_1^prefix=-11/40+(q_6-q_7),
p_2^prefix=141/1600+q_6^2+q_7^2.
```

Since `0<q_6-q_7<1/12`, one has `|p_1^prefix|<11/40` and
`p_2^prefix>141/1600`. The exact cross-residual identity `p_1^2=p_2`
for the full pair, together with (7), therefore implies

```
sqrt(141)/40 < sqrt(p_2)=|p_1| < 11/40+q_12.
```

Thus every such continuation satisfying the reciprocal identity obeys

```
q_12 > (sqrt(141)-11)/40,
a_12 < 2(sqrt(141)+11) < 46.                          (8)
```

This bound is independent of the number of later columns. Its extra
hypothesis is the critical pair sign order, which is required by the
full moments and retained by the admissible Euler-word construction.
It does not assert (8) for an arbitrary tail whose signs violate that
order, or for arbitrary choices of the first five magnitudes.

## Baseline spacing excludes every nonnegative tail shift

Retain the first five prescribed magnitudes, the critical sign order for
pair `(1,6)`, and suppose

```
A=a_6>=12,       B=a_7>=A+2,       a_12>=B+10.         (9)
```

Every prescribed barycentric baseline has all successive increments at
least two, so it satisfies (9). Any order-preserving nonnegative shift
after column eleven preserves these inequalities. The proof below only
uses (9), so it also applies to any other continuation satisfying these
same conditions.

The prefix signed sum is negative, with absolute value
`11/40-1/A+1/B`. Equation (7) and (9) imply

```
|p_1| <= U(A,B):=11/40-1/A+1/B+1/(B+10),
p_2 >= 141/1600+1/A^2+1/B^2.                          (10)
```

Here `U(A,B)>0`. The gap between the second bound and the square of
the first is strictly positive. To prove this, put `u=1/A`, `v=1/B`,
`w=1/(B+10)`, and `P=11/40`. Then

```
G(A,B):=141/1600+1/A^2+1/B^2-U(A,B)^2
      =141/1600+v^2-(P+v+w)^2+2u(P+v+w).
```

For fixed `B`, this increases strictly with `u`. Since `A<=B-2`,

```
G(A,B) >= G(B-2,B)
        = Q(B)/[80B(B-2)(B+10)^2],
Q(B)=B^4-26B^3-124B^2+6120B+28000.                    (11)
```

The denominator is positive because `B>=14`. For `B=14+t`, `0<=t<=6`,

```
Q(14+t)=56448-1664t-40t^2+30t^3+t^4
       >=56448-1664*6-40*36=45024>0.
```

For `B=20+t`, `t>=0`,

```
Q(20+t)=52800+1960t+716t^2+54t^3+t^4>0.
```

Thus (10) gives `p_2>U(A,B)^2>=p_1^2`, contradicting the exact
cross-residual identity `p_2=p_1^2`. This excludes all prescribed
barycentric baselines and all their order-preserving nonnegative tail
shifts from full moment solutions, without a large-shift threshold.

The hypotheses concern the explicit prefix, its early magnitudes and
spacing, and a necessary critical sign order. This does not exclude
arbitrary increasing magnitudes on that prefix, or other critical
templates, and gives no new endpoint lattice-circle bound.

## Scope: other prefix magnitudes escape the one-sided deficit test

The fixed early magnitudes matter. On the same eleven masks, take
`a_j=2^j`, `1<=j<=11`. For all sixteen cross pairs, the prefix signed
reciprocal sum now satisfies `p_1^2>p_2`. The exact smallest gap is
`695/1048576`, at pair `(0,5)`, and the smallest ratio `p_1^2/p_2` is
`1062961/1048901`, at `(3,7)`. The eight prefix sums `beta_i` have order
`4,5,0,1,2,3,6,7`, the order required of the full reciprocal sums.
These are prefix statements only; they do not establish the final
reciprocal order or the reciprocal identities after extension.

Arbitrary positive initial increments are compatible with the
first-moment cone. On an Euler word with later visits to all fourteen
axial heights `delta+e_i,delta-e_i`, assign any positive increments
`b_t` at all proper prefixes, retaining the prescribed early increments.
Put `v=sum_t b_t(delta-r_t)`. At a later visit to `delta+e_i`, add
`max(v_i,0)`, and at a later visit to `delta-e_i`, add `max(-v_i,0)`.
These corrections cancel `v` coordinate by coordinate. The resulting
increments remain positive, preserve the initial magnitudes, and
satisfy the first-moment equations by summation by parts. The positive
circulation used in the Euler construction supplies those later visits.

Thus the one-sided prefix deficit test used above, together with
first-moment feasibility, does not exclude arbitrary choices of initial
magnitudes. This fixture asserts neither full reciprocal feasibility nor
a full moment solution.

## Verification

The [exact checker](check_critical_moment_reciprocal_prefix_obstruction.py)
replays the prefix masks and heights, verifies the two disjoint supports,
checks the rectangle identity as a polynomial in arbitrary reciprocal
row sums, the initial critical tail pattern, and the exact constants
in both tail estimates, and the rational gap and polynomial expansions
in the spacing contradiction. It also checks all sixteen geometric-prefix
deficits, their minimum values, and the prefix reciprocal order.
No numerical failure of a nonlinear solver is used.
