# The integer five-star lift and its relation lattice

The full four-row Boolean profile has an exact, symmetric five-star lift to
Gaussian integers. Its ten star brackets are precisely the ten edge-block
norms times the original small residues (with the four imaginary row
coordinates supplying the edges incident to the anchor). This is the
integer analogue of the edge-bracket identity in
[the polynomial star pencil](four_row_star_pencil_reduction.md), but it does
not produce the polynomial pencil's constant five-term relation.

The primitive rank-three lattice of integer relations between the five
stars has covolume `exp(w+o(w))`. This gives the generic first-minimum bound
`exp(w/3+o(w))`; its ten edge congruences do not add an exponential
content saving. A stronger first-minimum bound would need arithmetic
control of the actual edge phases or of their coupled CRT classes. None is
claimed here.

## 1. Exact lift, including the correcting factors

Write `I={1,2,3,4}` and suppose

```text
P_a=K_a product_(T subset I, T nonempty, a in T) H_T=X_a+iY_a,
n_T=Norm(H_T),
Delta_ab=Im(conjugate(P_a)P_b)
        =t_ab product_(T containing a,b)n_T.
```

All fifteen `H_T` are odd and conjugate-primitive, their norms are
pairwise coprime, and every `t_ab` is nonzero. Assume
`log n_T=w+o(w)` and
`log Norm(K_a), log max(1,|Y_a|), log max(1,|t_ab|)=o(w)`.
For the endpoint application the `P_a` are conjugate-primitive; this
also makes `Y_a` nonzero once `Norm(P_a)>1`.

Relabel the fifteen blocks as the five vertices and ten edges of `K_5`:

```text
F_0=conjugate(H_I),                  F_a=H_{a},
F_0a=conjugate(H_(I\{a})),          F_ab=H_{ab}  (1<=a<b<=4).
```

Put `n_uv=Norm(F_uv)` for an edge, and define the actual Gaussian stars

```text
Q_0=F_0 product_(a=1)^4 F_0a,
Q_a=K_a F_a product_(b in {0,1,2,3,4}\{a}) F_ab.        (1)
```

These are Gaussian integers even when the `K_a` overlap core primes.
Each contains five distinct core blocks, so

```text
log |Q_u|=(5/2)w+o(w),          u=0,1,2,3,4.          (2)
```

Cancellation of the common `0a` edge gives the exact identity

```text
conjugate(Q_0) Q_a=n_0a P_a.                           (3)
```

Therefore, setting `D_uv=Im(conjugate(Q_u)Q_v)`, one gets

```text
D_0a=Y_a n_0a.                                         (4)
```

For `a,b>0`, equation (3) gives

```text
D_ab=(n_0a n_0b/Norm(Q_0)) Delta_ab.
```

The four shared original blocks of rows `a,b` are `H_ab`, the two
triples containing `a,b`, and `H_I`. Meanwhile `Norm(Q_0)` is the norm
of `H_I` times all four triple norms. Thus their incidence products
cancel exactly:

```text
(n_0a n_0b/Norm(Q_0))
       product_(T containing a,b)n_T=n_ab,
D_ab=t_ab n_ab.                                         (5)
```

Equations (4)--(5) are an `S_5`-symmetric edge-bracket profile with
`t_0a=Y_a`. The four correcting factors belong in `Q_a`; omitting them
would rotate `P_a` and lose the small bracket in (4).

## 2. Exact content and the generic lattice bound

Let `A` be the `2 by 5` integer matrix whose columns are
`(Re Q_u,Im Q_u)`. Every minor `D_uv` is nonzero. The ten edge norms
`n_uv` are pairwise coprime. Write

```text
h=gcd_(u<v)|D_uv|,
Lambda=ker(A:Z^5 -> Z^2).
```

For any two distinct edges `e,f`, the elementary valuation inequality

```text
gcd(t_e n_e,t_f n_f) | t_e t_f      when gcd(n_e,n_f)=1
```

shows that `h | |t_e t_f|`. Hence `log h=o(w)`. By the saturated-kernel
complementary-minor formula, already proved generally in
[the Gaussian dual-lattice note](auxiliary_form_route.md),

```text
covol(Lambda)=sqrt(sum_(u<v)D_uv^2)/h=exp(w+o(w)),
det(basis(Lambda) on the three columns complementary to {u,v})
        = +/-D_uv/h.                                   (6)
```

Minkowski's first theorem in rank three consequently supplies **one**
nonzero primitive integer relation

```text
sum_u lambda_u Q_u=0,
max_u |lambda_u| <= exp(w/3+o(w)).                       (7)
```

This is a first-minimum statement, not a bound on all three successive
minima. In the other direction, a relation supported on any three stars
is proportional to their three pair brackets. Those brackets have
pairwise-coprime edge norms of logarithmic height `w+o(w)`, and their
common divisor has logarithmic height `o(w)`. Its primitive coefficient
height is therefore `exp(w+o(w))`. A relation satisfying (7) must use at
least four stars for large `w`. This support observation gives no lower
bound for a relation supported on four or all five stars.

## 3. What the ten contact congruences do and do not give

For each edge `e={u,v}`, `F_e` divides both `Q_u,Q_v`. Any relation in
`Lambda` therefore obeys

```text
sum_(c notin {u,v}) lambda_c Q_c = 0 mod F_e.          (8)
```

Because `F_e` is odd and conjugate-primitive,
`Z[i]/(F_e)` is the cyclic ring `Z/n_e Z`. These are ten genuine,
coupled three-term congruences on the complementary labels. The
inequality `|lambda_c|<|F_e|` from (7) does not make an individual term
in (8) small: `Q_c` is of size `exp(5w/2+o(w))`, and its residue modulo
`F_e` is uncontrolled by the bracket-height assumptions.

The clean-support coefficient-minor rule in
[the Gaussian dual-lattice note](auxiliary_form_route.md) has the following
valuation-depth version, which permits overlap with the correcting
factors. Fix a split prime `p=pi bar(pi)`, put

```text
b_c=v_pi(Q_c),       b_min=min_c b_c,
d_pi(J)=max(0,min_(c notin J) b_c-b_min),
```

and let `C` have `r` integer relation vectors as rows, with `|J|=r`.
Then

```text
p^d_pi(J) divides det(C_J).                            (9)
```

Indeed, divide the functional `sum_c lambda_c Q_c=0` by
`pi^b_min`. Its coefficients are integral at `pi` and at least one is
a unit. If `d_pi(J)>0`, every minimum-attaining index is in `J`.
Modulo `pi^d_pi(J)`, the complementary terms vanish, and `C_J` kills
a column with a unit entry. Completing that column to an invertible
matrix shows its determinant is zero modulo `pi^d_pi(J)`.
Since `det(C_J)` is an ordinary integer, this is precisely the
ordinary divisibility in (9). The same argument at `bar(pi)` gives
the exponent `d_barpi(J)`; the two exponents combine by taking their
**maximum**, not their sum.

For a core block of norm divisible by `p^e`, choose `pi` in its
orientation and let `T` be its star support, of size one or two. With
`K_0=1` and `kappa_c^+=v_pi(K_c)`, `kappa_c^-=v_barpi(K_c)`, disjoint
rational-prime core supports give the exact valuations

```text
v_pi(Q_c)=e 1_(c in T)+kappa_c^+,
v_barpi(Q_c)=kappa_c^- .                               (10)
```

For `r=1,2`, the complement `J^c` has at least three labels, so it
contains a label outside `T`. Hence

```text
d_pi(J)<=max_c kappa_c^+,
d_barpi(J)<=max_c kappa_c^-,
max(d_pi(J),d_barpi(J))
    <=sum_c(kappa_c^++kappa_c^-).                       (11)
```

Summed over all core primes, the logarithmic modulus supplied by this
rule is at most `sum_c log Norm(K_c)=o(w)`. The same bound holds for
split primes supported only on corrections, by setting `e=0`.
Subtracting `b_min` separately at both orientations already accounts
for shared Gaussian factors and ordinary integer contents. No row
primitivity of the stars is required for (9)--(11).

For `r=3`, a leading core contribution is possible only when `J^c`
is the core edge `T`. In that case `d_pi(J)=e+O(max_c kappa_c^+)`;
otherwise (11) still applies. The maximal-minor identity `D_uv/h` in
(6) accounts for the resulting divisibility exactly, including all
correction and content valuations. Thus this localized minor rule
adds no leading exponential saving beyond the generic covolume.
The coupled residue classes in (8) remain a possible source of a gain.

It would be invalid to delete every prime appearing in a correction
and charge the entire deleted core mass to that correction: `p` may
occur once in `K_c` and to an arbitrarily large exponent in a core
block. The depth calculation above retains those exponents.

The polynomial star-pencil theorem applies to five *polynomials* whose
brackets have degree two. Here the `Q_u` are five numbers. Equations
(4)--(5) hold at one numerical configuration and do not give polynomial
degree reduction, a rank-four coefficient matrix, or a subpower
five-term relation. The exact CRT contact-lattice analysis in
[the shortest-lift note](contact_lattice_shortest_lift.md) likewise
isolates a least-representative height problem; local compatibility
alone does not bound the height of its required frame.

The finite standard-library checker
[check_four_row_integer_five_star_lift.py](check_four_row_integer_five_star_lift.py)
verifies the exact identities (1), (3)--(6) on disjoint split-prime
blocks with nontrivial correcting factors, and verifies the three-term
congruences (8).
It also checks (9)--(11) on actual integer relation lattices with
overlapping correction support, both prime orientations, large core
exponents, and nonprimitive stars.
It also audits the nearest-square threshold and both congruence signs
from [the one-coordinate divisor reduction](four_row_one_coordinate_divisor_reduction.md)
on an existing full-support three-row fixture. It is an identity check,
not a construction of four-row subpower residues.
