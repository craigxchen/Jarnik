# Integral pair defects and the shared-prime normalization budget

This note continues the maximal-order question left open in
[multiquadratic_nearsquare_audit.md](multiquadratic_nearsquare_audit.md).
Shared conductor factors do permit useful integral combinations of
near-square defects. The correct combination includes a small residue
term that the raw product omits. We prove that construction, compute
its exact trace-lattice index, and show that its improvement in
covolume matches the loss in smallness of the displayed basis.
For three rows, further integral combinations give exact identities
with the existing primitive residues. These results retain the actual
shared norm factors. They do not establish the uniform endpoint bound
or exclude other maximal-order constructions.

## 1. A valid integral normalization at every prime

Let

```text
h_i=x_i+i t_i,       q_i=N(h_i),       gcd_G(h_i,bar(h_i))=1,
E=Q(sqrt(q_1),...,sqrt(q_m)).
```

Choose the signs of the numerators so that `x_i>0`. For a pair define

```text
G_ij=N(gcd_G(h_i,h_j)),
h_ij=h_j bar(h_i)/G_ij=x_ij+i t_ij,
q_ij=N(h_ij)=q_i q_j/G_ij^2.                       (1)
```

The pair numerator is again a conjugate-primitive Gaussian integer.
The two exact coordinate identities are

```text
G_ij x_ij=x_i x_j+t_i t_j,
G_ij t_ij=x_i t_j-x_j t_i.                         (2)
```

The short-arc setting makes `x_ij>0` as well. Put

```text
delta_i=sqrt(q_i)-x_i,
psi_ij=sqrt(q_ij)-x_ij.
```

Then the following expression is an algebraic integer in `E`:

```text
psi_ij=(delta_i delta_j+x_i delta_j+x_j delta_i-t_i t_j)/G_ij.
                                                               (3)
```

Indeed its numerator is `sqrt(q_i q_j)-x_i x_j-t_i t_j`.
Equation (2) identifies its quotient as
`sqrt(q_ij)-x_ij`; both terms are integral. This proves integrality
at every prime and to the full valuation of the denominator. It
does not merely check reduction modulo the radical of `G_ij`.

In particular, the small rational correction `-t_i t_j/G_ij`
allows a valid construction even though `delta_i delta_j/G_ij`
alone need not be integral. At a shared conductor prime the latter
product can be a unit divided by that prime, as proved in the
previous note. There is no conflict: the mixed terms in (3) change
the local numerator.

The resulting defect satisfies exactly

```text
psi_ij=t_ij^2/(sqrt(q_ij)+x_ij),
Norm_(E/Q)(psi_ij)=(-t_ij^2)^(n/2),                (4)
```

where `n=[E:Q]` and the norm formula assumes `q_ij` is nonsquare.
Thus (3) constructs the primitive pair near-square defect already
attached to the actual reanchored pair. It does not produce a
smaller residue than `t_ij`.

## 2. Exact lattice improvement, with its archimedean cost

Assume that the squareclasses represented by

```text
1, q_i, q_i q_j (i<j)
```

are pairwise distinct. The available six-wise separation is more
than sufficient; separation of nonempty products of at most four
norms is enough here. Let `L=1+m+binom(m,2)`. Consider the integral
lattices

```text
Lambda_raw = span_Z(1, delta_i, delta_i delta_j : i<j),
Lambda_pair= span_Z(1, delta_i, psi_ij : i<j).
```

Integral triangular changes of basis give

```text
Lambda_raw = span_Z(1, sqrt(q_i), sqrt(q_i q_j)),
Lambda_pair= span_Z(1, sqrt(q_i), sqrt(q_ij)).       (5)
```

Because `sqrt(q_i q_j)=G_ij sqrt(q_ij)`, we obtain the exact index

```text
Lambda_raw subset Lambda_pair,
[Lambda_pair:Lambda_raw]=product_(i<j) G_ij.       (6)
```

Distinct quadratic characters are orthogonal for the trace inner
product. Consequently the exact Gram determinants are

```text
det Gram(Lambda_pair)=n^L product_i q_i product_(i<j) q_ij,
det Gram(Lambda_raw) =n^L product_i q_i product_(i<j) q_i q_j.
                                                               (7)
```

Their ratio is precisely `product G_ij^2`. This is a genuine
enlargement inside the maximal order, not the identical-lattice
change from radicals to raw products discussed in the earlier note.

For a tuple with the central pair moments,

```text
log|h_i|=W/4+o(W),       log|h_ij|=W/4+o(W),
log max(1,|t_i|,|t_ij|)=o(W),
```

equation (1) yields `log G_ij=W/4+o(W)`. The chosen positive real
values satisfy

```text
log delta_i=-W/4+o(W),
log(delta_i delta_j)=-W/2+o(W),
log psi_ij=-W/4+o(W).                             (8)
```

Each pair basis vector becomes larger by the same leading factor
by which the lattice covolume decreases. An exact comparison is

```text
psi_ij/(G_ij delta_i delta_j)
 = [t_ij^2/(t_i^2 t_j^2)]
   [(sqrt(q_i)+x_i)(sqrt(q_j)+x_j)]
   /[sqrt(q_i q_j)+x_i x_j+t_i t_j].              (9)
```

Its logarithm is `o(W)` in the stated setting: the second ratio
tends to two, and the nonzero residue factors have subpower height.
Thus, for this displayed normalization, the shared-prime index
does not supply an unmatched leading gain in a comparison between
covolume and the product of the small distinguished basis values.
Equations (6)--(9) concern these lattices and this comparison only;
they are not an impossibility theorem for other vectors or for the
whole maximal order.

## 3. Three balanced words in the same quadratic character

Here there is a useful further exact identity. We first state it
in the uncorrected independent-block model; the formula is not
asserted unchanged when row corrections `K_i` are present.

Group the independent oriented blocks by their restriction to
three rows, writing them as `H_abc` with `a,b,c in {0,1}`. Their
norms are `n_abc`. Each row numerator is the product of the blocks
whose corresponding bit is one. Define the primitive Gaussian
numerators of the three balanced words by

```text
g_1=H_100 H_010 bar(H_001) H_110^2 H_111,
g_2=H_100 bar(H_010) H_001 H_101^2 H_111,
g_3=bar(H_100) H_010 H_001 H_011^2 H_111.
```

They represent, respectively, the phase combinations
`phi_1+phi_2-phi_3`, `phi_1-phi_2+phi_3`, and
`-phi_1+phi_2+phi_3`. Set

```text
D_1=n_110,     D_2=n_101,     D_3=n_011,
Q=n_100 n_010 n_001 n_111,
X_a=Re(g_a),  Delta_a=D_a sqrt(Q)-X_a.
```

All prime powers are retained; `Q` need not be squarefree. Direct
factorization gives `N(g_a)=D_a^2 Q`, so all `Delta_a` are integral
and lie in the same quadratic subfield. Their exact rational
differences are

```text
D_2 Delta_1-D_1 Delta_2=-2 t_1 t_23,
D_3 Delta_1-D_1 Delta_3=-2 t_2 t_13,
D_3 Delta_2-D_2 Delta_3=-2 t_3 t_12.               (10)
```

For example,

```text
X_1/D_1-X_2/D_2
 =sqrt(Q)[cos(phi_1+phi_2-phi_3)-cos(phi_1-phi_2+phi_3)]
 =2 sqrt(Q) sin(phi_1) sin(phi_3-phi_2).
```

Since `sqrt(q_1) sqrt(q_23)=D_1 D_2 sqrt(Q)`, this becomes
`D_2 X_1-D_1 X_2=2 t_1 t_23`, proving the first line of (10).
The other two follow by the same exact coordinate calculation.
Eliminating the three defects gives

```text
D_3 t_1 t_23-D_2 t_2 t_13+D_1 t_3 t_12=0.        (11)
```

Here `G_12=n_111 D_1`, `G_13=n_111 D_2`, and
`G_23=n_111 D_3`; hence (11) is precisely the usual three-row
primitive residue cocycle after cancellation of `n_111`.

The endpoint exponents explain why (10) does not immediately
vanish. With inherited uniform weights `log n_abc=W/8+o(W)`,

```text
log D_a=W/8+o(W),     log Q=W/2+o(W),
log|g_a|=3W/8+o(W).
```

The original phase precision is `O(exp(-W/4))`. Thus
`0<Delta_a=|g_a|(1-cos(arg g_a))<=exp(-W/8+o(W))`
whenever the indicated word is nonreal. The products `D_b Delta_a`
have upper bound `exp(o(W))`, matching the actual nonzero integer
residue products on the right of (10). The normalization reaches
the integer scale without crossing it. If a word is real its
defect is zero, and (10) still holds; the displayed strict
positivity is not needed for the identities.

These formulas were checked using exact integer Gaussian arithmetic
on 1,000 independent seven-block configurations, with distinct split
prime supports and exponents from one to three. Every pair
normalization and all three lines of (10), equivalently their real
coordinate versions, passed. The algebraic derivation proves the
unrestricted statement.

## 4. Products of all pair defects retain a cut spectrum

The normalized pair defects also allow a calculation that keeps
the growing field degree and all graph incidences. Include anchor
zero, set `q_0=1`, and let the complete graph have `k=m+1` vertices.
For each edge `e={i,j}` use the actual primitive pair norm `q_e`,
positive real coordinate `x_e`, and nonzero residue `t_e`. Define

```text
psi_e=sqrt(q_e)-x_e,
lambda_e=log((sqrt(q_e)+x_e)/|t_e|),
Psi=product_e psi_e.
```

For an embedding `sigma` of `E`, put
`sigma(sqrt(q_i))=epsilon_i sqrt(q_i)` and `epsilon_0=1`.
Equation (1) shows that the edge character is exactly
`epsilon_i epsilon_j`. Therefore

```text
log|sigma(Psi)|
 =sum_e log|t_e|-sum_e lambda_e
    +2 sum_(e crosses the epsilon cut) lambda_e.  (12)
```

This holds even if the realized sign patterns form a proper
subgroup of the full binary cube. Four-wise norm separation makes
the edge characters nontrivial and pairwise distinct, so their
orthogonality gives the exact mean and variance

```text
(1/n)sum_sigma log|sigma(Psi)|=sum_e log|t_e|,
(1/n)sum_sigma(log|sigma(Psi)|-sum_e log|t_e|)^2
 =sum_e lambda_e^2.                              (13)
```

In particular,

```text
|Norm_(E/Q)(Psi)|=product_e |t_e|^n.              (14)
```

For equal leading weights, write `lambda_e=W/4+b_e` and let
`b=max_e |b_e|`. The complete graph cut bound gives

```text
max_sigma log|sigma(Psi)|
 <=sum_e log|t_e|+floor(k/2) W/4+binom(k,2)b.      (15)
```

The distinguished embedding has leading logarithm
`-binom(k,2)W/4`; the largest possible conjugate has only a linear
`kW/8` leading exponent when the errors in (15) are negligible.
Thus the cut structure genuinely improves a crude estimate that
would allow every edge factor to be large at once. Nevertheless
the exact norm is still (14). Character orthogonality, rather
than an omitted denominator, accounts for the compensation.
If every residue has modulus one, all these defects and their
product are units in the maximal order.

The identities leave open constructions involving carefully chosen
sums of several characters, divisibility of those sums at specified
prime ideals, or additional restrictions on their traces. Neither
the raw-product obstruction nor the valid pair normalization proves
that such a construction is impossible. The concrete progress here
is the integral formula (3), the exact index and matching size cost
(6)--(9), and the reduction of the first repeated-character
construction to the existing residue identity (10)--(11).

## 5. The rank-three matrix and a nonzero minor

The metric matrix suggested by the pair construction has the exact
form

```text
M_ij=sqrt(q_i)sqrt(q_j)-x_i x_j-t_i t_j
    =G_ij psi_ij,
M=r r^T-x x^T-t t^T,       r_i=sqrt(q_i).          (16)
```

Its diagonal is zero and its rank is at most three. Thus every
four-by-four determinant vanishes identically. Entrywise division
by `G_ij` is not generally a diagonal change of basis; the matrix
`(psi_ij)` need not retain that rank bound. In particular one must
retain the conductors when using a vanishing minor.

For four indices, put

```text
a=G_12 G_34 psi_12 psi_34,
b=G_13 G_24 psi_13 psi_24,
c=G_14 G_23 psi_14 psi_23.
```

The vanishing determinant is exactly

```text
a^2+b^2+c^2-2ab-2ac-2bc=0.                       (17)
```

There can be large common integer factors in its coefficients, and
these should be removed. For example, in the uncorrected uniform
four-row model, group the blocks by their restriction to these four
indices and write their norms as `n_T`, `T subset {1,2,3,4}`. The
common factor of the three matching coefficients is exactly

```text
F=product_(|T|=3) n_T * n_{1234}^2.              (18)
```

Indeed a cut meeting the four rows in three indices contributes once
to every matching; one meeting all four contributes twice; a cut of
size two contributes to just its matching. Pairwise independence of
the blocks proves the claimed exact gcd. With uniform weights,

```text
log F=3W/8+o(W),
log(G_12 G_34/F)=W/8+o(W),
```

and likewise for the other matchings. Thus (17) divided by `F^2`
is an integral algebraic identity with smaller coefficients. Each
of `a/F,b/F,c/F` has logarithmic size at most `-3W/8+o(W)` in
the distinguished embedding. Its vanishing remains the rank
identity (16), rather than a nonzero small integer. No additional
prime divisibility can be inferred from that vanishing alone.

The first nonzero principal minor gives a more direct check of the
norm budget. For any three indices let `P_3` be the zero-diagonal
matrix with off-diagonal entries `psi_12,psi_13,psi_23`. Then

```text
det(P_3)=2 psi_12 psi_13 psi_23,
|Norm_(E/Q)(det(P_3))|
 =2^n |t_12 t_13 t_23|^n.                        (19)
```

If the three edge characters are nontrivial and distinct, their
product is the trivial character, and expansion also gives

```text
Tr_(E/Q)(det(P_3))
 =2n [q_1 q_2 q_3/(G_12 G_13 G_23)-x_12 x_13 x_23].
                                                               (20)
```

The rational number `q_1 q_2 q_3/(G_12 G_13 G_23)` is an integer:
it is the product of the three integral radicals
`sqrt(q_12),sqrt(q_13),sqrt(q_23)`. Under the endpoint size bounds
the difference in (20) is positive and has upper bound
`exp(W/4+o(W))`; it is not forced below one. Formula (19) shows
that this nonzero minor has no uncharged conductor factor in its
norm. More elaborate maximal-order combinations of minors are not
settled by this calculation.
