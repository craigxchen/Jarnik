# Five-product capacity audit

## Status

The proposed universal five-product resultant does **not** exist in the
constant-coefficient determinant/Laurent class.

Five products are still the first arity at which four relative phases are
available, but that dimension count is not enough.  Signed reuse of one
conductor atom across the five products forces an exact finite obstruction:

```text
best universal reduced-determinant order = 5/3,
required order                            = 4.
```

The obstruction persists when the exponent-difference lattice has saturated
rank four and there is no common oriented factor.  More generally, a fixed
Laurent invariant with fourth-order angular vanishing cannot have one-use
denominator cost for all signed occurrence patterns.

This note gives exact certificates and replaces the failed universal lemma by
an observed-pattern capacity test.  It does not prove the endpoint bound.

---

## 1. Signed occurrence words

Let

```text
T_0,...,T_4
```

be five primitive signed Gaussian products.  Group the underlying conductor
factors into atoms `Gamma_a` with disjoint rational-prime support.  At product
`r`, an atom has one of three states:

```text
epsilon_a(r) = +1   if Gamma_a occurs,
               -1   if conjugate(Gamma_a) occurs,
                0   if the atom is absent.
```

Thus every atom has an occurrence word

```text
epsilon_a in {-1,0,1}^5.
```

Put `h_a=log Norm(Gamma_a)`.  When the reduced determinant between products
`i,j` is estimated by primitive integer spacing, this atom contributes

```text
c_epsilon(ij) h_a,
c_epsilon(ij) = |epsilon(i)-epsilon(j)|/2.          (1.1)
```

There are three cases:

```text
same state       -> cost 0;
present/absent   -> cost h_a/2;
opposite signs   -> cost h_a.
```

The last case is essential.  Treating two opposite orientations over one
rational prime as two independent atoms violates the one-use conductor model.

---

## 2. Exact determinant linear program

Give the reduced pair determinant on edge `ij` a nonnegative multiplicity
`x_ij`.  Multiplying the resulting inequalities gives angular vanishing order

```text
q = sum_(i<j) x_ij.                                 (2.1)
```

For a set `P` of occurrence words, the product respects one conductor use
exactly when

```text
sum_(i<j) c_epsilon(ij) x_ij <= 1
for every epsilon in P.                             (2.2)
```

Consequently define

```text
Cap_pair(P) = max sum_(i<j) x_ij

subject to x_ij >= 0 and (2.2).                     (2.3)
```

If `Cap_pair(P)>=4`, clearing rational multiplicities gives an actual integer
product with fourth-order angular decay and one-use conductor cost.  If it is
less than four, no multiplication or fractional averaging of reduced pair
determinants can prove the desired five-product lemma.

The exact dual is

```text
min sum_(epsilon in P) y_epsilon

subject to y_epsilon >= 0 and
sum_epsilon y_epsilon c_epsilon(ij) >= 1
for every usable edge ij.                           (2.4)
```

This dual formulation is the useful extraction diagnostic for an actual
selected five-tuple.

---

## 3. Universal signed capacity is `5/3`

Take all `3^5=243` signed occurrence words.  Put

```text
x_ij = 1/6
```

on every one of the ten edges.  For any signed word, the left side of (2.2) is
at most one, so this is feasible and has value

```text
10/6 = 5/3.                                        (3.1)
```

For the matching dual certificate, take the ten full-support sign words whose
`+1` class has size two, one representative of every `2 | 3` cut, and give
each weight `1/6`.  Every edge crosses exactly six of the ten cuts.  Hence

```text
sum_epsilon y_epsilon c_epsilon(ij)
  = 6*(1/6) = 1
```

for every edge, while the total dual mass is `5/3`.  Primal and dual equality
prove

```text
Cap_pair({-1,0,1}^5) = 5/3.                        (3.2)
```

The averaging half of this argument is formalized in

```text
GaussianChain/FiveProductCapacity.lean.
```

Literal integer multiplicities are worse: their optimum is one.

---

## 4. Two restricted models also fail

### Rooted binary words

If the five objects are rooted point products rather than signed transition
products, occurrence words lie in `{0,1}^5`.  The same ten balanced cuts, now
with crossing cost `1/2`, give

```text
Cap_pair({0,1}^5) = 10/3.                           (4.1)
```

The uniform primal weight is `x_ij=1/3`; the dual gives every balanced binary
cut weight `1/3`.  This is closer, but still strictly below four.  The literal
integer optimum is three.

### Consecutive chronological words

If the five selected transitions are consecutive, the nonzero signs of one
prime word must alternate.  There are 63 such words, including the zero word.
Their exact fractional capacity is

```text
Cap_pair(P_chronological) = 2.                      (4.2)
```

A primal optimum puts weight `1/2` on edges

```text
02, 03, 12, 13.
```

A four-word dual of total mass two proves optimality.  Thus chronology does
not repair the determinant product.

The exact rational certificates and literal integer searches are reproduced
by

```text
python experiments/five_product_capacity.py.
```

---

## 5. Rank four does not repair the lemma

The ten balanced full-support sign words already impose the dual bound `5/3`.
Adjoin the four one-hot words

```text
(0,1,0,0,0),
(0,0,1,0,0),
(0,0,0,1,0),
(0,0,0,0,1).
```

Relative to product zero, their four difference vectors are the standard
basis of `Z^4`.  Hence the occurrence-difference lattice is saturated of rank
four.  No nonunit atom has one fixed oriented state in all five products.
Nevertheless the ten balanced words remain, so the capacity is still at most
`5/3`.

Therefore neither of the proposed hypotheses

```text
rank four exponent differences,
no common oriented factor
```

controls the denominator multiplicity needed for fourth-order decay.

---

## 6. Newton-width obstruction for every fixed Laurent invariant

The failure is not limited to monomials in pair determinants.

Write the norm-one phase of product `r` as

```text
u_r = T_r/conjugate(T_r).
```

Let

```text
P(u_0,...,u_4) = sum_(alpha in A) c_alpha u^alpha
```

be a fixed Laurent polynomial.  For an atom with occurrence word `epsilon`,
the exponent of its phase in monomial `u^alpha` is `alpha dot epsilon`.
Clearing Gaussian denominators therefore costs

```text
Norm(Gamma)^(w_epsilon(P)/2),

w_epsilon(P)
  = max_(alpha in A) alpha dot epsilon
    - min_(alpha in A) alpha dot epsilon.           (6.1)
```

A universal one-use bound requires

```text
w_epsilon(P) <= 2
for every epsilon in {-1,0,1}^5.                    (6.2)
```

But

```text
max_epsilon w_epsilon(P)
  = max_(alpha,beta in A) ||alpha-beta||_1.         (6.3)
```

Indeed, the upper bound is immediate, while equality is obtained by choosing
each coordinate of `epsilon` to be the sign of
`alpha_i-beta_i` for a pair attaining the `l1` diameter.

Now impose common-rotation invariance, which is necessary for an invariant of
five rays whose common sector center is not fixed.  All exponent vectors in
`A` then have one common coordinate sum.  Under (6.2), the difference of any
two support exponents has coordinate sum zero and `l1` norm at most two, so it
is either zero or a root

```text
e_i-e_j.
```

After translating by one support monomial, pairwise compatibility forces all
nonzero roots in the support to share one positive index or one negative
index.  Thus, up to inversion, `P` has the form

```text
c_0 + sum_j c_j u_i/u_j.                            (6.4)
```

If (6.4) vanishes to order at least two on the diagonal, the first logarithmic
derivative in coordinate `j` forces `c_j=0` for every `j`, and the value on the
diagonal then forces `c_0=0`.  Hence a nonzero common-rotation-invariant
Laurent polynomial satisfying universal one-use width has diagonal
multiplicity at most one.

In particular it cannot have the required multiplicity four.

Without common-rotation homogeneity, an `l1`-diameter-two lattice support lies
in one of three elementary clique types: a translated cross-polytope, a unit
square, or the four-point even tetrahedron.  Direct Euler-derivative checks
give multiplicity at most two.  Thus dropping rotation homogeneity still does
not approach order four.

---

## 7. What the no-go theorem covers

The Newton-width argument covers:

1. every fixed Laurent invariant in the five output phases;
2. multihomogeneous polynomials in `T_r` and `conjugate(T_r)` after passing to
   phases;
3. constant-coefficient resultants of those expressions;
4. division by pattern-forced monomial content, which only translates the
   relevant exponent support.

It does not rule out:

1. coefficients that retain the oriented block phases as independent
   variables;
2. non-polynomial primitive-part or Gaussian-gcd operations;
3. a Diophantine approximation theorem that improves the elementary
   denominator-width lower bound;
4. an invariant chosen for the actual weighted set of occurrence words rather
   than uniformly for all words.

Those loopholes require genuinely new arithmetic information.  They are not
ordinary Plucker cancellation in another notation.

---

## 8. Correct replacement target

For an actual five-tuple, group atoms by occurrence word and let `h_epsilon`
be the total log norm in that class.  A fixed Laurent polynomial `P` has
weighted denominator cost

```text
D_h(P) = (1/2) sum_epsilon h_epsilon w_epsilon(P).  (8.1)
```

The precise polynomial target is now:

> Find a nonzero common-rotation-invariant Laurent polynomial `P` whose
> diagonal multiplicity is at least four and for which
>
> ```text
> D_h(P) <= H + O_C(1),                             (8.2)
> ```
>
> or prove that failure of (8.2) forces a structured Gaussian exceptional
> component admitting descent or Walsh isolation.

Unlike the rejected universal condition, (8.2) allows large Newton width on
absent or light patterns.  Its pair-determinant specialization is exactly the
finite linear program (2.3).

The support-only character/Bernoulli model can make all obstructing patterns
heavy, so no theorem based only on rank, return counts, or occurrence support
can guarantee (8.2).  The remaining gain must use the actual Gaussian values
or classify the failure of weighted capacity arithmetically.

---

## 9. Exact squarefree calibration

The obstruction is visible in an honest squarefree six-point configuration,
not only in an abstract pattern model.  Put

```text
Q = 5*13*17*29*41*61*73*89
  = 520699108865.
```

The following six primitive integer points all have norm `Q`:

```text
(707201,143408),
(706841,145172),
(706663,146036),
(706288,147839),
(706151,148492),
(705944,149473).
```

Their five primitive half-gap products are

```text
14429+18i,
1636+i,
2351+3i,
2163+i,
12955+9i.
```

The signed rational-prime occurrence rows, in the prime order displayed in
`Q`, are

```text
(-1,-1,-1,+1, 0, 0,+1,+1),
( 0, 0,+1,-1, 0,+1, 0,-1),
(+1,+1,-1, 0,+1,-1, 0, 0),
(-1,-1,+1,+1, 0, 0,-1, 0),
( 0,+1, 0,-1,-1,+1, 0,+1).
```

Their four differences from the first row have rank four; the minor on
columns `5,13,17,29` has determinant `-6`.  No rational-prime atom has one
fixed oriented state in all five rows.

For this actual pattern set, the pair-determinant capacity is

```text
13/6.
```

One optimal fractional product has

```text
x_02=1/2, x_03=1/2, x_04=1/3,
x_13=1/6, x_14=1/6, x_24=1/2,
```

with every one of the eight atom constraints tight.  Thus even this concrete
rank-four squarefree configuration is far below order four in the pair
determinant class.

An exact dual certificate, in the eight prime-column order above, is

```text
(1/3, 1/12, 7/12, 5/12, 0, 0, 1/12, 2/3).
```

Its total mass is `13/6`; its ten edge loads are

```text
23/12, 1, 1, 1, 4/3, 1, 1, 5/4, 1, 4/3.
```

This proves optimality over the rationals and is checked by the experiment
script.

Numerically, with `H=log Q`, the five product arguments have diameter

```text
s = 0.000813731233179208...,
H-4log(1/s) = -1.4770835217147....
```

Using the largest angle from the positive real ray instead gives

```text
H-4log(1/max arg T_i) = 0.3225009751293....
```

This does not disprove a coefficient-four inequality with an additive
constant.  It shows that the coefficient is already numerically tight at a
small squarefree example and that any proof must track the exact sector
convention and additive term.

---

## 10. Direct five-vector Subspace-Theorem audit

There is a natural alternative to constructing a resultant.  Over
`K=Q(i)`, put

```text
x = [1:u_1:u_2:u_3:u_4].
```

At the infinite place use the five forms

```text
X_0, X_1-X_0, ..., X_4-X_0,
```

and at every finite conductor place use the coordinate forms.  With
normalized absolute values, the product is exactly

```text
product_v product_(j=0)^4
  (|L_(j,v)(x)|_v / |x|_v)
  = H(x)^(-5) product_(i=1)^4 |u_i-1|.            (10.1)
```

If `W` is the active squarefree split-conductor height, then in the rooted
model

```text
log H(x)
  = (1/2) sum_(p : the rooted four-bit word at p is nonzero) log p
  <= W/2.                                         (10.2)
```

Endpoint proximity therefore gives

```text
H(x)^(-5) product_i |u_i-1|
  <= C^4 exp(-W) H(x)^(-5)
  <= C^4 H(x)^(-7).                               (10.3)
```

Thus the numerical exponent is not the obstacle.

For one fixed five-tuple, normalize every active prime by its proportion of
`log H(x)`.  There are only 16 four-bit valuation templates and nine distinct
global linear forms.  The fixed-weight quantitative theorem of
[Evertse--Ferretti](https://arxiv.org/abs/1008.2340) then gives a bounded number
of exceptional hyperplanes with constants independent of the number of
active primes.

The hyperplanes, however, depend on the complete normalized weight system.
That system changes with the selected five points: it records both the
four-bit word of every prime and that prime's normalized log weight.  Applying
the fixed-weight theorem separately therefore produces tuple-dependent
hyperplanes, and one tuple lying in one such hyperplane has no grid
consequence.  The quantitative product-inequality corollary that covers the
varying weights at once reintroduces an explicit exponential dependence on
the number of places.  This is the finite-maximum uniformity gap in this
five-vector normalization.

Aggregating the primes into 16 Gaussian block variables does not remove the
gap.  After deleting the invisible `0000` block, the character map

```text
(G_m)^15 -> (G_m)^4,
t |-> (product_b t_b^(b_1), ..., product_b t_b^(b_4))
```

has rank four and an 11-dimensional generic fiber.  A proper exceptional
component in block space can consequently dominate all four output phases.
For example, setting any one block `t_b=1` is a proper divisor, but the
remaining binary columns still span rank four, so its image is dominant.
Such a component supplies neither a relation among the five products nor a
grid bound.

The direct Subspace-Theorem route would become effective only after proving
one of the same two missing statements already isolated elsewhere:

1. a common exceptional hyperplane family uniform over every 16-template
   weight system; or
2. a classification showing that every block-space exceptional component
   dominating the phase space forces arithmetic descent.

So (10.1)--(10.3) give a strong exact inequality, but the route currently
repackages rather than resolves the uniform exceptional-set bridge.

---

## 11. Verdict

The proposed five-product lemma was worth testing because five is the smallest
arity with four relative phases.  The finite test rejects the proposed proof
mechanism before any large resultant calculation:

```text
four relative phases do not imply four units of one-use arithmetic decay.
```

The next viable target is not a universal five-product resultant.  It is the
weighted capacity/exceptional-structure dichotomy (8.2), or an atom-dependent
Gaussian divisibility theorem that lies outside the fixed Laurent/Plucker
algebra audited here.
