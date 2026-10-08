# The direct monomial-vector Subspace test pays a width height

The scalar packing obstruction for complete balanced cuts also survives a
specific simultaneous approximation upgrade. Put finitely many integral
phase characters into one projective point, approximate its coordinates
by their fixed torsion targets, and use only their first-order errors.
The exact projective height is a primewise **range**, rather than the
largest individual height or an uncharged common denominator. On the
complete balanced profile this range is already too expensive for the
Subspace theorem inequality, for every number of distinct characters.

This is an obstruction to that stated test. It does not exclude Subspace
theorem arguments using higher vanishing orders, additional finite-place
information, or special additive identities. No endpoint configuration is
constructed and no uniform lattice-point theorem is proved.

## 1. Exact projective height with all common primes retained

Use the fixed sign model of
[the coupled-product note](coupled_block_roth_packing.md). Its nonzero
Gaussian blocks are coprime to their conjugates, and therefore have
factorizations

```text
kappa_j = eta_j product_p pi_p^((e_jp)_+)
                         bar(pi_p)^((-e_jp)_+),
eta_j in {1,-1,i,-i},   e_jp in Z,   Norm(pi_p)=p.
```

All primes p here split in the Gaussian integers. For an integral vector v
define the rational Gaussian direction and its signed valuations by

```text
chi_v = product_j (kappa_j/bar(kappa_j))^v_j in Q(i),
E_p(v)=sum_j v_j e_jp.
```

The complex absolute value of every chi_v is one. Fix vectors
v_0,...,v_q and consider X=[chi_(v_0):...:chi_(v_q)]. Multiplication
of all coordinates by chi_(-v_0) lets us assume v_0=0. For each prime put

```text
m_p=min_(0<=t<=q) E_p(v_t),
M_p=max_(0<=t<=q) E_p(v_t),
rho_p=M_p-m_p.
```

In this normalization m_p<=0<=M_p. The Gaussian integer

```text
D=product_p pi_p^(-m_p) bar(pi_p)^M_p
```

clears every denominator. Apart from a Gaussian unit, coordinate t of
the integral vector D X is exactly

```text
product_p pi_p^(E_p(v_t)-m_p) bar(pi_p)^(M_p-E_p(v_t)).  (1)
```

At each Gaussian prime at least one coordinate has valuation zero. Thus
the coordinates in (1) have Gaussian gcd one. Every coordinate has
complex absolute value |D|. With the absolute multiplicative projective
height H on P^q(Q(i)), this proves the exact identity

```text
log H(X) = log |D|
         = (1/2) sum_p rho_p log p.                    (2)
```

For clarity, the normalized absolute value at the unique complex place
of Q(i) is the ordinary complex modulus; at the Gaussian prime pi_p it
is p^(-ord_(pi_p)/2). This gives precisely the factor one half in (2).
The scalar logarithmic norm h(v) in the coupled-product note is
**twice** log H([1:chi_v]).

Formula (2) includes arbitrary shared primes, opposite orientations,
nested powers, and complete character cancellation. It is invariant
under translating every v_t by one common vector. No disjointness or
nonunit assertion is used in its proof.

If the rational prime supports of the blocks are disjoint, (2) becomes

```text
log H(X)=(1/2) sum_j r_j log Norm(kappa_j),
r_j=max_t v_tj-min_t v_tj.                             (3)
```

Changing the orientation of a disjoint block does not alter r_j. Thus
the denominator cost is the coordinate range of the entire exponent
set. Replacing it by max_t log H([1:chi_(v_t)]) can lose a positive
power of the height.

## 2. What the direct Subspace test would require

Write V for the rational span of the full row differences. Suppose
v_t-v_0 belongs to V intersect Z^k for every t. The fixed rational
certificates in the coupled-product note eliminate the common arc
center and every moving parameter in ker(V). After passing to one of
finitely many branches, they give fixed algebraic numbers alpha_t of
absolute value one with

```text
|chi_(v_t-v_0)-alpha_t| <= K Delta,
Delta <= C R^(-1/2).                                  (4)
```

Here the alpha_t are roots of unity, and K depends only on the fixed
certificates. Independent Gaussian row units only change the finite
branch set. Neither a moving center nor a moving algebraic target is
being frozen in (4).

Let n=q+1. For the primitive Gaussian vector x=D(1,chi_(v_1),...,chi_(v_q))
use the fixed complex linear forms

```text
L_0(x)=x_0,      L_t(x)=x_t-alpha_t x_0  (1<=t<=q).
```

All |x_t| equal H=|D|, and (4) gives

```text
product_(t=0)^q |L_t(x)| <= H^(q+1) (K Delta)^q.       (5)
```

The direct first-order Subspace test therefore needs a positive exponent
gap of the form

```text
Delta^q < H^(-(q+1)-epsilon).                         (6)
```

One can verify this normalization using only the classical rational
Subspace theorem, without assuming the target coefficients lie in Q(i).
Take real and imaginary parts of each L_t. These are 2n independent real
linear forms with real algebraic coefficients in the 2n integral real
coordinates of x. Their product is at most
H^(2n)(K Delta)^(2q). The Euclidean norm of that integral vector is
sqrt(n)H. Thus (6), with its positive gap, makes that product a negative
power of the Euclidean norm. The theorem places these vectors in finitely
many proper rational subspaces; it does not by itself make the set finite.

The primary source used for this theorem is J.-H. Evertse and H. P.
Schlickewei, *A quantitative version of the Absolute Subspace Theorem*,
J. reine angew. Math. 548 (2002), 21--127, Introduction, equation (1.1)
and the preceding theorem statement
([authors' preprint, printed page 1](https://pub.math.leidenuniv.nl/~evertsejh/00-abssub.pdf)).
The normalized height convention also appears explicitly in J.-H.
Evertse and R. G. Ferretti, *A further improvement of the Quantitative
Subspace Theorem*, Introduction, equations (1.1)--(1.2)
([authors' preprint, printed pages 1--2](https://math.leidenuniv.nl/reports/files/2010-17.pdf)).

At an asymptotic height ratio log H/log R -> rho, inserting only the
endpoint error bound from (4) into (6) requires

```text
rho < q/(2(q+1)).                                     (7)
```

Failure of (7) means the supplied upper bounds do not establish the
Subspace inequality. It does not mean a particular actual error product
cannot be smaller. Extra smallness and the exceptional subspaces are
separate questions.

## 3. A range obstruction on every complete balanced profile

Let M>=4 be even, B=binom(M,M/2), and index the block coordinates by all
balanced cuts T. Put S_i(T)=2*1_(i in T)-1. Let

```text
Lambda=span_Q{S_i-S_0} intersect Z^B,
a_M=(M-2)/(2(M-1)).
```

The exact integral lattice theorem in
[the balanced phase minimum](coupled_balanced_phase_l1_minimum.md) says

```text
(1/B) sum_T |v(T)| >= a_M     for 0!=v in Lambda.       (8)
```

Take any q+1 pairwise distinct vectors v_0,...,v_q in an affine translate
of Lambda, with q>=1, and define their average range

```text
rho(v_0,...,v_q)
  =(1/B) sum_T [max_t v_t(T)-min_t v_t(T)].             (9)
```

**Range theorem.** Every such finite exponent set satisfies

```text
q=1:   rho >= a_M > 1/4;
q>=2:  rho >= (3/2)a_M >= 1/2 > q/(2(q+1)).           (10)
```

For q=1 the range is |v_1-v_0|, so (8) applies. For q>=2 retain any
three distinct vectors and translate the first to zero. Pointwise,

```text
max(0,a,b)-min(0,a,b)
       = (|a|+|b|+|a-b|)/2.                           (11)
```

Each of v_1-v_0, v_2-v_0, and v_1-v_2 is a nonzero member of Lambda.
Averaging (11), using (8), and observing that additional coordinates
cannot decrease a range proves (10). Since a_M>=1/3, the last
inequality is strict for every finite q. No bound on the sizes of the
integer exponents, no rank assumption on their span, and no finite
enumeration is needed.

The q>=2 lower bound is attained at M=6: set

```text
v_0=0,
v_1=(S_1+S_2)/2,
v_2=-(S_3+S_4)/2.
```

The remaining difference is -(S_5+S_6)/2. All three nonzero
differences have average absolute value 2/5, giving rho=3/5.
Attainment of the displayed lower bound is not claimed for other M.

Now consider the disjoint-support height benchmark

```text
log Norm(kappa_T)=w+o(w) uniformly in T,
log |d|=o(w),                w -> infinity.            (12)
```

This is equality of leading logarithmic weights; exact equality of
distinct coprime block norms is neither needed nor possible for nonunits.
For every fixed exponent set, (3) and the actual radius formula give

```text
log R=(B/2)w+o(w),
log H(X)=(B/2)rho(v_0,...,v_q)w+o(w).
```

Consequently (10) contradicts the necessary strict numerical condition
(7), for every direct monomial tuple of distinct characters. The gap
persists in a sufficiently small fixed relative weight window, depending
on M and q. The benchmark is a statement about this height functional;
no actual short-arc sequence satisfying (12) is asserted to exist.

Duplicating characters cannot increase the usable dimension: if two
exponent vectors coincide, their projective coordinates coincide
identically, giving a fixed linear relation. With optional fixed torsion
coordinate multipliers they are still proportional by a fixed algebraic
constant. Delete repetitions before the test. Accidental equalities at
particular values are additional exceptional relations, not extra
independent approximation coordinates.

## 4. Scope and a quantitative next target

The theorem covers arbitrary finite projective monomial lists whose
pairwise exponent differences belong to the fixed phase lattice. Thus it
allows a common nonconstant monomial factor and removes that factor
projectively. It covers all kernel dimensions for which the complete
balanced profile is defined. It is stronger than checking one previously
selected basis, but deliberately concerns only the first-order estimates
(4)--(5).

For shared supports, the arithmetic quantity that must improve is known
exactly: it is the primewise range sum in (2). One must prove, for actual
arc data and a useful fixed tuple, an estimate

```text
(1/2) sum_p [max_t E_p(v_t)-min_t E_p(v_t)] log p
       <= (q/(2(q+1))-epsilon) log R.                 (13)
```

Alternatively, first-order errors may be replaced by a proved extra
vanishing order. Either improvement must retain the actual primes,
coefficients, and fixed targets; establishing the Subspace inequality
would still leave its exceptional subspaces to control. Exterior powers
with nonmonomial minors and higher-order auxiliary forms are outside the
range theorem, because their smallness and their heights require further
information. Merely increasing the number of monomial coordinates does
not supply it.

Higher-order certificates have already been investigated in other
coordinates: the eight-point harmonic Gale pullbacks satisfy the
all-degree generic cut-order ceiling `sum_S nu_S(F)<=53 deg(F)` in
[the Newton-support note](gale_torus_newton_cut_orders.md), with a separate
[phase-coordinate and multiplicity audit](gale_torus_newton_independent_audit.md).
Those are additional scoped obstructions, so passing beyond the present
first-order test is not by itself an untried route or a promised gain.

Run [the standard-library checker](check_monomial_vector_subspace_height.py)
with `python3 docs/check_monomial_vector_subspace_height.py`. It verifies
primitive Gaussian denominator clearing, exact common-prime height
identities, translation invariance, the scalar range identity, and finite
balanced-cut examples. The all-exponent minimum (8) and the Subspace
theorem remain the explicitly identified prose inputs.
