# Multipoint continuation: rational normals and chord recurrences

## Status

This continuation does **not** prove the uniform endpoint bound. It gives
two elementary, fully explicit geometric reductions and a four-point
endpoint family that rules out a tempting recurrence argument. These
calculations do not invoke a Subspace Theorem or an exceptional-set
uniformity assertion.

The starting point is the audited frontier in
`docs/codex_uniform_bound_research.md`: its conductor determinant already
equals the weighted pair-distance/Plotkin calculation, and the surviving
complete hypersimplex profiles cannot be excluded by products of integral
row relations. The calculations below retain the ordering of the points
along the arc, rather than replacing every chord by the full arc length.

## 1. A precise rational-normal bound

Let an arc have radius \(R>0\), length \(L\leq \pi R\), midpoint radial
direction \(e\), and \(M\) lattice points. Let \(n\in\mathbf Z^2\setminus
\{0\}\) be primitive, put \(Q=|n|\), and let \(\varepsilon\) be the angle
between \(n\) and \(e\). Then

\[
 M\leq 2\left(1+
 \left\lfloor Q\left(\frac{L^2}{8R}
                   +L|\sin\varepsilon|\right)\right\rfloor\right).
 \tag{1}
\]

To prove this, write the arc as

\[
z(u)=R(e\cos u+e^\perp\sin u),
\qquad |u|\leq h=\frac{L}{2R}\leq\frac\pi2.
\]

The range of \(n\cdot z(u)\) has width at most

\[
QR\bigl(|\cos\varepsilon|(1-\cos h)
               +2|\sin\varepsilon|\sin h\bigr)
\leq Q\left(\frac{L^2}{8R}+L|\sin\varepsilon|\right).
\]

At lattice points this scalar product is an integer. An interval of width
\(w\) contains at most \(\lfloor w\rfloor+1\) integers, and each line
\(n\cdot z=k\) meets a circle in at most two points. This proves (1).

Consequently, on an endpoint arc \(L\leq C\sqrt R\),

\[
M\leq2\left(1+
 \left\lfloor Q\left(\frac{C^2}{8}
                    +C\sqrt R|\sin\varepsilon|\right)\right\rfloor
\right).
\tag{2}
\]

In particular, bounded-height rational midpoint normals, with
\(Q\leq Q_0\) and
\(\sqrt R|\sin\varepsilon|\leq A\), give the genuine uniform bound

\[
M\leq2\bigl(1+\lfloor Q_0(C^2/8+CA)\rfloor\bigr).
\tag{3}
\]

The restriction on \(Q\) is substantive. Ordinary rational approximation
does not provide bounded \(Q\) for arbitrary arc directions.

### The normal supplied by the actual endpoints

If the first and last lattice points in the selected arc are \(P,Q'\),
write their chord as

\[
d=Q'-P=g v,\qquad g=\gcd(|d_x|,|d_y|),
\qquad v\in\mathbf Z^2\text{ primitive}.
\]

For the shorter containing arc, \(n=v^\perp\), with the appropriate sign,
is exactly its midpoint radial normal. Thus \(\varepsilon=0\), and, if
\(L_0\) is that actual containing arc length,

\[
M\leq2\left(1+\left\lfloor
\frac{|v|L_0^2}{8R}\right\rfloor\right).
\tag{4}
\]

For \(M>2\), the weaker but convenient consequences are

\[
|v|\geq\frac{4(M-2)R}{L_0^2},
\qquad
g\leq\frac{|d|L_0^2}{4(M-2)R}.
\tag{5}
\]

On an endpoint arc these become

\[
|v|\geq\frac{4(M-2)}{C^2},
\qquad
g\leq\frac{C^3\sqrt R}{4(M-2)}.
\tag{6}
\]

Equation (4) applies to every consecutive subcluster. It retains both
the exact local span and the ordinary integer gcd of the extreme chord.
It should not be confused with the Gaussian gcd of the endpoint position
vectors. Those two gcds have different sizes and different roles in the
conductor formulation.

## 2. Exact recurrence of adjacent chords

Order distinct points counterclockwise along an arc of angular width less than \(\pi\), and
let

\[
d_j=z_{j+1}-z_j,
\qquad D_j=\det(d_{j-1},d_j)>0.
\]

The \(D_j\) are positive integers. For three consecutive chord vectors,
linear algebra gives the exact recurrence

\[
d_{j+1}
=\frac{\det(d_{j-1},d_{j+1})}{D_j}\,d_j
-\frac{D_{j+1}}{D_j}\,d_{j-1}.
\tag{7}
\]

If \(\delta_j>0\) is the angular gap and
\(\ell_j=2R\sin(\delta_j/2)\), then

\[
D_j=\ell_{j-1}\ell_j
       \sin\frac{\delta_{j-1}+\delta_j}{2},
\tag{8}
\]

and

\[
\det(d_{j-1},d_{j+1})
=\ell_{j-1}\ell_{j+1}
  \sin\frac{\delta_{j-1}+2\delta_j+\delta_{j+1}}2.
\tag{9}
\]

Thus a short total angular span bounds the sine factors but does not
bound the ratios of the adjacent chord lengths in (7). Even strengthening
integrality to unimodular consecutive primitive directions does not bound
the recurrence coefficients, as the next family proves.

## 3. Four endpoint points with an unbounded unimodular recurrence

Let \(t\geq4\) be an even integer, and put

\[
B=\frac{t^2}{2},\qquad R=\sqrt{B^2+1}.
\]

The four lattice points

\[
z_0=(t,B-1),\qquad z_1=(1,B),\qquad
z_2=(-1,B),\qquad z_3=(-t,B-1)
\tag{10}
\]

have common squared norm \(B^2+1\), since

\[
t^2+(B-1)^2=2B+(B-1)^2=B^2+1.
\]

They occur in the displayed order on the upper arc. Its exact length is

\[
L_t=2R\arctan\frac{t}{B-1}.
\tag{11}
\]

For \(B\geq8\),

\[
\frac{L_t}{\sqrt R}
\leq\frac{2t\sqrt R}{B-1}
\leq\frac{2\sqrt{2B(B+1)}}{B-1}
\leq4.
\tag{12}
\]

Here \(R\leq B+1\), and the last inequality follows from
\(B^2-5B+2\geq0\) for \(B\geq8\). Moreover,

\[
\frac{L_t}{\sqrt R}\longrightarrow2\sqrt2.
\tag{13}
\]

The primitive consecutive chord vectors are

\[
p_0=(1-t,1),\qquad p_1=(-1,0),\qquad p_2=(1-t,-1).
\tag{14}
\]

They satisfy

\[
\det(p_0,p_1)=\det(p_1,p_2)=1,
\qquad
p_2=2(t-1)p_1-p_0.
\tag{15}
\]

Consequently the positive recurrence coefficient \(2(t-1)\) tends to
infinity although all four points lie on a fixed-\(C\) endpoint arc and
both consecutive primitive determinants equal one. The actual chord
determinants are both equal to two. The middle chord has length two,
whereas each outer chord has length asymptotic to \(t\).

This is an explicit counterexample to each proposed intermediate claim
that endpoint arcs have uniformly comparable consecutive gaps, uniformly
bounded adjacent-chord recurrence coefficients, or uniformly bounded
coefficients after imposing consecutive primitive determinant one. It is
not a counterexample to the desired point-count bound: it has four points.

The family also explains the geometry of (1): the horizontal normal rows
\(y=B\) and \(y=B-1\) each contain two points. Its primitive radial normal
has height one, so the rational-normal method correctly supplies a
uniform bound here despite the unbounded recurrence coefficient.

## 4. What this route still needs

For a proof based on (1), one would need a new argument that sufficiently
many endpoint points force a primitive integer normal of uniformly
bounded height whose angular error is \(O_C(R^{-1/2})\), or a comparably
strong replacement. Such a statement is not obtained from Minkowski's
theorem or from (4): the primitive normal supplied by the extreme chord
can have height growing with the radius.

For a proof based on (7), a bounded-coefficient assertion is false even
under the strongest elementary determinant normalization. Any successful
recurrence method must use a genuinely many-point restriction that is
absent from the four-point family (10), rather than bounded coefficients
or bounded gap ratios alone.

Neither missing assertion has been proved in this continuation. The
uniform endpoint theorem remains open in the repository.

## Appendix A. Independent audit of the dependent-pair lemma

The argument in docs/uniformity_continuation.md, including its explicit
incidence bound, checks independently. In particular:

1. A multiplicative relation modulo \(\mu_4\) gives a rank-one relation
   between the two nonzero integer valuation rows. They can be written
   \(\mathbf d=a\mathbf c\), \(\mathbf d'=b\mathbf c\), with
   \(\gcd(a,b)=1\). Neither the coordinates of \(\mathbf c\) nor the two
   integers \(a,b\) need be positive. Orienting each Gaussian prime
   according to the sign of \(c_j\) gives one conjugate-coprime \(W\).
   The bound \(\max(|a|,|b|)|c_j|\leq e_j\) then proves
   \(|W|^{\max(|a|,|b|)}\leq Q\leq R\), without requiring a common
   orientation for other pairs in the cluster.
2. Bezout coefficients recover a ratio \(\eta W/\bar W\), with an
   arbitrary Gaussian unit \(\eta\). Its numerator
   \(\eta W-\bar W\) cannot vanish: if it did, conjugate-coprimality
   would force \(W\) to be a unit. Thus the lower bound \(1/|W|\)
   remains valid with units. The telescoping upper bound and the stated
   elementary estimate for \((\log R)R^{-1/12}\) give the claimed
   explicit cutoff for \(\max(|a|,|b|)\geq3\).
3. After the cutoff, one coefficient has absolute value one. The unit
   in \(v=\omega u^b\) is forced to be one by
   \(|\omega-1|\leq3C/\sqrt R<\sqrt2\). The arguments of all non-root
   ratios have the same sign. This excludes negative \(b\), and
   distinctness excludes \(b=1\), leaving precisely the square relations.
   For literally positive arguments one chooses the counterclockwise
   endpoint, or conjugates the configuration; either endpoint gives
   the same-sign property that the proof needs.
4. For \(v=u^2\), the pair-dependent conductor orientation gives
   \(D=G_{\mathbf c}/W^2\), with \(|D|\leq C^2\). Although the
   orientation varies, its modulus \(Q=|G_{\mathbf c}|\) is fixed
   throughout the original circle. Therefore fixing \(D\) fixes
   \(|W|^2=Q/|D|\). The bound \(|\zeta W-\bar W|\leq C\) restricts
   \(W\) to finitely many integer lines; each line meets this fixed
   circle in at most two points. Counting the four units and both
   ordered square relations gives exactly the valid coarse bound
   \[
   16(2\lfloor C\rfloor+1)(2\lfloor C^2\rfloor+1)^2.
   \]

No unbounded unit factor, conductor-support count, or moving-radius
parameter has been omitted in that count. The argument is unconditional,
but gives no bound for pairwise multiplicatively independent ratios.

## Appendix B. Exact defect introduced by bounded row restriction

Let \(M\) be even, and consider a weighted family of cuts of an \(M\)-row
set, each cut containing exactly \(M/2\) rows. Write \(w_a>0\) for the
cut weights and \(W=\sum_a w_a\). Thus every cut has zero defect in the
full \(M\)-row system.

Choose uniformly a subset \(T\) of \(2s\) rows, where \(2s\leq M\), and
write \(r_a=|T\cap S_a|\). The defect of the restriction is

\[
D_T=\sum_a w_a(r_a-s)^2.
\tag{B1}
\]

For each fixed cut, \(r_a\) is hypergeometric: one samples \(2s\) objects
without replacement from \(M\), of which \(M/2\) are marked. Consequently

\[
\mathbb E r_a=s,\qquad
\mathbb E(r_a-s)^2
=\frac{2s}{4}\frac{M-2s}{M-1}
=\frac{s(M-2s)}{2(M-1)}.
\tag{B2}
\]

For completeness, the variance follows directly by writing
\(r_a=\sum_{j=1}^{2s}X_j\): each indicator has variance \(1/4\), and
two distinct sampled indicators have covariance \(-1/[4(M-1)]\).
Linearity of expectation now gives, for arbitrary cut weights,

\[
\boxed{\mathbb E D_T
=\frac{s(M-2s)}{2(M-1)}W.}
\tag{B3}
\]

Compare this with the even-\(2s\) endpoint determinant allowance \(sW/2\),
using the same valuation-layer weight convention as Section 5 of
docs/codex_uniform_bound_research.md. The exact difference is

\[
\boxed{\frac{sW}{2}-\mathbb E D_T
=\frac{s(2s-1)}{2(M-1)}W.}
\tag{B4}
\]

For fixed \(s\), this margin tends to zero relative to \(W\) as
\(M\to\infty\). It supplies no extraction whose restricted defect
vanishes.

This is sharp as a statement about weighted balanced cuts. Take the
complete family of \(M/2\)-subsets with equal weights, or one
representative from each complementary pair with equal weights. The
square \((r_a-s)^2\) is unchanged by complementing a cut. Permutation
symmetry therefore makes \(D_T\) the same for every \(2s\)-subset \(T\),
and this common value equals (B3). If \(M>2s\), it is positive. Thus
there is **no** \(2s\)-row restriction with every cut balanced, even
though the original system has zero defect at every cut.

For the eight-point rational-ray argument, \(s=4\), so

\[
D_T=\frac{2(M-8)}{M-1}W
\]

in this equal-weight example, and its gap below the allowance \(2W\) is
\(14W/(M-1)\). An exact globally balanced system of arbitrarily many rows
therefore does not automatically contain the exact eight-row balanced
system required by that argument.

Conditioning on a common Gaussian-unit class cannot universally improve
this conclusion. Assign the same unit to every row of the preceding
complete cut family. Unit conditioning then leaves every row and every
defect unchanged. Pigeonholing the at most four unit classes is a valid
way to remove row units, but does not preserve, or create, exact balance
in a bounded row restriction.

These are obstructions to a proposed extraction from valuation balance
and unit labels alone. They do not construct an endpoint Gaussian
cluster realizing the complete cut system with the required simultaneous
small phases.

## Appendix C. The weaker product-unit condition and its selection cost

The unit-general phase identity behind the eight-row rational-ray argument
can be checked directly. Suppose the seven non-root products satisfy
\(\prod_i A_i=G^4\), put \(B_i=G/A_i\), and write

\[
u_i=\zeta_i A_i/\bar A_i,\qquad
\zeta_i=i^{e_i},\qquad S=\sum_i e_i.
\]

If the positive lifted argument of \(u_i\) is \(2\theta_i\), with
\(0\leq\theta_i\leq\Delta/2\), then
\(\arg A_i\equiv\theta_i-e_i\pi/4\pmod\pi\). The product identity implies,
for one common integer \(n\),

\[
\arg B_i\equiv
\frac14\sum_j\theta_j-\theta_i
+\frac{e_i\pi}{4}-\frac{S\pi}{16}+\frac{n\pi}{4}
\pmod\pi.
\tag{C1}
\]

Thus arbitrary row units can introduce rays at multiples of \(\pi/16\),
not merely \(\pi/8\). But the condition \(\prod_i\zeta_i=1\), equivalently
\(S\equiv0\pmod4\), makes every limiting ray a multiple of \(\pi/4\).
This product condition suffices for the rational Gaussian-ray argument;
individual equality of the seven units to one is unnecessary.

There is an elementary exact selection bound for this weaker condition:

> Any eleven Gaussian-unit labels contain eight whose product equals one.

To prove it, represent the labels by exponents modulo four. Among any
seven exponents one can make three disjoint pairs whose two entries
have the same parity: if the even and odd counts are \(a,b\), then
\(\lfloor a/2\rfloor+\lfloor b/2\rfloor=3\). Each pair sum is zero or
two modulo four. Two of the three sums agree, so their four entries
sum to zero modulo four. From eleven labels first remove such a
four-element zero-sum subset; seven labels remain and contain another.
Their union gives the desired eight.

The number eleven is sharp for unit labels alone. Seven labels equal
to one and three equal to \(i\) form a ten-element sequence in which
every eight-element subset contains one, two, or three labels equal
to \(i\), and hence has product different from one.

For eight absolute point units \(\epsilon_0,\ldots,\epsilon_7\), the
product of the seven relative units at any chosen root is

\[
\prod_{i\ne0}\frac{\epsilon_i}{\epsilon_0}
=\left(\prod_{i=0}^7\epsilon_i\right)\epsilon_0^{-8}
=\prod_{i=0}^7\epsilon_i.
\tag{C2}
\]

Therefore the product-one condition is independent of which of the
eight points is chosen as root. This removes the unit obstruction at
a bounded selection cost, but leaves the exact-balance obstruction in
Appendix B intact.

## Appendix D. A conditional finite local-test criterion

This appendix formulates a sufficient new arithmetic statement. It does
**not** prove that statement or assert that a suitable test function is
known. Symmetry improves the elementary sampling error from
\(O(M^{-1/2})\) to \(O(M^{-1})\).

### D.1 Normalization and weighted near-balance

Start with an even number \(M\geq8\) of points on an arc of length at
most \(\frac12\sqrt R\). Divide every point by their common Gaussian
gcd. This preserves cardinality and makes the endpoint constant no
larger: division by \(d\) changes it to \(1/(2\sqrt{|d|})\). Inert and
ramified factors disappear, and at each remaining split prime the
minimum and maximum valuations across the full cluster are \(0,e_p\).
With the resulting radius still called \(R\), put

\[
W=\sum_p e_p\log p=\log(R^2).
\]

Index the prime-exponent layers by \(a=(p,t)\), give each layer weight
\(w_a=\log p\), and let \(r_a\) be the number of the \(M\) points whose
chosen \(\pi\)-valuation is at least \(t\). The even determinant identity
gives

\[
\sum_a w_a(r_a-M/2)^2
\leq \frac M4W-M(M-1)\log2.
\tag{D1}
\]

Dropping the negative term and writing \(p_a=r_a/M\) yields

\[
\sum_a w_a(p_a-\tfrac12)^2\leq\frac{W}{4M}.
\tag{D2}
\]

The nonnegativity of the left side of (D1) also gives the useful bound

\[
W\geq4(M-1)\log2.
\tag{D3}
\]

These weights are prime-layer weights \(\log p\), not
\(\log|\pi|=\frac12\log p\). One can halve every weight consistently,
but then \(W=\log R\), and the additive constants in all local
inequalities must be halved as well.

### D.2 Sampling eight rows

Let \(f:\{0,\ldots,8\}\to\mathbb R\) satisfy \(f(r)=f(8-r)\), and define

\[
H=\max f-\min f,\qquad
\mu_f=2^{-8}\sum_{r=0}^8\binom8r f(r).
\]

Choose a uniform eight-element row subset \(T\), and set
\(r_{a,T}=|T\cap S_a|\). Then

\[
\boxed{\left|
\mathbb E_T\sum_a w_a f(r_{a,T})-\mu_fW
\right|\leq\frac{42H}{M}W.}
\tag{D4}
\]

In particular the right side is at most
\(84\|f\|_\infty W/M\).

Here is a direct proof. Couple sampling eight population indices with
replacement to sampling without replacement. The chance of a repeated
index in the former is at most \(\binom82/M=28/M\). For a fixed cut of
density \(p\), it follows that

\[
\left|\mathbb E f(\operatorname{Hyp}(M,pM,8))
-\mathbb E f(\operatorname{Bin}(8,p))\right|
\leq\frac{28H}{M}.
\tag{D5}
\]

The coupling can be realized by keeping the independent sample if its
indices are distinct, and otherwise replacing it by an independent
uniform ordered sample without replacement. Its latter marginal is
uniform without replacement in either branch.

Write \(b_f(p)=\mathbb E f(\operatorname{Bin}(8,p))\). Symmetry of \(f\)
gives \(b_f(p)=b_f(1-p)\), so \(b'_f(1/2)=0\). Differentiating its
Bernstein-polynomial expression gives

\[
b''_f(p)
=56\,\mathbb E\bigl[
f(Y+2)-2f(Y+1)+f(Y)\bigr],
\qquad Y\sim\operatorname{Bin}(6,p).
\]

The second difference has modulus at most \(2H\). Taylor's theorem
therefore gives

\[
|b_f(p)-\mu_f|\leq56H(p-\tfrac12)^2.
\tag{D6}
\]

Weighting (D5) and (D6), and using (D2), gives respectively
\(28HW/M\) and \(14HW/M\), proving (D4). No independence among
different conductor layers was used.

Without symmetry, the Bernoulli coupling gives the weaker bound
\[
\left|\mathbb E_T\sum_a w_a f(r_{a,T})-\mu_fW\right|
\leq H W\left(\frac{28}{M}+\frac4{\sqrt M}\right),
\]
using weighted Cauchy--Schwarz and (D2). Thus the proposed
\(O(\|f\|_\infty/\sqrt M)\) criterion also holds, but symmetry gives
the stronger explicit estimate above.

### D.3 The missing arithmetic inequality

Suppose one proves the following statement for one such fixed \(f\):
there exist \(a\in\mathbb R\), \(B\geq0\), and \(W_0\geq0\), independent
of the radius, conductor support, and selected points, such that every
eight-point endpoint subcluster, with its **inherited** prime layers,
satisfies

\[
\sum_a w_a f(r_{a,T})\leq aW+B
\quad\text{whenever }W\geq W_0,
\qquad
\epsilon:=\mu_f-a>0.
\tag{D7}
\]

The quantifiers in (D7) are the new, unproved input. Averaging it over
all eight-element subsets and applying (D4) gives

\[
\epsilon\leq\frac{42H}{M}+\frac BW.
\tag{D8}
\]

By (D3), for \(M\geq8\),

\[
\epsilon
\leq\frac{42H}{M}+\frac{B}{4(M-1)\log2}
\leq\frac1M\left(42H+\frac{B}{2\log2}\right).
\]

Consequently

\[
\boxed{M\leq
\frac{42H+B/(2\log2)}{\epsilon}.}
\tag{D9}
\]

If \(W<W_0\), (D3) instead gives
\(M\leq1+W_0/(4\log2)\). Removing one point when necessary handles odd
cardinalities. Thus (D7), if proved, supplies a uniform endpoint bound,
with an explicit constant and no uncontrolled bounded-height remainder.
For a general prescribed arc constant \(C\), partition into at most
\(\lceil2C\rceil\) arcs with constant \(1/2\).

### D.4 Inherited layers, signs, and units

Several normalization conditions in (D7) are necessary.

- A chosen eight-point subset may acquire new common Gaussian factors.
  The expectation in (D4) still uses the conductor of the **full**
  \(M\)-point cluster. Its layers with \(r_{a,T}=0\) or \(8\) must not
  silently be omitted, nor may its \(W\) be replaced by a smaller
  subset conductor.
- A sufficient way to lift a statement proved after normalizing each
  eight-point subset is \(f(0)=f(8)\leq a\), provided the same additive
  constant works for every endpoint constant at most \(1/2\). Indeed,
  the discarded constant layers then contribute at most \(a\) times
  their weight. When this endpoint-value condition fails, a separate
  quantitative lifting argument is required; the improved arc constant
  after common-factor division can matter.
- Reversing the orientation of a split prime complements the cut and
  replaces its sampled size by \(8-r\). Symmetry of \(f\) makes the
  statistic orientation-independent. Threshold layers handle arbitrary
  prime exponents; squarefreeness is not assumed. Signed valuation
  differences relative to a chosen root are unnecessary in (D1)--(D7).
- Units do not affect cuts or weights. If the new arithmetic proof
  requires a common Gaussian-unit class, one can first retain a largest
  such class, losing at most a factor four in cardinality, and then
  normalize and apply the argument to that cluster. The distribution
  in (D4) is the uniform sample from this retained cluster. One cannot
  instead condition individual eight-row samples on arbitrary unit
  labels and retain (D4) without proving a new sampling estimate.

### D.5 The existing quadratic test has no strict margin

For the known determinant test

\[
f(r)=(r-4)^2,
\]

the binomial mean is exactly

\[
\mu_f=\operatorname{Var}(\operatorname{Bin}(8,\tfrac12))=2.
\]

The eight-point Ramana inequality, with the original inherited
conductor, is

\[
\sum_a w_a(r_{a,T}-4)^2
\leq2W+56\log(1/2).
\tag{D10}
\]

Its coefficient is \(a=2=\mu_f\), so \(\epsilon=0\). Its favorable
additive constant does not create a positive binomial-mean gap.
The exact balanced eight-row phase exclusion also does not by itself
prove (D7): excluding the single profile in which every cut has size
four supplies no quantitative weighted inequality with
\(\mu_f>a\).

Accordingly (D7) is a precise finite arithmetic target for future work,
not a completed step toward the uniform theorem.

One especially concrete candidate is, for \(0<\delta\leq1\),

\[
f_\delta(r)=(r-4)^2+\delta\,\mathbf1_{\{r=4\}}.
\]

Its oscillation and binomial mean are

\[
H=16-\delta,\qquad
\mu_{f_\delta}=2+\frac{35}{128}\delta.
\]

Thus it would suffice to prove the following **currently unproved**
inequality, with one \(B\geq0\) independent of the conductor and points:

\[
\boxed{
D_T+\delta W_{T,4}\leq2W+B,\qquad
W_{T,4}=\sum_{a:\,r_{a,T}=4}w_a,
}
\tag{D11}
\]

for every eight-point endpoint subcluster, using its inherited layers.
Here \(D_T=\sum_a w_a(r_{a,T}-4)^2\). The common-unit case alone
would suffice by the initial factor-four restriction described above.
If (D11) were proved, (D9) would give, for an even normalized cluster
in its valid height range,

\[
M\leq
\frac{128}{35\delta}
\left(42(16-\delta)+\frac{B}{2\log2}\right).
\tag{D12}
\]

The added term asks for a fixed positive logarithmic contribution from
balanced layers beyond the ordinary determinant defect. The present
exact-balance phase exclusion does not establish this weighted
contribution, and (D11) is left as the explicit missing lemma.
