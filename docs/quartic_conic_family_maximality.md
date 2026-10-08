# Maximality of the entire nondegenerate rational-conic quartic family

This note extends the single-specialization result in
[quartic_extension_maximality.md](quartic_extension_maximality.md) to both
real branches of the conic construction in
[three_row_quartic_construction.md](three_row_quartic_construction.md).

**Theorem.** Let \(F_1,F_2,F_3\in\mathbb R[X]\) be any nondegenerate member
of the conic family below. If a polynomial \(F\in\mathbb R[X]\) satisfies

\[
 F-F_j\mid FF_j+1\qquad(j=1,2,3)
\tag{1}
\]

whenever \(F\ne F_j\), then \(F\) is one of \(F_1,F_2,F_3\). There is
no fourth compatible polynomial, of any degree. In particular this holds
for every nondegenerate rational parameter, with rational coefficients.

This is a theorem about this polynomial family and condition (1). It does
not exclude a fourth lattice direction whose primitive numerator has a
nonconstant imaginary part, a rational-function numerator, or a different
factor profile. It does not prove the uniform-in-radius count bound.

The exact algebraic certificate is
[check_quartic_conic_maximality.py](check_quartic_conic_maximality.py). It
uses only rational polynomial arithmetic from Python's standard library;
neither parameter sampling nor floating point arithmetic enters the proof.

## 1. Parameters and the nondegenerate domain

Write

\[
 q=v^2+8,\qquad s=\frac{6v}{q},\qquad
 r=\frac{3(8-v^2)}q,\qquad r^2+8s^2=9.
\tag{2}
\]

The four initial roots are

\[
 a=i,\quad b=3s+2i,\quad c=4s-i,\quad d=s-2i.
\tag{3}
\]

The remaining three roots on the minus branch are

\[
\begin{aligned}
e_-&=\frac{8(v^2-4)}{vq}-\frac{3v^2}{q}i,\\
f_-&=\frac{v(v^2+32)}{2q}-\frac{24}{q}i,\\
g_-&=\frac{16(v^2+2)}{vq}+\frac{3v^2}{q}i.
\end{aligned}
\tag{4}
\]

On the plus branch they are

\[
\begin{aligned}
e_+&=\frac{v(16-v^2)}{2q}-\frac{24}{q}i,\\
f_+&=\frac{16(v^2+2)}{vq}-\frac{3v^2}{q}i,\\
g_+&=\frac{v(v^2+32)}{2q}+\frac{24}{q}i.
\end{aligned}
\tag{5}
\]

Here the roots (a,b,c,d,e,f,g) correspond respectively to blocks
\(A,B,C,D,E,F,G\). Define \(H_1=(X-a)(X-b)(X-c)(X-d)\), and put

\[
\begin{aligned}
H_2&=H_1-\lambda N((X-a)(X-b)),\\
H_3&=H_1-\mu N((X-a)(X-c)),\\
\tau&=12s(1+s^2),\qquad F_j=\frac{\operatorname{Re}H_j}{\tau}.
\end{aligned}
\tag{6}
\]

With \(w=v^2\), the real coefficients are

\[
\begin{array}{c|c|c}
&\lambda&\mu\\ \hline
-&\displaystyle\frac{(w-16)(w+4)}{(w-4)(w+16)}
 &\displaystyle\frac{(5w+16)(w-16)(w+4)}
 {8(w-4)(w+2)(w+16)}\\[4pt]
+&\displaystyle\frac{(w-4)(w+16)}{(w-16)(w+4)}
 &\displaystyle\frac{(w-4)(w+16)(w+20)}
 {(w-16)(w+4)(w+32)}.
\end{array}
\tag{7}
\]

All three \(H_j\) have constant imaginary part \(\tau\); their root
sets are exactly

\[
 ABCD,\qquad ABEF,\qquad ACEG.
\tag{8}
\]

The term **nondegenerate** means that the construction gives three distinct
finite quartics, \(\tau\ne0\), and that its seven roots are distinct and
none is a conjugate of any of the seven roots. This is the same separation
needed for the seven independent oriented cut factors.

The following real parameter values cannot be nondegenerate:

* \(v=0\): \(\tau=0\), and some displayed roots are undefined.
* \(v^2=4\) or \(16\): one branch has a pole in \(\lambda\), while the
  other has \(\lambda=0\) and hence \(H_2=H_1\).
* \(v^2=8\): on either branch, \(f=\bar g\). In fact
  \(f=5v/4-3i/2\) and \(g=5v/4+3i/2\).

These are all the real zeros that will arise in the exceptional rank
calculation below. We do not use an assertion that other root-separation
conditions follow automatically from avoiding this list.

## 2. Why sixteen divisor directions cover every polynomial degree

Modulo \(F-F_j\), condition (1) is equivalent to

\[
 F-F_j\mid F_j^2+1.
\tag{9}
\]

The real polynomial \(F_j^2+1\) is squarefree of degree eight and is a
nonzero constant times the product of the four real irreducible norm
quadratics associated with the four roots in row \(j\). Consequently its
monic real divisors are precisely

\[
 P_U(X)=\prod_{H\in U}(X-h)(X-\bar h),
 \qquad U\subseteq\text{row }j.
\tag{10}
\]

There are sixteen such directions, of degrees \(0,2,4,6,8\). Every
possible \(F\) therefore lies on a line

\[
 F_j+\mathbb R P_U
\tag{11}
\]

through each anchor. This also proves \(\deg F\le8\); no arbitrary
degree truncation has been imposed.

For a pair of anchors, their difference is a nonzero real scalar times
\(P_{AB}\), \(P_{AC}\), or \(P_{AE}\), respectively. Thus lines with
directions \(P_U,P_V\) can meet only if the three coefficient vectors

\[
 P_U,\quad P_V,\quad P_{AB}
\tag{12}
\]

have rank at most two; replace \(AB\) by \(AC\) or \(AE\) for the other
pairs. This rank test is independent of the nonzero scalar in the anchor
difference.

The three anchors are noncollinear: their two differences from \(F_1\)
have the distinct norm directions \(P_{AB}\) and \(P_{AC}\). Therefore
any common completion must be captured by a pair of nonparallel lines.
It cannot evade the search by lying on coincident lines at every pair of
anchors.

## 3. The complete parameter-dependent rank certificate

The certificate clears denominators separately from each of the three
vectors in (12), and computes gcds of their \(3\times3\) coefficient
minors in \(\mathbb Q[v]\). A nonzero gcd can vanish only at a parameter
where every necessary minor vanishes. Only a subset of minors is needed
when their gcd is already constant.

Before comparison, the certificate removes powers of

\[
 v,\quad v^2-4,\quad v^2-16,
 \quad v^2+a\quad(a=2,4,8,16,32).
\tag{13}
\]

The real zeros of the first three factors are excluded degenerations;
the other factors have no real zeros. This removal is not an omission of
any permissible real parameter.

For each branch the computation covers all \(3\cdot16^2=768\) pairs of
directions. Exactly 108 direction triples have rank at most two identically
in \(v\): 36 per anchor pair. They are precisely:

1. One direction equals the known anchor-difference direction. These
   intersections give an existing anchor or its common pair pencil.
2. Both directions equal a proper common subproduct: \(1\), or either
   of the two common single norm factors. These are parallel distinct
   lines, because the anchor difference has degree four.
3. The following two direction pairs for each anchor pair:

\[
\begin{array}{c|c|c}
\text{anchors}&\text{existing third anchor}&\text{companion}\\ \hline
1,2&(AC,AE)&(BD,BF)\\
1,3&(AB,AE)&(CD,CG)\\
2,3&(AB,AC)&(EF,EG).
\end{array}
\tag{14}
\]

For every direction triple whose generic rank is three, the reduced gcd
is either constant or one of the following polynomials, up to a nonzero
scalar:

\[
\begin{gathered}
 v^2-8,\qquad 5v^2+16,\qquad v^2+20,\qquad
 v^4+52v^2+64,\\
 5v^6+276v^4+1152v^2+1024\quad\text{(minus branch)},\\
 v^6+72v^4+1104v^2+1280\quad\text{(plus branch)}.
\end{gathered}
\tag{15}
\]

Every polynomial in (15) except \(v^2-8\) is strictly positive for all
real \(v\). The exceptional \(v^2-8\) occurs only for the three direction
pairs \((F,G),(AF,AG),(AEF,AEG)\) at anchors 2 and 3; these are the
already excluded conjugate collision \(f=\bar g\).

Thus (14) is exhaustive at every nondegenerate real parameter, including
parameters that were absent from any finite sample. The coefficient
vectors of the nontrivial direction pairs in (14) remain independent:
their distinct root sets prevent proportionality under the nondegeneracy
assumption.

## 4. Excluding all three companions

Let \(J_{12}\) be the companion in the first row of (14). Since it differs
from \(F_1\) by a multiple of \(P_{BD}\) and from \(F_2\) by a multiple
of \(P_{BF}\), its Gaussian polynomial \(J_{12}+i\) has the three roots
\(b,d,f\). Its degree is at most four. None of these three roots occurs
in \(F_3+i\), whose root set is \(a,c,e,g\).

If \(J_{12}\) were compatible with \(F_3\), its difference from \(F_3\)
would divide \(F_3^2+1\). Each real norm factor of this difference would
give a shared root of \(J_{12}+i\) and \(F_3+i\). There can be at most
one such shared root, so

\[
 \deg(J_{12}-F_3)\le2.
\tag{16}
\]

Both the quartic and cubic coefficients of that difference must vanish.
The same argument applies to the other companions, whose Gaussian
polynomials contain the three roots \(c,d,g\) and \(e,f,g\).

The certificate computes those two coefficient differences as rational
functions of \(v\), reduces them, and takes the gcd of their numerators.
Up to a nonzero constant, the results are

\[
\begin{array}{c|c|c}
&J_{12}-F_3\text{ and }J_{13}-F_2&J_{23}-F_1\\ \hline
-&(v^2-16)(v^2+4)(v^2+8)^2
 &(v^2-16)(v^2+4)(v^2+8)\\
+&(v^2-4)(v^2+16)(v^2+8)^2
 &(v^2-4)(v^2+16)(v^2+8).
\end{array}
\tag{17}
\]

The only real zeros are the already excluded repeated-row or pole
parameters. Hence no companion can pass its remaining anchor. Together
with the exhaustive direction classification, this proves the theorem.

## 5. Verification and scope

Running

```text
python3 docs/check_quartic_conic_maximality.py
```

checks both branches, all 1,536 parameter-dependent direction pairs, the
generic incidence lists, positivity or degeneracy of every exceptional
rank gcd, and all six companion cancellation gcds. A preceding finite
search checked 94 distinct rational systems and suggested the stable
six-candidate pattern; that finite evidence is not used in the theorem.
The root orchestrator independently read the exhaustion argument and
certificate and ran both branches successfully. No Lean verification
of this theorem is claimed.

The conclusion strengthens maximality of one explicit quartic triple to
maximality of its entire nondegenerate conic family. It remains a scoped
construction obstruction. A route to the uniform count theorem must
address higher-degree families, nonconstant residues, or other arithmetic
mechanisms beyond this polynomial condition.
