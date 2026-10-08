# Direct elimination of all multiplicatively dependent endpoint pairs

## Status

The uniform endpoint arc bound remains unproved. The lemma below strengthens
Section 8 of `codex_uniform_bound_research.md`: its uniformly bounded incidence
conclusion holds for **all** multiplicatively dependent pairs, without invoking
Corvaja--Zannier and without assuming an a priori bound on the relation
exponents. The remaining pairs can therefore be required to be pairwise
multiplicatively independent, even modulo Gaussian roots of unity.

This does not bound the finite residual or its clique number. In particular,
the complete hypersimplex obstruction in the earlier note already has distinct
nonzero valuation rows that are pairwise linearly independent.

## 1. A direct exponent cutoff

Fix \(C>0\). Let \(z_0,z_1,z_2\in\mathbf Z[i]\) be distinct points on a
circle centered at the origin of radius \(R\), on one arc of length at most
\(C\sqrt R\), with \(z_0\) at an endpoint of the arc. Put

\[
u=z_1/z_0,\qquad v=z_2/z_0,\qquad \delta=C/\sqrt R.
\]

The two ratios have positive lifted arguments in an interval of length at
most \(\delta\), and

\[
|u-1|\leq\delta,\qquad |v-1|\leq\delta.
\tag{1}
\]

Write the split part of \(R^2\in\mathbf Z\) as
\(\prod_jp_j^{e_j}\), choose \(p_j=\pi_j\bar\pi_j\), and set

\[
g_j=\pi_j/\bar\pi_j,\qquad
Q=\prod_jp_j^{e_j/2}\leq R.
\]

Unique factorization gives

\[
u=\zeta\prod_jg_j^{d_j},\qquad
v=\xi\prod_jg_j^{d'_j},\qquad
|d_j|,|d'_j|\leq e_j,
\tag{2}
\]

with \(\zeta,\xi\in\mu_4\). Ramified and inert factors cancel in the
ratios, up to a Gaussian unit. If \(\delta<1\), a non-root ratio cannot
be a root of unity, so both valuation rows in (2) are nonzero.

More explicitly, if \(a_{0j},a_{1j},a_{2j}\in\{0,\ldots,e_j\}\)
are the exponents of \(\pi_j\) in the three original Gaussian points,
then \(d_j=a_{1j}-a_{0j}\) and \(d'_j=a_{2j}-a_{0j}\). Hence

\[
\max(0,d_j,d'_j)-\min(0,d_j,d'_j)\leq e_j.
\tag{2a}
\]

The two differences may have opposite signs at a prime when the root
exponent is interior. No common orientation for the entire divisor grid
is assumed. Below, only the primitive row \(\mathbf c\) of the one
dependent pair is oriented.

Suppose \(u,v\) are multiplicatively dependent modulo \(\mu_4\). Taking
valuations, and using that the integer kernel of a primitive linear equation
is saturated, gives

\[
d_j=a c_j,\qquad d'_j=b c_j,
\qquad (a,b)=1,\qquad ab\ne0,
\tag{3}
\]

for one nonzero integer vector \(\mathbf c\). Put

\[
M=\max(|a|,|b|),\qquad
W=\prod_{c_j>0}\pi_j^{c_j}
   \prod_{c_j<0}\bar\pi_j^{-c_j}.
\]

Then \(W\ne1\), \((W,\bar W)=1\), and

\[
|W|^M\leq Q\leq R,
\qquad
M\leq \frac{\log R}{\log\sqrt2}.
\tag{4}
\]

The first assertion follows term by term from
\(M|c_j|\leq e_j\). The second uses the weaker convenient bound
\(|W|\geq\sqrt2\).

Choose Bezout coefficients \(s,t\in\mathbf Z\) satisfying
\(sa+tb=1\) and \(|s|+|t|\leq2M\). Then

\[
w=u^s v^t=\eta W/\bar W\qquad(\eta\in\mu_4).
\]

For unit complex numbers, telescoping gives
\(|q^k-1|\leq |k|\,|q-1|\) for every integer \(k\). Hence (1) implies

\[
|w-1|\leq2M\delta.
\]

On the other hand, \(\eta W-\bar W\) is a nonzero Gaussian integer:
if it vanished, coprimality would make \(W\) a unit. Thus

\[
|w-1|\geq |W|^{-1}\geq R^{-1/M}.
\tag{5}
\]

If \(M\geq3\), (4)--(5) give

\[
R^{1/6}\leq2CM
\leq\frac{2C\log R}{\log\sqrt2}.
\tag{6}
\]

The elementary maximum of \((\log x)x^{-1/12}\), for \(x\geq1\), is
\(12/e\). Therefore (6) is impossible when

\[
R>\left(\frac{24C}{e\log\sqrt2}\right)^{12}.
\tag{7}
\]

Consequently every dependent pair has \(M\leq2\) beyond (7).

## 2. Removing the remaining units and signs

Define the explicit threshold

\[
R_0(C)=\max\left\{
4,\quad 81C^2,\quad
\left(\frac{24C}{e\log\sqrt2}\right)^{12}
\right\}.
\tag{8}
\]

Assume \(R>R_0(C)\), so \(\delta<1/9\). If \(M\leq2\), one of
\(|a|,|b|\) equals one. Suppose first that \(|a|=1\), and replace
\((a,b,\mathbf c)\) by their simultaneous sign reversal when necessary,
so \(a=1\). The valuation rows give

\[
v=\omega u^b,\qquad \omega\in\mu_4,
\qquad b\in\{-2,-1,1,2\}.
\]

Using (1),

\[
|\omega-1|=|v-u^b|
\leq |v-1|+|u^b-1|
\leq3\delta<\sqrt2.
\]

Distinct fourth roots of unity are separated by at least \(\sqrt2\), so
\(\omega=1\). The positive lifted arguments of \(u,v\) and
\(\delta<1/9\) rule out \(b<0\); equality \(b=1\) contradicts distinctness.
Thus \(v=u^2\). The case \(|b|=1\) gives \(u=v^2\).

We have proved the following unconditional statement:

> **Dependent-pair lemma.** For \(R>R_0(C)\), two distinct non-root
> ratios from an endpoint arc of length at most \(C\sqrt R\) that are
> multiplicatively dependent modulo roots of unity necessarily satisfy
> \(v=u^2\) or \(u=v^2\).

The one-sided arc hypothesis is used only at this final sign step. Without
it, inverse and inverse-square relations must also be retained.

## 3. An explicit uniform incidence bound

Consider all the non-root ratios from the same arc. For a relation
\(v=u^2\), write \(u=\zeta W/\bar W\) using the orientation of its
valuation vector. Define

\[
G_{\mathbf c}=\prod_{c_j>0}\pi_j^{e_j}
                  \prod_{c_j<0}\bar\pi_j^{e_j}
                  \prod_{c_j=0}\pi_j^{e_j}.
\]

Then \(W^2\mid G_{\mathbf c}\), \(|G_{\mathbf c}|=Q\), and
\(D=G_{\mathbf c}/W^2\in\mathbf Z[i]\setminus\{0\}\). Gaussian
integrality and endpoint proximity yield

\[
\frac{\sqrt R}{C}\leq|W|\leq\sqrt Q\leq\sqrt R,
\qquad
|D|\leq C^2,
\qquad
|\zeta W-\bar W|\leq C.
\tag{9}
\]

There are at most \((2\lfloor C^2\rfloor+1)^2\) possible \(D\)'s. For
fixed \(D\), the norm \(|W|^2=Q/|D|\) is fixed, even though the orientation
of \(G_{\mathbf c}\) may vary.

Writing \(W=x+iy\), the last inequality in (9), for the four choices of
\(\zeta\), bounds respectively one of \(2x,2y,\sqrt2(x+y),\sqrt2(x-y)\)
in absolute value by \(C\). Each bounded integer linear form has at most
\(2\lfloor C\rfloor+1\) values. Intersecting each resulting line with the
fixed circle gives at most two choices of \(W\). Accounting for all four
units gives at most \(8(2\lfloor C\rfloor+1)\) choices for each \(D\).

Thus there are at most

\[
B_{\mathrm{dep}}(C)
=16(2\lfloor C\rfloor+1)(2\lfloor C^2\rfloor+1)^2
\tag{10}
\]

ordered dependent pairs among all the non-root ratios, counting both
\(v=u^2\) and \(u=v^2\). This bound is deliberately coarse. It is independent
of \(R\), of the number and identities of the conductor primes, and of their
exponents.

For \(0<C<1\), (9) would force
\(\sqrt R<\sqrt R/C\leq |W|\leq\sqrt R\), a contradiction. Thus there
are **no** dependent pairs at all once \(R>R_0(C)\) in that range.
The nontrivial Gaussian numerator actually has modulus at least
\(\sqrt2\) in this unit-conjugate form, which could sharpen this constant;
the stated corollary uses only the lower bound already proved.

In particular, deleting at most \(B_{\mathrm{dep}}(C)\) vertices leaves
pairwise multiplicatively independent non-root ratios. Therefore an
unbounded endpoint cluster would necessarily yield unbounded clusters whose
distinct non-root ratios are pairwise multiplicatively independent modulo
\(\mu_4\). This is a direct, unconditional restriction on a hypothetical
counterexample. It supplies no further bound for those independent clusters.

## 4. Literature check

The primary problem statement remains Cilleruelo--Granville, *Close lattice
points on circles*, Section 12, Problem 2
([author PDF](https://matematicas.uam.es/~franciscojavier.cilleruelo/Papers/closelattice.pdf)).
The searches in this continuation did not identify a new theorem with the
uniform endpoint quantifiers needed here. This is a report on the search,
not a proof that no such theorem exists.

No support-independent assertion about the Corvaja--Zannier finite residual
is inferred from its qualitative theorem. The new lemma above is elementary
and does not use any of the previously audited fixed-weight substitutions.

## 5. Formalized ingredients

[GaussianChain/DependentEndpoint.lean](../GaussianChain/DependentEndpoint.lean)
proves Gaussian-integer separation, its unit-conjugate specialization, and
the analytic primitive-exponent cutoff. The ordinary-power version assumes

\[
R>\left(\frac{24C}{\log\sqrt2}\right)^{12}.
\]

This sufficient threshold is weaker than (7), because the formal proof uses
\(\log R\leq12R^{1/12}\) and omits the elementary improvement by \(1/e\).
The Lean theorem explicitly assumes the capacity and separation inequalities
derived in (4)--(5); it does not formalize the factorization of circle
ratios, Bezout construction, Gaussian-unit removal, or incidence counting.
It also proves the Gaussian-integrality input for the \(C<1\) corollary.

The new module compiles directly with the project's Lean toolchain. It has
no placeholders and makes no claim to have formalized the full dependent-pair
lemma or the still-unproved uniform endpoint theorem.
