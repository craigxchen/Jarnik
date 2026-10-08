# Actual cotangent cliques can saturate derivative and resultant concentration

This note audits a direct polynomial route to the integer lcm height target
in [integer_cotangent_lcm_height_target.md](integer_cotangent_lcm_height_target.md).
It gives an exact root-depth formula for the conductor and globally integral
positive cliques of every fixed size that saturate the elementary derivative
budget at arbitrarily deep split prime powers. Thus the proposed improvement
cannot come from reducing that derivative exponent universally, even after
imposing all actual pair divisibilities. No endpoint height improvement is
proved here.

This is separate from the generic polynomial contact obstruction: the
examples below are ordinary integer cotangent cliques, with every pair
quotient integral. They are not short-arc counterexamples to the uniform goal.

## 1. Derivative and resultant identities at actual integer roots

Let \(L>0\), let \(X_1,\ldots,X_k\) be distinct positive integers, and
assume
\[
 Q_{ij}=(X_iX_j+L^2)/(X_i-X_j)\in\mathbb Z.
\]
Write \(C_i=X_i^2+L^2\) and \(F(t)=\prod_i(t-X_i)\). Since
\(X_i-X_j\mid C_i\),
\[
 \boxed{F'(X_i)\mid C_i^{k-1}.} \tag{1}
\]
The other basic identity is
\[
 \boxed{F(iL)F(-iL)=\prod_i C_i,} \tag{2}
\]
which is the resultant with \(t^2+L^2\). The all-edge conductor does not
record these total root depths; it records their extremes at both Gaussian
orientations.

Specifically fix an odd split prime \(p=\pi\bar\pi\), and put
\[
 \alpha_i=v_\pi(X_i-iL),\qquad
 \gamma_i=v_\pi(X_i+iL),\qquad \tau=v_p(L).
\]
The difference between the two arguments is \(2iL\), of valuation \(\tau\).
If \(\alpha_i>\tau\), then \(\gamma_i=\tau\), and conversely. If neither
exceeds \(\tau\), they are equal. Thus for the phases
\(q_i=(X_i+iL)/(X_i-iL)\), including the anchor \(q_0=1\),
\[
 \boxed{v_p(N)=
 (\max_i\alpha_i-\tau)_++(\max_i\gamma_i-\tau)_+.} \tag{3}
\]
This is exactly the signed-valuation width of the primitive realization.
By contrast, the valuation of (2) is
\(\sum_i\alpha_i+\sum_i\gamma_i\). In principle the extremes in (3) can
be recovered from the Newton polygons of \(F(iL+t)\) and \(F(-iL+t)\),
but replacing the totals by extrema is the arithmetic information one must
actually retain. An identity involving only the resultant cannot silently
make that replacement.

## 2. A positive integral construction of arbitrary clique size

Fix any \(k\ge2\), choose a split prime \(p>k\), and put
\[
 L=\operatorname{lcm}(1,2,\ldots,k-1),\qquad P=p^e,
 \qquad e\ge1.
\]
Choose \(a>0\) satisfying
\[
 a^2+1\equiv0\pmod P,\qquad
 v_p((a+jP)^2+1)=e\quad(0\le j<k). \tag{4}
\]
Such an \(a\) always exists. A root modulo \(p^e\) lifts by the elementary
nonsingular square-root lifting rule. Starting with any representative
\(a_0\), replace it by \(a_0+cP\). Modulo \(p\), the normalized expression
\(((a_0+(c+j)P)^2+1)/P\) is a nonconstant affine function of \(c+j\).
Each of the \(k\) indices forbids one residue of \(c\); since \(p>k\),
some residue remains. A positive representative can be chosen.

Set
\[
 \boxed{X_j=L(a+jP),\qquad 0\le j<k.} \tag{5}
\]
These are distinct positive integers with \(X_j\ge L\). For every pair,
\[
 Q_{ij}=\frac{L}{i-j}
 \frac{(a+iP)(a+jP)+1}{P}\in\mathbb Z. \tag{6}
\]
The first factor is integral by the definition of \(L\); the second is
integral by \(a^2+1\equiv0\pmod P\). Thus all the required pair conditions
hold globally, without CRT assumptions about unspecified other primes.

Since \(p\nmid L\) and \(p>k\), every pair difference has valuation \(e\).
By (4), every \(C_i\) has valuation \(e\). Hence
\[
 \boxed{v_p(F'(X_i))=(k-1)e,\qquad v_p(C_i)=e.} \tag{7}
\]
The exponent \(k-1\) in (1) is therefore sharp for actual cliques, for
every fixed \(k\), at arbitrarily large depths.

All \(a+jP\) approach the same chosen square root of \(-1\) at this
prime. Consequently every anchor phase has the same nonzero signed
valuation \(e\) or \(-e\), while every finite-finite phase quotient has
signed valuation zero. Formula (3) gives
\[
 \boxed{v_p(N)=e,\qquad v_p(F(iL)F(-iL))=ke.} \tag{8}
\]
At the deep Gaussian orientation every summand of \(F^{(r)}(iL)\) is a
product of \(k-r\) factors of valuation \(e\), so
\[
 v_\pi(F^{(r)}(iL))\ge(k-r)e.
\]
Derivative gcds therefore retain this large shared depth as well; they do
not automatically turn the resultant's multiplicity into an extra conductor
factor.

For a small explicit case, take \(k=4\), \(p=13\), \(e=1\), \(a=5\).
Then \(L=6\) and \(X=(30,108,186,264)\) form such a clique. All six finite
pair quotients are integral. Each root derivative has 13-adic order three,
while the exact all-edge conductor has 13-adic order one.

## 3. Its actual height lies outside the endpoint regime

In (4), choose the permitted shift with \(1\le c\le p\), so that
\(P\le a\le(p+1)P\). Then
\[
 A/L=a\asymp_p P.
\]
The true all-edge conductor also has a precise order:
\[
 \boxed{N\asymp_{k,p}P^{k+1}.} \tag{9}
\]
To prove it, cancel the common integer scale \(L\) and write
\(H_j=a+jP+i\). Every reduced Gaussian denominator of
\(H_j/\bar H_j\) has norm comparable to \(P^2\); the possible self-conjugate
factor has norm at most two. All these denominators share exactly one
Gaussian orientation over \(p\) to depth \(e\), and none has the opposite
orientation. Remove that common factor. A gcd of two remaining denominators
divides the bounded integer \(i-j\): before removal it divides
\(P(i-j)\), and both orientations over \(p\) are absent afterward.
Their Gaussian lcm norm is therefore comparable, with constants depending
only on \(k,p\), to the product of the remaining norms. Restoring the one
common factor gives \(P\cdot P^k=P^{k+1}\). The Gaussian denominator lcm
norm equals the exact all-edge ordinary conductor, so this proves (9).

Consequently, for \(k\ge4\),
\[
 \frac{N}{(A/L)^4}\asymp_{k,p}P^{k-3}\longrightarrow\infty.
\]
The construction therefore **does not refute an endpoint-restricted or
low-conductor refinement** of the derivative budget. It rules out only the
corresponding unconditional refinement.

## 4. What low conductor already forces in the weighted budget

Suppose instead that an actual clique has
\(N=o((A/L)^4)\), and let \(m=k+1\) include its anchor. Its actual endpoint
constant tends to zero, so eventually the known pair-norm lower bounds at
constant two apply. With \(W=\log N\), the full cut deficit satisfies
\[
 D=\sum_{\text{prime layers}}(n-m/2)^2\log p\le mW/4.
\]
A layer in which all \(k\) finite rows share one deep isotropic orientation
and only the anchor lies on the other side has deficit
\((m-2)^2/4\). Hence the total logarithmic weight of such completely
concentrated layers is at most
\[
 \frac{m}{(m-2)^2}W.
\]
Thus low conductor does force these layers to occupy only an \(O(1/m)\)
fraction of the source norm. This is a genuine restriction missing from
the examples in Section 2, but it is the existing Plotkin-deficit budget.
It does not by itself turn a fixed-size exponent below four into the
endpoint exponent four. A further weighted derivative/resultant estimate
would have to use an additional relation between the remaining balanced
prime layers and the actual real root values. No such relation is proved
by (1)--(3) or the concentration example.

## 5. What this rules out, and what it does not

The construction rules out a universal replacement of (1) by a lower power
of \(C_i\), even for arbitrarily large fixed clique sizes. It also rules out
counting repeated resultant depth as if each occurrence supplied a distinct
prime-power contribution to the all-edge lcm. The examples retain actual
integer ordering and actual pair divisibility, so those hypotheses alone do
not repair either step.

It does not refute the desired height inequality. The other prime factors
of these constructed cliques can make their conductor much larger than the
endpoint threshold. An improvement could still exploit the interaction of
all prime-depth clusters with the real sizes of the roots. The derivative
and resultant identities above do not presently provide that interaction.
In particular the fixed-\(L\) nature of the example is used only to disprove
the proposed derivative refinement; it is not being promoted to a solution
of the general varying-\(L\) problem.

The [persistent exact checker](check_integer_cotangent_derivative_concentration.py)
uses 36 positive cliques with \(4\le k\le12\)
and \(1\le e\le4\). They checked every pair quotient, the actual all-edge
lcm against a Gaussian realization, and all root derivative valuations in
(7)--(8). All checks passed.
