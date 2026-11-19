### The five-point index inside a barrier-scale cluster

**PROVED.** Fix \(b>0\). Retain a single unit class, normalize, and shorten as in Part I to angular width \(\Delta\le e^{-W/4}/2\), where \(W=\log N\). Assume \(W\le b m\log m\). Of the ordered five-tuples of distinct rows, all but \(O_b(m^4)\) have valuation lattice *exactly* \(\mathbb Z^4\). The result holds for gap-admissible angle models as well. It requires the ordinary monomial gap and the common small arc. No arithmetic distribution of prime angles is used.

*Proof.* Put \(d_{ij}=\sum_p |a_{ip}-a_{jp}|\log p\). The arc and the monomial gap, in one unit class, give \(d_{ij}>W/2\) for every pair (Part I §2). For each \(p\), one closed half of \([0,e_p]\) contains the exponents of \(s\ge m/2\) rows. The total weight of the ranges of these rows is \(W'\le W-(e_p/2)\log p\). Each threshold layer separates at most \(s^2/4\) pairs. Their mean separation consequently is at most \(s W'/(2(s-1))\). Comparison with \(W/2\) gives \(e_p\log p\le4W/m\). Hence all primes below \(\sqrt m\) together contribute at most \(4W/\sqrt m<W/2\) once \(m>64\). Projection to columns of primes at least \(\sqrt m\) is injective, since agreeing there would leave a pair of separation below \(W/2\). On these retained columns all entries belong to \(\{0,\dots,D\}\) for \(D=\lceil8b\rceil\).

An affine \(j\)-space over \(\mathbb Q\) in the retained coordinates has at most \((D+1)^j\) rows: projection to some \(j\) coordinates is injective on the space, and each such coordinate has at most \(D+1\) values. Expose the five rows sequentially, without replacement. At the first failure of affine independence the span of previous rows has dimension at most three. A union bound shows that a proportion \(O_b(1/m)\) of tuples fails the test. For every remaining tuple the matrix of row differences has a nonzero four-column minor amongst the retained primes. Its magnitude is at most \((2D)^4\) by the Euclidean row bound for determinants. Every prime dividing its full valuation index thus belongs to the fixed range \(q\le(2D)^4\).

For each such \(q\), any residue class of full exponent rows modulo \(q\) contains at most \(q\) of the cluster points. Write its representative as \(c\). The angles of \(y_i=\exp(2i\langle(a_i-c)/q,\phi\rangle)\) have their \(q\)-th powers in a common arc of width \(\Delta\), by the unit-class restriction. In any one of the \(q\) inverse arcs, two points would give a nonzero integral exponent difference \(v=(a_i-a_j)/q\) and
\[
 |\sin\langle v,\phi\rangle|\le\Delta/(2q)< e^{-W/(2q)}\le e^{-w(v)/2}.
\]
Here \(w(v)\le W/q\); the strict inequality uses \(q\ge2\). This contradicts the monomial gap. Thus an affine \(j\)-space over \(\mathbb F_q\) meets at most \(q^{j+1}\) original rows. Sequential exposure as above bounds the fraction of dependent five-tuples over this field by \(O(q^4/m)\). Sum over the finite range of \(q\) already obtained. Any full-rank matrix losing no rank at these primes has index one by the determinantal-divisor formula. Tuple gcd reduction subtracts constants columnwise and therefore leaves all four row differences unchanged. This proves the claim.

**PROVED — scope.** In particular none of the rate-barrier profiles, including those with repeated exponents, has a positive fraction of high-index tuples. The estimate leaves room for a constant multiple of \(m^4\). It does not compare its constant with the degree of the exceptional set in the fixed-cover theorem, exclude special sparse collections of covered tuples, or settle the uniform-bound regime with \(W/(m\log m)\) unbounded.

**PROVED — balance occurs on almost all five-tuples.** Fix a tolerance \(0<\eta<1\). On every normalized pair-separated profile above, all but \(O_\eta(m^4)\) ordered five-tuples have each of their fifteen nonconstant unoriented cut weights in \([(1-\eta)W/16,(1+\eta)W/16]\). The omitted constant cut has the same estimate, so Gaussian-gcd normalization leaves the tuple log norm within \(\eta W/16\) of \(15W/16\). In particular almost every balanced tuple has index one in the barrier-scale regime.

*Proof.* Choose a threshold layer \(\ell=(p,h)\), \(1\le h\le e_p\), with probability \(\log p/W\). Let \(f_i(\ell)\in\{\pm1\}\) encode the side of row \(i\). The Gram matrix \(G_{ij}=\mathbb E_\ell f_i f_j\) has diagonal one and off-diagonal entries \(1-2d_{ij}/W<0\). Since \(\sum_{ij}G_{ij}\ge0\), this gives \(\sum_{i\ne j}|G_{ij}|\le m\) and \(\operatorname{tr}(G^2)\le2m\). Choose five independent uniform row indices \(I_j\). For an even nonempty subset \(A\subseteq[5]\), put \(T_A=\mathbb E_\ell\prod_{j\in A}f_{I_j}(\ell)\). With independent layers \(\ell,\ell'\) define \(q(\ell,\ell')=m^{-1}\sum_i f_i(\ell)f_i(\ell')\). Independence of row indices gives
\[
 \mathbb E_I T_A^2=\mathbb E_{\ell,\ell'}q^{|A|}
 \le\mathbb E_{\ell,\ell'}q^2=\operatorname{tr}(G^2)/m^2\le 2/m.
\]
For \(\epsilon\in\{\pm1\}^5\), the fraction of layers having pattern \(\epsilon\) or \(-\epsilon\) on this tuple equals
\[
 \frac{1}{16}\sum_{A\subseteq[5],\ |A|\ {\rm even}}
 T_A\prod_{j\in A}\epsilon_j,\qquad T_\varnothing=1,
\]
by expanding the two pattern indicators. Except with probability \(O_\eta(1/m)\) all fifteen nonconstant Walsh moments have absolute value \(\le\eta/15\), by the second-moment bound. The conclusion holds simultaneously for all sixteen patterns. Repeated indices cost \(O(1/m)\) and every nonconstant threshold survives gcd normalization, which removes precisely the constant layers. 
