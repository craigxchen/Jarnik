## IX.7. For rational prime witnesses, one unit on the coordinate rays already recovers the circle

**PROVED (rigidity).** Assume the ordinary monomial gaps (G). For each label \(p_j\) suppose the residues at exterior odd prime powers come from one \(r_j\in\mathbb Q(i)^*\) with \(r_j\bar r_j=p_j\), allowing an arbitrary finite exceptional set of primes for each parent and label. Suppose the strengthened e-all inequality holds, at every finite cutoff on the **same** data, just for
\[
                 v=n e_j,\qquad n=1,2,\ldots,\qquad \text{unit }-1 .
\]
Then \(r_j\) is an actual Gaussian prime factor and \(\phi_j\equiv\pm\arg r_j\bmod\pi\) for each \(j\). These conclusions realize the marks on the integer circle: one can choose which of the two conjugate primes to use at each label. If the same one-unit tests are added for \(v=n(e_i+e_j)\) on the edges of any spanning tree, all signs are equal relative to the given \(r_j\). The finite exceptional sets can depend arbitrarily on the instance. There is no height or uniformity assumption on them. Rational coherence outside such a set is a global hypothesis. General compatible local residues need not admit these lifts.

*Proof.* Enlarge a fixed finite set \(E\) to contain 2, all parent primes, all exceptional primes and all denominators. At a label write \(r_j/\bar r_j=u\alpha/\bar\alpha\), where \(u\) is a Gaussian unit and \(\alpha\) is an integral product on just one side of each conjugate prime pair. Put \(J=|\alpha|^2\) and include its primes in \(E\). At each split prime the difference of the valuations of \(r_j\) has the parity of their sum, which is one at \(p_j\) and zero elsewhere. In particular \(J\ge p_j\) and \(J\) is odd. The angle of \(\alpha\) divided by \(\pi\) is irrational, since \(\alpha,\bar\alpha\) have disjoint prime supports and \(J>1\).

Choose a positive multiple \(L\) of four such that \(\alpha^L\equiv1\pmod\ell\) for all \(\ell\in E\) with \(\ell\nmid J\), possible by finite unit groups. For \(n=Lk\), \(k\ge1\), the integer \(\Re(\alpha^n)\) is coprime to **every** prime of \(E\). At \(\ell\nmid J\) it is 1 modulo \(\ell\), including \(\ell=2\). At \(\ell\mid J\) the two components of \(\alpha^n\) in the split algebra are zero and a nonzero unit, so its trace \(2\Re(\alpha^n)\) is nonzero there modulo the rational prime. For each external prime power dividing \(\Re(\alpha^n)\), the reduction of \((r_j/\bar r_j)^n\) is \(-1\), since \(u^n=1\) and the conjugate denominator is a unit. Thus the *entire* absolute real part is charged, with e-all constant one for the even unit index \(s=2\). In particular
\[
 |\cos(n\phi_j)|\ge p_j^{-n/2}|\Re(\alpha^n)|\quad(L\mid n).
                                                                    \tag{IX.9}
\]
Along infinitely many multiples of \(L\), irrational rotation gives \(|\Re(\alpha^n)|\ge J^{n/2}/2\). This proves \(J\le p_j\), hence equality. The two valuations of \(r_j\) at \(p_j\) have sum one and difference \(\pm1\). Every other pair has sum and difference zero, and the norm condition rules out inert or ramified factors. Therefore \(r_j\) is an integral Gaussian prime of norm \(p_j\).

Here is the angular part for both kinds of word. Take \(z=r_j\), or \(z=r_i r_j\) with distinct labels on a tested edge. Write \(\beta=\arg z\), \(\varphi=\phi_j\) or \(\phi_i+\phi_j\), respectively. The weight of the unpowered word is now exactly \(\log|z|^2\). Both \(\beta/\pi\) and \(\varphi/\pi\) are irrational, by disjoint genuine factors and by (G) on multiples of the word. Repeating the clean-real-part argument using \(z\) gives some \(L\) such that
\[
 |\cos(kL\varphi)|\ge |\cos(kL\beta)|\qquad(k\ge1).
\]
Put \(a_k=\cos(2kL\varphi)-\cos(2kL\beta)\). Squaring the preceding inequality gives \(0\le a_k\le2\). Its Cesàro mean is zero by the finite geometric sum, so the mean of \(a_k^2\le2a_k\) is also zero. If \(L\varphi\not\equiv\pm L\beta\bmod\pi\), expanding this square would give mean one: each cosine square has mean one half, and both frequencies in the mixed product have zero mean. Thus
\[
                 \varphi\equiv \varepsilon\beta+\pi t/L\pmod\pi
\]
for a sign \(\varepsilon\) and an integer \(t\).

Finally use powers through \(n=1\). Let \(B\) be the \(E\)-part of \(|\Re z|\), which is finite and positive. There is a multiple \(L'\) of \(L\) such that for all positive \(n\equiv1\bmod L'\) the \(E\)-part of \(|\Re(z^n)|\) equals \(B\). At primes dividing \(|z|^2\) the real part never vanishes modulo the rational prime. At each remaining prime take the period of \(z\) modulo one power higher than its real-part valuation at \(n=1\). The exterior divisors now give \(|\cos(n\varphi)|\ge |\cos(n\beta)|/B\). Along this progression \(n\varphi\bmod\pi\) is dense. On a subsequence tending to \(\pi/2\), the right numerator tends to \(|\sin(\pi t/L)|\), since \(nt/L\equiv t/L\bmod1\). Therefore \(L\mid t\). This proves the single-label assertion and, when requested, the analogous angle congruence for each tested pair.

For an edge, subtracting the three congruences gives an integral relation modulo \(\pi\) between \(\phi_i\) and \(\phi_j\) unless the two single signs and the edge sign are identical. Such a relation is excluded by (G). A tree therefore synchronizes the signs. Even without the tree, choose \(r_j\) or \(\bar r_j\) according to the single-label sign. Substituting in \(2a-e\) now recovers the mark phases up to their Gaussian units: the additional \(\pi\)-multiples contribute the same sign to every mark since \(2a_{ij}-e_j\equiv e_j\bmod2\). This keeps the drifting common phase as well as all row windings. \(\square\)

**HEURISTIC (interpretation).** The all-height P3-char construction in IX.3 already satisfies the rational-lift hypothesis. Within that subclass, adding powers of *individual* labels at one unit has the full force of requiring genuine angles. For arbitrary prime-adic witnesses, or for finite cutoffs on rational lifts of unbounded height, I do not have this conclusion. Irrational rotations and the real cosine comparison are the archimedean steps.
