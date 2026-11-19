# Question B on five Gaussian lattice points

**PROVED — main result.** There is a proper closed subset of the similarity torus outside which $A_5(1/31)$ holds for every primitive five-tuple whose Gaussian valuation columns generate a lattice of index at least $7$ in $\mathbb Z^4$. The primes, their number, their exponents, the index, and the lattice itself may all vary. Index at least $5$ already gives an approximation exponent strictly above $10/3$. A separate elementary argument gives exponent $4$, with no exceptional set, when the four ratios are dependent modulo squares in the rational norm-one group. On asymptotically balanced tuples that latter exponent improves to $5$.

**PROVED — consequence for the open lemma.** In the one-unit-class formulation of $(L_5^\circ)$, it suffices to treat **three indices, $1,3,5$**. On the complementary family the required positive lower bound for $H_\Sigma/w$ is uniform. Index $1$ already supports simple balanced squarefree profiles. I obtain no improvement of the essential exponent for unrestricted rational points and no improvement of the growth bound for $M(C)$.

**PROVED — a separate obstruction on the same fourfold.** No nondegenerate, ray-generic member of the explicit reflection-symmetric family of round 24, §5.1, descends to the twisted torus over $\mathbb Q$. Two necessary square conditions contradict each other over $\mathbb Q_2$. This strengthens ineffective finiteness to emptiness on that family.

**HEURISTIC — choice of problem.** I chose B because a valuation lattice of fixed rank gives actual rational lifts of a point to finite covers. A lift changes the volume in the Roth-type theorem. A second argument controls the complexity of the cover when its degree grows. Both operations can be checked on a single five-tuple.

## 1. Heights, columns and conventions

**PROVED.** Let $z_1,\ldots,z_5$ be distinct, Gaussian-gcd-primitive, of norm $H=R^2$. Put $W=\log H$, $U=R^1_{\mathbb Q(i)/\mathbb Q}\mathbb G_m$, $y_j=z_j/z_5$ for $1\le j\le4$, $y_5=1$, and
$$
 \Delta=\max_{i,j}|y_i-y_j|=\operatorname{diam}(z)/R .
$$
Thus $y\in T(\mathbb Q)$ for $T=U^4$. A lower bound $\Delta\gg H^{-1/g}$ gives approximation exponent at least $g$. The arc condition of $A_5(\delta)$ gives $\Delta\le H^{-(1+2\delta)/4}$. This keeps the common phase free.

**PROVED.** By unique factorization, primitivity excludes ramified and inert primes from $H$. At a split $p=\pi_p\bar\pi_p$ with $p^{e_p}\parallel H$, write $a_{ip}=v_{\pi_p}(z_i)$. The conjugate valuation is $e_p-a_{ip}$. The minimum and maximum of the five $a_{ip}$ are $0,e_p$, by primitivity at the two primes. Set
$$
 b_p=(a_{1p}-a_{5p},\ldots,a_{4p}-a_{5p})^t,\qquad
 L=\sum_{p\mid H}\mathbb Z b_p,\qquad g_p=\pi_p/\bar\pi_p .
$$
Then $y_j=\epsilon_j\prod_p g_p^{b_{jp}}$, with $\epsilon_j\in\mu_4$. Write $I=[\mathbb Z^4:L]$ when $L$ has rank four. Changing $\pi_p$ to its conjugate negates a column, and changing the anchor acts by an integral unimodular matrix. Hence $I$ is intrinsic. For full rank it is the gcd of the absolute four-column determinants. Always $\|b_p\|_\infty\le e_p$ and $\sum e_p\log p=W$.

**PROVED.** Every $s\in U(\mathbb Q)$ can be written as $\epsilon A/\bar A$, with $A,\bar A$ coprime Gaussian integers and $A$ of odd norm. Indeed $s=\epsilon\prod_p g_p^{k_p}$, and one takes $\pi_p^{k_p}$ or $\bar\pi_p^{-k_p}$ in $A$ according to sign. Put $D(s)=|A|^2=\prod p^{|k_p|}$. Then $h(s)=\tfrac12\log D(s)$ for the absolute logarithmic Weil height. The complex absolute value of $s$ is $1$, and its reduced denominator has ideal norm $D(s)$, which proves the formula. If $s\ne1$, the nonzero Gaussian integer $\epsilon A-\bar A$ proves
$$
 |s-1|\ge D(s)^{-1/2}. \tag{1}
$$
In particular $d_{ij}:=\log D(y_i/y_j)=\sum_p |a_{ip}-a_{jp}|\log p$, and $\Delta\ge e^{-d_{ij}/2}$ for every distinct pair.

## 2. Uniformity over all valuation lattices of sufficiently large index

**PROVED — theorem.** For any fixed integer $d_0\ge1$ set $B_{d_0}=\tfrac45(70d_0)^{1/4}$. For every $0<g<B_{d_0}$ there are a proper closed $Z\subset T$ and $H_0$ such that every primitive five-tuple with $\operatorname{rank}L=4,\ I\ge d_0,\ H\ge H_0,\ y\notin Z$ satisfies $\Delta>H^{-1/g}$. Both $Z$ and $H_0$ depend only on $d_0,g$.

### 2.1 Bounded index and a ramified cover

**PROVED.** Fix a full-rank sublattice $L$, choose a matrix $P$ whose columns are an integral basis, and form the isogeny $f_P:U^4\to T$ given by that matrix on cocharacters. It has degree $I=|\det P|$. If $b_p=P h_p$, with $h_p\in\mathbb Z^4$, the point with coordinates $x_j=\prod_p g_p^{h_{jp}}$ is rational and $y=\epsilon f_P(x)$. There are only $4^4$ possible $\epsilon$. This constructs the lift with all denominator primes present.

**PROVED.** More generally fix an isogeny $f:T'\to T$ of degree $d$ and a translation $t\in T(\mathbb Q)$. On $t f(T'(\mathbb Q))$ the claimed bound holds for each $g<B_d$, with $Z,H_0$ allowed to depend on $t,f,g$.

*Proof.* Use the projective toric model $X_0$ given by the normal fan of the $A_4$ root polytope. Its ample polarization $B$ pulls back to the anticanonical bundle on the twisted Losev–Manin space. For the height and volume there is a direct check. The polytope is the convex hull of the roots $e_i-e_j$. Its support function at a valuation vector is $\max_i a_i-\min_i a_i$. The two conjugate places each have weight $1/2$, and the complex contribution is zero. Hence its standard logarithmic toric height is $W$. Also $B^4=70$. Indeed the polytope in $\sum x_i=0$ is given by $\sum x_i^+\le1$. Count lattice points in its $n$-th dilate by their $k$ strictly positive coordinates. The coordinate hyperplanes contribute $O(n^3)$, and the count off them is
$$
 \sum_{k=1}^4 {5\choose k}\sum_{s=1}^n
 {s-1\choose k-1}{s-1\choose 4-k}.
$$
Multiplication of the leading coefficient by $4!$ gives $\sum_{k=1}^4 {5\choose k}{3\choose k-1}=70$, the toric intersection number. The Galois involution $x\mapsto-x$ preserves the polytope and its polarization. [Round 22, Prop. 1.1](/home/user/Jarnik/docs/claude_review/round22/toric.md)

Normalize $X_0$ in the function field extension induced by $\psi=t f$. This yields a finite projective morphism $\bar\psi:X_f\to X_0$ of degree $d$. Over the torus it is $\psi$, finite étale in characteristic zero. The ample bundle $B_f=\bar\psi^*B$ has $B_f^4=70d$, and $H_{B_f}(x)\asymp H_B(\psi x)$ for these fixed data. Every target $q\in\bar\psi^{-1}(e)$ is algebraic and smooth, of multiplicity one.

Apply McKinnon–Roth, Corollary 4.2, in dimension four, and Theorem 5.1 at the single real place. The former gives $\beta_q(B_f)\ge (4/5)(70d)^{1/4}$. The latter gives, for $s>1/\beta_q(B_f)$, only finitely many rational solutions of $\operatorname{dist}(x,q)\le H_{B_f}(x)^{-s}$ outside a proper closed subset. These statements allow algebraic targets and apply to the ample bundle on $X_f$. [McKinnon–Roth, §§4–5](https://arxiv.org/abs/1306.2976) · [Full-text statement check, items 4–5](/home/user/Jarnik/docs/claude_review/papers_check/check.md)

Choose $1/B_d<s<1/g$. There are finitely many targets. By properness any real lift of $y$ sufficiently close to $e$ lies near one of them, and by étaleness its distance to that target is $O_{t,f}(\Delta)$. Thus a sufficiently large solution below the claimed downstairs threshold would satisfy the upstairs inequality. A finite morphism sends the finite union of proper exceptional sets to a proper closed set downstairs. Finite exceptions can be absorbed in $H_0$. This proves the fixed-cover assertion. For bounded $I$ there are finitely many sublattices by Hermite normal form, and the Gaussian lift uses only the finite list of unit translations. $\square$

### 2.2 Unbounded index and a logarithm estimate

**PROVED.** For all full-rank tuples one has, with absolute implied constants,
$$
 \frac{-\log\Delta}{W}
 \ \ll\ I^{-1/4}\log(2I)+\frac{\log(2+W)}{W}.
 \tag{2}
$$
In particular a sequence with $I\to\infty$ cannot have any fixed positive exponent of power saving in $\Delta$.

*Proof.* The dual lattice $L^*$ has determinant $1/I$. Minkowski's theorem applied to a cube gives $0\ne v\in L^*$ with $\|v\|_\infty\le 2I^{-1/4}$. Let $m$ be its order modulo $\mathbb Z^4$. Then $m\le I$ and $c=mv\in\mathbb Z^4$. Set
$$
 u=\prod_{p\mid H}g_p^{\langle v,b_p\rangle}\in U(\mathbb Q).
$$
The exponents are integers. At least one is nonzero since $L$ has rank four, hence $u$ is nontorsion. Moreover
$$
 2h(u)=\sum_p |\langle v,b_p\rangle|\log p\le4\|v\|_\infty W .
$$
The unit in $y^c=\epsilon^c u^m$ disappears after taking fourth powers. Telescoping products of complex numbers of modulus one gives $|u^{4m}-1|\le4\|c\|_1\Delta$.

For nontorsion $u\in\mathbb Q(i)$ on the unit circle, Baker–Wüstholz gives
$$
 -\log |u^n-1|\ll (1+h(u))\log(2n). \tag{3}
$$
Indeed take principal logarithms and choose the nearest integer $k$ to $n\arg(u)/(2\pi)$. The nonzero form $n\log u-2k\log(-1)$ has imaginary part in $[-\pi,\pi]$, where $|e^{i\theta}-1|\ge2|\theta|/\pi$. The logarithms have bounded imaginary parts, the field degree is at most two, the modified height of $-1$ is constant, and the integral coefficients are $O(n)$. This gives (3) from the linear-forms-in-logarithms theorem, with its one-logarithm case when $k=0$. The number of logarithms is independent of the number of primes. [Baker and Wüstholz, *Logarithmic forms and group varieties*](https://doi.org/10.1515/crll.1993.442.19)

Now $I\le24(W/\log 5)^4$ by expanding any nonzero four-column determinant, and $\|c\|_1\le8I$. Apply (3) with $n=4m$, substitute the height bound and absorb $O(\log(2I))$ into $O(\log(2+W))$. This proves (2). $\square$

**PROVED — conclusion.** Given $g$, choose a cutoff $I_*(g)$ so that the first term in (2), including its absolute constant, is $<1/(2g)$ for $I>I_*(g)$. Make $W$ large enough for the second term. All remaining indices in $[d_0,I_*]$ are covered by §2.1 with a finite union of exceptional sets. This proves the theorem with uniform quantifiers. In particular $B_5>17/5>10/3$, and $B_7>94/25>124/33=4/(1+2/31)$. The last comparisons use $4.7^4<490$ and cross multiplication. They give $A_5(1/31)$ on $I\ge7$.

**PROVED — interpretation.** The volume bound on the original variety gives only $(4/5)70^{1/4}$. Its amplification upstairs applies to rational points admitting the stated lifts. The cover can ramify along the boundary. This is compatible with the absence of nontrivial geometric étale covers of the projective model. The proper closed set for an arbitrary translated coset can depend on that coset. Uniformity in the index theorem comes from the Gaussian lift, the finite list of unit translates, and (2).

### 2.3 Deleting a fixed set of primes

**PROVED.** Fix once and for all a finite set $S_0$ of split primes. The theorem of §2 remains valid for full-rank tuples if $I$ is computed from the columns *outside* $S_0$, assigning index infinity to rank below four there. The bounds and closed set may also depend on $S_0$. For example, deleting a fixed prime carrying column $e_1$ can leave the lattice with basis $7e_1,e_2,e_3,e_4$, even though the original index is one. These columns are realized by singleton allocations at distinct split primes.

*Proof.* For finite bounded outside index, resolve the residue of $b_p$ modulo the outside lattice separately for each $p\in S_0$, taking $b_p=0$ when absent. The residues have finitely many representatives. Their $g_p$-products give only finitely many extra translates of the cover, since these primes themselves are fixed. The rest of each column lifts integrally.

For finite unbounded outside index take $v,m,c$ from the outside lattice as in §2.2. For $p\in S_0$ write $\langle c,b_p\rangle=m q_p+r_p$, $|r_p|\le m/2$. Include $g_p^{q_p}$ in $u$. Then $h(u)\le 2\|v\|_\infty W+O_{S_0}(1)$, and $u$ is nontorsion by full rank of the outside columns. In place of (3), use logarithms of $u,-1$ and the fixed $g_p$ in the expression for $(y^c)^4=u^{4m}\prod_{p\in S_0}g_p^{4r_p}$. The winding number and the other logarithm coefficients are $O_{S_0}(m)$, and this power is nonidentity by full rank. Baker–Wüstholz gives exactly the bound used in (2) with constant depending on $S_0$.

If the outside rank is $r<4$, choose $0\ne c\in\mathbb Z^4$ annihilating those columns with $\|c\|_\infty\ll(1+W)^r$. For $r>0$ take the signed maximal minors on $r+1$ coordinate positions containing a nonzero $r$-minor of $r$ independent columns, and extend by zeros. For $r=0$ use a coordinate vector. Now $(y^c)^4$ is a nonidentity product of only the fixed $g_p$, with coefficients $O((1+W)^4)$. The same logarithm estimate and telescoping give $-\log\Delta\ll_{S_0}\log(2+W)$. This proves the extension. The deleted set is fixed and the tuple may contain arbitrarily many new primes. $\square$

## 3. Relations modulo squares

**PROVED — theorem.** If $y_1,\ldots,y_4$ have linearly dependent images in $U(\mathbb Q)/U(\mathbb Q)^2$, then
$$
 \Delta\ge 2^{-1/2}H^{-1/4}. \tag{4}
$$
This holds with arbitrary Gaussian units and no exceptional set. In one unit class the hypothesis is exactly failure of the columns $b_p$ to span $\mathbb F_2^4$ modulo two.

*Proof.* A nonzero parity relation on the four anchored ratios has a unique even support on the five marks, obtained by including the anchor if the original support had odd size. Its size is two or four. For any integral lift of its coefficients of sum zero the corresponding character value is a square, since two lifts differ by even coefficients. If $s=t^2\ne1$, both $t-1$ and $t+1$ are nonzero. Their larger absolute value is at least $\sqrt2$ and each is at least $D(t)^{-1/2}$ by (1). Thus
$$
 |s-1|\ge\sqrt2 D(s)^{-1/4}. \tag{5}
$$
For a two-mark relation take $s=y_i/y_j$ and use $D(s)\le H$.

For four marks consider the three partitions $ij|kl$ into pairs and their words $s=y_i y_j/(y_k y_l)$. Each satisfies $|s-1|\le2\Delta$. At a prime the average of the three numbers $|a_i+a_j-a_k-a_l|$ is at most $e_p$. Indeed the sum is convex separately in each entry of $[0,e_p]^4$, so it suffices to test cube vertices. For $0,1,2,3,4$ high entries the sums in units of $e_p$ are $0,3,2,3,0$. After summing over primes, some word has $\log D(s)\le W$. If it is nonidentity, (5) settles (4).

If a word equals $1$, there is an exact affine dependence on the exponent rows. Here is the bound for any such dependence. Normalize its positive and negative coefficients separately to probability vectors $P,Q$ on disjoint sets. Sample $i\sim P,j\sim Q$ independently. At each prime $X=a_{ip}/e_p,\ Y=a_{jp}/e_p$ have the same mean $\mu$. Since $|X-Y|\le X+Y-2XY$, their expected distance is at most $2\mu(1-\mu)\le1/2$. Some pair therefore has $d_{ij}\le W/2$. By (1) it gives $\Delta\ge H^{-1/4}$. Finally, in one unit class all the anchored unit factors are one. A product of the $g_p$ is a square precisely when all its exponents are even, proving the mod-two criterion. $\square$

**PROVED — balanced version.** Suppose each of the fifteen cut weights lies in $[(1-\eta)w,(1+\eta)w]$, with $\eta<1/31$. Under the theorem's squareclass hypothesis,
$$
 \Delta\ge 2^{-1/2}e^{-3(1+\eta)w}. \tag{6}
$$
For one unit class this gives $H_\Sigma\ge(1-7\eta)w-O(1)$. In particular $(L_5^\circ)$ holds quantitatively on this family.

*Proof.* The cut at level $t$ of $p$ partitions the marks according as $a_{ip}\ge t$ or $a_{ip}<t$, and contributes $\log p$ to its unoriented cut weight $\log n_A$. The weights sum to $W$. Each pair is separated by eight cuts, so $d_{ij}\ge8(1-\eta)w>15(1+\eta)w/2\ge W/2$. The common-mean argument excludes all affine dependences. Thus all the words in the proof above are nonidentity. For a coefficient vector $c$ of sum zero, decomposing each column into its nested threshold layers gives
$$
 \log D(y^c)\le\sum_A\left|\sum_{i\in A}c_i\right|\log n_A ,
$$
over the fifteen cuts, with either choice of side. For the coefficients $1,1,-1,-1,0$ these absolute sums total $12$. Orient away from the mark with zero coefficient. The eight odd subsets give $1$ each and the two same-sign pairs give $2$ each. For a two-mark word the total is $8$. Equation (5) now gives (6), and $W\ge15(1-\eta)w$ gives exponent $5(1-\eta)/(1+\eta)$. If $ij$ realizes $\Delta$, then $2|t_{ij}|=\Delta e^{d_{ij}/2}$. Its logarithm yields the claimed bound for $H_\Sigma$, since in one unit class the reduced ratio $y_a/y_b=A/\bar A$ gives $|t_{ab}|=|\operatorname{Im} A|\ge1$. $\square$

**PROVED — reduction to three indices.** Fix $g=94/25$ and use §2 for $I\ge7$. For a balanced tuple outside its $Z$, the diameter pair gives
$$
 H_\Sigma\ge\left[4(1-\eta)-15(1+\eta)/g\right]w-O(1).
$$
The bracket is $(1-751\eta)/94$, hence at least $1/188$ for $\eta\le1/1502$. Balance excludes rank below four by the common-mean argument. Every even $I$ gives rank deficiency modulo two, hence (6) in one unit class. Taking the smaller of the two positive constants and increasing $w_0$ leaves exactly $I=1,3,5$ to treat. For arbitrary unit classes the torsion coordinate in $U(\mathbb Q)/U(\mathbb Q)^2$ must be retained in (4).

**PROVED — the family includes full-rank, balanced, Zariski-dense actual points.** Orient cuts away from mark five. For each nonempty $A\subseteq[4]$, choose a different split prime $p_A$, and give it the two allocations $e_A$ on $A$ and zero off $A$, using common unit $1$. Choose $e_A\log p_A\in[w,w+\sqrt w]$, with $e_A$ odd for the three subsets $12,13,14$ and even otherwise. For all sufficiently large $w$ this is possible. Every prime has both allocations present, so the Gaussian tuple is primitive. The four singleton columns give rank four over $\mathbb Q$ and distinguish its points. The three odd columns span the even-parity hyperplane of $\mathbb F_2^4$. Its unique nonzero relation has support four, which excludes every pair ratio from the square class of one. Every cut is asymptotically balanced and every prime contributes to just one cut.

Hold the nonsingleton exponents fixed at each scale and vary the four even singleton exponents independently across their intervals. The anchored ratios form a product of four sets of distinct elements, up to coordinatewise translation, since the $g_p$ are nontorsion. The cardinality in each factor tends to infinity. A nonzero Laurent polynomial with denominators cleared cannot vanish on the whole grid when each factor has more elements than its degree in that variable, by successive univariate interpolation. This proves Zariski density. The construction describes the domain of the bound and makes no assertion of rapid archimedean approximation.

## 4. What happens at index one

**PROVED — a squarefree subfamily.** Write $w_E=\log n_E$ for the weight of the cut with two-mark side $E$. For every such $E$,
$$
 \min_{i<j} d_{ij}\le (4W+2w_E)/7. \tag{7}
$$
Consequently $\Delta\ge H^{-2/7}$ if even one of the ten two-mark cuts is absent. In particular this gives exponent $7/2>10/3$ on all squarefree primitive norms with at most nine prime factors. With at most six prime factors the pointwise exponent is at least $4$. The primes may vary freely.

*Proof.* Fix $E=\{1,2\}$. Average over the seven pairs meeting $E$. The cut $E$ crosses six of those edges and every other cut crosses at most four. To verify the last count, the incident degrees in the seven-edge graph are $4,4,2,2,2$, and for a two-mark side its crossing number is the sum of its degrees minus twice its internal edge count. Summing cut weights proves (7). Each prime of a squarefree norm supplies just one cut. For the six-cut assertion, if two disjoint two-mark cuts $E,F$ are absent, average over the four pairs between $E$ and $F$. Only $E,F$ can cross all four, and every other cut crosses at most two (this graph is a four-cycle). Otherwise the absent two-mark cuts pairwise intersect. With at most six cuts present, at least four two-mark cuts are absent. Four or more pairwise intersecting edges in $K_5$ have a common vertex, since any pairwise intersecting nonstar is contained in a triangle. Thus all four edges at one vertex are absent, and the six cuts present must be the other six edges. Each present cut crosses just two of the four edges incident to that vertex. This gives $\min d_{ij}\le W/2$ in both cases, and (1) concludes. $\square$

**PROVED — the angle-model obstruction is sharp even for the unrestricted five-point exponent.** For all sufficiently large $X$, choose ten distinct primes $p_E\equiv1\bmod4$ in $[X,2X]$, indexed by pairs $E$ of marks. Such primes exist by the prime number theorem in arithmetic progressions. Give each prime the single cut with side $E$. Then $I=1$, and there are real angles satisfying **all** monomial gaps of the record with $\Delta\asymp X^{-3}\asymp H^{-3/10}$. Their similarity classes can avoid any prescribed proper closed $Z$. This is a five-point angle model with angles chosen in the model. The proof makes no assertion about the arguments of the primes in the Gaussian factorization.

*Proof.* In $\mathbb Z^5/\mathbb Z{\bf1}$ the columns are $e_i+e_j$. Their differences generate the image of the sum-zero lattice, since $(e_i+e_k)-(e_j+e_k)=e_i-e_j$. The quotient by that image is $\mathbb Z/5$, in which a column has class $2$. Thus the column lattice is the whole lattice. Write $B_0$ for the resulting anchored $4\times10$ matrix. The map on angle tori $\phi\mapsto\theta=B_0\phi\bmod 2\pi$ has connected kernel, by Smith normal form. The model point ratios are $\exp(2i\theta_j)$.

Put $s=X^{-3}$. Sample the four $\theta_j$ independently and uniformly from intervals of length $s$ centered at $10js$, $j=1,\ldots,4$, and put $\theta_5=0$. Conditional on $\theta$, sample $\phi$ with Haar probability on the fibre. For a nonzero integer vector $v$ on the ten primes put $w(v)=\sum_E |v_E|\log p_E$. We must avoid the intervals of radii $\arcsin(e^{-w(v)/2})$ about $\frac\pi2\mathbb Z$, and of radii $\arcsin(e^{-w(v)/2}/\sqrt2)$ about $\frac\pi4+\frac\pi2\mathbb Z$, for the angle $v\cdot\phi$.

If $v$ is outside the integral row lattice, it is a nontrivial character of the connected kernel and hence uniform on the angle circle in every fibre. The conditional probability of any such forbidden event, summed over these $v$, is
$$
 O\left(\sum_{0\ne v\in\mathbb Z^{10}} X^{-\|v\|_1/2}\right)=O(X^{-1/2}).
$$
If $v$ is in the row lattice, write $v_{ij}=c_i+c_j$, with $\sum c_i=0$, $c_i\in\mathbb Z$. The angle in question is $\sum_{j=1}^4 c_j\theta_j$. Write $M(c)=\sum_{i<j}|c_i+c_j|$. For a root $c=e_i-e_j$, $M(c)=6$, and the absolute angle is between $9s$ and $41s$, so both sets of gaps hold for large $X$. Every other nonzero $c$ has $M(c)\ge8$. Indeed $2M(c)\ge3\|c\|_1$, by summing $\left|\sum_{j\ne i}(c_i+c_j)\right|=3|c_i|$. At $\|c\|_1=4$ the patterns up to signs and permutation are $(2,-2,0,0,0),(2,-1,-1,0,0),(1,1,-1,-1,0)$, giving $12,12,8$. All larger norms are at least six.

For each of these remaining characters, condition on three coordinates, keeping one with nonzero coefficient random. An interval of length $|c_j|s$ meets $O(1+|c_j|s)$ periodic target intervals, giving probability $O((1+s^{-1})X^{-M(c)/2})$ for its forbidden events. The number of $c$ at a given $M(c)$ grows at most polynomially by the norm estimate. The sum is $O(X^{-1})$. Hence with positive probability all the gaps hold simultaneously. The distribution of $\theta$ has a density on an open four-box. The zero set of any nonzero Laurent polynomial in the $\exp(2i\theta_j)$ has measure zero there, by successive univariate restriction. Thus $Z$ can be avoided. The box forces $\Delta\asymp s$, and $H=\prod p_E\asymp X^{10}$. $\square$

**PROVED — relevance to the two thresholds.** This example tests $10/3$ on arbitrary profiles. For the balanced index-one obstruction to the cover argument, assign instead one distinct split prime from $[X,2X]$ to each of all fifteen cuts, at exponent one. The four singleton columns generate $\mathbb Z^4$ in anchored coordinates, and the cut weights differ by at most $\log2$. The actual Gaussian points for those allocations are primitive. Their true archimedean distances are unspecified. These profiles explain why a bound on $I$ or on the support of the cuts does not by itself give the unrestricted theorem.

## 5. A local obstruction to the round-24 family

**PROVED — statement and scope.** Consider the explicit degree-two projection and reflection-symmetric normal form of round 24, §5.1, including its nondegenerate real family $-3<m<-1$. Assume its thirty distinct boundary places occupy the thirty labelled rays, once each. No such member of this parametrized family defines a curve over $\mathbb Q$ on the twisted torus. The assertion includes its permutations of the marks. It concerns the actual geometric family, and imposes no assumption that a numerical computation identified a whole component of the moduli space. [The normal form](/home/user/Jarnik/docs/claude_review/round24/search.md:280) · [The earlier descent and finiteness argument](/home/user/Jarnik/docs/claude_review/round24/referee_search.md:200)

**PROVED — the descent supplies two square conditions.** Write
$$
 \begin{split}
 c&=m^3-11m^2+39m-13,\qquad D=m^4+70m^2-48m+9,\\
 u&=-\frac{(m+3)c}{16(m-3)^2(m-1)^2},\qquad
 v=-\frac{(m+1)c}{4(m-3)^3(m-1)},\\
 G(Y)&=(Y^2-Y+u)^2+vY .
 \end{split}
$$
If the nondegenerate member descends, $u$ is a nonzero rational square **and** $G$ splits into distinct rational linear factors.

*Proof.* The normal form uses four roots $b_r$ of $b^4+v b^2-v^2 b+uv^2=0$, with elementary symmetric functions $e_1=0,e_2=v,e_3=v^2,e_4=uv^2$. It has $a_rb_r=t\ne0,\ t^2=\prod b_r$, and four distinct $b_r^2$. The labelled boundary places project on the auxiliary $X$-line to the sixteen values $\pm a_r\pm i b_r$. The partners for opposite rays have opposite signs of $i b_r$ with the same sign of $a_r$.

Here is why a common affine rescaling suffices even for a complex member that purportedly descends. Forget the lifted coordinate $u_4$ of the parametrization. The three remaining ratios on the $X$-line have fourteen nonzero primitive valuation vectors, each once, namely the signed indicators of nonempty subsets of $[3]$. This makes the map from that line to the normalization of the projected image birational. Above the image of any one such place, its primitive vector forces ramification index one, and any other preimage would have a positive multiple of the same vector. Uniqueness and the degree formula now give degree one. The three projected ratios are quotients of monic quartics. Contact four at $X=\infty$ makes the numerator minus denominator a nonzero constant, so this is the unique point above the projected identity. The projected normalization thus descends with a rational point at this infinity. Take a rational affine coordinate $x$ there, with $X=\lambda(x-x_0)$.

The same primitivity and uniqueness argument identifies the normalization of the original curve with its parametrizing conic. Every boundary place is individually rational over $\mathbb Q(i)$, since Galois fixing $i$ fixes its ray. Complex conjugation pairs precisely opposite rays. The set of sixteen distinct projected values is Galois-stable, so their average $x_0$ is rational. Taking the half-sum and the difference divided by $2i$ of the conjugate pair $x_0+(a_r\pm ib_r)/\lambda$ gives $a_r/\lambda,b_r/\lambda\in\mathbb Q$. In particular $t/\lambda^2\in\mathbb Q$, and
$$
 e_2/\lambda^2=-\tfrac12\sum_r (b_r/\lambda)^2\ne0,\qquad
 u=e_4/e_2^2=\left(\frac{t/\lambda^2}{e_2/\lambda^2}\right)^2 .
$$
The invariant $v=e_3^2/e_2^3$ is rational. It is nonzero because $e_3=0$ would make the $b$-quartic even, repeating a square. Also $Y_r=-b_r^2/e_2$ is rational. In the displayed normalization, the $b$-quartic gives $Y_r^2-Y_r+u=b_r$, so these four distinct numbers are exactly the roots of $G$.

The parameter itself is rational on this chart. Direct substitution gives
$$
 \frac vu=\frac{4(m^2-1)}{m^2-9},\qquad
 u(m)-u(-m)=-\frac{256m^3}{(m^2-9)^2(m^2-1)^2}.
$$
At regular finite parameters the first gives $m^2\in\mathbb Q$. Unless $m=0$, the second excludes conjugacy of $m,-m$ over $\mathbb Q$. The exceptional $m=-3,-1$ and the zeros of $c$ give $uv=0$, and $1,3$ are poles. At $m=\infty$, $u=-1/16$ is impossible for the twisted descent. This covers the projective parameter line whenever the asserted nondegenerate normal form is present. $\square$

**PROVED — contradiction over $\mathbb Q_2$.** The quartic discriminant identity is
$$
 \operatorname{disc}G
  =u\,(2D)\left[\frac{2(m+1)c}{(m-3)^5(m-1)^3}\right]^2.
$$
It follows by expanding the monic discriminant, an identity checked symbolically in the accompanying script. The Vandermonde formula makes the discriminant a square when all four roots are rational. Hence $A=-(m+3)c$ and $B=2D$ must both be squares. If $v_2(m)<0$, then $v_2(B)=1+4v_2(m)$ is odd. For $m\in2\mathbb Z_2$, $D$ is odd. For odd $m$, substituting $m=r+8k$ into $A,B$ gives the table.

| $r\bmod 8$ | $A\bmod256$ | $B\bmod256$ |
|---:|---:|---:|
| 1 | 192 | 64 |
| 3 | 64 | 128 |
| 5 | 0 | 192 |
| 7 | 128 | 0 |

The residue $128$ has odd valuation, and $192$ has valuation six with odd part $3\bmod4$. Every row excludes simultaneous squares. This proves emptiness on the stated family. It does not require that the putative branches generate a particular real quadratic field. $\square$

**PROVED — split comparison.** The same nondegenerate normal form also has no split rational member. In the split form the sixteen projected values in the rational coordinate are individually rational, giving $a_r/\lambda,ib_r/\lambda\in\mathbb Q$. Thus $it/\lambda^2$ is rational, $-u$ is a nonzero square, and the $Y_r$ are still distinct rational numbers. Consequently $-A,-B$ must be squares. The parameter recovery above applies, including $e_2/\lambda^2\ne0$ since the $ib_r/\lambda$ are real. At $m\equiv1,2\bmod3$, $-B\equiv2\bmod3$. At $m\equiv0,3\bmod9$, $-A\equiv6,3\bmod9$, respectively. At $m\equiv6\bmod9$, $-B\equiv54\bmod81$, of valuation three. For $v_3(m)<0$, or at $m=\infty$, $u\equiv v\equiv-1\bmod3$. The reduction of the monic integral $G$ is then $Y^4+Y^3+2Y^2+Y+1$, with only the root $1$ in $\mathbb F_3$, and differs from $(Y-1)^4$. A split monic integral polynomial over $\mathbb Q_3$ has integral roots, so this excludes those parameters too. The two forms have different local obstructions to descent in this family. $\square$

**COMPUTED.** An independent finite check verifies the pair-cut lattice index, the cut certificates, the character costs and the numerical thresholds using exact arithmetic. It also verifies both coordinates of (8) on a genuine ten-prime tuple. [check_profiles.py](check_profiles.py)

**COMPUTED.** The second script checks the discriminant and parameter identities in $\mathbb Q(m)$, all odd residues modulo $256$, and the split control at $3$. Run it with Python and SymPy. [check_component.py](check_component.py)

## 6. Where the argument stops

**PROVED — exact loss.** At $L=\mathbb Z^4$ the lift furnished by the valuation columns has degree one, so §2.1 has no volume amplification. Also $L^*=\mathbb Z^4$, so §2.2 gives no fractional character of small norm from the index. The ten-pair-cut example above shows that all monomial gaps and the exact value of the index together still permit exponent $10/3$ in the angle model. On fifteen balanced, simple, squarefree cuts the singleton columns also prevent the mod-two gain. The local obstruction in §5 concerns one parametrized family at $15/4$. A sequence of sporadic tuples need not lie on any member of that family.

**EVIDENCE — relation to the record.** The record already has affine rigidity at small arc constant, Gaussian character gaps, and the volume bound on the degree-one model. I found no statement there that treats the unbounded index by a dual lattice logarithm estimate and bounded indices simultaneously by ramified covers, no five-mark squareclass bound of (6), and no local emptiness result for the round-24 family. These results specialize known tools. I see no conflict with the polynomial or fake-circle barriers as stated with their hypotheses. [Affine rigidity citation in round 4](/home/user/Jarnik/docs/claude_review/round4/outside.md:420) · [The original cover discussion](/home/user/Jarnik/docs/claude_review/papers_check/check.md:53) · [The previous rational-parameter search](/home/user/Jarnik/docs/claude_review/round24/referee_search.md:216)

**PROVED — an additive identity for the index-one test.** On an actual one-unit-class tuple with ten simple pair cuts and no singleton cuts, orient every prime toward its two-mark side. For the pair $ij$ let $\beta_{ij}=\beta_{ji}$ be its Gaussian block, chosen in that orientation, with norm $n_{ij}$. Put $P_i=\prod_{j\ne i}\beta_{ij}$. Then $z_i$ is a common scalar times $P_i/\bar P_i$. Define signed $t_{ij}=\det(P_i,P_j)/n_{ij}$, viewing Gaussian integers as plane vectors. This is an integer because $\beta_{ij}$ divides both vectors in the Gaussian ring. Its absolute value is precisely the reduced chord cofactor: $D(y_i/y_j)=|P_i|^2|P_j|^2/n_{ij}^2$ by disjoint prime supports, and $|y_i/y_j-1|=2|\det(P_i,P_j)|/(|P_i||P_j|)$. All $t_{ij}$ are nonzero for distinct points. The exact identity is
$$
 \sum_{i=1}^5\frac{P_i}{\prod_{j\ne i}t_{ij}}=0. \tag{8}
$$
Indeed for every homogeneous cubic $F$ one has $\sum_i F(P_i)/\prod_{j\ne i}\det(P_i,P_j)=0$. In a real determinant-one coordinate frame with all second coordinates $Y_i$ nonzero, put $s_i=X_i/Y_i$. This is the coefficient of degree four in the Lagrange interpolation formula for the cubic $F(s,1)$, divided by the common $\prod Y_i$. Apply it in the original plane with cubic $F(P)=|P|^2 P$, component by component, and use $|P_i|^2=\prod_{j\ne i} n_{ij}$. This proves (8). Clearing denominators by $\prod_{a<b}|t_{ab}|$ gives a relation with integral coefficients of absolute value at most $\exp(H_\Sigma)$. The coefficients move with the tuple, and the $P_i$ can all be large. Equation (8) alone supplies no lower bound for $H_\Sigma$.

**HEURISTIC — next check.** I would ask for an independent referee of §2, especially the passage from the real inverse branches to a single closed exceptional set downstairs. For the next mathematical attack I would start with $I=1$, one varying squarefree prime per two-mark cut and no singleton cuts. I would check whether (8), together with the shared Gaussian divisors of the $P_i$, forces a power-size coefficient outside a fixed proper locus. Even a uniform gain over $10/3$ for those actual Gaussian arguments would distinguish them from the explicit angle models of §4. For the balanced lemma, the corresponding first test has fifteen simple cuts. Neither test follows from controlling the number of prime factors.
