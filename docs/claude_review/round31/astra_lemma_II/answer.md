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

# Part II — The parent norm, covers at index one, and the ten pair cuts

**PROVED — principal arithmetic result.** Take the one-unit-class ten-pair model of (8), with pairwise coprime block norms. If the norm vector belongs to a fixed rational subspace $U$ of dimension at most two, equation (8) and the shared Gaussian divisors give
$$
 \max_{ij}|t_{ij}|\gg_U (\max_{ij}n_{ij})^{1/6}.
$$
The constant may depend ineffectively on the subspace $U$. This gives an approximation exponent at least $38/11>10/3$ for all its genuine simple ten-pair tuples, and $40/11$ when the norms are comparable. An elementary argument already gives exponent $4$ on a generic plane with comparable norms, including an explicit family of squarefree index-one profiles. The harder planes lead to a fixed genus-five curve and a real five-star lemma. See II.6–II.7 for the proofs.

**PROVED — an independent short-interval result.** In the simple ten-pair family there is a pointwise gain over $10/3$ if the six block norms separating even one pair lie in $[X,X+cX^{1/6}]$, for an absolute $c>0$. That pair's scalar cofactor has size $\gg\sqrt X$. If the other four block norms are comparable to $X$, the resulting lower bound is $\Delta\gg H^{-1/4}$, with no exceptional set. A second test below uses the residues of the ten norms modulo a large integer, without requiring a short interval. These restrictions on the norms are additional hypotheses. The general case of comparable norms with no further structure stays open.

**PROVED — scope of the cover argument.** For index $7$ the height comparison in the question goes the other way. The primitive subtuple norm generally needs to be *smaller by a power* than the parent norm to apply §2. Large clusters do have this property on almost every subtuple, but at the growth scale almost every subtuple, including in the angle models, has index one. There is an interior-branched cover that lifts index-one points with bounded cofactors. Its target point upstairs is singular.

**EVIDENCE — status of the polynomial test.** I do not know if the full function-field system has a nondegenerate solution, even with linear blocks. Below are an exact reconstruction from four base blocks, an obstruction valid in every degree for norm polynomials spanning a pencil, and a finite-field exclusion for linear blocks. The linear-block moment problem was already stated explicitly in the record. [Previous five-point system and numerical search](/home/user/Jarnik/docs/claude_review/growth/cotangent.md:245)

## II.1. Which height enters the grid argument?

**PROVED.** Use a primitive parent of norm $N$, write $W=\log N$, and let its $m$ points occupy an angular arc of width at most $C N^{-1/4}$. Angular width bounds chordal width. For any five distinct marks $S$, write $G_S$ for their Gaussian gcd and $H_S=N/|G_S|^2$. Division by $G_S$ leaves all ratios and valuation-column differences unchanged. At each split $p$ its effect is to replace $e_p$ by $\max_{S}a_{ip}-\min_{S}a_{ip}$. In particular (1) gives $H_S\ge C^{-2}N^{1/2}$. The sufficient comparison for §2 at a given $g$ is
$$
 H_S\le C^{-g}N^{g/4},\qquad\text{equivalently}\qquad
 |G_S|^2\ge C^g N^{1-g/4}. \tag{II.1}
$$
Take $g=94/25$ at index at least $7$. This requires the exponent $47/50$ on the primitive norm, whereas a bounded-factor comparison gives exponent $1$. The §2 theorem therefore does not imply the proposed fixed-set inclusion at $H_S\asymp N$. At index at least $9$, $70\cdot9>5^4$ and one can use $4<g<B_9$. Then every $H_S\le N$ qualifies for large $N$. For squareclass dependence, (4) gives the necessary condition $|G_S|^2\le4C^4$. In particular every such tuple with larger gcd is excluded, without a closed set.

**PROVED — bounded gcd really can occur with dependence.** Let $u_j=t+c_j+i$ for $c_j=0,1,3,7$ and take five of the six two-element sets $A\subset[4]$. Put $z_A=\prod_{j\in A}u_j\prod_{j\notin A}\bar u_j$ for positive integral $t$. Their norms agree. All five monic quartics have the same cubic coefficient, so their differences have degree at most two, and their common radius is asymptotic to $t^4$. Their quadratic coefficients have distinct imaginary parts, since the six sums $\sum_{j\in A}c_j$ are distinct. Thus they occupy a critical arc with an absolute fixed $C$. At least two complementary pairs survive and give $z_A z_{A^c}=z_B z_{B^c}$, an exact squareclass relation. The eight linear factors are distinct and no one divides all five polynomials. Bézout over $\mathbb Q(i)[t]$, followed by clearing denominators, bounds the Gaussian gcd at every specialization. Divide by it and take the five points as the parent. Now $H_S=N\to\infty$ and the relation persists. This example concerns squareclass dependence. Its valuation rows are affinely dependent. [Exact polynomial checks](check_part_ii.py)

**PROVED — almost every gcd supplies the missing power.** For every $\epsilon>0$, all but $O_{C,\epsilon}(m^4)$ ordered five-tuples satisfy
$$
 \left|\frac{\log H_S}{W}-\frac{15}{16}\right|\le\epsilon . \tag{II.2}
$$
This also holds for the angle model with the axial monomial gap as stated in the expert note, even with differing units. Here and below we may discard bounded $N$, since the number of distinct lattice points, or distinct allocation rows with their units, is then bounded.

*Proof.* Expand each prime into $e_p$ threshold layers $l=(p,k)$ of weight $w_l=\log p$. Write $\sigma_{il}=2{\bf1}_{a_{ip}\ge k}-1$ and $V_{il}=\sqrt{w_l/W}\,\sigma_{il}$. Write $d_{ij}=\sum_p|a_{ip}-a_{jp}|\log p$. Split the arc into a bounded number of arcs of width at most $\tfrac12N^{-1/4}$, and separate the four unit classes, giving a total of $L=O_C(1)$. Within a piece $d_{ij}>W/2$ for distinct rows. Indeed the half phase difference is within $\tfrac14N^{-1/4}$ of $\pi\mathbb Z$ after fixing the units, and the axial gap is at least $\arcsin(e^{-d_{ij}/2})$. For genuine points (1) suffices. A repeated allocation in a piece would give the very same point. Thus the Gram matrix $A$ of the rows in one piece has diagonal one and nonpositive off-diagonals. The symmetric entrywise nonnegative matrix $\operatorname{Id}-A$ has largest eigenvalue at most one by $A\ge0$, and its spectral radius equals its largest eigenvalue. Hence $\|A\|\le2$. Stacking the pieces gives $\|V\|^2\le2L$.

Choose five row indices independently and uniformly. The fraction $\eta$ of layer weight constant on all five is
$$
 \eta=2^{-4}\sum_{J\subset[5],\,|J|\ {\rm even}} f_J,\qquad
 f_J=\sum_l(w_l/W)\prod_{j\in J}\sigma_{i_jl}.
$$
This follows by adding the two indicators for all plus or all minus at a layer. For $|J|=2,4$, independence gives
$$
 \mathbb E f_J^2=\sum_{l,l'}\frac{w_lw_{l'}}{W^2}
 \left(\frac1m\sum_i\sigma_{il}\sigma_{il'}\right)^{|J|}
 \le \frac{\operatorname{tr}((VV^t)^2)}{m^2}\le\frac{2L}{m}.
$$
The trace bound uses $\operatorname{tr}(VV^t)=m$. The empty term is one, and Chebyshev on the finitely many remaining terms proves the assertion since $\log H_S/W=1-\eta$. This proof keeps pairs drawn from different pieces. $\square$

**PROVED — corrected grid conclusion.** Choose $0<\epsilon<47/50-15/16=1/400$. Then (II.2) implies (II.1) at $g=94/25$ for large $N$. By §2, all the qualifying index-at-least-seven tuples lie on a single proper closed $Z$. In split coordinates clear the Laurent denominators in any nonzero equation of $Z$ and then substitute $z_i/z_5$. Clearing the resulting power of $z_5$ gives a nonzero polynomial in five variables of some fixed total degree $D$ (specialize $z_5=1$ to check nonvanishing). It has at most $D m^4$ zeros on the grid, by Schwartz–Zippel. Thus the **total** number of such indices is at most $(D+O_C(1))m^4$, including the atypical gcds. Counting just those satisfying (II.1) costs only $D m^4$. At index at least nine the latter count applies to all five-tuples. For squareclass dependence the total is $O_C(m^4)$ directly by (II.2) and the gcd bound. A fixed positive proportion in any of these counts bounds $m$ uniformly in the proportion and $C$. Having just more than $D m^4$ tuples of unspecified height at index seven is not that criterion.

## II.2. The index in a growth obstruction

**PROVED — arbitrary fixed arc constant, all unit classes.** Suppose $W\le A m\log m$, as on a putative obstruction to improving the known growth order. All but $O_{A,C}(m^4)$ ordered five-tuples have full rank and index **exactly one**. This holds for genuine clusters and for gap-admissible angle models. It rules out a positive proportion at index $\ge5$, or even index, in any such family including versions with prime powers. The parallel audit now also gives $O_A(m^3)$ on its normalized short one-unit arc, with no parity exceptions there. The following proof gives the conclusion on the original grid without imposing that normalization on the five chosen marks. [Parallel audit, Part II.4](../free/answer.md)

*Proof.* For any $s$ rows lying in a single piece of the preceding proof, write $W_B$ for the weighted sum of their coordinate ranges. Each threshold splits at most $s^2/4$ unordered pairs, so $\binom s2 W/2< s^2 W_B/4$ if $s>1$. At any given prime at least half the rows lie in one of the two closed halves of $[0,e_p]$. Among these take the largest piece, so $s\ge m/(2L)$ and $W_B\le W-e_p\log p/2$. It follows that $e_p\log p\le4LW/m$ for $m$ large. All primes $\le\sqrt m$ together contribute at most $4LW/\sqrt m=o(W)$. At larger primes $e_p\le8LA$.

Discard the $O_C(m^4)$ tuples with $W_S=\log H_S>31W/32$, by (II.2). Also discard those with any pair satisfying $|\langle V_i,V_j\rangle|>1/128$. The trace bound above makes their number $O_C(m^4)$. Each remaining pair has $d_{ij}=(1-\langle V_i,V_j\rangle)W/2$. After deleting the small primes this is still greater than $31W/64$ for large $m$. The weighted range of the retained five rows on the large primes is at most $31W/32$. These five rows are affinely independent there. Otherwise separate a real affine relation into probability distributions on disjoint supports with the same barycentre. Independent samples from the two sides in a coordinate rescaled to $[0,1]$ give $\mathbb E|X-Y|\le\mathbb E(X+Y-2XY)\le1/2$. Some cross pair then has distance at most half the weighted range, a contradiction. In particular some four-column minor is nonzero and has absolute value at most $24\lceil8LA\rceil^4=:D_0$. Every prime divisor of the index is at most $D_0$.

For fixed prime $\ell$, every residue class of allocation rows modulo $\ell$ contains $O_C(\ell)$ rows. In one piece write $2\langle a_i,\phi\rangle=\theta+\delta_i+2\pi n_i$ with $\delta_i$ in an interval of length $\tfrac12 N^{-1/4}$. Two matching rows in the same winding class $n_i\bmod\ell$ give $v=(a_i-a_j)/\ell$ integral of weight $\le W/\ell$ and $\operatorname{dist}(\langle v,\phi\rangle,\pi\mathbb Z)\le N^{-1/4}/(4\ell)$. The gap excludes $v\ne0$ for every $\ell\ge2$. Equality of the rows in a piece is also excluded. There are at most $L\ell$ choices. The proof for genuine points may take their actual prime angles. In an affinely dependent ordered five-tuple modulo $\ell$, the first dependent row lies in a previous affine span of dimension at most three, hence with at most $\ell^3$ residues. Fix the other four slots to obtain $O_C(\ell^4m^4)$ tuples. Sum over the finitely many $\ell\le D_0$. This also gives the even-index bound without the growth hypothesis. $\square$

**PROVED — two qualifications.** For genuine clusters without any growth hypothesis, combine the index-at-least-seven count with the fixed-prime counts for $2,3,5$ and (II.2) and the typical pair-distance argument for real rank (using all primes). The same $O_C(m^4)$ upper bound holds for all failures of index one, using §2. For a squarefree parent every anchored column can separately be signed into $\{0,1\}^4$. A binary $4\times4$ matrix embeds in a $5\times5$ sign matrix with determinant $16$ times its own, by bordering with ones and subtracting rows. Hadamard gives $|\det|\le5^{5/2}/16<4$. Thus *every* full-rank squarefree index is in $\{1,2,3\}$, without any arc assumption. For a squarefree example of even index, use anchor row zero and four binary columns with supports $12,13,14,23$ among the first four rows. Their determinant has absolute value two. Four distinct split primes realize these columns, with no critical-arc assertion. Prime powers can create large formal indices (doubling all columns of an index-one example gives index $16$). The gap excludes positive density of these in growth-scale clusters. This argument does not exclude an isolated high-index tuple with prime powers.

## II.3. Which covers can lift index-one points?

**PROVED — boundary branching.** Over $\mathbb C$, every connected finite étale cover of the torus is an isogeny after choosing a point above its identity. Indeed the finite-index subgroup of $\pi_1((\mathbb C^*)^4)=\mathbb Z^4$ determines a monomial isogeny, and the equivalence for finite covers identifies them algebraically. Over $\mathbb Q$, if the source has a rational point, translate its image to the identity. The pointed group law descends because the lift of multiplication taking the two chosen origins to the origin is unique, and the same holds for inverse. Consequently a geometrically connected cover of the compactification étale above $T$ gives only a translated isogeny there. On the cocharacter lattices conjugation acts by $-1$ also on the finite-index source lattice. Outside the support of the translating point, rational liftability requires each prime column to belong to the image sublattice. For the identity translate this excludes $I=1$ in every nontrivial degree. A fixed nonunit translate allows exceptions at its finitely many primes, as in §2.3. The assertion about divisibility by a ramification index applies at a good integral model and on the particular sheet met by a lifted section. A cover can have unramified sheets above a branch divisor.

**PROVED — no fixed finite list of covers lifts the unrestricted ten-prime profile.** More precisely, the rational five-point ratios obtained from ten distinct split prime norms, one per pair, are not thin in the four-torus. This holds even with all ten norms in $[X^2/2,X^2]$ as $X$ varies. Thus allowing branch loci in the torus, including the specified exceptional hypersurfaces, cannot provide a universal lift from the squarefree profile alone.

*Proof.* Put $K=\mathbb Q(i)$. Start with ten independent variables $B_{ij}$ in $\operatorname{Res}_{K/\mathbb Q}\mathbb G_m$, hence with an open subset of $\mathbb A^{20}$. The map to $T$ first sends each $B$ to $B/\bar B$ and then applies the ten pair columns in the anchored coordinates. The first homomorphism has geometrically connected kernel in each factor. The second has geometrically connected kernel because its column lattice is all of $\mathbb Z^4$. To see saturation intrinsically, use $\mathbb Z^5/\mathbb Z(1,1,1,1,1)$. Differences of two incident edges give all $e_i-e_j$, and an edge has sum two, generating the remaining quotient $\mathbb Z/5$. Thus the composite is geometrically surjective with geometrically connected torus fibres (over a split field the lattice map splits). Its function field extension is regular. Every fixed generically finite dominant geometrically integral cover of degree $>1$ therefore pulls back generically to an integral cover of the same degree, also geometrically integral. Removing the exceptional loci costs only proper closed subsets. The inverse image of any specified thin set is thin in this affine parameter space.

The quantitative Hilbert irreducibility theorem bounds the number of integral points of any fixed thin subset of $\mathbb A^n(\mathbb Q)$ in the box of side $2X$ by $O(X^{n-1/2}\log X)$ with an implied constant that can depend on the thin set. On the other hand the prime number theorem for $1\bmod4$ supplies $\gg X^2/\log X$ distinct primes in the indicated interval. For each, choose one Gaussian generator, which lies in the box. Ten ordered distinct choices give $\gg X^{20}/(\log X)^{10}$ parameter points, exceeding the bound. They give primitive genuine tuples by taking $z_i=(\prod_e\bar B_e)P_i/\bar P_i$. The five marks are distinct since the supports distinguish them. This proves the assertion. No angle constraint on these points has been imposed. [S. D. Cohen, *The distribution of Galois groups and Hilbert's irreducibility theorem*, Proc. London Math. Soc. 43 (1981), 227–250](https://doi.org/10.1112/plms/s3-43.2.227) $\square$

**PROVED — a cubic cover that does lift.** In homogeneous torus coordinates $u_i$ modulo a common multiplier, set
$$
 J(u)=\frac{1}{(2i)^{10}}\,
       \frac{\prod_{i<j}(u_i-u_j)}{(\prod_i u_i)^2}.
 \tag{II.3}
$$
This is defined over $\mathbb Q$ on the norm-one twist because simultaneous inversion leaves it unchanged, since there are ten sign changes and every index occurs four times. Evaluate at $u_i=P_i/\bar P_i$. Then $u_i-u_j=-2i\,t_{ij}n_{ij}/(\bar P_i\bar P_j)$. Each edge occurs twice in $\prod_i|P_i|^2$ and $H=\prod_e n_e$. Cancellation in (II.3) gives the exact identity
$$
 J(u)=\frac{\prod_{i<j}t_{ij}}{H^3}. \tag{II.4}
$$
For any fixed rational cube class $b$, the normalization in $v^3=J/b$ is a degree-three cover branched in the torus along the ten diagonals, each with simple branch divisor downstairs. If $\prod t_{ij}=b r^3$ with $r\in\mathbb Q$, the point has the rational lift $v=r/H$. A fixed bound on all $|t_{ij}|$ needs only finitely many such covers. The formula includes squarefree tuples of index one and uses arithmetic of the cofactors beyond contact at the boundary.

**PROVED — why the smooth-target calculation does not repeat.** At the identity the branch equation has order ten in four local parameters. The hypersurface $v^3-J/b=0$ has multiplicity three there and is already normal near the identity by the hypersurface criterion. It is smooth off $v=0$ and off the singular locus of its reduced branch divisor within $v=0$, which has codimension at least two upstairs. There is one point above the identity and $|v|=O_b(\Delta^{10/3})$. The elementary McKinnon–Roth volume/multiplicity bound on the pullback of $B$ has the same ratio $(3\cdot70)/3=70$ as downstairs. This leaves open a finer calculation of the approximation constant at this singular point. I have no lift to a specified cover branched on the fifteen exceptional hypersurfaces for the shrinking-arc points.

**PROVED — uniqueness for a limited divisibility calculation.** Over the split torus an invertible rational function away from the diagonals is, up to a constant, a character times powers of the differences. In five homogeneous coordinates write it as $\prod u_i^{a_i}\prod_{i<j}(u_i-u_j)^{v_{ij}}$, with $\sum a_i+\sum v_{ij}=0$ and integral exponents. This follows by unique factorization of Laurent polynomials in four anchored coordinates. After substitution, factor out the indicated powers of the signed cofactors and constants. The exponents of $\beta_{ij},\bar\beta_{ij}$ are then respectively $v_{ij}+a_i+a_j$ and $v_{ij}-s_i-s_j-a_i-a_j$, where $s_i=\sum_{j\ne i}v_{ij}$. To force a $k$-th power by cancellation in each block formally, both must vanish modulo $k$. Write $A=\sum a_i$. The first condition gives $v_{ij}=-a_i-a_j$ and $s_i=-3a_i-A$ modulo $k$. The second gives $a_i+a_j+2A=0$ for every pair. All $a_i=h$, and $12h=0$. Homogeneity gives $15h=0$. Thus $3h=0$ and $v_{ij}=h$. The exponent pattern is generated modulo powers by (II.3) (with $a_i=-2,v_{ij}=1$). In particular the cubic is the only primitive nontrivial possibility in this formal test. This puts no restriction on lifts forced by further additive equations or specializations.

## II.4. A power-size cofactor on an integer subfamily

**PROVED — six norms in a short interval.** There are absolute $c_0,c_1>0$ and $X_0$ such that if integers $x_\nu\in[X,X+c_0X^{1/6}]$, $1\le\nu\le6$, with $X\ge X_0$, satisfy
$$
 \prod_{\nu=1}^6 x_\nu=s^2+t^2,\qquad s,t\in\mathbb Z,\quad t\ne0 ,
$$
then $|t|\ge c_1\sqrt X$.

*Proof.* Take $X$ to be the minimum integer and write $x_\nu=X+k_\nu$, $0\le k_\nu\le K$, one $k_\nu=0$. Replacing a real left endpoint by the minimum only relaxes the window. Use $K\ge1$, enlarging to one if needed. Put $F(Y)=\prod(Y+k_\nu)$ and let $e_j$ be the elementary symmetric functions of the six $k_\nu$. The cubic
$$
 Q(Y)=Y^3+\tfrac12e_1Y^2+
    (\tfrac12e_2-\tfrac18e_1^2)Y+
    \tfrac12e_3-\tfrac14 e_1e_2+\tfrac1{16}e_1^3
$$
has $16Q\in\mathbb Z[Y]$, and $F-Q^2$ has degree at most two. Its coefficients in degrees $2,1,0$ are $O(K^4),O(K^5),O(K^6)$ by homogeneity. Assume $|t|\le c_1\sqrt X$ and take $s\ge0$. Then $Q(X)\asymp X^3$, and
$$
 |s-Q(X)|\le\frac{|F(X)-Q(X)^2|+t^2}{s+Q(X)}
 \ll \frac{K^4}{X}+\frac{K^5}{X^2}+\frac{K^6+t^2}{X^3}<\frac1{16}
$$
for large $X$. The denominator bound gives equality $s=Q(X)$. Since $X\mid F(X)$, we now have $X\mid (16Q(0))^2+(16t)^2$. This integer is positive and $\ll K^6+256t^2$. Take $c_0,c_1$ sufficiently small that it is $<X$, a contradiction. $\square$

**PROVED — application to the actual five points.** Write $[A,B]=\operatorname{Im}(\bar A B)$, $m_{ij}=[P_i,P_j]=t_{ij}n_{ij}$, and $s_{ij}=\operatorname{Re}(\bar P_i P_j)/n_{ij}$. Both quotients are integers. Then
$$
 s_{ij}^2+t_{ij}^2=D_{ij}:=\prod_{k\ne i,j}n_{ik}n_{jk},\qquad
 |u_i-u_j|=2|t_{ij}|/\sqrt{D_{ij}} . \tag{II.5}
$$
The six factors are exactly the six cuts separating the pair. On the stated window, $\sqrt{D_{ij}}\asymp X^3$, so $|u_i-u_j|\gg X^{-5/2}$. With the other four block norms comparable to $X$, $H\asymp X^{10}$. This proves the exponent $4$ on this subfamily. The short interval is additive. Equal logarithmic weights up to $o(\log X)$, or even a common dyadic interval, are much weaker.

**PROVED — a congruence which does use the shared blocks.** For distinct $i,j,k$ write $l,m$ for the other two marks. The two-dimensional vector identity $m_{jk}P_i-m_{ik}P_j=-m_{ij}P_k$ gives, by taking squared norms and cancelling $n_{ij}n_{ik}n_{jk}$,
$$
 2 t_{ik}t_{jk}s_{ij}
 =t_{jk}^2n_{jk}n_{il}n_{im}
 +t_{ik}^2n_{ik}n_{jl}n_{jm}
 -t_{ij}^2n_{ij}n_{kl}n_{km}=:A_{ijk}(n). \tag{II.6}
$$
Take $b$ to be any positive divisor of $D_{ij}$. If each of the ten $n_e$ has a signed residue $r_e\bmod b$ of size at most $K$, and $T=\max|t_{ab}|$, then (II.5) and (II.6) give
$$
 b\mid A_{ijk}(r)^2+4(t_{ij}t_{ik}t_{jk})^2,\qquad
 b\le9 T^4K^6+4T^6 . \tag{II.7}
$$
The last inequality uses positivity, since all $t_{ab}\ne0$. For example, if $K\le b^\delta$, $\delta<1/6$, it forces $T\gg b^{\min(1/6,(1-6\delta)/4)}$. This applies to separated block sizes close to arbitrary multiples of $b$ as well. A triangle suffices, so it is the plane-vector relation underlying (8) that supplies this extra congruence. Neither window nor congruence condition follows from the cut balance in the growth problem.

## II.5. The polynomial equations

**PROVED — normalization and specialization.** Use coefficientwise conjugation on $\mathbb C[x]$ with $x$ real. In the simple pure profile require the ten $\beta_e$ of the same positive degree $d$ and their conjugates to have pairwise disjoint root support across the twenty polynomials. Repetition within one is allowed. Require $[P_i,P_j]/n_{ij}=t_{ij}$ to be nonzero real constants. Nonconstant common support between any two of the twenty factors changes the pure simple profile or cancels some of its contact. One can assume all blocks monic. In fact the vanishing of the leading term of $[P_i,P_j]$ forces all five leading coefficients of the stars to have the same argument modulo $\pi$. Dividing each block by its leading coefficient thus scales each reduced determinant by a nonzero *real* constant. Conversely monic solutions suffice for the existence problem.

For a solution over $\mathbb Q(i)[x]$, clear the denominators by a common positive integer multiplier of the blocks. Specialize $x$ to large positive integers and set $z_i=(\prod_e\bar\beta_e)P_i/\bar P_i$, which are Gaussian integers. The primitive norm is $\asymp x^{20d}$. Indeed all pairwise polynomial gcds among the twenty factors are trivial, so all pairwise gcds of distinct factors' values divide fixed nonzero Gaussian integers by Bézout. At a Gaussian prime at most one factor value can have valuation above a fixed bound supported on the finitely many resultant primes, and each factor is absent from some $z_i$. This bounds the common gcd uniformly. Likewise reduction of any pair can only reduce the absolute scalar cofactor by bounded factors, because the unreduced six-factor word already has nonzero constant imaginary part. Equation (II.5) gives $\Delta\asymp x^{-6d}$. Hence the primitive tuples attain exponent $10/3$ on this curve, up to constants. A calculation over $\mathbb C$ alone gives no such rational specialization. Also one curve outside today's known exceptional components can itself be put into a fixed exceptional set. To determine an **essential** approximation constant, one needs such approximants outside *every* fixed proper set.

**PROVED — what the literal polynomial wording permits.** With no coprimality assumed, $t_{ij}\ne0$ already forbids a real root in any block, overlap between a block and its own conjugate, and any common root of two different oriented blocks. In the last case choose $i\in e\setminus f$, $j\in f\setminus e$. Then the reduced word $B_{ij}$ and its conjugate would both vanish at that root, although they differ by the nonzero constant $2it_{ij}$. The same argument excludes conjugate overlap of *incident* blocks by taking a shared vertex and a vertex outside their union. It excludes a block's self-overlap by taking any pair it separates. Opposite-conjugate overlap on **disjoint** edges need not be excluded by this argument. Such a prime belongs to two layers, outside the private, pairwise coprime norm family used for II.6–II.7. The no-pencil proof below works even for the literal equations. I have no unrestricted existence decision under either version.

**PROVED — equivalent moments.** The reduced word $B_{ij}=\bar P_i P_j/n_{ij}$ is monic of degree $6d$ and $\operatorname{Im}B_{ij}=t_{ij}$. Writing $r$ for block roots counted with multiplicities, the condition that its nonconstant coefficients be real is equivalent to
$$
 \operatorname{Im}\sum_{j\ne i}\ \sum_{r\in\beta_{ij}} r^k
 \quad\hbox{independent of }i,\qquad 1\le k<6d . \tag{II.8}
$$
Indeed Newton's triangular identities say that the first $6d-1$ coefficients below the leading term are real exactly when the corresponding power sums are real. In the reduced word the shared edge cancels from the difference of the two star sums. Nonzero $t_{ij}$ remains a separate condition. The real unknowns here are $20d$ root coordinates, and the $4(6d-1)$ scalar equations can have dependent or degenerate components. Counting them does not decide existence.

**PROVED — no pencil of norm polynomials in any degree.** In a solution the ten $n_e=\beta_e\bar\beta_e$ span a real vector space of dimension at least three. In particular, for linear blocks they must span the full space of real quadratics.

*Proof.* Fix one triangle. Eliminating $s_{ij}$ from (II.5) and (II.6) yields a homogeneous sextic
$$
 4t_{ik}^2 t_{jk}^2\prod_{h\ne i,j}n_{ih}n_{jh}
        - A_{ijk}(n)^2 = 4(t_{ij}t_{ik}t_{jk})^2>0 . \tag{II.9}
$$
It cannot specialize to this constant on a one-dimensional span of nonconstant monic polynomials. If the span has dimension two with real basis $p,q$, the left is a nonzero binary sextic in $p,q$. Over $\mathbb C$ it factors into linear forms. Every factor evaluated at $p(x),q(x)$ must be a nonzero constant, since their product is a nonzero constant in $\mathbb C[x]$. Taking real or imaginary parts shows that $1$ belongs to the real span. All $n_e$ are monic of the same degree, so $n_e=G+c_e$ for one nonconstant monic real $G$ and real constants $c_e$. Now (II.6) expresses $s_{ij}=S(G)$ with $S$ a real polynomial. Substitution in (II.5) gives the polynomial identity
$$
 \prod_{h\ne i,j}(Y+c_{ih})(Y+c_{jh})=S(Y)^2+t_{ij}^2 .
$$
At a root of a factor on the left the right is strictly positive. This contradiction proves the claim. More generally, if all $n_e=F_e(G)$ for one nonconstant real $G$ and real polynomials $F_e$, (II.5) and (II.6) show that *none* of the $F_e$ can have a real zero. $\square$

**PROVED — reconstruction from four blocks with exact saturation conditions.** Choose monic $\gamma_i=\beta_{i5}$, $1\le i\le4$, degree $d$, disjoint as above. Set $U_i=\prod_{h\le4,h\ne i}\gamma_h$. There is a unique complex polynomial $D_i$ of degree $<3d$ satisfying $[U_i,D_i]=1$. In fact the real-linear map $D\mapsto[U_i,D]$ has target the real polynomials of degree $<6d$, of the same real dimension as its domain. Its kernel is zero since $\bar U_iD=U_i\bar D$ and $(U_i,\bar U_i)=1$ imply $U_i\mid D$. This also constructs $D_i$ over $\mathbb Q(i)$ for rational Gaussian input.

Choose four nonzero real constants $h_i$, pairwise distinct, and set
$$
 C_i=U_i+h_iD_i .
$$
The entire problem now is to arrange that each of the six monic gcds $E_{ij}=(C_i,C_j)$ has degree $d$ and that the ten blocks $\gamma_i,E_{ij}$ and their conjugates satisfy the original disjointness condition. These conditions are necessary and sufficient. Indeed $P_5=\gamma_i U_i$, $P_i=\gamma_i C_i$ in the reconstruction. In any solution $[U_i,C_i]=t_{5i}$ and $C_i-U_i$ has degree $<3d$, proving necessity of the formula with $h_i=t_{5i}$. Conversely the disjoint gcds meeting $i$ exhaust the monic $C_i$ of degree $3d$. For real $x\to+\infty$,
$$
 \operatorname{Im}(C_i/U_i)=h_i/|U_i|^2,\qquad
 \arg(P_i/P_5)=h_i x^{-6d}+O(x^{-6d-1}).
$$
Thus $[P_i,P_j]/n_{ij}$, a polynomial by the shared factor, tends to $h_j-h_i$. It is identically that constant. In any solution the same computation also proves $h_i\ne h_j$ since $t_{ij}\ne0$. For $d=1$ this leaves four complex base roots, four real $h_i$, and six cubic-resultant equations over $\mathbb C$, with the stated gcd, reality and noncollision tests. This is a lower-dimensional exact elimination problem. A zero resultant alone admits extraneous components. Nor does an answer at $d=1$ settle arbitrary $d$.

**PROVED — a real-geometric exclusion for linear blocks.** The ten distinct nonreal roots of a degree-one solution cannot all lie on any generalized circle (a circle or a line) in the complex plane. Write $G_{ij}=B_{ij}=R_{ij}+it_{ij}$ with $R_{ij}$ real monic. Its six roots use each edge adjacent to $ij$ once, three with conjugation and three without. Let $M$ be the adjacency matrix of the line graph of $K_5$, indexed by its ten edges. It is invertible and $6M^{-1}$ is integral. In fact the latter has entries $1$ on the diagonal and on adjacent pairs and $-2$ on disjoint pairs, as direct multiplication shows.

*Proof.* If a nonhorizontal line contains the roots, let $a$ be its intersection with the real axis. At $a$ the three factors of either orientation in $G_{ij}$ have cancelling arguments modulo $\pi$, so $G_{ij}(a)$ is real, impossible. For a horizontal line of nonzero height $b$ the six roots of one $G_{ij}$ have heights $b,b,b,-b,-b,-b$. Reality of power sums of degrees $2,3,4$ implies equality of the power sums of degrees $1,2,3$ of the abscissas in the two triples because in $\operatorname{Im}(x+ib)^k$ the leading term is $kb x^{k-1}$ and the lower terms have already cancelled. Newton gives equality of the abscissa multisets, contrary to disjointness among the original blocks.

If a circle meets the real axis at distinct $a,b$, the numbers $w_e=(a-r_e)/(b-r_e)$, for the ten roots $r_e$, have common argument modulo $\pi$. Thus $G_{ij}(a)/G_{ij}(b)$ is real. Both evaluations have the same nonzero imaginary part, so their ratio is one. Taking logarithms of absolute values gives $M(\log|w_e|)_e=0$. All roots then lie on the perpendicular bisector of $a,b$, a line already excluded. At tangency $a$, the numbers $1/(a-r_e)$ have constant imaginary part. Hence $G'_{ij}(a)/G_{ij}(a)$ is real. Since $G'_{ij}(a)$ is real and $G_{ij}(a)$ is nonreal, it is zero. Taking real parts and inverting $M$ gives $\operatorname{Re}(1/(a-r_e))=0$ for every edge, again a line.

Finally suppose the circle has centre $u+iv$ and radius $\rho<|v|$. Set $a=u+i\sqrt{v^2-\rho^2}$ and $b=\bar a$. The circle is a locus of constant $|(a-r)/(b-r)|=k>0$, as follows by squaring (neither limiting point lies on the circle). For a conjugate root the ratio is $1/\bar w_e$. Thus $|G_{ij}(a)|=|G_{ij}(b)|$ by the three factors of each orientation. With $G_{ij}=R_{ij}+it_{ij}$ this forces $\operatorname{Im}R_{ij}(a)=0$, so $G_{ij}(a)=G_{ij}(b)$. Taking arguments gives $M(\arg w_e)_e=0\bmod 2\pi$. All $w_e/k$ are sixth roots of unity because $6M^{-1}$ is integral. A Möbius transformation cannot take ten distinct roots into these six values. $\square$

**COMPUTED — local obstruction for linear blocks.** There is no solution of the first five moment equations with ten roots in $\mathbb F_{49}$ whose twenty roots and conjugates are distinct. Consequently any rational Gaussian degree-one solution has bad reduction at $7$ after *every* real rational affine change of variable for which the roots are integral at $7$. This is a restriction on reduction, not a $7$-adic nonsolvability claim. Identify $\mathbb F_{49}$ with $\mathbb F_7[i]$. Enumerate a possible four-root star by choosing four of the 21 nonreal conjugate pairs and one of the two signs for each. Store its five sums of imaginary parts of powers $1,\ldots,5$. There are 95,760 stars in 16,660 classes, maximum size 14. Within each class the script tries to select five stars with intersection size one for each pair, no element in three stars, and consistent conjugate choices. That describes exactly an unlabelled copy of the required edge incidence because the ten pair intersections consume all twenty star incidences. The count is zero. The script also gives zero if conjugate collisions on disjoint edges are allowed while oriented roots remain distinct and each star avoids its conjugates. [Complete finite-field enumeration](tenpair_modular.py)

**PROVED and COMPUTED — an exact slice of the linear-block moment test.** If the first two *complex* star power sums are constant across the five vertices, the ten edge roots take at most four different values. Thus this slice of (II.8) is empty with nonzero cofactors, even allowing conjugate overlaps on disjoint edges.

*Proof.* Subtract the common first star sum divided by four from every root. The five star sums are now zero, and the square star sums are still equal. The first condition gives a five-dimensional complex vector space $V$. Indeed the unsigned vertex-edge incidence matrix has rank five, since a dependency would have $a_i+a_j=0$ for all pairs. The four homogeneous quadratic equalities cut out a subset of $\mathbb P(V)=\mathbb P^4$. There are six projective points obtained by putting $+1$ on a Hamilton cycle and $-1$ on its complementary cycle (12 cycles, identified in pairs). A further ten points come from isolating one vertex with its four values zero and putting $1,\omega,\omega^2$ on the three perfect matchings of the remaining $K_4$, where $\omega$ is a primitive cube root of unity. For each isolated vertex there are two assignments modulo scaling. Both constructions satisfy the equations. At a representative of each orbit a $9\times9$ Jacobian minor of the five linear and four quadratic equations, omitting a factor two from each quadratic row, is respectively $32$ and $16\omega+8\ne0$. The script specifies the columns and computes the determinants exactly. Permutation and scaling preserve rank. Thus all sixteen projective points are isolated. The sum of the degrees of *all* irreducible components of the common zero locus is at most $2^4$ by successive intersection. Intersect successively with the four quadrics, retaining contained components and applying projective Bézout to every other component. The sixteen points exhaust that bound. The affine origin is also allowed before restoring the shift, and no case gives ten distinct oriented roots. [Exact Jacobian certificates](tenpair_two_moments.py) $\square$

**PROVED and COMPUTED — one determinant in the ten pole residues.** For disjoint linear blocks put $q_i=P_i/\bar P_i$, $D=\prod_e\bar\beta_e$ and $F=(D(q_i-q_5))_{i\le4}$. Every component has degree at most four, since
$$
 D(q_i-q_j)=-2it_{ij}\,\beta_{ij}\!\prod_{f\cap\{i,j\}=\varnothing}\bar\beta_f .
$$
The ten quadratic monomials of $F$ lie in the nine-dimensional space of degree-at-most-eight polynomials. Thus some nonzero homogeneous quadric vanishes on the entire image. Put $v_i$ equal to the four coordinate vectors for $i\le4$ and $v_5=-(1,1,1,1)^t$. At the root of $\bar\beta_{ij}$ just the two terms with indices $i,j$ survive in $F=\sum_i v_iDq_i$. Both are nonzero, so $F$ there represents $v_i+h_{ij}v_j$, where $h_{ij}$ is the ratio of the two residues of $q_j,q_i$. Evaluating all ten quadratic monomials at these ten vectors gives a square matrix whose determinant must vanish. With every $h_{ij}=1$ its determinant is $-3$ in lexicographic order. Consequently at least one $h_{ij}$ stays outside a fixed complex neighbourhood of one. At every nonarchimedean place of residue characteristic other than three they cannot all be integral and congruent to one. For example, if at each pole all six other factors $\beta_f/\bar\beta_f$ occurring in $q_j/q_i$ are units congruent to one, we obtain a contradiction even if the twenty roots had collisions on reduction. This test, and the good-reduction test at seven, still allow singular residue configurations. [Exact quadric determinant](tenpair_quadric.py)

**PROVED — one fixed rational plane contains approximants outside every fixed closed set.** Fix ten distinct positive offsets, for example $1,\ldots,10$. Multiply all offsets by a single integer $B$ divisible by $2$ and by every prime dividing an original offset or a nonzero difference of their squares. Denote the new offsets by $c_e$. The plane of norm vectors spanned by $\mathbf1$ and $(c_e^2)_e$ is now fixed for this whole assertion. We can arrange $I=1$, all $n_e=k^2+c_e^2$ in a bounded additive window, and $u\to e$ outside any prescribed proper closed subset $Z$ of the four-torus.

First, for $k\equiv1\bmod B$ the norms are pairwise coprime and odd, and $\gcd(k,c_e)=1$. Indeed a prime shared by two norms must divide a difference of squared offsets, hence divides $B$, while all the norms equal one modulo $B$. We allow an independent rational rotation of each Gaussian block $k+ic_e$. Put $\gamma_e=m_e+i$, where each positive $m_e$ is an even multiple of $B$ and of all the earlier $|\gamma_f|^2$ in some fixed ordering. Thus the $|\gamma_e|^2$ are odd and pairwise coprime, avoid the primes of $B$, and each $\gamma_e$ is coprime to its conjugate. Use the new block
$$
 \tilde\beta_e=(k+ic_e)\gamma_e/\bar\gamma_e .
$$
For every Gaussian prime $\lambda\mid\bar\gamma_e$ require $v_\lambda(k+ic_e)=v_\lambda(\bar\gamma_e)$. This is an ordinary congruence with a choice of nonzero next digit for $k$: at a split prime over $p$, the rational integers map isomorphically modulo $p^a$ to $\mathbb Z[i]/(\lambda^a)$ for all $a$. The conjugate valuation is zero since $p\nmid2c_e$. These conditions make each $\tilde\beta_e$ integral and coprime to its conjugate by transferring the entire selected prime power to the opposite orientation. All norms stay fixed as functions of $k$. Also choose one new split prime per edge, outside all preceding moduli, and impose valuation exactly one in $k^2+c_e^2$ there using a simple root modulo that prime. The Chinese remainder theorem realizes all conditions together with $k\equiv1\bmod B$. The dedicated primitive columns give $I=1$.

As $k$ goes to infinity along this progression, the torus point tends to the pair-column image of $\big((\gamma_e/\bar\gamma_e)^2\big)_e$. Given any neighbourhood of the identity, take every $m_e$ sufficiently large to put this limit in that neighbourhood. The possible limits with this restriction are Zariski dense. Each admissible prefix of the $m_e$ has infinitely many possible next entries, so induction on polynomial coefficients gives Zariski density of these parameters in $\mathbb A^{10}$. The map to $T$ is dominant by the rational parametrization of the norm-one torus, squaring, and the surjectivity of the pair-column map in II.3. We may therefore choose the limit outside $Z$. Then all sufficiently large $k$ also miss $Z$. Shrinking the neighbourhoods proves the assertion with unbounded height on this *one* plane. This construction gives approaching points, without asserting sharpness of any of the lower bounds or squarefreeness of these block norms.

## II.6. Finitely many bounded-cofactor tuples on each fixed rational plane of norm vectors

**PROVED.** Fix a rational subspace \(U\subset\mathbb Q^{10}\) of dimension at most two and a constant \(T\). In the genuine ten-pair-cut model of Part I, equation (8), only finitely many private pairwise coprime block-norm vectors \(n=(n_{ij})\in U\) have \(0<|t_{ij}|\le T\) for all pairs. This allows the norms to grow at unrelated rates a priori. The finiteness assertion is individual to \(U,T\).

### From the cofactor bounds to a conic

**PROVED.** Write \(P_i=\prod_{j\ne i}\beta_{ij}\), \(n_{ij}=|\beta_{ij}|^2\), \(t_{ij}=\det(P_i,P_j)/n_{ij}\), and \(s_{ij}=\operatorname{Re}(\bar P_iP_j)/n_{ij}\). Both quotients are integral. The norm gives
\[
 s_{ij}^{\,2}+t_{ij}^{\,2}=\prod_{a\notin\{i,j\}}n_{ia}n_{ja}.
 \tag{II.10}
\]
For distinct \(i,j,k\) and complementary indices \(l,m\), define
\[
 A_{ijk}=t_{jk}^{\,2}n_{jk}n_{il}n_{im}
 +t_{ik}^{\,2}n_{ik}n_{jl}n_{jm}
 -t_{ij}^{\,2}n_{ij}n_{kl}n_{km}.
\]
Put \(m_{ij}=t_{ij}n_{ij}\). Taking squared absolute values in \(m_{jk}P_i-m_{ik}P_j=-m_{ij}P_k\), using \(|P_i|^2=\prod_{a\ne i}n_{ia}\), and dividing by \(n_{ij}n_{ik}n_{jk}\) gives \(2t_{ik}t_{jk}s_{ij}=A_{ijk}\). In particular \(s_{ij}\), with a fixed choice of \(k\), is a real homogeneous cubic in the \(n\)'s once the \(t\)'s are fixed.

Fix the signed \(t\)'s along a putative infinite subsequence. Dimension one gives finiteness by (II.10) and the homogeneous cubic, since their difference is a homogeneous sextic taking a nonzero constant value. In dimension two choose a lattice basis of \(U\cap\mathbb Z^{10}\). Write \(n_e=L_e(p,q)\) for rational integer linear forms and integer \(p,q\). Substituting for fixed \(i,j,k\) gives one homogeneous rational binary sextic \(F(p,q)=c\ne0\), with \(c=t_{ij}^2\) and \(F=\prod_{a\notin\{i,j\}}L_{ia}L_{ja}-s_{ij}(L)^2\). There are infinitely many projective parameter slopes since \(F\) takes a fixed nonzero value.

After clearing denominators each primitive integer irreducible factor of \(F\) has a fixed nonzero integer value on an infinite subsequence, by divisibility. Any irreducible factor of degree at least three is excluded by Thue's theorem on binary forms. Two independent rational linear factors fix \(p,q\). A linear factor with fixed value confines the parameters to an affine line. On this line all equations (II.10), with their cubic substitutions, hold identically by the infinite subsequence. At least one \(L_{ab}\) is nonconstant since the map into \(U\) is injective on parameters. Evaluating (II.10) for \(i=a\) and \(j\notin\{a,b\}\) at its real root is impossible since \(t_{ij}\ne0\). Finally, two nonproportional irreducible quadratic factors with fixed nonzero values would give a nonzero homogeneous quadratic vanishing on infinitely many slopes by taking a linear combination. A definite quadratic with fixed nonzero value has finitely many integer representations. The only remaining possibility is an indefinite nonsplit rational quadratic \(Q(p,q)=d\ne0\).

**PROVED.** Diagonalize the real conic as \(u v=1\) by a real linear change of the two coordinates. The point \((p,q)\) on the conic is a real linear combination of \(u\) and \(u^{-1}\). Pass to one of its two real branches and take \(u>0\). Every \(n_e\) now has the form \(a_eu+b_e/u\) with \(a_eb_e\ne0\), since the isotropic lines of a nonsplit rational quadratic have irrational slopes and \(L_e\) is rational and nonzero. All the norm identities (II.10) and cubic reconstruction identities extend over the conic by equality at infinitely many parameters. Every \(n_e\) stays positive on \(u>0\). Indeed it is positive at actual tuples, and (II.10) for a pair whose product includes \(n_e\) precludes a real zero. Hence \(a_e,b_e>0\). Distinct forms \(L_e,L_f\) are nonproportional. If they were proportional, write \(n_e/n_f=A/B\) in lowest terms. Coprimality forces \(n_e=A,n_f=B\) at all sample points, impossible for a nonconstant Laurent polynomial taking a fixed value infinitely often. Consequently the ten ratios \(b_e/a_e\) are distinct.

Define \(V_{ij}=(s_{ij}+it_{ij})/(s_{ij}-it_{ij})\), a rational function of \(u\). On every sample point it equals \((P_j/\bar P_j)/(P_i/\bar P_i)\), so \(V_{ij}V_{jk}=V_{ik}\) identically. At zero and infinity, (II.10) shows that \(s_{ij}\) has a pole of order three. Thus \(V_{ij}=1+O(u^3)\) at zero and \(1+O(u^{-3})\) at infinity. The only remaining possible zeros and poles of \(V_{ij}\) are the two roots \(\pm i\sqrt{b_e/a_e}\) for any of the six edges occurring in (II.10). At each such root exactly one of \(s_{ij}\pm it_{ij}\) vanishes simply, since the six edge forms occurring in (II.10) have disjoint simple roots and \(2it_{ij}\ne0\). Real conjugation exchanges the zero and the pole. At a fixed edge root the valuations of \(V_{ij}\) form a cocycle on the five indices. The valuation vanishes for pairs \(ij\) on the same side of the cut given by that edge and is \(+1\) or \(-1\) for pairs on opposite sides. A cocycle with these properties is a signed difference of the indicators of membership in the edge. We can thus choose \(y_e=\pm\sqrt{b_e/a_e}\) independently for each edge so that
\[
 V_{ij}=\prod_e\left(\frac{u-iy_e}{u+iy_e}\right)^{
                  {\bf1}_{j\in e}-{\bf1}_{i\in e}} .
\]
The divisor determines the right side up to its constant, which is fixed by \(V_{ij}(\infty)=1\). Comparing linear terms at infinity and logarithmic derivatives at zero shows that both \(\sum_{e\ni i}y_e\) and \(\sum_{e\ni i}1/y_e\) are independent of \(i\). The following real lemma excludes such a configuration.

### A real five-star lemma

**PROVED.** Suppose \(y_{ij}=y_{ji}\) for the ten edges of \(K_5\) are nonzero reals distinct up to sign. Their star sums and their reciprocal star sums cannot both be constant.

*Proof.* Suppose the star sums are \(s\) and the reciprocal sums \(r\). Then the five monic quartics of the star values are
\[
 P_i(z)=\prod_{j\ne i}(z-y_{ij})
       =z^4-sz^3+c_i z^2-r d_i z+d_i.
\]
If \(r=s=0\), each \(P_i\) is even and its roots cannot be distinct up to sign. Otherwise consider the real rational parametric plane curve
\[
 C(z)=[1-rz:z^4-sz^3:z^2].
\]
The degree of the parameterization as a linear system is four. The homogenized forms have no base point. On the open third-coordinate chart, two coordinates are \(Y=z^2-sz\) and \(X=z^{-2}-rz^{-1}\). Since \(Y\) has degree two, failure of birationality requires its involution \(z\mapsto s-z\) to fix \(X\). The only pole of \(X\) is zero, forcing \(s=0\), then \(r=0\) by invariance. Thus \(C\) is an irreducible plane quartic with normalization \(\mathbb P^1\). At the distinct points \(C(0)\) and \(C(\infty)\) the normalization has only one preimage, as the third coordinate can vanish only at these two parameters. Its differential vanishes at each because both vanishing local coordinates have order at least two. These are two distinct real unibranch singularities.

The five lines \(w_1+c_iw_2+d_iw_0=0\), in coordinates \([w_0:w_1:w_2]\), are distinct. Equality would identify two star quartics with four common roots although their star-value sets intersect in exactly one element. The intersection of lines \(i,j\) is \(C(y_{ij})\). No three lines concur, since that would give a point with three distinct finite nonzero preimages while \(Y=z^2-sz\) has degree two. A nonsingular real conic is tangent to all five lines. In the dual plane the five points with no three collinear impose five independent conditions on conics (four points are covered by two lines avoiding the fifth). The resulting conic cannot be reducible since two lines cannot cover five such points, even over \(\mathbb C\). The dual conic is real with real tangencies, hence is split over \(\mathbb R\).

Identify this plane with \({\rm Sym}^2\mathbb P^1\) over \(\mathbb R\) by assigning a point its unordered pair of tangents to the conic, counted with multiplicity over \(\mathbb C\). In coordinates on the dual conic this is the familiar identification with binary quadratics up to scale. The five lines become divisors of pairs containing one of five distinct real parameters \(h_1,\ldots,h_5\). Choose all \(h_i\) finite. Pull back a real equation of \(C\) to a symmetric polynomial \(D(x,y)\) of bidegree at most \((4,4)\) on the product. Since no one of the lines is contained in \(C\), writing \(f(x)=\prod_i(x-h_i)\) gives \(D(h_i,y)=k_i f(y)/(y-h_i)\) with \(k_i\in\mathbb R^\times\). Let \(g\) be the interpolation polynomial of degree at most four with \(g(h_i)=k_i\). Interpolation in the first coordinate shows
\[
 D(x,y)=\frac{f(x)g(y)-f(y)g(x)}{x-y}.
\]
Indeed both sides have degree at most four in each variable and agree at the five distinct first coordinates \(h_i\). Therefore \(C\) parametrizes unordered pairs in a common fiber of the real degree-five map \(h=f/g\), including a double point on the diagonal at every ramified preimage. The fiber of \(h\) at zero consists of five distinct real points.

We now classify the possible singularities of this plane curve over \(\mathbb C\). At an unordered pair of distinct points in a fiber with local degrees \(m,n\), the map taking ordered to unordered pairs is locally an isomorphism, and the completed analytic local equation is \(X^m-Y^n=0\), after a change of local parameters. This follows by taking a local coordinate on the target and its analytic \(m\)th and \(n\)th roots at the respective sources. The only singular unibranch possibility with \(m+n\le5\) is \(\{m,n\}=\{2,3\}\), an ordinary cusp of delta invariant one. At a diagonal point of local degree \(m\), the equation is \((X^m-Y^m)/(X-Y)\) in the symmetric parameters \(v=XY,w=X+Y\). For \(m=2,3\) the equations \(w=0,w^2-v=0\) are smooth. For \(m=4\) one gets \(w(w^2-2v)=0\). For \(m=5\) one gets \(w^4-3w^2v+v^2=0\), with two branches. The simple diagonal case does not belong to \(C\). It follows that the two unibranch points each force a different ramified fiber of type \((3,2)\). Both critical points in each of these fibers are real because the unordered pair is real and the local degrees distinguish its two members.

The two fibers account for six of the eight ramification units of a degree-five rational map. The genus formula for an irreducible rational plane quartic gives total delta invariant three. After the two cusps there must be another singularity. With ramification budget two, the sole choice that gives any singularity in the local equations above is one fiber of type \((2,2,1)\), giving an ordinary node. The other choices are one triple critical point or two simple critical points over distinct values, neither of which gives a singularity.

Apply real Möbius changes to put the two cusp values at zero and infinity, the triple zero at \(x=0\), the double zero at \(x=\infty\), and the triple pole at \(x=1\). The double pole is \(b\in\mathbb R\setminus\{0,1\}\). The new function, up to real scaling, is
\[
 H(x)=\frac{x^3}{(x-1)^3(x-b)^2}.
\]
The two remaining critical points are distinct roots of \(q(x)=2x^2+x-3b\). Their \(H\)-values coincide. Neither is a pole. Thus the coefficient vectors of the linear remainders of \(x^3\) and \((x-1)^3(x-b)^2\) modulo \(q\) are proportional. Their constant/linear determinant is \(b^2(9b+16)(24b+1)/16\), by polynomial division (see the [exact critical-value calculation](tenpair_belyi_check.py)). The choice \(b=-1/24\) makes the two roots equal. Hence \(b=-16/9\), giving discriminant \(1+24b=-125/3<0\). The two remaining critical points are nonreal.

The only real critical values contributed by *real* source points are now zero and infinity. The number of real inverse images is constant on each interval of \(\mathbb P^1(\mathbb R)\) complementary to these two, since the restriction to the real circle is proper and a local diffeomorphism there. On either side, sufficiently close to zero there is one real inverse near the triple zero and either none or two near the double zero, by their real local degrees. Thus every regular real fiber has at most three real points. This contradicts the original fiber with five distinct real points, a property invariant under both real Möbius changes. This proves the lemma and the fixed-pencil finiteness theorem. \(\square\)

**PROVED — scope.** The preceding argument treats bounded cofactors on each fixed rational plane. The upgrade to growing cofactors is in II.7. The polynomial obstruction for a two-dimensional real span of norms is shorter. A homogeneous binary sextic constant on a nonconstant polynomial parameterization forces a nonzero linear combination of the two parameter polynomials to be constant, by factorization over \(\mathbb C\). With real independent input take real and imaginary parts to obtain a real such combination. On the resulting affine line at least one norm form has a real root in the line parameter. Equation (1) and cubic reconstruction give the same contradiction. This polynomial argument is independent of the real five-star lemma.

## II.7. A power lower bound on each fixed rational plane of norm vectors

**PROVED.** Let \(U\subset\mathbb Q^{10}\) be a fixed rational subspace of dimension at most two. For every one-unit-class ten-pair tuple with pairwise coprime block norms \(n_e\), \(n\in U\), and nonzero signed cofactors \(t_e\) as in Part I, equation (8), one has
\[
       T=\max_e|t_e|\ \gg_U\ \big(\max_e n_e\big)^{1/6}.
\]
If the ten coordinate forms restricted to $U$ are pairwise nonproportional, the cofactor exponent is $1/4$. The implied constant may be ineffective. The proof uses Faltings' theorem for a curve of genus five, Roth's theorem on approximation to a fixed algebraic number, and the real five-star lemma of II.6.

**PROVED — integral parameters.** In dimension one coprimality leaves at most one positive integral vector on the fixed ray. For dimension two choose once and for all a lattice basis of \(U\cap\mathbb Z^{10}\), so \(n_e=L_e(p,q)\) with fixed integral linear forms. Put \(B=\max(|p|,|q|)\). The basis and its left inverse give \(B\asymp_U \max n_e\). Pairwise coprimality implies \(\gcd(p,q)=1\). Each form is nonzero whenever the family has a sample.

If two forms are proportional, their positive values have a fixed rational ratio, and coprimality fixes both values. The integer parameters therefore lie on a fixed affine line and have the form \((p_0,q_0)+v(p_1,q_1)\) for \(v\in\mathbb Z\). The origin and direction of this parameterization can be fixed from \(U\) alone. One edge \(ab\) has nonconstant norm on the line. Fix \(i=a\), \(j\notin\{a,b\}\), and a third index \(k\ne i,j\). Use the integer polynomial
\[
 4t_{ik}^2t_{jk}^2\prod_{h\notin\{i,j\}} n_{ih}n_{jh}
 -A_{ijk}(n)^2-4t_{ik}^2t_{jk}^2t_{ij}^2
 \tag{P1}
\]
from (II.6). This polynomial vanishes at the sample parameter. Its coefficients as a polynomial in \(v\) are \(O_U(T^6)\) and it is nonzero because evaluating at the real zero of \(n_{ab}\) gives a strictly negative number. Its leading nonzero coefficient has absolute value at least one. Cauchy's root bound gives \(B\ll_U 1+|v|\ll_U T^6\). Henceforth all ten forms are pairwise nonproportional.

**PROVED — forcing a line in the Grassmannian.** Write \(M_{ij}(p,q)=t_{ij}L_{ij}(p,q)\), with skew-symmetric extension in the two indices. For each quartet its Plücker relation is a homogeneous quadratic binary form \(E\), with integral coefficients of absolute value \(O_U(T^2)\). At the sample parameters it vanishes since the entries are the minors of \(P_1,\ldots,P_5\). If any \(E\) is a nonzero polynomial then \(B\ll_U T^2\). To see this, remove any powers of the two coordinate variables, apply the rational-root theorem to the resulting polynomial, and use primitivity. On a coordinate axis \(B=1\). We may assume all relations vanish identically in the two parameters.

The coefficient bivectors of \(p,q\) in \(M\) are independent since the forms span the dual parameter plane. Both bivectors are decomposable and their sum is too. Their wedge therefore vanishes, and their rational two-dimensional support spaces meet in a rational line. Choose a generator \(u\in\mathbb Q^5\) of this line. There are rational vectors \(v,w\) such that \(M=u\wedge(pv+qw)\). None of the five \(u_i\) vanishes, since then all the incident entries \(M_{ij}\) would be proportional binary forms. With \(A_i=(v_i/u_i,w_i/u_i)\in\mathbb Q^2\),
\[
       t_{ij}L_{ij}(p,q)=u_i u_j\langle A_j-A_i,(p,q)\rangle .
\]
Thus the ten pair directions of the five affine points are fixed by \(U\) and are all different. In particular no three of these points are collinear. The configuration is uniquely determined up to translation and nonzero homothety. After matching the first two points of any two such configurations by such a transformation, all other points agree since each is the intersection of the two specified nonparallel lines from the first two points. Fix a rational representative \(A_i^0\) of the shape from one solution, if any. With the nonzero rational antisymmetric constants \(c_{ij}\) given by \(\langle A_j^0-A_i^0,(p,q)\rangle=c_{ij}L_{ij}(p,q)\), one has
\[
       t_{ij}=\lambda u_i u_j c_{ij},\qquad \lambda\in\mathbb Q^\times.
 \tag{P2}
\]
**PROVED — finite projective cofactor patterns.** Put \(\ell_i=\langle A_i^0,(p,q)\rangle\) and \(D_i=\prod_{j\ne i}c_{ij}\). The five columns \(W_i=(u_i,u_i\ell_i)\) have all their minors equal to \(1/\lambda\) times the minors of the \(P_i\) at the sample, so the two rank-two matrices differ by one invertible real row transformation (match two independent columns, then take minors against them). Equation (8) of Part I transfers to the \(W_i\). As \(\prod_{j\ne i}t_{ij}=\lambda^4 u_i^3D_i\prod_j u_j\), its two coordinates give
\[
       \sum_i\frac1{u_i^2D_i}=0,\qquad
       \left\langle\sum_i\frac{A_i^0}{u_i^2D_i},(p,q)\right\rangle=0.
 \tag{P3}
\]
For each $i$ choose a perfect matching $\{ab,cd\}$ on the other four vertices and put $m_i=(t_{ab}/c_{ab})(t_{cd}/c_{cd})$. Equation (P2) gives $m_i=\lambda^2(\prod_j u_j)/u_i$, with the *same* numerator for all five $i$. The $m_i$ have bounded denominators depending only on $U$ and size $O_U(T^2)$. The second equation of (P3) becomes $\langle\sum A_i^0m_i^2/D_i,(p,q)\rangle=0$. After clearing fixed denominators its two integral coefficients are $O_U(T^4)$. If either is nonzero then primitivity gives $B\ll_U T^4$. Otherwise both vanish.

The remaining rational projective vectors \([1/u_1:\cdots:1/u_5]\) have their coordinatewise squares on the fixed projective line \(\Lambda\subset\mathbb P^4\) cut out by the three rows with columns \((1,A_i^0)/D_i\). Each triple of columns is independent by noncollinearity. The five coordinate-zero intersections of \(\Lambda\) are therefore five distinct points, each a simple zero. Two coordinates cannot vanish simultaneously, since the remaining three columns are independent. No coordinate vanishes identically, since a second coordinate would then vanish there at some projective point. Over \(\mathbb C\) the cover of the line defined by taking the four coordinate ratios' square roots is connected of degree \(16\). Indeed for \(i=1,\ldots,4\), the ratio of coordinate \(i\) to coordinate five has odd valuation at the \(i\)-th coordinate-zero point, and the other three ratios have valuation zero there. These four squareclasses are independent. The only branch points are the five coordinate zeros, each of inertia order two (units over \(\mathbb C[[z]]\) have square roots). Riemann–Hurwitz gives \(2g-2=-32+5\cdot8=8\), hence \(g=5\). All our projective points have nonzero coordinates, on the smooth étale locus, and thus lift to rational points of the smooth normalization defined over \(\mathbb Q\). Faltings' theorem gives finitely many such projective vectors. By (P2) this leaves finitely many projective cofactor vectors \([t_{ij}]\), depending only on \(U\).

**PROVED — size for each fixed pattern.** Write \(t_{ij}=\lambda t_{ij}^0\) for one of these rational patterns. For a fixed ordered pair and a third index define the real cubic \(s_{ij}(n)=A_{ijk}(t^0,n)/(2t^0_{ik}t^0_{jk})\). It is independent of \(\lambda\) and at a sample gives the true scalar product divided by the common norm. Then the rational binary sextic
\[
   F(p,q)=\left(\prod_{h\notin\{i,j\}}L_{ih}L_{jh}
                         -s_{ij}(L)^2\right)/(t_{ij}^0)^2
\]
is fixed, nonzero, and takes the value \(\lambda^2\). We claim \(|F(p,q)|\gg_F B\) at primitive parameters with \(F(p,q)\ne0\) unless \(F\) is a constant times a sixth power of a rational linear form or a cube of a nonsplit indefinite rational quadratic. Indeed away from the projective real roots, \(|F|\gg B^6\). Near a particular root of multiplicity \(m\), every other factor stays bounded away from zero on the normalized direction, so \(|F|\) is comparable to \(B^6\) times the \(m\)-th power of the projective distance to that root. For rational roots, primitivity and integrality of a rational linear form give distance \(\gg_F B^{-1}\), enough for \(m\le5\). For irrational algebraic real roots, Roth gives \(\gg_F B^{-2-1/4}\) in a fixed bounded projective chart, enough for \(m\le2\). Multiplicity three forces an irreducible quadratic with both conjugate roots of multiplicity three. Higher multiplicity is impossible for an irrational root. These cases give \(|\lambda|\gg_U B^{1/2}\), sufficient here.

For the two exceptional sextic shapes, infinitely many distinct positive integral samples for a fixed projective pattern would supply infinitely many distinct parameter rays by primitivity. For each sample rescale \((p,q)\), and therefore every norm, by \(k=|\lambda|^{-1/3}>0\), and every Gaussian block by \(\sqrt{k}\). All \(P_i\) change by \(k^2\), so every cofactor changes by \(k^3\). The new *real* samples have fixed cofactors \(t^0\) up to a common sign. Retain an infinite subsequence with one sign. Their distinct rays remain distinct, and the rescaled parameters lie on \(F=1\). In the sixth-power case they lie on one fixed affine line after passage to one of at most two lines. The nonconstant-norm root argument using (P1) excludes infinitely many points there. In the indefinite quadratic-cube case the parameters lie on a fixed real hyperbola with irrational asymptote directions in the rational parameter plane. The conic proof of II.6 applies. Infinitely many samples extend the rational identities to the curve. Every norm is positive on the chosen branch, the ten forms are pairwise nonproportional, and no nonzero rational form vanishes on either isotropic direction. The five-star lemma then applies without an integrality hypothesis. The resulting star-sum contradiction excludes an infinite set. Thus for the finitely many exceptional patterns only finitely many samples need absorption into the constant. All cases give the claimed exponent. \(\square\)

**PROVED — a nonempty balanced instance.** Assign the ten ranks \(r=0,\ldots,9\) in any fixed order to the edges and take \(L_r(p,q)=p+rq\). There are arbitrarily large samples with the ten block norms distinct comparable primes \(1\bmod4\). The Green–Tao–Ziegler theorem for affine-linear forms of finite complexity applies to \(4x+4ry+1\) on boxes \(X\le x\le2X,\ X/20\le y\le X/10\). No two forms have proportional linear parts, and the local factors are positive at every prime, witnessed by \(x=y=0\) modulo that prime. Each prime supplies a Gaussian block. For comparable norms on any fixed plane and large $X=\max n_e$, no two coordinate forms can be proportional, since proportionality fixes their two norms by coprimality. The largest cofactor is then $\gg_U X^{1/4}$ and its pair distance is $\gg_U X^{-11/4}$. Hence $\Delta\gg_U H^{-11/40}$ for $H\asymp X^{10}$, allowing dependence on the comparison constants.

**PROVED — exponent with no comparability assumption.** First suppose the ten forms are pairwise nonproportional and let $Y=\min n_e$. All other nine norms are $\asymp_U B$, since on the compact projective circle no two normalized independent forms can be simultaneously arbitrarily small. Thus $H\asymp_U B^9Y$. The pair with maximum cofactor gives $\Delta\gg_U B^{-11/4}$ by (II.5). Choose also any pair separated by the cut of norm $Y$. Its nonzero integral cofactor gives $\Delta\gg_U B^{-5/2}Y^{-1/2}$. Raise these inequalities to powers $8/19$ and $11/19$, respectively, and multiply. The result is $\Delta\gg_U (B^9Y)^{-11/38}$.

If two forms are proportional use the fixed affine line from the start of the proof. Let $h$ be the number of nonconstant norms on that line, so $1\le h\le8$. For unbounded positive samples these $h$ norms are all $\asymp_U B$ and all others are fixed, hence $H\asymp_U B^h$. Some pair is separated by at most $4h/7$ of the nonconstant cuts. To check this, the ten separation counts form $M{\bf1}_S$, where $M$ is the edge adjacency matrix of II.5 and $|S|=h$. Their average is $3h/5$. For $h=1,2,3,4,6,7,8$ rounding down gives at most $4h/7$. For $h=5$ not all counts can equal three because $M$ is invertible and $M{\bf1}=6{\bf1}$, so this would imply $\mathbf1_S=\mathbf1/2$. That pair gives $\Delta\gg_U H^{-2/7}$, which also implies the asserted bound $\Delta\gg_U H^{-11/38}$. Bounded samples and a one-dimensional $U$ are absorbed in the constant. No threshold here is uniform over rational planes.

**PROVED — generic quartet, including moving small coefficients.** Only the six norms on a single quartet need have the form $L_{ij}(p,q)$, where the integral forms have coefficients bounded by $K\ge1$ and $\gcd(p,q)=1$. Suppose the three binary quadratics $L_{ij}L_{kl},L_{ik}L_{jl},L_{il}L_{jk}$ are linearly independent. Then, writing $X_4$ for the largest of these six **norms**, one has
$$
 \max_{ab\ {\rm in\ quartet}} |t_{ab}|\ \gg\ X_4^{1/2}K^{-3/2}
$$
with an absolute implied constant even when the forms vary. Indeed the Plücker quadratic in the proof above is nonzero by independence, has coefficients of size at most $6K^2 T_4^2$, and vanishes at a primitive pair. The rational-root theorem gives $\max(|p|,|q|)\le \max(1,6K^2 T_4^2)$. Since $X_4\le2K\max(|p|,|q|)$ this proves the assertion. For example assign ranks $0,1,2,3,4,5$ of $L_r=p+rq$ to edges $12,34,13,24,14,23$ in that order, with $6,\ldots,9$ assigned to the remaining edges. The coefficient triples of the three products are $(1,1,0),(1,5,6),(1,9,20)$, of determinant $32$. The progression of split primes above gives nonempty squarefree families with precisely these assignments and index one. For *all* sufficiently large such tuples with comparable norms, $T\gg X^{1/2}$ and $\Delta\gg H^{-1/4}$, including when the progression has diameter comparable to its smallest norm. Neither squarefreeness nor primality is necessary for the inequality. The use of the prime-pattern theorem is solely for the squarefree nonvacuity example. This is an arithmetic progression of **block norms**, a different condition from one on the cotangent coordinates of the five marks.

**COMPUTED — an explicit squarefree test case.** The ten terms $34913+2100r$, $0\le r\le9$, are distinct primes $1\bmod4$. The script certifies primality by trial division, constructs Gaussian blocks and checks that all ten reduced cofactors are nonzero integers. [Exact integer checks](check_part_ii.py)

**COMPUTED and PROVED — arithmetic progressions in every assignment.** The entire ten-term progression has the same exponent $4$ under arbitrary relabelling of its edges, assuming comparability. Indeed if all five Plücker forms vanished identically, the five affine points from (P2) would have ten distinct unoriented directions with slopes $0,1,\ldots,9$ because the coefficient vectors of the ten forms are $(1,r)$. Such a set of five real points does not exist. Translate and rescale the ends of the unique horizontal segment to $(0,0),(1,0)$. Each of the remaining points has two anchor slopes $a,b$ among $1,\ldots,9$, which determine it as $(b/(b-a),ab/(b-a))$. The script tests every triple with six distinct anchor slopes and asks whether its three internal slopes are the remaining distinct integers. There are 10,080 triples and zero solutions. Thus for every assignment some nonzero integral Plücker form has coefficients $O(T^2)$ at the sample. Pairwise coprimality of the ten norms implies $\gcd(p,q)=1$, and the rational-root theorem gives $T\gg\sqrt{\max n_e}$ with an absolute constant for these ten forms. The earlier specific assignment has an entirely symbolic certificate. [Exact direction enumeration](tenpair_progression_directions.py)

**PROVED — a moving progression with gaps.** If, on the independent quartet of the generic-quartet test, $n_e=p+r_e q$ with $|r_e|\le K$, $K\ge1$, and the origin of this progression has $p\asymp X=\max n_e$, the same rational-root calculation gives $T\gg\sqrt X/K$. This covers $K=X^{o(1)}$ without a template-dependent height threshold. The independence of the three complementary products must still be checked. The existence of the long progression of primes used above follows from the finite-complexity linear forms theorem with positive local factors. [Green and Tao, *Linear equations in primes*](https://doi.org/10.4007/annals.2010.171.1753)

**EVIDENCE — what remains at index one.** The growth reduction only needs a parent-dependent equation of small degree and tolerates atypical tuples. The fixed-plane result still fixes a rational subspace of the **norms**, while (II.2) controls their logarithmic sizes. I have no implication from small cofactor product to such a plane of controlled height, or to a quartet with the moving-progression condition. Generic quadratic norm polynomials of linear Gaussian blocks already span three dimensions. The real five-star lemma excludes the two-ended Pell alternative in II.6. A one-ended polynomial curve in the three-dimensional span need not yield reciprocal star sums. This is the step at which the polynomial decision and the general arithmetic bound both remain open. [Growth hypothesis and triangle compatibility](../free/answer.md)

**HEURISTIC — referee request.** I would prioritize independent checks of two steps in II.7, namely the transfer of (8) to (P3), including the perfect matching on the four complementary vertices, and the exclusion of the Pell conic by the real quartic's two cusps. For the polynomial problem my next exact test would be the six saturated resultants in II.5 over the rational Gaussian coefficient field, separating genuine components from components supported on collisions. The prime-seven computation is a check on any proposed rational answer, not a proof of emptiness.

### Twelve-line summary

1. **PROVED.** For index seven the sufficient subtuple-height cutoff is $H_S\le C^{-94/25}N^{47/50}$, while $H_S\asymp N$ does not supply it.
2. **PROVED.** Almost every five-tuple has $\log H_S/\log N$ near $15/16$, and at index nine the fixed-set argument works at every height.
3. **PROVED.** At the growth scale only $O_{A,C}(m^4)$ tuples fail index one, also in admissible angle models. Full squarefree index is at most three.
4. **PROVED.** Pure boundary covers reduce to translated isogenies. Unrestricted squarefree ten-pair profiles escape any finite list of nontrivial covers.
5. **PROVED.** Fixing the cube class of $\prod t_{ij}$ gives a rational lift to a cubic cover branched on the diagonals. Its target has multiplicity three.
6. **PROVED.** On every fixed rational plane of the ten pairwise coprime block norms, $T=\max|t_{ij}|\gg_U(\max n_{ij})^{1/6}$.
7. **PROVED.** The approximation exponent is at least $38/11$ there and $40/11$ for comparable norms. The constants may depend ineffectively on the plane.
8. **PROVED.** Six separating norms in $[X,X+cX^{1/6}]$ force their own cofactor to be $\gg\sqrt X$, without a fixed-plane hypothesis.
9. **COMPUTED and PROVED.** Ten comparable norms in arithmetic progression give exponent four in every edge assignment, including squarefree index-one examples.
10. **PROVED.** Polynomial norm pencils are impossible in all degrees. Linear blocks also fail on every generalized circle and on the two-complex-moment slice.
11. **EVIDENCE.** The unrestricted polynomial system, a power gain at unrestricted index one, and any improvement in cluster growth remain open.
12. **HEURISTIC.** Please referee the genus-five step from (8) and the real two-cusp exclusion, then test the six-resultant system with collisions removed.
