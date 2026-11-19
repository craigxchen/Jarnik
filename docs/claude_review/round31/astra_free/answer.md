# Short arcs of circles — a divisor count and an audit of the barriers

**EVIDENCE — verdict.** I do not have an improvement of the integer growth rate or a uniform bound. For growth, I would work on the primitive divisors counted in (1) below. A counterexample to growth would give quadratically many *different rational integers* of this kind at one parent norm. One can pose the counting problem without specifying an arc, prime angles, a configuration of points, or an exceptional algebraic set. I would try to refute this stronger counting hypothesis before attempting to prove it. An exact test revealed an extra fact: at a squarefree parent, any three candidate divisors with square product already close up as angles. In profiles of five or more vertices, separate primitive representations on every edge force a genuine circle. My main negative result about analytic inputs concerns a different proposed opening. The pointwise, all-scale, multiplicative spacing condition in round 20 actually determines the prime angles up to the ambiguities of a genuine cluster.

**PROVED using COMPUTED finite cases.** The 21 rank-five cases left in round 30, §7, all fail. In its squarefree real function-field model at degree eight and critical width, the maximum affine rank is four. The maximum size is eight and there is one maximizing sign code up to the stated symmetries. The last section gives the proof and the independently regenerated census. There is no step in this argument transferring the upper bound to an integer circle.

Labels throughout distinguish a proof (including a corollary of a named standard theorem), an exact computation with retained code, evidence or research judgment, and a heuristic or open hypothesis. All logarithms are natural.

## 1. What I think the fake-circle test establishes

**PROVED — logical scope.** Let a specified list of axioms admit models for a fixed arc constant with unbounded sizes \(m\), distinct positions, and \(W=\log N\le b m\log m\), for one fixed \(b\). Those axioms cannot imply the desired little-oh estimate. Indeed \(R=\exp(W/2)\) and \(W/\log W=O_b(m)\) when \(W\to\infty\), by monotonicity of \(u/\log u\). This argument applies to the conjunction of the axioms satisfied by the *same models*. From models for each of two lists separately, no corresponding implication for their conjunction follows. Likewise a construction for every fixed residue exponent \(A\), with \(b=b(A)\), does not address a test choosing \(A\) after \(b\).

**EVIDENCE — F1.** The skew double core is a convincing obstruction to using cut inequalities, monomial gaps, and their rational-character consequences. The bounded-modulus constructions and the analytic grafts have additional hypotheses. Some combinations really were constructed together, including the robust prime ensemble and character-level residues. The overview itself assigns different meanings to several of its barrier rows. In particular the single-bump limitation for restriction estimates is not an adversarial model satisfying every possible restriction input coupled to arithmetic. The small-height regime and twisted residue assertions have the qualifications recorded there. I would keep those distinctions in any statement of F1. There is also a literal slip in the expert note's opening: a uniform bound conditional on Vojta immediately implies a little-oh bound conditional on Vojta. The meaningful distinction is the lack of a known intermediate growth argument from a weaker standard conjecture. [Overview §§5.1–5.5](/home/user/Jarnik/docs/claude_review/overview/overview.md) · [Round 8 construction](/home/user/Jarnik/docs/claude_review/round8/lw.md) · [Round 10 analytic referee](/home/user/Jarnik/docs/claude_review/round10/referee_analytic.md)

**PROVED — the missing pair input.** There is a difference between bounding a reduced chord from below and requiring its scaled length to be an integer. In one unit class the latter says that an explicitly specified divisor \(d\) of \(N\) has the form \(s^2+t^2\), with coprime rational integers \(s,t\), and that the actual angular displacement is \(2\arg(s+it)\). Section 2 discards the displacements but retains a count of distinct such divisors. Section 4 shows that enforcing the exact displacements on a sparse set of edges already reconstructs a genuine circle. Thus even some conditions involving only *pairs of vertices* go beyond the fake axioms.

**EVIDENCE — F2.** I accept the stated one-section barriers in their specified spaces. For the five-point anticanonical problem the recorded maximal order is \(\lfloor 10m/3\rfloor\) in degree \(m\). The desired essential exponent exceeds \(15/4\), and the balanced single-section calculation runs exactly into the boundary. Bihomogeneous scaling cancels the same Gaussian gcd on the two sides of the Liouville comparison. But the unqualified claim that every determinant is critical is already contradicted by route 5 in the record. Choosing the box after the profile gives **supercritical forced vanishing** in its growth regime. What is missing is nonvanishing at the actual point for a form in the height-limited space. Small model boxes and the failed interpolation criterion do not settle that zero estimate. I give an independent check of the four-point counterexample, including a unit convention, below. [Round 23 width calculation](/home/user/Jarnik/docs/claude_review/round23/width.md) · [Overview §§6.2, 7.2](/home/user/Jarnik/docs/claude_review/overview/overview.md) · [Round 28 test](/home/user/Jarnik/docs/claude_review/round28/README.md)

**EVIDENCE — which openings I would fund.** P3-all* has the right quantifier for growth in its version (T′), and the all-monomial demands are genuinely stronger than assigning coherent residues to standalone characters. Its one-angle rigidity problem can still end with the original arithmetic problem. The low-height window \(P_+(N)\le m^{2-\delta}\) looks worthwhile for a restricted divisor count. Haar scarcity controls a probability, and pinned profiles have no free Haar coordinates. Even the absence of nonarithmetic window fakes would leave the genuine clusters. For the sharp sieve I distinguish prime-only support, a nonmultiplicative modification, and completely multiplicative modification on all primitive ideals. Section 3 classifies the last of these for a precise leading coefficient. I would put a worst-case divisor count ahead of trying to infer growth from sieve rigidity. I found no theorem establishing the necessity of Vojta-strength input. [P3 status in overview §6.3](/home/user/Jarnik/docs/claude_review/overview/overview.md) · [Round 21](/home/user/Jarnik/docs/claude_review/round21/README.md)

## 2. A counterexample to growth would give quadratically many primitive near-square divisors

**PROVED — reduction.** For an odd integer \(N\) supported on split primes put \(W=\log N\). For fixed \(B>0\), let \(D_B(N)\) count the positive divisors \(d\) satisfying
\[
\begin{split}
 16\sqrt N\le d\le\sqrt N\,W^B,\qquad
 s=\lfloor\sqrt d\rfloor,\qquad d-s^2=t^2,\\
 t\in\mathbb Z_{>0},\qquad (s,t)=1,\qquad 16t^2\sqrt N\le d .
\end{split} \tag{1}
\]
If the little-oh growth bound fails for any fixed arc constant, then for some fixed \(B,b>0\) there are unbounded normalized norms \(N\) with
\[
 D_B(N)\ \ge\ b\left(\frac{W}{\log W}\right)^2. \tag{2}
\]
Moreover the divisors witnessing this inequality have distinct short-arc edges realizing them on *one* circle. The definition of \(D_B\) makes no such compatibility demand.

**PROVED — proof.** Subdivide the original arc into a bounded number of pieces of length at most \(\tfrac12\sqrt{R_{\rm orig}}\), keep a most populated piece, and retain one Gaussian unit class for a fixed choice of split prime generators. Divide all retained points by their common Gaussian gcd \(g\). The radius becomes \(R=R_{\rm orig}/|g|\), and the arc length becomes at most \(\sqrt R/(2\sqrt{|g|})\). The residual norm \(N=R^2\) has only odd split factors since ramified and inert factors were common. The common unit class persists up to a common unit. Write \(m\) for the number retained and \(a_{ip}\) for the exponent of \(\pi_p\) in row \(i\), with \(p^{e_p}\Vert N\).

For an edge, orient the arguments in increasing order within the chosen arc and put
\[
 d_{ij}=\prod_p p^{|a_{ip}-a_{jp}|}.
\]
The ratio of its endpoints is \(A/\bar A\), with \(A\in\mathbb Z[i]\) of norm \(d_{ij}\) and \((A,\bar A)=1\) in \(\mathbb Z[i]\). Indeed take \(\pi_p\) or \(\bar\pi_p\) with exponent the absolute difference, according to its sign. Choose the overall sign of \(A=s+it\) so its argument is half the positive angular difference. Then \(s,t>0\), \((s,t)=1\), and
\[
  \frac{2t}{\sqrt{d_{ij}}}
  =\left|\frac{z_i}{z_j}-1\right|
  \le \frac{1}{2N^{1/4}}.
\]
This proves the lower and angular cutoffs in (1). The archimedean argument uses only differences of positions.

Set \(L_{ij}=\log(d_{ij}/\sqrt N)>0\). At each level \((p,h)\), let \(S\) be the number of rows with \(a_{ip}<h\). That level separates exactly \(S(m-S)\le m^2/4\) unordered pairs. Summing with weight \(\log p\) gives
\[
 \sum_{i<j} L_{ij}\le \frac{m^2W}{4}-\binom m2\frac W2=\frac{mW}{4}.
\]
At least half the pairs have \(L_{ij}\le W/(m-1)\). Failure of the original rate, after the bounded losses, gives \(m\gg_C W/\log W\) along a subsequence, with the implicit constant also depending on the size of the failure. This follows from \(W\le 2\log R_{\rm orig}\), monotonicity, and \(W\to\infty\) since \(m\to\infty\). Choose fixed \(B\) so that \(W/(m-1)\le B\log W\). In fact the floor-square property holds on *every* edge, even without this upper cutoff. Since \(d\mid N\),
\[
 0<\sqrt d-s=\frac{t^2}{\sqrt d+s}
   \le \frac{\sqrt d}{16\sqrt N}\le\frac1{16}.
\]
Thus \(s=\lfloor\sqrt d\rfloor\); there is only one admissible positive representation per divisor.

It remains to count divisors, since counting pairs with multiplicity could invalidate (2). At a retained \(d\), both \(s\) and positive \(t\) are determined by \(d\). The oriented pair therefore has the form \((g'\bar A,g'A)\) for some \(g'\in\mathbb Z[i]\), by coprimality. Its first endpoint confines \(g'\) to one arc on the circle of radius \(\sqrt{N/d}\). That arc has length at most
\[
 \tfrac12 N^{-1/4}\sqrt{N/d}\le\tfrac18 .
\]
It contains at most one Gaussian integer. Thus all the selected unordered edges give different \(d\)'s, and \(D_B(N)\ge\tfrac12\binom m2\). This proves (2). The small-arc distinct-ratio phenomenon is already present in the record's multiplicative Sidon observation. The use here is to count distinct *ordinary divisors* after taking floors. [Round 4 referee](/home/user/Jarnik/docs/claude_review/round4/referee_outside.md)

**HEURISTIC — open target.** Prove \(D_B(N)=o_B((W/\log W)^2)\) for every fixed \(B\). One may restrict the target to \(P_+(N)\le W^A\) for each fixed \(A\) by the record's near-extremal smoothness result. Even in the more restrictive window, one needs a worst-case estimate. Round 8, §8.3, counts smooth Gaussian integers of the common norm \(N^2\) after multiplying points by a conjugate anchor. Here the norms vary over divisors of \(N\), multiplicities from common multipliers are gone, and the relevant offsets are polynomial in \(\log N\). The ordinary-divisor question can be tested or disproved independently of the existence of a common arc. I did not find this count in the cited sections. [Round 8, §§3, 8.3](/home/user/Jarnik/docs/claude_review/round8/stability.md)

**PROVED — an additional compatibility test at the same parent.** A failure sequence also supplies, for some fixed larger \(B\), \(\gg(W/\log W)^3\) different unordered triples of divisors in (1) with square product. Increase the Markov cutoff on \(L_{ij}\) until fewer than \(1/6\) of the edges are discarded. Each missing edge spoils \(m-2\) triangles, so at least half of all triangles survive. The product \(d_{ij}d_{jk}d_{ki}\) is a square because the sum of the three absolute exponent differences is even at every prime. In fact each of the three quotients of the product of two by the third is an integer square. Global injectivity of the surviving edge labels gives distinct divisor triples. Thus if a large count of individual divisors refutes the first target, there is still a purely integer compatibility test before imposing any angles. Subcubic growth for this triangle count alone would suffice.

**PROVED — an average estimate and its limit.** For \(X\) large,
\[
 \sum_{X\le N\le 2X}D_B(N)\ll_B X^{3/4}\log\log X,
\]
where one may include all integer parents in the upper bound. The floor and window conditions put \(s\) between constant multiples of \(X^{1/4}\) and \(X^{1/4}(\log X)^{B/2}\). At each \(s\) the angular cutoff gives at most \((s+1)/(4X^{1/4})\) positive choices of \(t\). Each has at most \(2X/(s^2+t^2)\) parent multiples. The resulting harmonic sum over this range of \(s\) is \(O_B(\log\log X)\). An infinite sparse sequence of parents in (2) is consistent with this mean. For pointwise repulsion, two distinct primitive candidates \(s+it,s'+it'\) have \(\mathrm N\gcd(s+it,s'+it')\mid st'-s't\), and similarly for the conjugate and the plus sign. This follows by factoring out the common Gaussian multiplier in the determinant. The bounds on these determinants are only \(O( N^{1/4}W^B)\). I do not obtain a worst-case count from them.

**PROVED — one tempting weakening fails at precisely the requested order.** Omit the coprimality in (1). Take \(L\) to be the product of the primes \(q\equiv1\pmod4\) in \([y/2,y]\), \(x=32L\), and \(N=L^2(x^2+1)\). All prime factors are split and odd. For every product \(k=q_1q_2\) of two distinct primes of \(L\), take
\[
 d=k^2(x^2+1)=(kx)^2+k^2.
\]
Here \(d\mid N\), \(\lfloor\sqrt d\rfloor=kx\) for large \(y\), and \(d/\sqrt N=k^2\sqrt{x^2+1}/L\ge32 k^2\). The prime number theorem in the progression \(1\bmod4\) gives \(W\asymp y\) and \(\gg(W/\log W)^2\) distinct \(k\)'s. All pass the upper cutoff with, for example, \(B=5\). Every one has \((s,t)=k\). This obstruction is to the unqualified nonprimitive count over all split parents. I make no assertion of smoothness for \(x^2+1\).

**COMPUTED — finite tests.** I exhaustively counted (1) for odd split-supported parents through \(10^7\) at \(B=4\). There were 18,180 parents with one candidate, 113 with two, and two with three. For the first maximizing parent \(N=2,859,545\), the triples \((d,s,t)\) are \((2,859,545,1691,8)\), \((219,965,469,2)\), and \((77,285,278,1)\). To check completeness, write \(N=dm\). The angular condition is exactly \(256t^4m\le d\), so \(m\le\sqrt{10^7/256}\). The other maximizing parent is \(8,424,065\). Over the entire census **zero** unordered triples of candidates have square product. Thus even the first three candidates at a single parent need not fit together as triangle edges. The integer loops cover every \(s,t,m\) in these ranges. A primitive odd sum of two squares has only split prime divisors, so only \(m\) needs the prime-factor test. The upper cutoff is automatic in this computation. Indeed \(N\ge256\), and \(\sqrt N<(\log N)^4\) throughout \([256,10^7]\), as seen by checking the endpoints and the single interior critical point of the ratio. [Exact divisor census](divisor_census.py)

**COMPUTED — test with and without the gcd.** On a genuine four-point specialization from the family below, the script keeps all four primary points after gcd normalization. The common norm has 417 digits, and the directed arc has length below \(N^{1/4}/2\), certified with rational inequalities. Its six separation divisors are distinct and their offsets are \(1284505600,84934656,16777216,16777216,1048576,26214400\). All four triangular products are squares. Thus even two equal offsets need not duplicate a divisor. At \(y=200\) the nonprimitive construction gives 45 distinct divisors of a 90-digit parent, all rejected by the gcd test. These checks exercise the arithmetic of the reduction. Their sizes say little about the uniform asymptotic in (2). [Exact edge and stress tests](near_square_check.py)

**EVIDENCE — decision.** The first arithmetic test I would pursue is whether a family of *primitive* candidates can violate the count while coming from unrelated parts of the same parent circle. That would kill this approach cleanly, even if the circle conjecture is true. I would use exact division of \(d-\lfloor\sqrt d\rfloor^2\), coprimality, and the angular cutoff in the search. A direct count of all smooth numbers near all squares loses the common parent that the theorem needs.

**PROVED — a pointwise restriction.** Within any interval of divisors of multiplicative width \(H\), suppose \(t\le T\). Each squarefree core \(D>HT^4\) contributes at most one primitive divisor. Write \(d_i=Dk_i^2=s_i^2+t_i^2\). The errors \(k_i\sqrt D-s_i=t_i^2/(k_i\sqrt D+s_i)\) are positive. Subtraction therefore gives
\[
 |k_j s_i-k_i s_j|\le \max(k_j/k_i,k_i/k_j)T^2/\sqrt D<1.
\]
The determinant vanishes. Both \(s_i/k_i\) and \(t_i/k_i\) are then independent of \(i\) (the latter follows by squaring and positivity). Two proportional positive primitive integer pairs agree. This excludes repetitions from one large core in a dyadic interval, but gives no saving for squarefree parents, where all the cores already differ.

### The small test that changed my view of triangle compatibility

**COMPUTED.** I repeated the divisor census over all \(2^{22}=4,194,304\) squarefree parents supported on the first 22 split primes, from 5 through 229. I *removed* the upper cutoff \(d\le\sqrt N W^B\) and retained every other condition in (1). The maximum candidate count was six. There were 404 incidences of a parent and a square-product triangle, comprising 339 different divisor triples. **All 339 satisfied an exact signed angle relation.** For each \(d=s^2+t^2\) use \(A=s+it\). The script checks the imaginary part of each possible signed product by multiplying Gaussian integers, with conjugation for a minus sign. It also finds no four-vertex binary clique with six distinct divisor labels in this pool. The completeness check enumerates every squarefree divisor of the pool, then every squarefree complementary multiplier \(m\) allowed by \(256t^4m\le d\). [Smooth-parent census and exact phase test](smooth_divisor_census.py)

**PROVED — triangle theorem.** Let three distances \(d_1,d_2,d_3\) come from the three edges of an exponent profile in \(\prod_{p\mid N}[0,e_p]\), with \(d_{ij}=\prod p^{|a_{ip}-a_{jp}|}\). Suppose that *separately* each \(d_j\) admits any positive primitive representation \(s_j^2+t_j^2\) with \(16t_j^2\sqrt N\le d_j\), and set \(\psi_j=\arg(s_j+it_j)\). Then one \(\psi_j\) is the sum of the other two. No phases at vertices, coherent choices at primes, or floor-square conditions are assumed. For squarefree \(N\), the sole profile condition on three divisors is that their product be a square: at a prime it occurs in zero or two of them, and three binary vertex exponents realize either pattern.

*Proof.* Write \(A_j=s_j+it_j\). For the four choices of \(\epsilon_2,\epsilon_3\in\{\pm1\}\), form the Gaussian product \(P_\epsilon\) of the three factors, conjugating the last two as dictated by their signs. Put \(\epsilon_1=1\). For fixed \(\pi_p\mid p\), the signed valuations \(k_{jp}=v_{\pi_p}(A_j)-v_{\bar\pi_p}(A_j)\) have absolute values \(v_p(d_j)\), by primitivity and odd norm. For
\[
 h_{\epsilon p}=\sum_j\epsilon_j k_{jp}
\]
all values are even. Indeed the three local distances between points on a line are \(a,b,a+b\) for some \(a,b\ge0\). More importantly,
\[
 \frac14\sum_\epsilon |h_{\epsilon p}|=a+b\le e_p .
\]
The fixed signs of the \(k_{jp}\) do not affect this average. One can average over all eight signs and divide out the simultaneous reversal. Conditional on the sign of the largest term, the remaining symmetric sum has absolute value at most that term, giving the displayed equality.

By unique factorization,
\[
 P_\epsilon=r_\epsilon\,\eta_\epsilon\, B_\epsilon^2,\qquad
 r_\epsilon\in\mathbb Z_{>0},\quad \eta_\epsilon\in\{1,i,-1,-i\},\quad
 Q_\epsilon=\mathrm N B_\epsilon=\prod_p p^{|h_{\epsilon p}|/2}.
\]
To obtain this identity remove the common rational prime powers, then halve the remaining even exponents. The mean of \(\log Q_\epsilon\) is at most \(W/2\). Choose one sign pattern with \(Q_\epsilon\le\sqrt N\). Write \(u=\sum_j\epsilon_j\psi_j\). The absolute value \(\sqrt{Q_\epsilon}|\sin(u/2)|\) is an absolute coordinate of a Gaussian integer rotated by a multiple of \(\pi/4\). It is zero or at least \(1/\sqrt2\). But \(\psi_j\le\arcsin(1/(4N^{1/4}))\), and \(\arcsin x\le\pi x/2\) on \([0,1]\). Thus
\[
 |\sin(u/2)|\le\frac{3\pi}{16N^{1/4}}
      <\frac{1}{\sqrt2\,N^{1/4}}\le\frac{1}{\sqrt{2Q_\epsilon}} .
\]
It must vanish. The bound on \(|u|\) makes \(u=0\), and positivity gives the claimed equality. This explains the exact outcome on all 339 triples.

**PROVED — five vertices recover a circle without assigning vertex angles.** Take \(m\ge5\) rows of a profile and a primitive near-real representation as above for every edge. Then there exist Gaussian points of norm \(N\) whose angular distances are precisely \(2\psi_{ij}\) and which lie in an arc of length
\[
 2\sqrt N\,\arcsin(1/(4N^{1/4})) .
\]
Up to column reversals and shifts, their exponent rows realize the given profile. This also holds for fewer vertices if the largest edge length occurs just once.

*Proof.* Regard \(\psi_{ij}>0\) as distances. By the triangle theorem every triangle is degenerate, with equality in its largest-side triangle inequality. Take diameter endpoints \(A,B\), with distance \(D\), and put \(x_i=\psi_{Ai}\). Then \(\psi_{Bi}=D-x_i\). For two interior vertices, the two endpoint triangles imply either \(\psi_{ij}=|x_i-x_j|\) or \(\psi_{ij}=D\) with \(x_i+x_j=D\): if the difference fails, the two sums must agree. In the second case there can be no third interior vertex \(k\). Two sides of length \(D\) in a degenerate triangle would force the third side to be zero, so both distances from \(k\) to \(i,j\) would have to be the ordinary coordinate differences. But
\[
 |x_k-x_i|+|x_k-x_j|<D
 \quad (0<x_i,x_j,x_k<D,\ x_i+x_j=D),
\]
contradicting the triangle inequality. Hence for \(m\ge5\), or for a unique diameter, all distances are ordinary differences on \([0,D]\).

Put \(w_i=e^{2ix_i}\). Every ratio \(w_i/w_j\) is \(A_{ij}/\bar A_{ij}\) or its inverse, and so is in \(\mathbb Q(i)\). Its absolute signed valuation at \(p\) is exactly \(|a_{ip}-a_{jp}|\). Since \(w_A=1\), all \(w_i\) belong to \(\mathbb Q(i)\) and their valuations vanish outside the specified primes. At each \(p\), shift the valuations \(h_i=v_{\pi_p}(w_i)\) by their minimum. They fit in \([0,e_p]\), since their pairwise differences in absolute value equal those of the given column. In fact a map between finite subsets of a line preserving every distance is a translation or a reflection followed by translation, as is seen using the two extreme points. Form the Gaussian products with these shifted exponents and complements to \(e_p\); norm one of the ratios gives the correct conjugate valuations. Restoring the residual Gaussian units as in §4 realizes the \(w_i\) up to common rotation at norm \(N\). The diameter bound gives the stated arc length.

**PROVED — even primitivity can be removed from the complete-graph test.** With arbitrary \(A_j=s_j+it_j\) at the same norms, its signed valuation has magnitude at most the local exponent \(v_p(d_j)\) and has the same parity. In the triangle proof the mean of the absolute signed sum over signs is nondecreasing in each nonnegative magnitude, since \((|x+u|+|x-u|)/2=\max(|x|,u)\) for \(u\ge0\). Thus its upper bound \(e_p\) and even parity survive. The square-factor proof and line embedding are unchanged. For the rational ratios obtained on the line we now have \(|h_i-h_j|\le |a_{ip}-a_{jp}|\), sufficient to fit their ranges in \([0,e_p]\) and reach the very same norm \(N\). The stronger recovery of the original column uses primitivity. Together with the automatic floor identity, this gives an integer-only test on a complete profile for \(m\ge5\): on every edge check that \(d-\lfloor\sqrt d\rfloor^2=t^2>0\) and \(16t^2\sqrt N\le d\). Any passing profile produces a genuine cluster. Conversely any sufficiently shortened genuine one-unit cluster passes with its own profile. In particular the coprimality dropped here remains essential for the individual-divisor counting proposal, as its explicit counterexample shows.

**COMPUTED — test with prime powers.** A further exhaustive check used the \(3^{13}=1,594,323\) parents having exponents 0, 1, or 2 at the 13 primes listed in the script. It allowed imprimitive edges. It found 53 parent/triple incidences with the required local exponent pattern, including 52 at nonsquarefree parents and seven with an imprimitive factor. They gave 42 distinct divisor triples. All passed the exact Gaussian phase test. The enumeration includes every complement of an eligible divisor within the parent exponent limits and the angular cutoff. [Prime-power test](prime_power_triangles.py)

**COMPUTED and PROVED — four vertices really are exceptional for this lift.** The smooth census gives \(N=131490079981\) and
\[
 d_1=6451601=2540^2+1,\quad d_2=28879877=5374^2+1,\quad
 d_3=92794693=9633^2+2^2.
\]
The support of \(d_3\) is the symmetric difference of the first two supports and their union is the support of \(N\). Take the four binary rows \(0,v_1,v_2,v_1+v_2\) over \(\mathbb F_2\), where \(v_i\) is the prime-support indicator of \(d_i\). Every edge passes the primitive test. Each of the three distinct lengths occurs on two opposite edges, so the positive diameter occurs twice. Four distinct points on a real interval have exactly one diametral pair. This profile with its assigned edge lengths cannot lift to the claimed four points on a short arc. The example illustrates an exception to the reconstruction theorem. The four-point counterexample to Question B is a different, Zariski-density assertion proved in §6. [Exact divisor data](smooth_divisor_census.py)

**EVIDENCE — consequence for the search.** In a squarefree parent, take the binary group with one coordinate per prime. A candidate divisor is a step in that group, given by its prime support. Count the allowed steps for the two-point growth target, or count triples of allowed steps with zero binary sum for its triangle refinement. If one imposes the primitive angular cutoff on *every* edge of a clique, the lifting theorem shows precisely why this is already as difficult as a genuine short arc for five or more vertices, up to the fixed arc constant and point units. In the reduction only a positive proportion of the edges have to meet the stronger polylogarithmic upper window in a putative growth counterexample. The asymptotic saving I would look for is in the number of these individual steps or their square-product triples, before imposing the full clique.

## 3. The sharp multiplicative sieve on primitive Gaussian ideals

**PROVED — rigidity theorem.** Let \(S\) be a finite set of split odd primes. On primitive Gaussian ideals (their generators have no rational integer factor \(>1\)), change the angles at \(\pi_p,\bar\pi_p\) for \(p\in S\) to \(\phi_p,-\phi_p\), keep all other angles genuine, and extend completely additively modulo \(\pi/2\). Write \(\psi_d\) for the resulting angle of an \(S\)-smooth primitive ideal \(d\), and \(h=\mathrm N d\). For an angle write \(\|\cdot\|_u\) for distance to \((\pi/2)\mathbb Z\). Assume there exists \(c>0\) such that for every such \(d\ne1\) and all primitive \(\alpha,\beta\) of norms coprime to \(\prod_{p\in S}p\),
\[
 \sin\|\psi_d+\arg\beta-\arg\alpha\|_u
 \ \ge\ \frac{c}{\max(h\,\mathrm N\beta,\mathrm N\alpha)}. \tag{3}
\]
Then for every \(p\in S\),
\[
 \phi_p=\pm\arg\pi_p\pmod{\pi/4}.
 \tag{4}
\]
Conversely these choices satisfy (3) with an absolute \(c>0\). A purported circle cluster whose angles satisfy (3) is a genuine norm-\(N\) cluster up to a common rotation, choices of conjugate prime, and individual Gaussian units. This resolves the pointwise multiplicative case of the round-20 heuristic. It requires all scales and multiplication across modified and unmodified primes. [Round 20 killer, §5.1, equation (**)](/home/user/Jarnik/docs/claude_review/round20/kill_composite.md)

**PROVED — first step, irrational directions.** Call a direction rational if it contains a nonzero integer vector. Fix \(d\) and let \(U\) be rotation by \(\psi_d\). Set \(q=\prod_{p\in S}p\), \(G=\mathrm{SL}_2(\mathbb R)\), and let \(\Gamma=\Gamma(q)\) be the principal congruence subgroup. If the direction is irrational then the diagonal \(G\)-orbit of \((\Gamma,U\Gamma)\) in \((G/\Gamma)^2\) is dense. Here is the deduction from Ratner's orbit-closure theorem. The diagonal group is connected and generated by one-parameter unipotent subgroups. Its orbit closure is a homogeneous orbit of a connected intermediate group with lattice stabilizer. A Lie algebra in \(\mathfrak{sl}_2\oplus\mathfrak{sl}_2\) containing the diagonal is either the diagonal or the whole algebra, since the quotient by the diagonal is the irreducible adjoint module. In the diagonal case \(\Gamma\cap U\Gamma U^{-1}\) is a lattice, hence finite-index in both lattices. Commensurate subgroups have the same parabolic fixed points, since some power of each parabolic belongs to their intersection. The parabolic fixed points of \(\Gamma(q)\) are exactly \(\mathbb P^1(\mathbb Q)\) by conjugating powers of the integral elementary unipotent at infinity. Consequently \(U\) must send the horizontal rational line to a rational line. This contradicts irrationality. The only imported homogeneous-dynamics result is Ratner's orbit-closure theorem. M. Ratner, *Raghunathan's topological conjecture and distributions of unipotent flows*, Duke Math. J. 63 (1991), 235–280.

Given \(\epsilon>0\) small, density provides \(g\in G\), \(\gamma_0,\gamma_1\in\Gamma\) with \(g\gamma_0\) and \(gU\gamma_1\) within \(\epsilon\) of the identity in operator norm. Use the *same* basis column \(e_j\) for both matrices, choosing \(j\) such that \(|g^{-1}e_j|\ge\|g^{-1}\|/\sqrt2\). Put \(\alpha=\gamma_0e_j,\ \beta=\gamma_1e_j\). These vectors are primitive and each has norm squared \(1\bmod q\). After applying \(g\), the two vectors \(\alpha,U\beta\) have determinant \(O(\epsilon)\), which is also their original determinant. Before applying \(g\) their lengths are both between \((1/\sqrt2-\epsilon)\|g^{-1}\|\) and \((1+\epsilon)\|g^{-1}\|\). Therefore the left side of (3) times its height denominator is \(O(h\epsilon)\). Since \(h\) is fixed and \(\epsilon\) arbitrary, (3) forces a rational direction. This argument needs no rate in Ratner's theorem.

**PROVED — second step, rational directions have a height cost.** Write a rational \(\psi_d=\arg b\) modulo units with \(b\in\mathbb Z[i]\) primitive, \(Q=\mathrm N b\). Then the infimum over \(\alpha,\beta\) of the left side of (3) times the height denominator is at most \(\sqrt{h/Q}\). Start with the integer vector \(v=b\) and an integer vector \(w\) with \(\det(v,w)=1\). Replace \(w\) by \(w+kv\) so that \(\mathrm Nw\) avoids \(S\). At a split \(p\in S\) there are at most two bad residues of \(k\), since the isotropic locus is two lines and the affine line \(w+kv\) is contained in neither by the determinant condition. Apply CRT. Now \(w\) is primitive.

For a rational covector \(\ell\) with \(\ell(v)\ne0\), rescale to an integral covector and choose a nonzero integral vector \(a\) in its kernel. The transformations
\[
 g_n=I+nqQ\,a\ell
\]
have determinant one. If \(B_b\) denotes the integer matrix of complex multiplication by \(b\), both \(g_n\) and \(B_b^{-1}g_n B_b\) belong to \(\Gamma(q)\), since \(\det B_b=Q\). Take \(\alpha_n=g_nw\) and \(\beta_n=B_b^{-1}g_n v\), primitive vectors of norms coprime to \(q\), the second obtained from the first standard basis vector. The determinant of \(b\beta_n,\alpha_n\) is one, and \(|\alpha_n|/|b\beta_n|\) tends to \(|\ell(w)/\ell(v)|\) provided the latter is nonzero. Since \(v,w\) form an integer basis, these limits can approach any positive real. If \(r=|\alpha_n|/|b\beta_n|\), the normalized sine for the two vectors before reduction modulo units is precisely \(\max(r,h/(Qr))\). Its infimum along these limits is \(\sqrt{h/Q}\). Distance modulo units can only decrease it. Thus \(Q\le h/c^2\).

Apply this fact to even powers of the ideal \(\pi_p\). Write a primitive rational direction at this prime as \(\arg b_0+t\pi/4\), where \(b_0\) has odd norm \(Q_0\) and \((b_0,\bar b_0)=1\). This decomposition follows by removing the possible single factor \(1+i\) from a primitive representative. The primitive norm for direction \(2k\phi_p\) is \(Q_0^{2k}\), since the ramified direction becomes a unit and the odd part remains primitive under powers. Our bound gives \(Q_0^{2k}\le p^{2k}/c^2\) for every \(k\), hence \(Q_0\le p\).

**PROVED — identifying the primes.** For each \(p\in S\), take the signed exponent vector of the odd primitive \(b_{0,p}\) over split primes, recording the difference of valuations at conjugate primes. Form the square matrix keeping only the rows with rational prime in \(S\). This matrix has full rank. Otherwise an integer relation among its columns gives a nonunit primitive \(S\)-smooth ideal \(d\), taking conjugate prime powers for negative coefficients. Its angle is a rational direction with zero signed valuation at every prime of \(S\). Its primitive representative \(\alpha\) has norm coprime to \(q\), so \(\beta=1\) gives an exact collision in (3). Order \(S\) increasingly. The bound \(Q_{0,p}\le p\) makes the matrix triangular. Each diagonal entry must be nonzero, which uses the entire norm budget \(p\). Thus \(Q_{0,p}=p\), with exactly the conjugation choices in (4).

Conversely a signed prime substitution in (4) conjugates some fixed split prime factors and adds a multiple of \(\pi/4\). On all primitive directions, record the free signed valuations at the split primes and the ramified exponent modulo two. The induced transformation changes signs of some free coordinates and adds a parity function of them to the ramified exponent. It is injective. Its resulting primitive representative has norm at most twice the original norm, since the ramified exponent is zero or one. Two distinct primitive ideals of norm at most \(X\) therefore give distinct primitive rational directions of norm at most \(2X\). Their determinant, after rotating by the minimizing Gaussian unit, is a nonzero integer. This gives sine distance at least \(1/(2X)\) and proves the converse. Finally, in the angle of a row \(\sum_p(2a_{ip}-e_p)\phi_p\), shifts \(t_p\pi/4\) add a common phase plus row multiples of \(\pi/2\). This recovers the arcs including their freely drifting common direction.

**PROVED — the operator bound and the significance of its constant.** Suppose the dual large-sieve norm on *all* primitive ideals up to \(X\), for frequencies \(-K,\ldots,K\) and angles \(4\psi\), has upper bound \(A X+B_0 K\) for all \(X,K\ge1\) with fixed \(B_0<4\). Then its point set has unit-reduced angular spacing \(\gg_{A,B_0}1/X\). Choose \(0<b<1\) with \(2(1+b)>B_0\), and \(K=\lceil(A+1)X/(2(1+b)-B_0)\rceil\). If two angles differ modulo units by at most \(\arccos(b)/(4K)\), the off-diagonal Dirichlet sum of their Gram matrix is at least \((2K+1)b\). Its larger eigenvalue is at least \((2K+1)(1+b)>AX+B_0K\). This proves spacing and hence (3) with a changed constant. The genuine separated-points sieve has leading term \(2K\). Consequently a fully multiplicative graft satisfying this precise all-scale sieve gives genuine circle clusters. The proof permits a very large finite height of the first violation. It says nothing about a sieve checked only to a prescribed height or only at prime ideals.

**PROVED — the order-of-magnitude bound differs.** There is an explicit irrational graft satisfying both families of monomial gaps for one prime and an \(O(X+K)\) sieve on *all* primitive ideals. At \(p=13\) set \(\theta=\arg(2+3i)\) and \(\phi=\theta/2\). The slope is \((\sqrt{13}-2)/3\). For odd \(n>0\), the imaginary part of \((2+3i)^n\) is a nonzero multiple of 3, and the real part a nonzero multiple of 2. Nonvanishing follows from unique factorization, since a real or purely imaginary power would be associate to its conjugate. Using \(|\sin 2x|\le2\sin\|x\|_u\) and \(|\cos 2x|\le2\sin\|x-\pi/4\|_u\) proves the axial gap \(13^{-n/2}\) and diagonal gap \(13^{-n/2}/\sqrt2\). For even \(n=2k\), the genuine lattice gaps at \(k\theta\) are stronger, by integrality of the real and imaginary components and of their sum and difference. Negative indices follow by symmetry.

Write the signed 13-exponent of a primitive ideal as \(2h+r\) with \(r=0,1\). Within each fixed \(r\) its angle is a translate by \(r\phi\) of a genuine primitive direction, obtained by replacing the exponent by \(h\). This map is injective on each class and the norm only decreases since \(|h|\le|2h+r|\). The genuine constant bounds each of the two sums separately, so Cauchy–Schwarz gives twice that constant on the union. Yet by the irrational-direction proof, (3) fails for every \(c>0\). This is a one-prime example distinguishing sieve statements. It makes no assertion about many-point fakes at the barrier scale.

**COMPUTED — a test on three fibres.** Set \(p=29\) and \(\phi=\arg(1+i\sqrt2)\). The primitive 29-free multipliers \(533+8422i,1613+994i,139+326i\) on the fibres of exponent \(0,1,2\) give heights \(71{,}214{,}173,104{,}104{,}345,105{,}627{,}077\). With \(X\) the largest, the three unit-reduced angles fit in an arc of width \(<0.0216/X\). The certificate uses integer arithmetic in \(\mathbb Z[\sqrt2]\). The same script certifies the \(\sqrt{h/Q}\) cost within a factor \(\sqrt{1.02}\) for six rational cases, including \(5\mapsto3+4i\) at powers 2, 4, 8, where the error compounds despite the absence of an exact collision with a 5-free direction. [Exact certificates](graft_check.py) · [Optional numerical search for the three multipliers](route_graft_search.py)

**PROVED — this test angle passes the monomial tests.** Write \((1+i\sqrt2)^n=a+i\sqrt2 b\) with integer \(a,b\). Neither is zero for \(n>0\), since \((1+i\sqrt2)/(1-i\sqrt2)\) has trace \(-2/3\) and hence is no root of unity. The axial gap is at least \(3^{-n/2}\). Also \(|a\pm\sqrt2 b|\ge1/(|a|+\sqrt2|b|)\ge1/(\sqrt2\,3^{n/2})\) by taking the quadratic norm. The sine of the distance to the diagonal directions is therefore at least \(1/(2\cdot 3^n)\ge29^{-n/2}/\sqrt2\). Thus the small three-fibre spacing survives the monomial tests imposed by the record. The exact all-scale rigidity theorem comes from the two-lattice argument, independently of this computation.

## 4. How little exact chord data already reconstructs a circle

**PROVED — scalar reconstruction.** Let \(m\ge3\), let the phases \(\theta_i\) be distinct modulo \(2\pi\), and give them any exponent profile with \(0\le a_{ip}\le e_p\) at the split primes of \(N\). Define \(d_{ij}\) as in §2 and put
\[
 t_{ij}=\sqrt{d_{ij}}\sin((\theta_i-\theta_j)/2),\qquad
 s_{ij}=\sqrt{d_{ij}}\cos((\theta_i-\theta_j)/2).
\]
Choose a set \(U\) of at least two vertices such that at each prime its exponent values straddle \(e_p/2\), inclusively. If \(t_{ij}\in\mathbb Z\) on every edge incident with \(U\), then some Gaussian points of norm exactly \(N\) have all the prescribed phase differences. This conclusion allows point units and does not identify the original artificial prime angles. The data do not even need to arise from a common angle vector.

**PROVED — proof.** On a triangle of tested edges with anchor 0, write \(s_i=s_{i0},t_i=t_{i0}\). Sine subtraction gives
\[
 t_i s_j-t_j s_i=t_{ij}\sqrt{d_{i0}d_{j0}/d_{ij}}=k\in\mathbb Z\setminus\{0\}.
\]
At each prime the radicand exponent is nonnegative and even by the triangle inequality and parity for differences of integers. Also \(s_i^2=d_{i0}-t_i^2\in\mathbb Z\). Squaring \(t_i s_j=k+t_j s_i\) shows \(2kt_js_i\in\mathbb Z\), so \(s_i\) is rational and hence integral. Repeat with exchanged anchors. Two centers give triangles through every other vertex. Thus every \(w_i=e^{i(\theta_i-\theta_0)}\) lies in \(\mathbb Q(i)\), since on an anchor edge its ratio is \((s+it)/(s-it)\). On every tested edge \(\operatorname{Im}(w_i/w_j)=2s_{ij}t_{ij}/d_{ij}\) gives integral \(s_{ij}\) by the same rational-square-root argument.

Write \(h_i=v_{\pi_p}(w_i)\) at a fixed split prime. On a tested edge,
\[
 |h_i-h_j|\le |a_{ip}-a_{jp}|,
\]
because the numerator \(s_{ij}+it_{ij}\) has norm \(d_{ij}\). For two arbitrary vertices write their exponents as \(a\le b\). A center in \([a,b]\) bounds the valuation change by \(b-a\). If both are below all centers, using the minimal center exponent \(c\le e_p/2\) costs at most \(2c\le e_p\). Above all centers the analogous bound is \(2(e_p-c')\le e_p\) with the maximal center exponent \(c'\). In the remaining case \([a,b]\) lies between consecutive center values \(c<c'\). The two path lengths bound the change by \(\min(a+b-2c,2c'-a-b)\le c'-c\le e_p\). Hence the range of all \(h_i\) at this prime is at most \(e_p\). There are no other prime factors in any \(w_i\), by the anchor formulas. Multiply the factors \(\pi_p^{h_i-\min h}\bar\pi_p^{e_p-h_i+\min h}\) over \(p\), taking the minimum separately at each prime. Their norms are \(N\); the ratios to the anchor have exactly the valuations of \(w_i\) at both conjugate primes because each \(w_i\) has norm one. The remaining factors are units. Choosing these units gives the assertion.

**PROVED — at the barrier scale few edges suffice.** For a normalized profile assume \(\log d_{ij}\ge W/2\) on all pairs. This holds on the short one-unit arc above by the pair gap alone. If \(S\) rows lie below a given prime threshold, summing the pair distances as in §2 gives
\[
 \sum_{p,h}(2S-m)^2\log p\le m W .
\]
In fact the left side is \(m^2W-4\sum_{i<j}\log d_{ij}\). Call a prime bad if either of its vertex sets \(a_{ip}\le e_p/2\), \(a_{ip}\ge e_p/2\) has size less than \(m/4\). Each bad prime has at least one threshold with \(|2S-m|\ge m/2\). There are at most \(4W/(m\log 5)\) bad primes. Include vertices attaining the minimum and maximum of each bad prime (normalization gives 0 and \(e_p\)). For \(r\) total primes, \(\lceil4\log(4r)\rceil\) independent uniform vertices hit both half-sets at every good prime with positive probability, by the union bound. Add two distinct vertices if necessary. Their incident edges number \(O(m\log(2r)+W+m)\), in particular \(O(m\log m)\) when \(W=O(m\log m)\). Every witness passing these sparse tests produces a genuine cluster. Conversely every genuine one-unit cluster supplies a witness on its own exponent profile. This compares existence questions without identifying the prime angles of a prescribed fake.

**COMPUTED — what can go wrong with two arbitrarily chosen centers.** At \(N=5\cdot13\cdot17\), rows \(000,010,100,111\) and the rational phases in the script have integer scaled sines of magnitudes \(2,1,12,1,2\) on edges \(01,02,03,12,13\). The missing edge has magnitude \(7/5\). The valuation range at a prime above 5 is two, whereas \(v_5(N)=1\). Thus rationality from triangles alone does not settle the common denominator. An independent exact triangle check tested 12,501 positive magnitude triples, found 264 survivors of the nonzero squared subtraction condition, and found an integral real component on every survivor. [Scalar chord certificate](scalar_chord_check.py) · [Earlier converse caveat in round 8](/home/user/Jarnik/docs/claude_review/round8/referee_stability.md)

## 5. How much of a five-point theorem growth actually needs

**HEURISTIC — an alternative open statement, allowing moving exceptions.** There is no canonical weakest conjecture across different methods. For a question about five-point geometry I would ask for the following count before asking for \(\alpha^{\rm ess}>15/4\). For each split parent \(N\) and fixed \(B\), form the set of ordered similarity classes of distinct, one-unit, primitive five-tuples with primitive norm dividing \(N\). Keep only classes satisfying **all** of these conditions.

- All angular displacements are at most \(\tfrac12 N^{-1/4}\).
- For each of the fifteen nontrivial unoriented cuts, the sum \(w_A\) of the prime-level weights \(\log p\) giving that cut satisfies \(|w_A-W/16|\le B\sqrt{W\log W}\).
- The ten scalar cofactors of §4 satisfy \(H=\sum_{i<j}\log|t_{ij}|\le B\log W\).

Call this set \(\mathcal B_B(N)\). Use four complex ratio coordinates \((z_2/z_1,\ldots,z_5/z_1)\). The hypothesis is that **for each fixed \(B\)**, there is a nonzero polynomial \(F_{N,B}\) in these four variables of total degree \(o_B(W/\log W)\) such that
\[
 \#\{u\in\mathcal B_B(N):F_{N,B}(u)\ne0\}
       =o_B((W/\log W)^5).
 \tag{5}
\]
The little-oh bounds must hold over all relevant parents, not as pointwise limits for each prime set. The possibility of choosing one polynomial per parent and leaving some points uncovered matters. This asks for no uniform exceptional set from a fixed-\(S\) theorem.

**PROVED — (5) suffices for growth.** Here is a self-contained extraction to fix the quantifiers and the gcds. As in §2 take a failure sequence with \(m\gg W/\log W\). Index threshold layers by \(l\), let \(w_l=\log p\) on its prime, and put \(v_{il}=\sqrt{w_l/W}\,\sigma_{il}\), with sign \(\sigma_{il}=\pm1\) according to the side of the threshold. The Gram matrix \(G=(\langle v_i,v_j\rangle)\) has diagonal one and off-diagonal entries \(1-2\log d_{ij}/W\le0\). Its operator norm is at most two. Indeed \(I-G\) is entrywise nonnegative, its largest eigenvalue is at most one by positivity of \(G\), and taking absolute values in the Rayleigh quotient bounds the absolute values of its other eigenvalues by its largest.

Take five independent uniform row indices \(i_1,\ldots,i_5\). For a nonempty subset \(J\) of slots set
\[
 f_J=W^{-1}\sum_l w_l\prod_{j\in J}\sigma_{i_jl}.
\]
For \(|J|=1\), its second moment is at most \(2/m\): apply the sign matrix with entries \(v_{il}\) to the unit vector \((\sqrt{w_l/W})_l\) and use its squared operator norm bound of two. For \(|J|\ge2\), expanding and using independence gives a weighted sum of powers \((m^{-1}\sum_i\sigma_{il}\sigma_{il'})^{|J|}\). Replacing these terms by squares is an upper bound. Summing the squares yields \(\operatorname{tr}(G^2)/m^2\le2/m\). Thus with probability at least \(3/4\) all 31 coefficients have absolute value \(\le K/\sqrt m\) for a fixed sufficiently large \(K\), by Markov's inequality and a union bound. Expanding each pattern indicator as a product of five factors \((1\pm\sigma)/2\) gives total weight \(W/32+O(W/\sqrt m)\) for each oriented pattern. Complementary patterns sum to \(W/16+O(W/\sqrt m)\). Dividing a tuple by its Gaussian gcd removes precisely the layers constant on it, retaining all fifteen nontrivial cut weights. Its primitive norm divides \(N\), its angles change by a common rotation, and its pair factors and their scalar cofactors are preserved.

The probability of a repeated row is at most \(10/m\). For ordered distinct rows the §2 estimate gives \(\mathbb E\sum_{i<j}L_{ij}\le5W/(m-1)\) on their ten edges, and \(0\le\log |t_{ij}|\le L_{ij}/2\). Markov's inequality bounds \(H\) by \(O(W/m)\) with probability at least \(3/4\) conditional on distinctness. Together these estimates give \(\gg m^5\) ordered five-tuples satisfying the three conditions for one fixed \(B\) on the failure sequence.

These tuples give just as many *distinct classes* on the chosen arc. Their first directed ratio fixes its primitive Gaussian factor up to sign: two Gaussian representatives with the same ratio differ by a real rational scalar, hence by sign if both primitive. At the fixed parent norm the possible common multipliers lie in an arc of length at most \(1/8\), as in §2, also for a negatively directed first pair. That fixes the first two points, and the other ratios then fix the rest.

Finally homogenize \(F_{N,B}\) by multiplying its evaluation on the four ratios by \(z_1^{\deg F}\). This gives a nonzero polynomial of total degree \(\deg F\) on the five original coordinates. Its zero set meets the Cartesian power of any \(m\)-element set in at most \((\deg F)m^4\) points. For completeness, induction on the variables proves this bound: take the leading coefficient for the last variable of degree \(j\), bound its exceptional choices by \((\deg F-j)m^3\) in the four previous variables, and otherwise count at most \(j\) last coordinates. Thus the zeros cover \(o(m^5)\) of the selected tuples. The remaining classes contradict (5).

**PROVED — comparison with the recorded statements.** The fixed-set \((L_5^\circ)\) implies (5): a fixed polynomial cuts out a hypersurface containing its exceptional set, normalized balance tends to exact balance, and \(B\log W\) is eventually below the asserted linear lower bound at scale \(W/16\). The same extraction also explains why \((L_5^\circ)\) gives *uniformity*: for any sequence with \(m\to\infty\), its selected tuples have relative balance error \(O(m^{-1/2})\) and \(H/W=O(1/m)\), and a fixed hypersurface misses almost all of them. For growth alone, even a superlogarithmic cofactor lower bound along the shrinking balance window would suffice. Formula (5) allows violations of that lower bound. I see no implication in the reverse direction. The elementary two-point count in §2 and (5) are different sufficient hypotheses. A failure of the count in (1) need not refute Question B or Vojta. [Round 20 five-point formulation](/home/user/Jarnik/docs/claude_review/round20/gauss.md) · [Earlier superlogarithmic target](/home/user/Jarnik/docs/claude_review/threshold/L8.md)

## 6. Independent check of the four-point obstruction

**PROVED (four points, including a choice of units).** Fix primary Gaussian prime generators \(\pi_p\equiv1\bmod(1+i)^3\) above split odd rational primes. For any proper closed \(Z\subset T_4\), there exist primitive four-point clusters outside \(Z\), in one unit class for these generators, with arbitrarily large radius \(R\), diameter \(O_Z(R^{3/7})\), all seven cut weights \(2\log n+O_Z(1)\), and \(H_\Sigma=O_Z(1)\), where \(R\asymp_Z n^7\). Hence \(A_4(\delta)\) fails for \(\delta<1/14\), and the primary-generator version of \((L_4^\circ)\) fails. Here fixing the generators before all parameters is essential. The proof makes no claim about an arbitrary preassigned nonprimary convention for the stronger one-unit-class formulation.

*Proof.* Take \(a=x,\ b=x(1-x)/(m^2+x),\ c=1-a-b,\ u=mb\). Then \(abc=u^2\) and \(a+b+c=1\). As \((a,b,c)\) runs through \((x,y,z),(x,z,y),(y,z,x)\), with \(y=x(1-x)/(m^2+x)\), \(z=1-x-y\), define the three monic quartics \(P_1,P_2,P_3\) whose roots are

\[
 i,\quad u+i(a-1),\quad u+i(b-1),\quad u(1+c)/c+ic .
\]

Each quartic has constant imaginary part on the real line. In fact the imaginary parts of the sums of its roots to powers \(1,2,3\) are \(0,0,3(1+c)(u^2/c-ab)=0\). The first equality uses \(a+b+c=1\), the second uses \(a+b-2+1+c=0\), and the third follows from
\(1+(a-1)^3+(b-1)^3+c^3=3(1+c)ab\). Newton's identities finish the assertion.

**COMPUTED (dominance certificate used in the proof).** The logarithmic Jacobian of the three rational functions \(F_j=P_j/\overline P_j\) in \((x,m,s)\), at \((-3,3,7)\), reduces modulo \(101\) with \(i\mapsto10\) to

\[
 \begin{pmatrix}15&12&93\\42&52&27\\81&73&35\end{pmatrix},
 \qquad \det=54\bmod 101.
\]

The check sums \(d(s-r)/(s-r)-d(s-\bar r)/(s-\bar r)\) over the four roots, verifies every denominator is a unit modulo 101, and certifies that the seven roots at this point are distinct up to conjugation. [four_point_rank.py](four_point_rank.py)

**PROVED (completion).** The nonzero Jacobian gives a dominant map from rational parameter three-space to the three-dimensional torus of similarity classes. For any specified proper closed \(Z\), its inverse image is contained in a proper closed subset of parameter space where the map is defined. A polynomial vanishing there is not identically zero. The condition on \((x,m)\) that all its coefficients in \(s\) vanish is a proper closed subset. Thus one can choose rational \((x,m)\) such that all but finitely many \(s\) map outside \(Z\), all seven roots remain nonreal and distinct up to conjugation, and the four phases are generically distinct (include all pair diagonals in \(Z\)). Fix that choice. Let \(d>0\) be an integer with every \(d r\in4\mathbb Z[i]\). Take \(n\equiv1\bmod4\), \(s=n/d\), and \(H_r=n-dr\). Form \(z_0=\prod_r\bar H_r\) and \(z_j\) by replacing \(\bar H_r\) with \(H_r\) exactly at the four roots of \(P_j\). All four Gaussian integers have modulus \(\asymp n^7\); their ratios to \(z_0\) are \(F_j\). All are primary. After division by the primary associate of their Gaussian gcd they remain primary, and after this normalization their remaining prime factors are split odd primes only. Thus all normalized points have the same unit for the fixed primary generators. Because \(P_j\) has constant imaginary part, \(|F_j-1|=O(n^{-4})\).

For the height assertions, the 14 linear forms \(n-dr,n-d\bar r\) have fixed pairwise nonzero differences. Let \(D\) be a positive integer divisible in \(\mathbb Z[i]\) by all these differences. At each Gaussian prime \(\lambda\) at most one form has valuation greater than \(v_\lambda(D)\). Stripping \(\min(v_\lambda(D),v_\lambda(H))\) copies of \(\lambda\) from every form changes each form by a divisor of \(D\), leaving pairwise coprime forms and preserving conjugation at the level of ideals. For the stripped forms, each prime contributes exactly the cut defined by the corresponding root, and the four products have trivial common gcd. Restoring the bounded factors changes the valuation lists by bounded amounts only at the finitely many primes over \(D\), hence moves at most \(O_D(1)\) threshold levels there, and gives a uniformly bounded common gcd. Explicitly, perturbing four exponents by at most \(b\) can change the cut at threshold \(t\) only within distance \(b\) of one of four original exponents. Each of the seven distinct unoriented cuts corresponds to exactly one root (its incidence vector among \(P_1,P_2,P_3\)), and its contact norm is therefore \(\asymp n^2\). The normalized radius satisfies \(R\asymp n^7\), the diameter is \(O(n^3)\), and each pair is separated by four cuts of contact norm \(\asymp n^2\). Its scalar cofactor is bounded by \(O(n^3 n^4/n^7)=O(1)\). These estimates yield the claimed balance, \(H_\Sigma\), and exponent. All \(O\)-constants may depend on the chosen rational parameter line. Since \(3/7<1/2-\delta\), this dependence poses no difficulty for the counterexample with arbitrarily large \(R\) and fixed \(Z\). \(\square\)

**EVIDENCE (scope of audit).** The original round-20 family had examples in more than one unit class. The arithmetic progression above supplies a common class with the primary convention. The geometric statement \(A_4(\delta)\) itself is independent of any choice of Gaussian associates. [Round 20, §4.1](/home/user/Jarnik/docs/claude_review/round20/gauss.md) · [Killer, §1.2](/home/user/Jarnik/docs/claude_review/round20/kill_gauss.md)

## 7. The squarefree degree-eight function-field test

**PROVED using the finite COMPUTED certificates below.** Let `t` be real under conjugation, let $\rho_j=x_j+i y_j$ with $y_j>0$ for $1\le j\le8$, and assume the sixteen numbers $\rho_j,\bar\rho_j$ are distinct. A sign row $s\in\{\pm1\}^8$ specifies the monic polynomial $z_s(t)$ having root $\rho_j$ for $s_j=1$ and root $\bar\rho_j$ otherwise. All such polynomials have the same squarefree norm. Let $\mathcal C$ be any collection satisfying $\deg(z_s-z_{s'})\le4$ for distinct $s,s'\in\mathcal C$. Then $\mathcal C$ has affine rank at most four and cardinality at most eight. Eight occurs. Every eight-row example has the same sign code up to coordinate permutation, row permutation and global sign. In particular this decides the rank question at $D=8,c=0$ for the real model with ordinary complex conjugation. The extremizer assertion classifies sign rows and makes no uniqueness claim for their realizing roots.

### Reduction to three moments

**PROVED.** Newton's identities identify equality of the three leading nonconstant coefficients below the monic term with equality of the first three power sums of the roots. On subtracting two power sums, the real parts cancel coordinate by coordinate. Consequently the cluster condition is precisely
$$
 \sum_j (s_j-s'_j)\operatorname{Im}(\rho_j^k)=0,\qquad k=1,2,3 .
$$
Write $D$ for a matrix of independent row differences, in this section only. Replacing sign rows by their $0,1$ versions scales the differences and leaves their span unchanged.

**PROVED.** Every nonzero $m\in\mathbb Z^8$ in the rational row span of $D$ has positive mass and negative mass both at least two. To see this, form
$$
 P(t)=\prod_{m_j>0}(t-\rho_j)^{m_j}\prod_{m_j<0}(t-\bar\rho_j)^{-m_j}.
$$
The log expansion of $P/\bar P$ at infinity has its first three inverse powers zero by the moment equations. Hence $\deg\operatorname{Im}P\le \|m\|_1-4$. This imaginary part is nonzero since the root multiset lacks the conjugate of every root it contains. Along the real axis the continuous argument of $P$ goes from $-\pi\sum_jm_j$ to zero, with the argument chosen for each monic factor to tend to zero at positive infinity. The intermediate multiples of $\pi$ force at least $|\sum_j m_j|-1$ distinct real zeros of $\operatorname{Im}P$. Thus $\|m\|_1-|\sum m_j|\ge3$. The left side is twice the smaller mass, which is integral. Applying this to two distinct binary rows shows that in each direction at least two bits change. This is the round-30 saturated mass test. [Round 30, Lemma 1.2](/home/user/Jarnik/docs/claude_review/round30/prove.md)

### Independent enumeration

**PROVED.** A cluster of affine rank at least five contains six affinely independent rows. To enumerate possible six-row subsets, start with the empty set, append each of the 256 possible rows to each representative, and require the preceding pair condition. Encode a coordinate by its column of bits across the current rows. Sort the eight columns to take the quotient by coordinate permutations. Minimize their packed integer encoding over row permutations and global complement. This is a complete invariant under these operations. Every admissible set is reached by induction, since removing a row preserves pair admissibility and all extensions of each representative are tried. At six rows keep affine rank five. Reject a remaining matrix if a nonzero vector in $[-2,2]^8\cap\mathbb Z^8$ of positive mass at most one annihilates a basis for $\ker D$. Such a vector belongs to the rational row span of $D$. Each rejection is justified by the mass lemma regardless of whether this box detects every violation.

**COMPUTED.** For one through six rows the pair-admissible orbit counts were $5,14,26,50,65,72$. Of the 72 terminal orbits, 64 had affine rank five and 21 passed the finite mass test. The canonical keys of the 21 classes agree as a set with those in the archived census. All subsequent class indices refer to its zero-based order. The standalone script embeds those representatives and checks their exhaustive coverage. [Independent census](exact_census.py) · [Round 30 census log](/home/user/Jarnik/docs/claude_review/round30/judge_checks/sym8_D8.log)

### The nonlinear test

**PROVED.** For a surviving six-row matrix let $H\in\mathbb Q^{8\times3}$ consist of a column basis of $\ker D$. Denote row $j$ by $h_j^T$. The moment equations give real vectors $u,v,w$ with
$$
 y_j=h_j^Tu,\qquad x_jy_j=h_j^Tv,\qquad 3x_j^2y_j-y_j^3=h_j^Tw .
$$
It follows that
$$
 (h_j^Tu)^4=h_j^T Q h_j,\qquad Q=3vv^T-\tfrac12(uw^T+wu^T).
$$
The restriction of $Q$ to the linear plane $u^\perp$ has rank at most one.

**PROVED.** Two proportional rows $h_j=\lambda h_k$ preclude a realization. A zero row contradicts positivity of $y$. Otherwise positivity forces $\lambda>0$, and the first two moments give $y_j=\lambda y_k$, $x_j=x_k$. The third gives $\lambda^3=\lambda$, so the roots coincide.

**PROVED.** Four distinct affine-collinear rows also preclude a realization. Parametrize the line by $h(t)$ affine in $t$ and set $Y(t)=h(t)^Tu$. If $Y'\ne0$, the third derivative of $Y(t)^4-h(t)^TQh(t)$ equals $24(Y')^3Y(t)$, of constant nonzero sign between the extreme parameter values since $Y$ is positive at both. Four roots contradict Rolle's theorem. Hence $Y$ is constant. Now set $V(t)=h(t)^Tv$, $W(t)=h(t)^Tw$. The quadratic polynomial $3V(t)^2-YW(t)-Y^4$ vanishes at four distinct parameters. Its quadratic coefficient forces $V'=0$, and all four $x=V/Y$ coincide. These two geometric obstructions apply for arbitrary kernel dimension whenever the first three moments agree.

**PROVED.** For the remaining cases take the reduced-echelon nullspace basis from the script. It has a row $(0,0,1)$, so the corresponding $y_j>0$ allows a common positive dilation of the roots to put $u=(a,b,1)^T$. Evaluate the six ternary quadratic monomials at the eight $h_j$. In these cases the resulting rational design matrix has rank six, as checked in the script. Interpolating the six coefficients of $Q$ from its values $(h_j^Tu)^4$ at six independent rows makes its entries polynomials in $a,b$. There remain two interpolation equations. Adjoin
$$
 \det(FQF^T)=0,\qquad
 F=\begin{pmatrix}1&0&-a\\0&1&-b\end{pmatrix},
$$
since the rows of $F$ span $u^\perp$. Compute a Gröbner basis for these three polynomials over $\mathbb Q$ in variables $b,a$.

**PROVED.** A unit ideal excludes every realization. There is another exact rejection when the product of necessarily positive factors has zero remainder. The factors used are the $y_j=h_j^Tu$ and, for $j<k$,
$$
 B_{jk}=(y_j-y_k)^2y_j^2y_k^2+
 \tfrac13(y_kh_j-y_jh_k)^T Q(y_kh_j-y_jh_k).
$$
Indeed $y_kh_j-y_jh_k\perp u$, so the final summand equals $(y_jy_k)^2(x_j-x_k)^2$. Thus $B_{jk}=(y_jy_k)^2|\rho_j-\rho_k|^2>0$ on any allowed configuration. Reducing after each multiplication preserves the product's residue class.

**COMPUTED.** All 21 survivors were excluded by the following certificates. [Exact Gröbner and geometric checks](exact_plane.py)

| Certificate | Archived classes |
|---|---|
| Proportional rows | 3, 5 |
| Four collinear rows | 8, 9, 11, 12, 13, 15, 16, 17, 19, 20 |
| Unit ideal | 0, 1, 4 |
| Zero remainder for a positive product | 2, 6, 7, 10, 14, 18 |

**PROVED, individual check.** In class 0 the rows of $H$ are $(-1,1,3)/2,(1,-1,3)/2,(1,1,1)/2,(1,-1,1),(-1,1,1),e_1,e_2,e_3$. Interpolating at the first six linearly independent quadratic evaluation rows gives, with $d=a-b$, the two residual equations
$$
 d(3d^2-5)=0,\qquad
 3d^4+12d^3-30d^2-20d-45=0 .
$$
At zero the second left side equals $-45$, and at $d^2=5/3$ it equals $-260/3$. These two equations already exclude this class. Combining the exhaustive exclusions proves the rank bound.

### Size and equality

**PROVED.** An affine $r$-plane admits an injective projection onto $r$ coordinate axes. Thus it contains at most $2^r$ binary rows. An eight-row set of rank at most three must have rank three and project to the entire binary three-cube. Every output coordinate on this cube is an affine function taking values in $\{0,1\}$, hence is constant or one input bit or its complement. Indeed varying a nonzero coefficient forces its value to be $1$ or $-1$, and the range on the cube is the sum of the absolute coefficients. Each of the three input directions needs at least four changing output coordinates by the pair condition. They are disjoint, requiring at least twelve coordinates. Consequently any eight-row example must have rank four.

**PROVED.** To bound size in rank four, choose five affinely independent rows as a basis. For every one of the 65 pair-admissible five-row orbits, retain the rank-four cases, intersect its rational affine span with all 256 binary rows, and test every extension of the basis satisfying pair admissibility. Every possible cluster of rank four is included, since it has such a basis.

**COMPUTED.** The maximum was eight. The computation produced three representatives after the five-row bases were put in canonical form. In the byte encoding, bit $j$ is one when sign $j$ is positive. The representatives were
$$
\begin{aligned}
 &(15,51,60,85,170,195,204,240),\\
 &(15,51,85,102,153,170,204,240),\\
 &(15,51,85,105,150,170,204,240).
\end{aligned}
$$
The coordinate permutation sending $0,1,\ldots,7$ to $0,4,1,5,2,6,3,7$ takes the first tuple to the second. The last tuple has identical rows 0 and 7 in the nullspace basis $H$, forcing a repeated root. This proves uniqueness of the oriented sign code, conditional only on the displayed exact finite computations. [Affine basis extensions](exact_maximum.py)

**PROVED, attainment.** One can check the known cube construction on a single rational example. Put $q=(-9/4,-2/3,1/4,8/3)$, $\epsilon=(1,-1,-1,1)$, $r_j=\epsilon_j/q_j+iq_j$. The eight roots $\{r_j,-r_j\}_{j=0}^3$ are distinct up to conjugacy because the four $|q_j|$ are nonzero and distinct and $\operatorname{Re}r_j\ne0$. In the degree-eight polynomial with these roots, conjugate the roots in any one of the eight subsets
$$
 \varnothing,\quad \text{all},\quad
 \{\pm r_a,\pm r_b\} ((a,b)=(0,1),(2,3),(0,2),(1,3)),\quad
 \{r_0,r_1,r_2,r_3\},\quad\{-r_0,-r_1,-r_2,-r_3\}.
$$
Each subset has zero sum of imaginary parts of the first three powers. For $r_j$ these three values are $q_j,2\epsilon_j,3/q_j-q_j^3$. For the choice of all four $r_j$ the sums vanish by direct substitution of the given fractions. For the paired sets odd moments cancel and $\epsilon_a+\epsilon_b=0$. The full set is a union of two such sets. The resulting eight monic polynomials are distinct and satisfy the required degree condition by Newton's identities. Selecting upper-half-plane representatives for the conjugate pairs puts the construction in the stated model. [Round 29, Theorem C](/home/user/Jarnik/docs/claude_review/round29/ff.md)

**COMPUTED, reproduction.** Run `python3 exact_plane.py` and `python3 exact_maximum.py` with SymPy and NumPy installed. The first imports `exact_census.py`, independently regenerates the 21 classes, then verifies each nonlinear exclusion. The second checks the maximum within rank four and the equality cases. Canonical encodings have at most 48 bits. The scripts check coefficient bounds before their fixed-width integer matrix multiplications. Gröbner reductions use exact rational arithmetic. Neither invocation needs the archived files. [Census](exact_census.py) · [Rank certificates](exact_plane.py) · [Size certificates](exact_maximum.py)

**EVIDENCE, relevance.** The realization constraint that worked here is the shared third moment of exact conjugate-root coordinates. It is available at the function-field place at infinity where the first three coefficients agree. The proof supplies no counterpart that forces three exact moment equations from a short arc at one real radius. It therefore settles the finite function-field test and gives no improved growth bound for integer circles.
