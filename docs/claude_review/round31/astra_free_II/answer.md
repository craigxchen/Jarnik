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


# Part II — Minimum vectors, specialization, grid determinants and the five-point index

**PROVED with the COMPUTED certificates below.** The sharp degree-eight divisor count is **12 up to sign**. A different rational template has nine such divisors. Its largest common-moment cluster has five marks, but only three can have every edge on the minimum shell. Seven divisors survive *simultaneously* on an infinite sequence of integer parents, with no four-clique using their associated near-real ratios. A six-row grid minor also answers the nonvanishing question affirmatively, but its forced exponent on a three-three cut is zero. The strongest asymptotic result here concerns the valuation index. In a normalized cluster of size \(m\) and log norm \(W\le b m\log m\), only \(O_b(m^3)\) ordered five-tuples can have index greater than one. This is an obstruction to using the cover on those clusters. On genuine arcs, the Baker–Wüstholz estimate from the parallel report removes the growth-scale assumption. A positive proportion of the marks can be kept with index one on every five-tuple. An **index-one-only** five-point hypothesis therefore suffices even for uniformity.

**HEURISTIC — research judgment.** I would carry forward the common-parent triangle problem and the index-one five-point problem. I have no saving in the lattice-point growth exponent and no family with a growing number of eligible primitive divisors. The finite examples establish that many separate divisor labels can coexist relative to the largest polynomial clique, and that the leading constant at specialization requires an arithmetic check of its own.

## II.1. The exact degree-eight count

**PROVED.** For eight distinct roots \(\rho_i=x_i+iy_i\) in the upper half-plane, write
\[
 h_i=(y_i,x_i y_i,3x_i^2 y_i-y_i^3).
\]
The greatest possible number of ternary moment-kernel vectors of weight four **up to sign** is **12**, or 24 oriented vectors. Equality holds precisely for configurations supporting the eight-row cube sign code of Part I. This is uniqueness of the code, with room for more than one choice of roots. Each vector up to sign gives a different support and thus a different squarefree real norm divisor.

**PROVED.** The three moment equations are exactly \(\sum_i m_i h_i=0\). The balance lemma of Part I for total weight four forces positive and negative mass two each. The ternary vectors therefore correspond to four distinct indices with \(h_i+h_j=h_k+h_l\). These points are distinct since \(y>0\) and \(x=(xy)/y\). Two different pairings on one four-support cannot both give equal sums, by subtraction. If repeated factors are allowed, one must also count every midpoint \(2h_i=h_j+h_k\). A bound that the minimum weight is *at least* four does not guarantee that any vector of weight four exists.

**PROVED.** No four distinct moment points lie on an affine line. A line with height \(Y\) varying has its other coordinates \(V,W\) affine in \(Y\). Its intersections solve \(Y^4+YW-3V^2=0\), a monic quartic with no cubic term. Four positive roots would have positive sum, contradicting Vieta. On a line of fixed positive height the surface is the nondegenerate parabola \(3V^2=YW+Y^4\), with no three collinear points.

**PROVED — exhaustive reduction.** Choose a generic linear functional and label the moment points in increasing order. For \(i<j<k<l\) the only possible equality of pair sums on those four points is \(h_i+h_l=h_j+h_k\). There are 70 candidate rows \(e_i+e_l-e_j-e_k\). For repeated-factor minima add the 56 rows \(e_i+e_k-2e_j\), \(i<j<k\). Given a span \(E\) of these relations, the formal point for label \(i\) is the image of \(e_i\) in \(\mathbb Q^8/E\). Any actual moment configuration satisfying those relations is a linear image of the formal points. Discard a span that identifies two formal points or puts four on a line. These exclusions persist in every extension. The formal affine dimension is \(7-\dim E\) because the rows all sum to zero. Thus no span of rank six or above can occur.

**PROVED — enumeration certificate.** Start at the zero span. At each node reduce all candidate rows modulo a basis of the span. Count the zero images. Group nonzero images by proportionality and extend by the first index in each group if that index exceeds the last basis index used. Every desired span is visited: take its lexicographically greedy independent basis. If a candidate for the next selected basis vector is proportional modulo the current span to an earlier-indexed row, that row also belongs to the target span, in contradiction to the greedy choice. At rank four the zero images and a single proportionality group are exactly the rows in the corresponding rank-five extension, so one may skip extensions containing fewer than the target count. This argument does not require any particular chamber of the ordering to be feasible. Extra chambers only enlarge the census. The calculations modulo 769 are exact for these ranks. Every candidate row has Euclidean norm at most \(\sqrt 6\), and differences of two point labels used to test affine dependences have norm \(\sqrt 2\). Nonzero minors of the stacked rows have size at most seven (all rows have sum zero), and absolute determinant at most \((\sqrt 6)^7<769\) by Hadamard. Ranks of the quotient incidence tests and collinearity tests therefore agree with those in characteristic zero.

**COMPUTED — ternary case.** With a threshold of 12 the numbers of processed nodes at ranks 0 through 5 were \(1,70,1799,20640,94190,50\); only two rank-four nodes and the 50 retained rank-five nodes reached threshold. Affine classification of the rank-five formal point sets gives just three classes. Orient an affine triangle of maximum area, map it to \((0,0),(1,0),(0,1)\), sort the eight rational vectors, and minimize over such choices. This is a complete affine invariant. In the order of first appearance the class counts are 28, 18, 4. Representatives, given as zero-based lexicographic indices among the 70 quadruples, are
\[
 (0,5,12,14,16),\quad(0,5,12,14,18),\quad(2,7,14,23,27).
\]
The first two give zero remainders for positive products in the exact three-moment test of Part I, §7 (last factors correspond respectively to pairs \((0,3)\) and \((1,2)\)). For the third, the rows of its nullspace matrix are
\[
 (2,1,-2),(2,0,-1),(1,1,-1),(0,2,-1),(2,-1,0),e_1,e_2,e_3.
\]
The height parameters can be normalized to \((a,b,1)\). The two quadratic interpolation consistency equations from §7, up to a nonzero factor, are
\[
 (a-b)^2(b-1)(2a+b-1)=0,\quad
 (b-1)(2a+b-1)(2a^2-2ab-2a+b^2+1)=0 .
\]
Here \(2a-1>0\), \(b>0\). If \(b\ne1\) the first equation gives \(a=b\) and the last factor of the second becomes \((a-1)^2\), a contradiction. When \(b=1\), rows 3, 6, 7 are collinear at the same height one, impossible by the fixed-height parabola. All three classes are excluded.

**PROVED using the COMPUTED census.** For both rank-four survivors the formal points, in the nullspace basis, are exactly \((a,b,c,1-a-b-c)\) for \(a,b,c\in\{0,1\}\), differing only in order. Thus the moment points are an affine image of the vertices of a three-cube. Its six four-subsets formed by taking any two opposite pairs and its two parity tetrahedra have the same cardinality and centroid. Their eight indicator rows give eight polynomials with equal first three moments, hence a critical cluster. Conversely the Part I unique eight-row code gives precisely twelve different four-supports of its pair differences, each occurring twice as an edge, and four other edges of weight eight. The rational cube example in Part I realizes twelve up-to-sign primitive divisors. The 28 unordered pairs therefore give twelve distinct degree-four divisors.

**COMPUTED — reproduction.** Compile with `c++ -O3 -std=c++17 ii_min_flat.cpp -o ii_min_flat`. Run `./ii_min_flat 12 > ii_min_cases.txt`, then `python3 -B ii_min_verify.py`. Python uses SymPy as in Part I. [Flat enumeration](ii_min_flat.cpp) · [Affine classes](ii_min_class.py) · [Plane interpolation](ii_min_plane.py) · [Checks](ii_min_verify.py)

**PROVED using a second COMPUTED census — allowing repeated factors.** The maximum number of *all integral kernel vectors of weight four*, including midpoints, is again 12 up to sign, with the same equality condition. On the 126-row run the numbers of processed nodes at ranks 0 through 5 are 1, 126, 5,887, 124,292, 1,086,867, 3,375. Ten rank-four nodes and the 3,375 rank-five nodes reach threshold. After affine classification there are 59 planar possibilities. Computing their actual number of midpoint and four-index relations without presuming that the displayed ordering is feasible gives
\[
\begin{array}{c|rrrrrrr}
\#\text{ relations}&12&13&14&15&17&18&19\\
\#\text{ shapes}&30&15&9&1&2&1&1 .
\end{array}
\]
The positive-product test of §7 excludes 56 of them. The three exceptions are excluded as follows, with row indices starting at zero and the same normalized height parameters as above.

* For \((0,8,19,23,33)\), the two excess interpolation equations have factors \((a-1)(a-b)^2(2a+b-1)\) and \((a-1)(a+b-2)^2(2a+b-1)\), with \(2a+b-1=y_4+b>0\). They force \(a=1\), and the midpoint triple \(4,5,7\) then has constant height.
* For \((6,13,20,24,27)\), an interpolation equation gives \((a-1)^2(b-1)(a+b-1)=0\) with \(2(a+b-1)=y_0+1>0\). At \(a=1\) the midpoint triple \(3,5,7\) has constant height, and at \(b=1\) the triple \(4,6,7\) does.
* For \((3,14,21,29,35)\), the first equation gives \((a-b)(a-5b+2)(a-2b+1)^2=0\), with \(a-5b+2=-y_0-1<0\). At \(a=b\) the second gives \((b-1)^3(2b-1)=0\). At \(a=2b-1\) it gives \((b-1)^3(3b-1)=0\). Since \(2b-1=y_2>0\), both force \(b=1\), at which the midpoint triple \(2,6,7\) has constant height.

Each constant-height triple contradicts the parabola argument. The ten rank-four records consist of the two cube records above and eight containing seven distinct points of the form \(C,C\pm U,C\pm V,C\pm(U+V)\). The latter configuration cannot occur on the moment surface. Translating the real part of every root and dilating all roots by a positive constant act linearly on moment triples and can put \(C=(1,0,-1)\). If \(C\pm(u,v,w)\) are on the surface, substitution into \(3V^2=YW+Y^4\), subtraction and addition give \(w=-3u-4u^3\) and \(v^2=u^2(1-u^2)\). For displacements \(U,V,U+V\), with height components \(u,u',u+u'\), the cubic identities force \(uu'(u+u')=0\). The displacement with zero height then vanishes in all three coordinates, a contradiction. Thus only the cube records survive for the extended count as well. Run `./ii_min_flat 12 all > ii_min_allcases.txt` followed by `python3 -B ii_min_all_verify.py`. [Full weight-four verification](ii_min_all_verify.py)

**PROVED using a COMPUTED census — the first gap.** An admissible configuration with eleven ternary directions already has twelve. Repeating the 70-row enumeration at threshold eleven produces five rank-four records, each of whose formal quotients carries twelve parallelograms. At rank five, only the classes represented by row-index sets `(0,5,12,17,30)` and `(1,10,17,24,29)` have exactly eleven parallelograms. The positive-product certificate from §7 rules out each. Compute the extra records with `./ii_min_flat 11 > ii_min_gap.txt` and verify them with `python3 -B ii_min_verify.py --gap`.

**PROVED using a COMPUTED instance — a different question about shortest vectors.** Take eight roots of the form \(x+i\), with \(x\in\{1,5,8,10,12,16,22,23\}\). Equal moment sums are precisely equal-cardinality multisets with the same sum and sum of squares. Among the multisets of cardinality \(k\), the number of unordered pairs of distinct multisets with matching signatures is \(0,0,0,0,0,0,13\) for \(k=1,\ldots,7\). The exhaustive script groups multisets by their two integer moments. The first collisions have disjoint supports because a common element would give a smaller collision. Thus this kernel has minimum 14 and thirteen shortest directions up to sign. The sharp twelve count applies to the requested degree-four shell, or to actual shortest vectors conditional on the lattice attaining weight four. [Multiset certificate](ii_min_longer.py)

## II.2. Unrelated polynomial divisors that survive integer specialization

**PROVED using the COMPUTED certificates below.** In degree eight there is a rational root configuration with seven signed ternary directions of weight four up to reversal and with maximum cluster size four. Five of its seven directions simultaneously give divisors counted by \(D_B(N)\) along one infinite sequence of odd split-supported integer parents, for *each* fixed \(B>0\). For all sufficiently large specializations no four lattice points realize all five selected labels as near-real ratios. This is a fixed lower bound and gives no violation of the growth conjecture for \(D_B\).

**COMPUTED — configuration.** Number the roots from 0 through 7. Their pairs \((x,y)\), denoting \(x+iy\), are
\[
 (-12,4),(-8,6),(-3,16),(0,10),(0,12),(3,16),(8,6),(12,4).
\]
Put \(u_j=(y_j,x_jy_j,3x_j^2y_j-y_j^3)\). A weight-four ternary solution must have two positive and two negative entries by the mass bound. Thus it comes exactly from an equality \(u_i+u_j=u_k+u_l\) with distinct indices. All pair-sum coincidences, abbreviated by their indices, are
\[
 04=13,\qquad 05=27,\qquad 06=17,\qquad 15=26=34,\qquad 36=47.
\]
Count all three pairwise equalities on the length-three chain. For a cluster of sign rows, the corresponding subsets must have one and the same sum of the \(u_j\), and this is sufficient by Newton's identities. Enumerating the 256 subsets gives fibre-size frequencies \(120,32,8,12\) for sizes one through four. The elementary checks use only integer arithmetic. [Certificate](ii_template.py)

**PROVED — the specialization cost.** Choose the five equalities \(04=13,06=17,15=34,26=34,36=47\). For any one of them take \(B_v(n)\) to be the product of \(n-\rho_j\) for the two indices on the left and \(n-\bar\rho_j\) for the two on the right. The equality of the three moments makes its first three non-leading coefficients real, so \(b_v=|\Im B_v(n)|\) is constant. In the listed order the five constants are \(4800,9600,13200,13200,4800\), evaluated at \(n=0\). Put
\[
 f_j(n)=(n-x_j)^2+y_j^2,\quad N_0(n)=\prod_j f_j(n),\quad
 g_v=\gcd(\Re B_v(n),b_v),\quad d_v=\mathrm N B_v(n)/g_v^2,\quad t_v=b_v/g_v.
\]
For the natural monic parent \(N_0\), the ratio \(d_v/(16t_v^2\sqrt{N_0})\) tends to \(1/(16b_v^2)\), so passage to integer polynomials alone loses the desired cutoff. Take instead \(N=\operatorname{lcm}_v d_v\) and \(G=N_0/N\). Then the ratio is multiplied by \(\sqrt G\).

**COMPUTED — a simultaneous CRT certificate.** For each prime in the first row below impose \(n\equiv n_p\pmod{p^{k_p}}\) with \(k_p=1+v_p(N_0(n_p))\). The last row is the valuation of \(G(n_p)\).
\[
\begin{array}{c|rrrrrrrr}
p&2&3&5&13&37&41&53&509\\
n_p&2&0&10&5&12&17&0&238\\
k_p&16&3&12&5&3&3&3&3\\
v_p(G)&15&2&9&1&1&1&1&1
\end{array}
\]
All five \(d_v(n_2)\) are odd. The CRT representative printed by the script has 43 digits and gives a 666-digit odd parent with five distinct candidates. At that parent it verifies the upper cutoff even at \(B=5\) and obtains \(t_v=(2,12,11,11,2)\). [Certificate](ii_template.py)

**PROVED — passage to the progression.** Each displayed valuation persists on the progression: \(p^{k_p}\) fixes the valuations of \(N_0,\mathrm N B_v\) since the latter divides the former, and of the gcd \(g_v\) since \(g_v^2\mid\mathrm N B_v\). Thus \(G\) is divisible by \(L=306439728192000000000>256\max_v b_v^4\). Each \(d_v\) is an odd primitive sum of squares, hence \(N\) is odd and split-supported. Uniformly on the progression,
\[
 \frac{d_v}{16t_v^2\sqrt N}
 =\frac{\sqrt G}{16b_v^2}\frac{\prod_{j\in\mathrm{supp}(v)} f_j}{\sqrt{N_0}}>1
\]
eventually because the last quotient tends to one. The \(G\)'s are bounded above independently of \(n\). Indeed each of the eight \(f_j\)'s occurs on some selected support, and each such \(f_j(n)\) divides \((\operatorname{lcm}_v b_v)^2 N(n)\). The quadratics \(f_j\) are pairwise coprime over \(\mathbb Q[t]\), so pairwise gcds of their integer values divide fixed nonzero resultants. In particular \(\prod f_j/\operatorname{lcm}_j f_j\) is bounded (it divides the product of the pairwise gcds). This gives \(N\asymp n^{16}\), \(d_v/\sqrt N=O(1)\) and the upper cutoff for every fixed \(B>0\) eventually. Distinct supports yield five distinct \(d_v\)'s for large \(n\), since there are finitely many possibilities for \(g_v\mid b_v\) and the different products of the monic, distinct irreducible \(f_j\)'s cannot be proportional as polynomials. The floor condition follows either from Part I or because \(t_v\) stays bounded while the reduced real part tends to infinity.

**PROVED using a finite check — meaning of unrelated.** No polynomial subset-sum fibre supplies more than three of the selected labels. More intrinsically, of their ten distinct triples exactly one has a signed sum of direction vectors zero. In zero-based order of *all seven* equalities (break the chain into \(15=26,15=34,26=34\)), it is \(v_0-v_2+v_6=0\). Gaussian unique factorization of the monic linear factors shows that any of the other three-direction products \(\prod (B_v/\bar B_v)^{\pm1}\) is a nonconstant rational function. Hence these other triangles still fail for all sufficiently large integer \(n\). Five distinct edges among four points would contain two triangles. For the near-real edge labels used in (1), each such triangle imposes signed phase closure, since \(B_v/\bar B_v\) is fixed up to inversion from its near-square divisor. This includes any realization as separation divisors of four points in the shortened one-unit formulation of Part I: their positive small offset is determined by the floor of the square root. Without the one-unit convention the assertion concerns the specified ratios up to inversion. The construction does not bound the size of integer clusters involving further divisors.

**EVIDENCE — growth.** This template supplies a fixed finite number while \((\log N/\log\log N)^2\) tends to infinity. Even an increasing sequence of polynomial examples would need a lower bound on the number of distinct eligible directions at one specialization in terms of that specialization's actual \(\log N\), including denominators, CRT modulus and angle cutoff. The existence statement of the bridge theorem by itself provides none of those uniform estimates.

**COMPUTED with the same PROVED progression argument — larger count.** The eight roots \((-23,4),(-11,10),(-8,16),(2,16),(5,10),(5,12),(13,6),(17,4)\) give nine ternary directions and fibre-size frequencies \(104,32,16,8\) for sizes \(1,2,3,5\) respectively. Seven of their directions simultaneously specialize to \(D_B(N)\) divisors on an infinite progression. No polynomial fibre accounts for all seven. The check at a 752-digit parent gives reduced offsets \(123,24,3,84,12,11,2\). In the script's order choose directions \(0,1,2,3,5,7,8\). The local tuples \((p,n_p,k_p,v_p(G))\) are \((2,3,16,15),(5,314,13,10),(41,12,3,1),(53,16,3,1),(349,131,3,1),(541,242,3,1),(953,544,3,1)\), giving \(L=125119585606720000000000>256(98400)^4\). There are also two midpoint relations in each configuration, taking the full weight-four counts up to sign to nine and eleven respectively. They are verified by testing all distinguished midpoints and all pairs of distinct endpoints. This improves the exhibited finite divisor count. The previous example gives the stronger four-vertex incompatibility. [Second certificate](ii_template.py)

**COMPUTED — the clique question with every edge on the minimum shell.** The largest such cliques are of sizes **4, 3, 3** in the cube, the seven-direction configuration and the nine-direction configuration, respectively. Requiring the edges also to have distinct divisor labels gives the same maxima. The full common-moment fibres have maxima **8, 4, 5**. Restricting to just the five or seven directions selected for integer specialization still gives maximum three in the last two configurations. The checker enumerates every subset of every fibre and tests each pair difference and its support. This also exhausts cliques among ternary directions anywhere in the integer lattice: if all pair differences are ternary, every coordinate has range at most one, so an integral translation places the clique in the binary cube. Thus the proposed Sidon description needs both an edge-length condition and a distinctness condition at the polynomial endpoint. [Minimum-clique certificate](ii_cliques.py)

**PROVED — the seven selected integer ratios themselves have no four-clique.** This holds for all sufficiently large specializations on the second progression. Associate to each selected label the rational function \(r_v=B_v/\bar B_v\), specialized at the integer parameter. Were four vertices to have all their pair ratios among these seven functions up to inversion for infinitely many specializations, fix the six signed edge labels along an infinite subsequence. The triangle identities would hold identically as rational functions. Normalizing one vertex to one, the other three now have well-defined exponent vectors in the free abelian group on the basic ratios \((t-\rho_j)/(t-\bar\rho_j)\). Freeness follows by taking valuations at the distinct upper roots. Every pair difference is one of the allowed nonzero ternary directions up to sign, yielding a four-clique in the just-enumerated polynomial graph. Its maximum is three. The finite-choices argument also shows that any polynomial triangle persists on specialization outside finitely many parameters. This controls edges drawn from these seven labels at the specified parent; additional eligible divisors of that parent are possible.

### Two obstructions to obtaining growing examples by symmetry or composition

**PROVED.** For any even \(D=2m\), suppose the \(D\) distinct upper roots lie on a common ellipse with real center and axes parallel to the real and imaginary axes. There are at most two ternary directions of weight \(m\), up to reversal. Indeed, parametrize points as \(x_0+Aw+B/w\), \(|w|=1\), with real \(A=(a+b)/2\), \(B=(a-b)/2\) for positive semiaxes \(a,b\). The imaginary parts of the powers through degree \(m-1\) are triangular real combinations of \(\Im w^j\) with nonzero diagonal \(A^j-B^j\). Thus the moment equations for any signed root selection are equivalent to the unit-circle equations. For a selection of \(m\) roots upstairs on that unit circle, let \(F(z)=z^m+\sum_{j=0}^{m-1}a_jz^j\) be the corresponding polynomial. Its interior coefficients are real by Newton's identities, and \(a_0\) is nonreal because a real polynomial would contain every selected conjugate. The unit-root identity \(a_j=a_0\bar a_{m-j}\) forces every interior coefficient to vanish. Consequently every selection is a whole fibre \(w^m=\eta\). Each pair of conjugate selections uses \(2m\) roots, distinct such pairs are disjoint, and the full upper and lower configuration has \(4m\) roots.

**PROVED.** Given a fixed rational Gaussian template with distinct nonreal roots \(\rho_j,\bar\rho_j\), there is one \(a\in\mathbb Z\) such that for *every* positive integer \(k\), substitution \(t\mapsto t^k+a\) introduces no new half-degree constant-offset divisors over \(\mathbb Q(i)[t]\). Choose distinct odd split rational primes \(p_j\) avoiding the denominators and a Gaussian prime \(\pi_j\mid p_j\) for each upper root. The map \(\mathbb Z/p_j^2\mathbb Z\to\mathbb Z[i]/\pi_j^2\) is an isomorphism by equal cardinalities and contraction; rational denominators prime to \(p_j\) are units there. Choose a residue \(a\pmod {p_j^2}\) for which \(v_{\pi_j}(a-\rho_j)=1\), and combine the conditions by CRT. At the conjugate root the analogous valuation holds at \(\bar\pi_j\). Every polynomial \(t^k+a-\rho_j\) and its conjugate is now Eisenstein in the relevant discrete valuation ring. The substituted polynomials are distinct monic irreducibles, so every monic divisor is exactly a union of whole fibres and equals \(F(t^k+a)\) for a divisor \(F\) downstairs. The degree is multiplied by \(k\), and \(\Im F(t^k+a)\) is constant precisely when \(\Im F\) is constant. This gives an obstruction for rational-polynomial divisors under these substitutions. The roots upstairs need not be rational Gaussian and the proof leaves split substitutions with partial fibres open.

## II.3. A grid determinant with automatic nonvanishing and no pair factors

**PROVED — normalization.** Write \(z_1,\ldots,z_m\) for the circle points and let \(A=\{\alpha\in\mathbb Z_{\ge0}^5:|\alpha|\le d\}\), \(n=\binom{d+5}5\), \(V=\sum_A|\alpha|=5\binom{d+5}6\). For a design \(G\) of \(n\) grid rows write \(F_G=\det((z_{g_1},\ldots,z_{g_5})^\alpha)\), a polynomial in the \(m\) *original* coordinates. Every term is homogeneous of degree \(V\). A common shift acts on the columns by a unitriangular change of monomial basis, so if \(F_G\ne0\) symbolically its diagonal vanishing order is exactly \(V\). At a split prime use entry costs \(c_{g\alpha}=\sum_h \alpha_h a_{g_h,p}\). If \(L_p,U_p\) are the minimum and maximum assignment costs at \(\pi_p\), its raw contribution to the forced norm exponent is \(e_p V+L_p-U_p\). This accounts for both conjugate places. Set \(\kappa_{\rm raw}=\sum_{p\mid N}(e_pV+L_p-U_p)\log p/(WV)\). The archimedean single-determinant bound on a critical arc has norm exponent \(V/2\) in \(N\), up to its coefficient/size factor. If equal monomials are collected after inserting the actual design, the new contribution is at least as large. These exponents describe forced monomial divisibility. They need not equal the determinant's actual valuation.

**PROVED — finite-size normalization.** For an ordinary univariate Vandermonde on \(n\) marks and a binary cut with \(s\) marks on one side, the collected cost divided by \(V\) is \((\binom{s}{2}+\binom{n-s}{2})/\binom n2\). An exactly half-half cut at even \(n\) therefore gives \((n-2)/(2(n-1))\), tending to \(1/2\); independently fair bits have expectation exactly \(1/2\).

**PROVED — an entirely non-pair design.** For \(d=1\) take six points \(z_0,\ldots,z_5\) in strict circular order in angular width \(<\pi/3\). Use the six grid rows \((z_i,z_{i+1},\ldots,z_{i+4})\), indices modulo 6. Every row contains five distinct marks. Up to a sign its determinant in \(1,x_1,\ldots,x_5\) is
\[
 \prod_{k=1}^5 \left(\sum_{j=0}^5 \xi^{kj} z_j\right),\qquad \xi=e^{2\pi i/6}.
\]
Indeed the six shifted columns form a circulant up to permutation and sum to the constant column times \(\sum_j z_j\). Diagonalize the circulant by the discrete Fourier basis and cancel the trivial character as a polynomial identity. None of the five factors is a pair difference or vanishes. For \(h=\xi^k\ne1\), summation by parts gives
\[
 (1-h)\sum_{j=0}^5h^j z_j=\sum_{j=1}^5(z_j-z_{j-1})(h^j-1).
\]
The nonzero \(h^j-1\) occupy a sector of width at most \(\pi-2\pi/6\). The chord increments of consecutively ordered points occupy a sector of width strictly less than \(\pi/3\) by the half-angle formula. The nonzero summands thus lie in a common open half-plane, and there is at least one.

**PROVED — exact cost.** A cyclic diagonal assigns each of the five variable columns to a distinct row at the *same* chosen mark \(z_j\); the remaining row gets the constant column. Thus \(L_p=5\min_j a_{jp}\), \(U_p=5\max_j a_{jp}\), and
\[
 \kappa_{\rm raw}=1-\frac1W\sum_p\big(\max_j a_{jp}-\min_j a_{jp}\big)\log p .
\]
Expanding the Fourier product also gives a nonzero fifth power of every single variable, so symbolic cancellations do not change the width. On exactly three-three cuts of these six marks, \(\kappa=0\). On six independent fair binary entries \(\kappa=1/32\) in expectation. On any shortened actual profile the valuation range sum is at least the separation \(\log d_{ij}\ge W/2+\log16\) of any two marks, hence \(\kappa\le1/2-\log(16)/W\). A determinant can thus have automatic nonvanishing with no pair factors and have less forced divisibility than the Vandermonde.

**PROVED — every total degree.** For each coordinate \(h\) choose nodes \(s_{h,0},\ldots,s_{h,d}\) in the arc, disjoint across coordinate lists if distinct marks per row are required; choose two fresh nodes \(u,v\). Start with the simplex rows \((s_{1,b_1},\ldots,s_{5,b_5})\) for \(b\in A\). Replace the rows indexed by \((d,0,0,0,0),(d-1,1,0,0,0)\) by \((u,v,s_{3,0},s_{4,0},s_{5,0})\) and its swap in the first two slots. In the unitriangular Newton basis
\[
 N_\alpha(x)=\prod_{h=1}^5\prod_{j<\alpha_h}(x_h-s_{h,j}),
\]
the columns indexed by the two removed rows vanish on all surviving rows because the two indices are maximal. The remaining simplex block is triangular in a linear extension of the coordinate order. Its determinant is the product of its diagonal entries, all products of pair differences. The final \(2\times2\) block has determinant, up to sign,
\[
 P(u)P(v)(u-v)(u+v-s_{1,d-1}-s_{2,0}),\qquad P(t)=\prod_{j<d-1}(t-s_{1,j}).
\]
The four-mark factor is nonzero: two chords of the same circle with the same nonzero midpoint coincide, and no antipodes lie in the short arc. Hence there are precisely \(V-1\) pair factors, one non-pair factor, and automatic nonvanishing. Under independent fair cuts the collected \(\kappa\) is \(((V-1)/2+1/8)/V\) because a primitive linear factor contributes one exactly when all its marks have the same bit. For uniform exactly half-sized cuts on all \(m\) marks replace \(1/2,1/8\) by \(2\binom{m/2}{2}/\binom m2,2\binom{m/2}{4}/\binom m4\). For an actual profile every one of the linear factors contributes at most \(1/2-\log(16)/W\) after normalization, by the two-mark comparison above.

**COMPUTED — actual assignments for a balanced model.** On 32 model rows indexed by \(\mathbb F_2^5\), take the 31 nonconstant Walsh cuts with *equal model weights*. Every column has 16 zeros and 16 ones. For the cyclic example choose the zero and five basis vectors as labels; every exponent is zero. For the two-shell-row modification index node \(s_{h,j}\) by \(5j+h-1\), \(u,v\) by \(5(d+1),5(d+1)+1\), using binary encoding of the five bits. The script checks all assignment optima using integer dual inequalities. It also evaluates the cyclic determinant on the six rational unit-circle parameters \(j/24\) after a common integral scaling, and verifies the nonzero Fourier product by a polynomial resultant. The collected column counts monochromatic supports of the linear factors. These are model calculations, with no assumption that 31 actual split primes have identical logarithms. [Exact cost checker](ii_grid_check.py)

| \(d\) | rows \(n\) | \(V\) | raw \(\kappa\) | collected \(\kappa\) |
|---:|---:|---:|---:|---:|
| 1 | 6 | 5 | 61/155 | 67/155 |
| 2 | 21 | 35 | 464/1085 | 513/1085 |
| 3 | 56 | 140 | 127/310 | 522/1085 |
| 4 | 126 | 420 | 1259/3255 | 524/1085 |

**PROVED — obstruction to an unspecified Gram minor.** For any \(d<m\), the original lower simplex design always gives an invertible minor, even using the same node list in each coordinate if needed. Its determinant is exactly the product of \(V\) pair differences by the Newton argument. Thus the minimum of \(\kappa\) over all nonsingular minors on the full grid is \(\le1/2-\log(16)/W\) on every normalized profile having the pair bound in Part I. If five distinct marks per grid row are wanted, disjoint coordinate lists give the same obstruction once there are at least five times the list length, in particular for degree sublinear in the number of marks. This holds also with split prime powers, since a pair factor's two valuations forced by the monomials sum to \(e_p-|a_{ip}-a_{jp}|\). Purposely choosing the displayed non-pair minor does not rescue this exponent. This argument does not classify designs tailored to a sparse prime profile.

**PROVED — a separate cap covering every total-degree minor on fair cuts.** At a uniform independent fair cut \(S\) put \(m_S=\min_{\beta\in\mathrm{supp} F_G}\sum_{i\in S}\beta_i\). At the cube vertex with zeros on \(S\) and ones off \(S\), the expansion of \(F_G\) has multiplicity at least \(m_S+m_{S^c}\): it lies in the \(m_S\)-th power of the zero-coordinate ideal and, by translation invariance, in the \(m_{S^c}\)-th power of the ideal of the other differences \(z_i-1\). These ideals involve disjoint variables in the Taylor expansion. The average multiplicity on the Boolean cube of a nonzero polynomial of degree \(V\) is at most \(V/2\). For a proof write the polynomial as a degree-\(j\) polynomial in its last variable, with leading coefficient \(H\). At a fixed previous tuple where \(H\) has multiplicity \(r\), apply a Hasse derivative of order \(r\) in the previous variables that makes the leading coefficient nonzero. The resulting univariate degree-\(j\) polynomial has total multiplicity at its two test points at most \(j\), and each original multiplicity is bounded by \(r\) plus the corresponding univariate multiplicity. Induction gives an average at most \((V-j)/2+j/2\). This proves \(\mathbb E(m_S+m_{S^c})\le V/2\), which covers the collected and hence the raw assignment lower bound. Independent fair *rows* at each cut are a different quantifier from a specified sparse, balanced family of cuts.

**PROVED — Laurent bookkeeping.** If the columns are Laurent monomials let \(T=\sum_{\alpha,h}\alpha_h\), \(Q=\sum_{\alpha,h}(-\alpha_h)_+\). Realizing negative entries by Gaussian conjugates multiplies *every term* of the evaluated Laurent determinant by \(N^Q\). The squared archimedean scale is \(N^{T+2Q}\Delta^{2v}\), up to the coefficient factor, where \(v\) is the vanishing order at the all-one torus point and \(\Delta\) bounds the angular width. The forced norm exponent at \(p\) is \(e_p(T+2Q)-\mathrm{width}_{a_p}(\mathrm{supp}F)\) after symbolic collection. Thus the missing inequality still compares the width sum against \(vW/2\). Under independent fair cuts the recorded mean-width theorem gives \(2\mathbb E\,\mathrm{width}(F)\ge v\). For proofs that adapt the Laurent space to sparse profiles these results leave the possibility open. [Mean-width theorem](/home/user/Jarnik/docs/claude_review/threshold/construct.md)

## II.4. How many five-tuples actually admit a nontrivial cover?

**PROVED — normalized counting theorem.** Write \(a_i=(a_{ip})_p\in\prod_p[0,e_p]\cap\mathbb Z^\omega\), \(w(v)=\sum_p|v_p|\log p\), and \(W=\sum_p e_p\log p\). Use the one-unit reduction of Part I, with the \(m\) phases \(2\langle a_i,\phi\rangle\), up to a common rotation, in an arc of angular width \(\Delta\le \tfrac12 e^{-W/4}\). Assume the axial monomial gap
\[
 \operatorname{dist}(\langle v,\phi\rangle,\tfrac\pi2\mathbb Z)
 \ \ge\ \arcsin(e^{-w(v)/2})\quad (0\ne v\in\mathbb Z^\omega).
\]
If \(W\le b m\log m\), then all but \(O_b(m^3)\) ordered five-tuples of distinct rows have index **one**. For each fixed odd prime \(q\), only \(O_q(m^3)\) tuples lose rank modulo \(q\), **without** the hypothesis on \(W\). No five-tuple loses rank modulo two. Here the index is that of the subgroup of \(\mathbb Z^4\) generated by the prime columns of the four differences from one anchor. Dividing a tuple by its own Gaussian gcd leaves those differences unchanged.

**PROVED — real affine independence.** The gap on \(a_i-a_j\) gives \(w(a_i-a_j)>W/2\) for \(i\ne j\). Indeed its sine in absolute value is at most \(\Delta/2\), while the gap gives a lower bound \(e^{-w(a_i-a_j)/2}\). If the full set of rows were affinely dependent, separate the two sides of a dependence into disjoint probability distributions with the same barycenter. Draw rows independently from the two sides. At a prime, their normalized coordinates \(X,Y\in[0,1]\) have the same mean \(\mu\), and \(\mathbb E|X-Y|\le\mathbb E(X+Y-2XY)=2\mu(1-\mu)\le1/2\). Summing with weights \(e_p\log p\) gives a cross pair at distance at most \(W/2\), a contradiction.

**PROVED — a bounded minor for every five-tuple at the growth scale.** At each prime choose \(s\ge m/2\) rows in one closed half of its exponent interval. The weighted sum of the coordinate ranges on these rows is at most \(W-e_p\log p/2\). Each integer threshold in a range separates at most \(s^2/4\) unordered pairs, so the mean pair distance is at most \(s(W-e_p\log p/2)/(2(s-1))\). Consequently \(e_p\log p\le4W/m\). The primes below \(\sqrt m\) together have weight at most \(4W/\sqrt m\). At every other prime \(e_p\le D=\lceil8b\rceil\).

Project any five rows onto those other primes. If their four anchored differences had rank \(r<4\), there would be a nonzero integer relation among the projected rows, of coefficient sum zero and \(\ell^1\)-norm bounded solely in terms of \(D\). Indeed, take \(r\) independent differences, express one remaining difference by Cramer's rule on \(r\) independent columns, clear the determinant, and include the anchor. The case \(r=0\) just uses two rows. Each minor has entries at most \(D\) in absolute value. For the same relation on the full rows the residual vector \(v\) is nonzero by real affine independence. Yet \(w(v)=O_b(W/\sqrt m)\); its character \(e^{2i\langle v,\phi\rangle}\) has argument \(O_b(\Delta)\) modulo \(2\pi\), since the coefficients sum to zero. This violates the monomial gap for large \(m\). Here \(W\to\infty\), already because real affine independence requires at least \(m-1\) columns. Thus *every* tuple has a nonzero minor of absolute value at most \((2D)^4\) on the retained primes, by Hadamard. Every prime divisor of its index is at most this constant. Bounded \(m\) is absorbed into the claimed counting constant.

**PROVED — two free rows in one residue class.** Fix an odd prime \(q\). A rank loss is a nonzero \(c\in\mathbb F_q^5\) with \(\sum c_i=0\) and \(\sum c_i a_{j_i}=0\pmod q\), coordinatewise. Fix \(c\), and fix the rows in three slots, leaving two slots with nonzero coefficients free. Write the two centered integer representatives as \(h,k\), and \(L=|h|+|k|\le q-1\). All sums \(s=h a_x+k a_y\) lie in a single residue class \(r\pmod q\). The points
\[
 y_s=\exp\big(2i\langle(s-r)/q,\phi\rangle\big)
\]
have their \(q\)-th powers in one arc of width at most \(L\Delta\). The center can drift with the cluster, the residue, and the chosen slots. Within one of the \(q\) inverse arcs, distinct sums have \(v=(s-s')/q\ne0\) and \(|\sin\langle v,\phi\rangle|\le L\Delta/(2q)\). The gap gives \(w(s-s')>qW/2\). On the other hand all these sums lie in an integer box of total weighted width at most \(LW\). For \(n\ge2\) distinct sums there, the same threshold count bounds their mean distance by \(nLW/(2(n-1))\). Hence \(n\le q/(q-L)\) on each branch. A specified sum \(s\) has at most two ordered representations in the two slots: a collision is a real affine relation, so the signed measures \(h\delta_x+k\delta_y\) agree on the affinely independent rows. The slots use distinct vertices. There are \(O_q(1)\) choices for the two remaining rows. Sum over the three fixed slots and the finitely many coefficient vectors to obtain \(O_q(m^3)\).

**PROVED — parity.** A nonzero affine dependence on five marks over \(\mathbb F_2\) uses two or four marks. For two marks their difference, divided by two, has weight at most \(W/2\). For four, split them into two plus and two minus signs. Averaged over the three pairings, the absolute signed sum on any binary threshold layer is at most one: the possible types have averages \(0,1,2/3\) according as the minority has size zero, one or two. The triangle inequality over the layers gives one choice whose half-sum \(v\) has weight at most \(W/2\). In either case \(v\) is integral by the parity relation and nonzero by affine independence. The full signed sum of the row phases is within \(2\Delta\) of \(2\pi\mathbb Z\), hence \(\langle v,\phi\rangle\) is within \(\Delta/2\) of \((\pi/2)\mathbb Z\). This contradicts the gap. The fixed-prime estimate and the bounded minor now prove the normalized counting theorem.

**PROVED — effect on the proposed grid test.** For any fixed positive degree, the total number of non-unit-index tuples is eventually smaller than that degree times \(m^4\) on these growth-scale profiles. This includes profiles in the angle model after the same one-unit shortening. Large powers at small primes do not rescue the proposed selection. This does not preclude a test on a specially chosen smaller grid with a different zero count.

**PROVED — the record's skew double core.** For this particular near-extremal fake, *every proper restriction* of the marks generates the entire anchored lattice. Let \(M=q+1\), and label the marks by a border \(*\) and the elements \(a\in\mathbb F_q\). In the two Paley columns indexed by \(a\), the binary columns have difference \(e_a-e_*\): the sign flip changes just the border and \(a\), with opposite signs. Project onto a proper subset. If the border is absent, these differences supply all the retained coordinate vectors. If the border is present, use an omitted \(a\) to get \(-e_*\), and then get every retained \(e_a\). Taking differences from any anchor maps this lattice onto the full anchored lattice. At all \(M\) marks, the index is instead \(M/2\). Indeed the differences generate the sum-zero lattice in \(\mathbb Z^M\), each binary column has sum \(M/2\), and anchoring means quotienting by \(\mathbb Z\mathbf1\), of coordinate sum \(M\). The quotient therefore has order \(\gcd(M,M/2)\). This example shows exactly how a growing index on the whole profile can disappear on every five-mark restriction. [Double-core columns](/home/user/Jarnik/docs/claude_review/round8/lw.md:360)

**COMPUTED.** Integral Hermite normal form gives whole-profile indices four and six at \(M=8,12\), respectively. All 56 and 792 five-subsets have index one. The Walsh cuts used in II.3 are a separate equal-weight assignment test. They were never required to satisfy the simultaneous phase gaps. [Paley index certificate](ii_paley.py)

**PROVED using Baker–Wüstholz — removal of the growth hypothesis for genuine points.** There is an **absolute** bound on all five-point indices occurring on any normalized genuine arc above. This uses the two-logarithm estimate in §2.2 of the parallel report. Here is the deduction including the change of height. Write \(H\mid N\) for the gcd-primitive norm of the five points, \(W_T=\log H\), and \(\delta_T\) for their chordal diameter after division by the radius. Then \(0<\delta_T\le\Delta\le\tfrac12 H^{-1/4}\). For a rank-four tuple of index \(I\) the logarithm estimate gives
\[
 \frac{-\log\delta_T}{W_T}
 \ll I^{-1/4}\log(2I)+\frac{\log(2+W_T)}{W_T}.
\]
The implicit constant is absolute. The left side is at least \(1/4\). For large \(W_T\) this bounds \(I\) absolutely, and for bounded \(W_T\) the elementary bound on four-column minors does the same.

For completeness, let \(b_p\in\mathbb Z^4\) be the prime columns. Minkowski gives a nonzero dual vector \(v\) with \(\|v\|_\infty\le2 I^{-1/4}\). Take its order \(n\le I\) modulo \(\mathbb Z^4\) and set \(c=nv\). The number \(u=\prod_p(\pi_p/\bar\pi_p)^{v\cdot b_p}\) is rational Gaussian, has absolute value one, and is nontorsion by full rank and unique factorization. Its absolute logarithmic height satisfies \(2h(u)\le4\|v\|_\infty W_T\), since \(2h(u)\) is the sum of absolute signed valuations times \(\log p\). Telescoping the four anchored ratios, with Gaussian units killed by a fourth power, gives \(|u^{4n}-1|\le4\|c\|_1\delta_T\). The Baker–Wüstholz bound \(-\log|u^{4n}-1|\ll(1+h(u))\log(2n)\) uses only \(\log u\) and \(\log(-1)\), with the winding integer \(O(n)\). Finally \(I\ll(1+W_T)^4\) by any nonzero minor and \(\|c\|_1\le8I\). Substitution proves the displayed estimate. Thus on genuine clusters the fixed-prime argument gives \(O(m^3)\) exceptional five-tuples **at any growth rate**, with an absolute constant after shortening. This use of the parallel result needs its logarithm estimate alone. [Parallel report, §2.2](../lemma/answer.md) · [Baker–Wüstholz](https://doi.org/10.1515/crll.1993.442.19)

**PROVED — a sufficient five-point hypothesis at index one alone.** For the uniform arc bound it suffices to establish \((L_5^\circ)\) **on index one**. Concretely, suppose there exist \(\eta,\gamma>0\) and a fixed proper closed subset of the four-ratio torus such that, outside it, on all sufficiently large one-unit primitive index-one tuples with fifteen cut weights in \([(1-\eta)w,(1+\eta)w]\), the sum of the ten logarithmic scalar cofactors from Part I is at least \(\gamma w\). If clusters of unbounded size existed, shorten and normalize them. The Gram and threshold calculation of Part I §5 gives \(\gg m^5\) ordered tuples with relative cut errors \(O(m^{-1/2})\), \(w\asymp W\), and sum of logarithmic scalar cofactors \(O(W/m)\). Removing \(O(m^3)\) tuples leaves \(\gg m^5\) at index one. A fixed hypersurface in the ratios contains only \(O(m^4)\) grid tuples by homogenization and Schwartz–Zippel. An outside tuple contradicts the hypothesis. For the growth saving alone, it likewise suffices to impose (5) of Part I just on its index-one set, by the same removal in that proof. These are sufficient hypotheses for clusters. Individual primitive five-tuples at indices three and five still matter to the unrestricted five-point theorem.

**PROVED — calibration of the cover exponent.** In the notation of the parallel report, \(A_5(1/31)\) applies when \(\delta_T\le H^{-33/124}\). The cluster bound \(\delta_T\le\tfrac12 N^{-1/4}\) suffices if \(H\le N^{31/33}\). The typical primitive log norm \((15/16+o(1))W\) has enough slack since \(15/16<31/33\). A bare assertion that the primitive norm is close to \(N\) does not give this inequality. The index bound above explains why the fixed closed set supplied by the high-index cover cannot catch a positive fraction even of the tuples whose heights qualify. [Parallel report, §§1–3](../lemma/answer.md)

### A positive proportion of the points can be kept with index one on every five-tuple

**PROVED — stronger selection theorem.** Suppose every prime divisor of every five-point index is at most \(Q\). This holds with an absolute \(Q\) for genuine clusters above and with \(Q=Q_b\) for their growth-scale angle models. After deleting at most \(K_Q\) points, *every* sufficiently balanced five-tuple has index one. Moreover one can keep at least \(m/145-K_Q\) points for which **every five-tuple**, with any cut weights, has index one. Balance here means relative balance of the fifteen nonconstant cut weights, and the tolerance \(1/200\) suffices. The conclusion controls an entire subcluster chosen once.

*Proof.* Choose the threshold layer \(\ell=(p,h)\) with probability \(\log p/W\), put \(B_{i\ell}=\mathbf1_{a_{ip}\ge h}\), \(f_{i\ell}=2B_{i\ell}-1\), and \(G_{ij}=\mathbb E_\ell f_{i\ell}f_{j\ell}=1-2w(a_i-a_j)/W\). Call a tuple eligible for this proof if all its off-diagonal \(|G_{ij}|\le1/12\). Pair separation gives \(G_{ij}<0\). In a relatively balanced tuple every pair crosses eight of the fifteen nonconstant cuts, so \(w(a_i-a_j)/W\le8(1+\eta)/(15(1-\eta))\); at \(\eta=1/200\) this is smaller than \(13/24\). Thus every such balanced tuple is eligible.

Fix an odd \(q\le Q\) and one nontrivial affine relation vector in \(\mathbb F_q^5\). Multiply by a nonzero scalar and choose centered integral representatives \(c_i\) so that \(S_2=\sum c_i^2\le5q(q+1)/12\). This follows by averaging over scalars, since the centered square of a uniform nonzero residue has mean \(q(q+1)/12\). Among eligible tuples \(T\) with this relation, the sums \(b_T=\sum_i c_i a_{T_i}\) are multiples of \(q\), coordinatewise. Write \(L=\sum_i|c_i|\). As in the two-slot proof, the possible phases of \(b_T/q\) lie in \(q\) inverse arcs of width \(L\Delta/q\) in the doubled-angle variable. Within each branch two distinct sums have separation at least \(qW/2-O(q)\), by the monomial gap and \(L/q\le5/2\). If a branch has \(k\) distinct sums, sample two independent uniform representatives, one representative tuple per sum. Their mean separation divided by \(W\) is at least \((1-1/k)(q/2-O(q/W))\).

There is also an upper bound on this mean. For a fixed layer put \(Z_\ell=\sum_i c_i B_{T_i,\ell}\) with randomness in the chosen representative. By the triangle inequality over the layers, the normalized mean separation is at most the layer average of \(\mathbb E|Z_\ell-Z'_\ell|\). Flip all bits on a layer together with probability one-half, and pool the resulting laws over layers. The Gini mean difference of this pooled law bounds the preceding average from above, because it is concave in the law: \(\mathbb E|Z-Z'|=2\int F(t)(1-F(t))\,dt\). For any real law it satisfies \(\mathbb E|Z-Z'|\le(2/\sqrt3)\sqrt{\operatorname{Var}(Z)}\). Indeed if \(P\) is its quantile function the left side is \(2\int_0^1(2u-1)(P(u)-\mathbb EZ)\,du\), by integration on the two half-squares, and Cauchy–Schwarz gives the bound.

For each representative tuple the orientation flip gives the same mean \(\sum c_i/2\). Its variance over the layers and the flip is \(c^tG_Tc/4\le (1+4/12)S_2/4\), since \(\sum_{i\ne j}|c_ic_j|\le4S_2\). Pooling over representatives keeps this bound. The squared Gini mean difference is therefore at most
\[
 \frac49 S_2\le\frac5{27}q(q+1)\le\frac{20}{81}q^2<\frac{q^2}{4}.
\]
Comparison with the lower bound forces \(k=O(1)\) for large \(W\). For each sum the active row set is unique: subtract two exact expressions for it, whose coefficient totals agree, and use real affine independence. Deleting the active rows of representatives for these boundedly many sums in all branches, for all \(q,c\), removes every eligible tuple of index larger than one. Bounded \(W\) is absorbed in \(K_Q\), and parity was already excluded.

Finally \(G\) is positive semidefinite, and \(I-G\) is symmetric entrywise nonnegative. Its largest eigenvalue is at most one, and bounds the absolute value of every eigenvalue since \(|x^t(I-G)x|\le |x|^t(I-G)|x|\). Thus \(\sum_{j\ne i}G_{ij}^2\le1\) at each mark. The graph of pairs with \(|G_{ij}|>1/12\) has maximum degree at most 144 and has an independent set of size at least \(m/145\) by greedy selection. All its five-tuples are eligible. Removing the same exceptional marks proves the second assertion.

**COMPUTED — a check on the margin.** For \(q=3,5,7\) one can test the Gini estimate directly on all projective affine-relation coefficient vectors. For a mixture of the fair five-bit law and the uniform law on nonconstant five-bit patterns, maximize over the mixture, minimize over nonzero scalar multiples, then maximize over relation classes. The exact values are \(315/256,239/128,633/256\), respectively, all below \(q/2\). The calculation enumerates four, eight and thirteen classes modulo coordinate permutations and scaling and tests each rational quadratic on its full parameter interval. [Finite Gini certificate](ii_index_gini.py)

## II.5. What the fake divisor count proves, and a common-parent arithmetic test

**PROVED — the fake analogue even distinguishes exponent divisors.** I agree with the objection to using the marginal angle count. In the squarefree case the signed differences are ternary. With prime powers the corresponding box is \(-e_p\le v_p\le e_p\). On the normalized one-unit arc the averaging in Part I gives \(\gg m^2\) edges with \(w(v)\le W/2+B\log W\) on a growth-failure sequence. Their \(\langle v,\phi\rangle\) are within \(\Delta/2\) of \(\pi\mathbb Z\), hence within a factor \(O(W^{B/2})\) of their Liouville lower bound. Distinct edges among these even yield distinct ordinary labels \(\prod p^{|v_p|}\) in the fake. For if \(|v_p|=|v'_p|\) at every prime and \(v\ne\pm v'\), then \((v+v')/2,(v-v')/2\) are nonzero integral vectors of total weight \(w(v)\le W\). One has weight at most \(W/2\). Its angle is within \(\Delta/2\) of \((\pi/2)\mathbb Z\), contradicting the gap. If \(v=\pm v'\), the two edges coincide by real affine independence. This uses the edge phases near \(\pi\mathbb Z\) obtained after the one-unit restriction. Thus distinct labels or a square-product-triangle count for those *angle-defined* labels do not escape F1. The extra assertion in (1) is the literal integer identity \(d-\lfloor\sqrt d\rfloor^2=t^2\) with its primitivity and size conditions.

**PROVED — a third input for compatible triangles.** Consider three genuine arithmetic divisor labels satisfying the hypotheses of the triangle theorem in Part I. There are positive integers \(a,b,c\) with labels \(ab,bc,ca\) and \(abc\mid N\). Indeed the three quotients \((d_{12}d_{31})/d_{23}\), etc., are integer squares, and the exponent of \(abc\) at each prime is the local range of the three profile entries. Write \(w,u,v\) for the small positive integer offsets of \(ab,bc,ca\) respectively, and \(ab=s^2+w^2\). Then
\[
 \boxed{\quad c=\frac{a u^2+b v^2\ \pm\ 2suv}{w^2},\qquad
             c\mid \frac{N}{ab}.\quad}                                      \tag{II.1}
\]
In particular each way of splitting one edge \(d=ab\) prescribes a positive binary quadratic form \([a,\pm2s,b]\) of discriminant \(-4w^2\). The value at the two other *integer* offsets must be \(w^2\) times a divisor of the fixed complement \(N/d\). On edges in (1) those offsets are at most \(W^{B/2}/4\).

To prove (II.1), use the triangle theorem to put three directions on the line of angles so the absolute differences are the three angles of the near-real Gaussian factors. Choose planar vectors \(X,Y,Z\) in these directions of lengths \(\sqrt a,\sqrt b,\sqrt c\). The absolute determinants of \((X,Y),(Y,Z),(Z,X)\) are respectively \(w,u,v\), by the sine formula. The signed vectors \(\det(Y,Z)X,\det(Z,X)Y,\det(X,Y)Z\) sum to zero by the two-dimensional determinant identity. Their squared lengths are \(A=a u^2,B=b v^2,C=c w^2\), and the absolute determinant of the first two is \(uvw\). Consequently
\[
 (A+B-C)^2=4AB-4U,\qquad U=(uvw)^2.
\]
Since \(AB-U=s^2u^2v^2\), this gives (II.1), including the two possible signs.

**PROVED — an exact integer parametrization.** The form condition is reversible for finding the other two representations. If \(d=ab=s^2+w^2\) is eligible and an integral \(c>0\) satisfies (II.1) with sign \(\epsilon\), then
\[
 ac-v^2=\left(\frac{a u+\epsilon s v}{w}\right)^2,\qquad
 bc-u^2=\left(\frac{s u+\epsilon b v}{w}\right)^2 .
\]
Expansion proves both identities, and a rational number whose square is integral is integral. Thus a single quadratic-form evaluation supplies the two real coordinates as well. The exact local profile condition is \(abc\mid N\) and \(\gcd(a,b,c)=1\). To see the second condition, let \(A_p,B_p,C_p\) be the valuations of \(a,b,c\). The three local distances are their three pair sums. The largest equals the sum of the other two exactly when \(\min(A_p,B_p,C_p)=0\); their local range is then \(A_p+B_p+C_p\). Such three distances can be realized on the integer interval \([0,e_p]\) by placing the middle mark between the extremes. For squarefree \(N\) the divisibility of \(abc\) already ensures this.

Consequently one can enumerate the precise triangle target in Part I by ranging over eligible \(d\), ordered factorizations \(ab=d\), positive \(u,v\le W^{B/2}/4\), and the two signs in (II.1). Require integral \(c\mid N/d\), \(\gcd(a,b,c)=1\), and impose the primitive and size cutoffs (1) on the resulting two representations. If desired require the three labels to be distinct. A triple of distinct labels gives at most six ordered records: ordering the three labels recovers \(a,b,c\) uniquely by square roots of quotients, its offsets by floors, and its sign by (II.1). The preceding triangle theorem proves exhaustiveness. A little-oh bound of order \((W/\log W)^3\) for this common-parent count would suffice for the original growth saving.

**PROVED — what an integral descent does and where it loses the parent.** Hold \(U\) fixed and reorder \(A\le B\le C\). The other quadratic root is \(C'=2(A+B)-C=((A-B)^2+4U)/C\), a positive integer. If \(C>A+B\), replacing \(C\) by \(C'\) strictly decreases the sum of the three entries. When the descent stops, \(0\le A+B-C\le A\), so \(4AB-4U\le A^2\le AB\), whence \(3AB\le4U\) and \(C\le A+B\le8U/3\). This gives finitely many seeds for fixed offsets. The obstruction to counting from them is exact: the Vieta move preserves this equation while losing the fixed-parent divisibility.

**COMPUTED.** At the triangle \(N=131490079981\) in Part I take \((a,b,c)=(1417,20381,4553)\), \((u,v,w)=(2,1,1)\). The squared lengths are \((5668,20381,4553)\) and \(U=4\). Reflecting the largest entry gives 61 with \(N\bmod61=40\). Here the reflected entry had offset multiplier \(v^2=1\), so it genuinely fails the complementary-parent condition. [Exact Heron and descent check](ii_triangle.py)

**HEURISTIC — what I would count.** I see no pointwise estimate for \(D_B(N)\) from the pair gcd inequality or an average over thin rectangles. Formula (II.1) instead suggests counting completions of one eligible \(d\) against its own complement, summed over \(d\), its factorizations and its two small offsets. It uses simultaneous integral coordinates and retains the prime information that a descent over all parents discards. I have no bound of the required size for that sum. I would test this before trying to prove the individual-divisor hypothesis: the seven-divisor example already separates the two tasks. A separate obstruction to a growing fully split template is that even one degree-\(k\) direction in degree \(2k\) supplies two disjoint conjugate \(k\)-element sets in \(\mathbb Q(i)\) with equal power sums through \(k-1\), an ideal Prouhet–Tarry–Escott solution. This follows directly from the moment equations. Keeping the specialization heights under control would be a further requirement.

## Twelve-line summary

1. **PROVED.** Fakes can match quadratically many distinct integer divisor labels with near-Liouville angles; integrality of the actual offsets remains essential.
2. **PROVED with COMPUTED certificates.** The degree-eight weight-four maximum is 12 up to sign, attained only by the cube code.
3. **COMPUTED.** Its 28 edges consist of four weight-eight edges and two copies of each of twelve weight-four directions.
4. **COMPUTED.** A rational configuration has nine primitive directions, a largest common-moment fibre of five, and a largest minimum-edge clique of three.
5. **PROVED.** Seven of its labels pass the integer divisor test simultaneously on infinitely many parents, with no four-clique on those seven ratios.
6. **EVIDENCE.** This gives a fixed lower bound; a growing counterexample to the divisor hypothesis and a growth saving for arcs remain open.
7. **PROVED.** A cyclic six-mark determinant has automatic nonvanishing, no pair factors and forced criticality constant zero on three-three cuts.
8. **PROVED.** An unspecified nonzero Gram minor cannot force supercritical divisibility; on independent fair cuts every total-degree minor has mean constant at most one-half.
9. **PROVED.** At most \(O_b(m^3)\) five-tuples have higher index on a normalized growth-scale profile; every proper restriction of the Paley double core has index one.
10. **PROVED using Baker–Wüstholz.** Every genuine shortened one-unit cluster has an all-index-one five-tuple subcluster of size at least \(m/145-O(1)\).
11. **PROVED.** The index-one case of \((L_5^\circ)\) suffices for uniformity; the index-one version of (5) suffices for a growth saving.
12. **HEURISTIC.** I would referee the branch-packing proof and test (II.1) for many small quadratic-form values dividing one fixed complementary parent.
