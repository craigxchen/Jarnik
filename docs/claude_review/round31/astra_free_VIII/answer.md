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

# Part III — Recovering the common cuts from small integer Plücker coordinates

**PROVED — main result.** For any sufficiently large normalized profile, the parent and the exact separation divisors determine the directions of *all* realizable integer stars up to one rational projectivity. This survives a missing tenth of the constraints at each mark (III.7). A two-against-three cut on five marks also forces independently chosen small quartet solutions to give compatible cross ratios (III.8). The restriction lies in finding those integer solutions and factoring them consistently. A density theorem gives a uniform bound for almost all parents in an explicit prime-product family (III.9). At five marks I construct a full primitive integer model with all prescribed minor divisibilities at the ten-pair exponent 10/3 (III.11). The worst-parent growth bound remains open.

**PROVED.** My choice is (c), using four columns of the integer star matrix at a time. There is a constraint here that the usual angle fakes fail. At one fixed parent, assign the same primes to the same set of balanced binary cuts in different orders. In the Paley double core with \(r=2m\) primes and \(m\) ordinary marks, at most a fraction \(1/\lfloor r/64\rfloor!\) of the assignments can satisfy the small-integer Plücker tests below. One can even choose the assignment so that **no four marks pass**. The excluded matrices may have arbitrary projective lines modulo the primes. No Gaussian factorization or fixed rational plane is imposed on them. On the narrow prime intervals used below at polynomial height, *every* assignment still admits an angle fake. A second use of the shared layers proves that a common rank-\(r_0\) unit group can account for only \(O(r_0+1)\) extensions of three anchors with polynomially small incident cofactors, even for rational cross ratios.

**EVIDENCE — verdict on growth.** The saving is in the number of prime assignments at a specified parent. A growth counterexample needs only one of those assignments. The natural way I tried to use the other \(m-4\) marks to amplify a single quartet congruence recovers exactly the original quartet modulus. I do not get a saving in the bound for \(m\), or in the constant \(2/3\). I also see no reduction from the growth problem to a necessary Vojta-strength assertion. The residual question is whether one exceptional assignment can support all the quartet identities simultaneously.

## III.1. A small-integer test involving three disjoint parts of one parent

**PROVED — notation and exact test.** Use the shortened one-unit convention of Part I, angular width \(\Delta\le \tfrac12 N^{-1/4}\), and exponents \(0\le a_{ip}\le e_p\). Orient the split generators once. Put
\[
 n_i=\prod_p p^{a_{ip}},\quad P_i=\prod_p\pi_p^{a_{ip}},\quad
 g_{ij}=\prod_p p^{\min(a_{ip},a_{jp})},\qquad
 d_{ij}=n_i n_j/g_{ij}^2 .
\]
Here \(P_i\) is a primitive integer plane vector. Its norm squared is \(n_i\). The doubled arguments of the stars agree with the point arguments up to one common phase. Since \(g_{ij}\) divides \(\bar P_iP_j\), the signed number \(t_{ij}=[P_i,P_j]/g_{ij}\) is an integer with \(0<|t_{ij}|\le\sqrt{d_{ij}}\,\Delta/2\). This formulation does not require the \(n_i\) to be mutually comparable.

For four distinct labels use the three matchings \(12|34,13|24,14|23\), in that order. Take the three products \(g_{12}g_{34},g_{13}g_{24},g_{14}g_{23}\) and divide them by their common gcd. Write \(A,B,C\) for the results. They are pairwise coprime. At a prime with four entries sorted as \(h_1\le h_2\le h_3\le h_4\), the pairing of the two lowest and two highest gives exponent \(h_1+h_3\) before division; the other two give \(h_1+h_2\). Thus the only surviving exponent is \(h_3-h_2\), on the matching induced by its two-two threshold cuts. In particular \(ABC\mid N\), and the three numbers depend only on the *unoriented* cuts. The two-dimensional determinant identity becomes
\[
 A x - B y+C z=0,\qquad
 (x,y,z)=(t_{12}t_{34},t_{13}t_{24},t_{14}t_{23}).                     \tag{III.1}
\]
This condition is necessary for **any** integer star matrix with the indicated factorizations of its minors. In the tests below one may choose three independent signed coefficients afresh on each quartet, subject to a size bound and nonvanishing. No consistency of coefficients between quartets is assumed. The underlying three-term identity is the ordinary Plücker relation; the result below concerns recovery and counting across different assignments of its shared factors.

**PROVED — a test using only integer arithmetic.** The congruence \(A x+C z=0\pmod B\) puts \((x,z)\) in the planar lattice spanned by \((B,0)\) and \((-CA^{-1}\bmod B,1)\) for \(B>1\). If its squared shortest-vector length exceeds \(2U^2\), there is no solution of (III.1) with \(|x|,|y|,|z|\le U\). Try each matching as modulus. For measured separation divisors, the same \(A,B,C\) can be computed without factoring anything. Let \(D_1=d_{12}d_{34},D_2=d_{13}d_{24},D_3=d_{14}d_{23}\) and \(L=\operatorname{lcm}(D_1,D_2,D_3)\). Then \((A,B,C)=(\sqrt{L/D_1},\sqrt{L/D_2},\sqrt{L/D_3})\), since the three \(D_j\) are one common integer divided by \(A^2,B^2,C^2\) and the three factors are coprime.

**PROVED — recovery from a common modulus.** Suppose two possible assignments give \((A,B,C)\) and \((A',B',C')\), with signed nonzero triples of coefficients bounded by \(U\) passing (III.1). Put
\[
 Q=\gcd(B,B'),\qquad u/v=(AC')/(A'C)\quad\text{in lowest positive terms}.
\]
Then \(Q\mid u x z'-v x' z\). Indeed reduce both equations modulo \(Q\), cross-multiply, and cancel the common factors in \(AC',A'C\), which are coprime to \(Q\). If
\[
                  2U^2\max(u,v)<Q,                                      \tag{III.2}
\]
the congruence is an equality of integers. Consequently every prime at which \(A/C\) and \(A'/C'\) have different valuations divides \(xx'zz'\). In particular if all four coefficients are coprime to the parent, \(A/C=A'/C'\), and coprimality gives \(A=A',C=C'\). Applying this reasoning with another matching as modulus also recovers \(B\). This argument works for prime powers with their actual valuations.

## III.2. Applying recovery on many quartets at one parent

**PROVED — the fixed-parent theorem.** Let \(m\equiv3\pmod4\) be a sufficiently large prime. Use the ordinary \(m\) rows of the round-8 Paley skew double core of order \(m+1\), deleting its border row. There are \(r=2m\) cut columns, distinguishable even up to taking complements. Fix \(r\) distinct split primes in \([P,P\exp(1/(10r))]\), and let \(N\) be their product. Assume \(P\to\infty\). Assign the primes bijectively to the columns. Call an assignment admissible if every four-row restriction admits some solution to (III.1) in signed nonzero integers of absolute value at most \(U=\lfloor P/4\rfloor\). Then:

1. For a uniform random assignment, the probability that **even one** quartet passes is \(O(m^4 2^{-c m})\) for an absolute \(c>0\).
2. Two admissible assignments differ on more than \(\lfloor r/64\rfloor\) primes. In particular at most \(r!/\lfloor r/64\rfloor!\) assignments are admissible.
3. Fix which unordered pair of primes goes to each original and modified pair of Paley columns. Of the \(2^m\) ways of ordering these pairs, at most **one** can be admissible.

Actual normalized clusters on this profile must be admissible. Each pair is separated by \(m+1\) of the \(2m\) cuts, whence \(d_{ij}/\sqrt N\le P\exp((m+1)/(10r))\). Thus \(t_{ij}^2\le d_{ij}/(16\sqrt N)<P/8\) and \(|t_{ij}t_{kl}|\le U\) for large \(P\). More generally the theorem excludes the relaxed integer-matrix model whenever it supplies the stated small cofactors. It makes no assumption about the residues defining its projective lines.

**PROVED — the needed cut counts.** In the ordinary rows, the two columns indexed by \(a\in\mathbb F_m\) have signs \(\chi(a-i)\) at \(i\ne a\), and opposite signs \(\pm1\) at \(i=a\); here \(\chi\) is the quadratic character. A row pair is separated at both exceptional columns' indices once each, and elsewhere at both copies whenever the two characters disagree. The identity \(\sum_a\chi((a-i)(a-j))=-1\) gives the claimed total \(m+1\). This identity follows by counting \(v^2=u(u-1)\) using \((2u-1-2v)(2u-1+2v)=1\). Distinct base columns differing in their field index disagree on \(m/2+O(1)\) ordinary rows, and the two columns of one base index differ at just that index. Each column contains \(m/2+O(1)\) bits of each kind.

For every quartet of distinct ordinary marks, each of its three two-two matching classes has \(r/8+O(\sqrt m)\) columns. Here is the standard Hasse-bound calculation. Away from its four marked indices the indicator of an unoriented sign pattern \(\epsilon\) is
\[
 \tfrac18\left(1+\sum_{i<j}\epsilon_i\epsilon_j\chi((a-i)(a-j))
             +\Big(\prod_i\epsilon_i\Big)\chi\Big(\prod_i(a-i)\Big)\right).
\]
The quartic sum has absolute value at most \(2\sqrt m+1\). In fact the smooth completion of \(v^2=\prod_i(a-i)\) is a genus-one curve with two rational points at infinity (four simple branch points), so Hasse's theorem gives this bound after counting its affine points. The four exceptional indices change the count by \(O(1)\), and there are two copies away from them. Consequently each matching has at least \(r/16\) primes for sufficiently large \(m\), uniformly over the quartet. The definition of the double core and its angle realization are in the record. [Round 8, §§5–6](/home/user/Jarnik/docs/claude_review/round8/lw.md)

**PROVED — one quartet already removes most assignments.** Fix its marks. Choose \(s=\lfloor r/64\rfloor\) disjoint pairs of column slots, taking one slot of each pair from class \(A\) and one from class \(C\). Partition the assignments into orbits of the \(2^s\) independent interchanges within these pairs. Within an orbit \(B=B'\ge P^{r/16}\); the reduced numerator and denominator of \(AC'/(A'C)\) are at most \((2P)^{2s}\). Inequality (III.2) holds for large \(m\). All coefficients have magnitude less than the least parent prime. Any two passing assignments in the same orbit must therefore have the same \(A\), which determines all \(s\) choices by unique factorization. The passing fraction is at most \(2^{-s}\). Sum over the quartets for assertion 1.

**PROVED — all quartets encode the marked prime cuts.** Suppose two assignments differ on \(h\le r/64\) primes. On each quartet and with any matching as modulus, \(Q\ge P^{r/16-h}\). Only the moved primes can occur in \(u,v\), to exponents at most two, so \(\max(u,v)\le(2P)^{2h}\). Since \(r/16-3h\ge r/64\) and \(U\le P\), (III.2) again holds. If both assignments pass on that quartet, they give the same three matching products. If a prime moved from cut \(S\) to cut \(S'\), orient the cuts to agree on at least half the marks and choose a mark where they disagree. On their agreement set they have \(m/4-O(1)\) or more common zeros and as many common ones, since each cut is itself balanced up to \(O(1)\) marks. Choose three agreeing marks so that this quartet is two-two for \(S\). It is one-three for \(S'\), contradicting equality of matching products. Thus both full assignments cannot pass. Fixing the locations of \(r-\lfloor r/64\rfloor\) primes now leaves at most one admissible completion, proving assertion 2. In fact the distinguishing quartet can be chosen in \(\gg m^3\) ways through the indicated mark, so the separation conclusion tolerates \(o(m^3)\) local exceptions through each mark in each assignment.

For assertion 3, any quartet distinguishes the two slots at just the four base indices of its marks. All other allowed pair interchanges move primes between slots of the same quartet type. Hence on that quartet at most eight primes can change type, even if every pair has been interchanged globally. The preceding proof uses only this local number in (III.2). Agreement on the matching products of every quartet again recovers every cut, proving uniqueness.

**PROVED — comparison with angle fakes.** Such fixed parents exist with \(W=O(m\log m)\). For example take a sufficiently large fixed \(D\), use the prime number theorem in the progression \(1\bmod4\) on \([m^D,2m^D]\), and partition its log interval into \(O(m)\) bins of width at most \(1/(10r)\). One bin contains the required primes. The proof of round 8, Theorem 6.1, gives a gap-admissible fake for *each assignment*. Its margins use the incidence matrix, the minimum prime weight and the log interval width, allowing its phases to vary with the assignment. Apply it before deleting the single exceptional winding row to obtain the one-unit subarc here. Thus an arbitrary angle assignment at this parent can pass every monomial gap while exact integer divisibility of the star minors fails on every quartet. In this last statement choose one of the profile assignments supplied by assertion 1 and then a good angle realization for it.

**HEURISTIC — what is still missing.** This is a constraint on the simultaneous assignment of primes, using the multiplicity of marks to distinguish even nearly identical cuts. Most of the permutation freedom at a fixed parent is irrelevant to the original conjecture: a large actual cluster needs one labeling, and does not automatically furnish another labeling by exchanging primes. No density restriction in this theorem eliminates that last labeling. In particular the result is compatible with the existence of an integral fake on an exceptional parent or assignment.

### Exact tests, including every prime assignment at a single parent

**COMPUTED.** I assigned successive split primes above \(10^5\) to the ordinary Paley rows using a fixed permutation (seed 2026). I allowed the same upper bound \(U=\lceil\max_{ij}d_{ij}/(16\lfloor\sqrt N\rfloor)\rceil\) on all three coefficients and tested the shortest vector for each of the three possible moduli. The table counts all quartets, with no sampling. [Exact test](iii_quartet.py)

| Ordinary marks | Parent digits | \(U\) | Quartets tested | Quartets surviving |
|---:|---:|---:|---:|---:|
| 7 | 71 | 6,285 | 35 | 0 |
| 11 | 111 | 6,327 | 330 | 0 |
| 19 | 191 | 6,399 | 3,876 | 0 |

**COMPUTED — every assignment of the first parent.** Take exactly the fourteen split primes from 100049 through 100333. Use \(U=6294\), large enough for *every* assignment since the largest possible separation divisor is bounded by the product of the eight largest primes. On the seven ordinary marks, 21 quartets have two primes per matching class and 14 quartets have three. For two per class, among all 45,045 unordered triples of disjoint two-prime products, 156 admit a nonzero-coefficient relation. For three per class, among all **560,560** unordered triples of disjoint three-prime products, **none** admit a relation. Thus **no permutation of this parent** passes the full seven-mark test. As a positive control, the genuine balanced four-cycle below passes (III.1), with the signed products computed from its six actual Gaussian ratios. [All-assignment certificate](iii_quartet_exhaust.py) · [Positive control](iii_quartet.py)

**PROVED — certificate method.** To enumerate coefficients exactly, lift \((B,0),(r,1)\), \(r=-CA^{-1}\bmod B\), to the kernel basis \((B,-A,0),(r,(-Ar-C)/B,1)\). Two-dimensional Gauss reduction works in its Euclidean span. In a reduced basis \(v,w\), \(\|v\|\le\|w\|\) and \(2|v\cdot w|\le\|v\|^2\), so \(v\) is shortest: for a combination with second coefficient of magnitude at least two use its component perpendicular to \(v\), and for magnitude one use nearest-integer reduction. Any vector in the coefficient cube has squared length at most \(3U^2\). The diagonal entries of the inverse Gram matrix, multiplied by \(3U^2\), bound the squares of its two basis coefficients by Cauchy–Schwarz. The script tests this entire rectangle, dropping nonprimitive vectors since dividing by their coefficient gcd retains all required nonzero entries.

## III.3. Why the other marks do not automatically enlarge this modulus

**PROVED — eliminating a common mark.** Fix four target marks \(i,j,x,y\) and use another mark \(k\) to compare the two quartets \(ijkx,ijky\). Write \(B(ij|kx)\) for the matching product attached to the separation \(ij|kx\). On \(Q_k=\gcd(B(ij|kx),B(ij|ky))\), the two Plücker congruences allow division by their crossed-edge cofactors for genuine data. At a prime in the corresponding two-two cut, every crossed edge has nonzero valuation difference; its primitive cofactor is prime to that prime. The quotients of the crossed matching products are likewise local units. Cancelling the common terms involving \(k\) gives exactly the direct \(ijxy\) congruence, on this smaller modulus. Indeed the raw cross quotient \(g_{ik}g_{jx}/(g_{jk}g_{ix})\), divided by its version with \(y\) in place of \(x\), equals \(g_{jx}g_{iy}/(g_{ix}g_{jy})\), independent of \(k\). The signed cofactor quotients telescope in the same way.

In fact
\[
                 \operatorname{lcm}_{k}Q_k\mid B(ij|xy).                 \tag{III.3}
\]
For squarefree parents, a prime in \(Q_k\) puts \(i,j\) on one side and \(k,x,y\) on the other. For prime powers, fix \(p\). A positive exponent for \(B(ij|kx)\) counts integer thresholds separating \(i,j\) from \(k,x\). For fixed \(i,j,k\), the possible thresholds for \(x\) and \(y\) are nested intervals beginning at the same end of the gap between \(k\) and the interval spanned by \(i,j\). Their common length is the minimum of their lengths. These thresholds all separate \(i,j\) from \(x,y\). Taking a maximum over \(k\) proves (III.3). Equality holds, for example, at a squarefree prime if another mark on the \(x,y\) side is available. Many witnesses for this congruence can therefore do no more than restore its direct modulus.

**COMPUTED and PROVED — the height loss for a single substitute mark.** For squarefree five-mark data with fair cut weights, classify respectively the restrictions \(1234,1235\) by \(A,B,C,0\), using consistent names for the three matchings of the shared anchors and writing zero for a different cut size. Among the 32 oriented patterns, \((0,0)\) occurs fourteen times; each \((J,0),(0,J),(J,J)\) occurs twice for \(J=A,B,C\). Thus a common \(B\) has log size \(W/16\), whereas *each* of the reduced numerator and denominator of \(AC'/(A'C)\) has log size \(W/8\). Here \(W\) is the ambient log weight, including the constant cut. The calculation is direct enumeration of five bits and is stable up to \(o(W)\) under asymptotic balance. It explains why the same small-height comparison that recovers a few exchanged prime labels does not recover a new mark in a typical five-tuple. Using all substitute marks reaches (III.3). It gives no excess beyond the direct four-point constraint. [Five-bit check](iii_quartet.py)

## III.4. The rank actually needed for a unit-equation count with fixed anchors

**PROVED — linear rank cost for one common system.** Fix three marks \(i,j,k\) of a genuine normalized cluster. For another mark \(x\), put
\[
 F_x=\frac{[P_x,P_i][P_j,P_k]}{[P_x,P_k][P_j,P_i]}.
\]
Then \(1-F_x=[P_x,P_j][P_i,P_k]/([P_x,P_k][P_i,P_j])\). Both are nonzero rational numbers. Suppose \(s\) distinct extensions have \(|t_{xi}|,|t_{xj}|,|t_{xk}|\le m^B\), and all \((F_x,1-F_x)\) belong to one coset of a rank-\(r_0\) multiplicative group. For fixed \(B\), uniformly in the parent and in the three anchors,
\[
                         s\le (16/3+o_B(1))r_0+O(1),\qquad m\to\infty.   \tag{III.4}
\]
No polynomial bound on the three anchor-to-anchor cofactors is required. This is an actual shared-atom obstruction to rank compression for the rational equation \(F+(1-F)=1\). It allows a different group and different coefficients absorbed as a fixed translate for *each* triple of anchors.

*Proof.* Let \(w(v)=\sum_p|v_p|\log p\) and \(W=\log N\). Recall from II.4 that \(w(a_x-a_y)>W/2\), that \(e_p\log p\le4W/m\), and hence that primes \(p<\sqrt m\) have total weight at most \(4W/\sqrt m\). The threshold sign matrix with rows \(f_x(p,h)=2\mathbf1_{a_{xp}\ge h}-1\), using mass \(\log p/W\) per layer, has operator norm at most \(\sqrt2\): its row Gram is positive semidefinite with unit diagonal and negative off-diagonal, and the symmetric nonnegative matrix given by subtracting this Gram from the identity has spectral norm at most one, as in II.4.

At each remaining prime let \(u_p,v_p\) be the smallest and largest exponents of the anchors and \(L_p=v_p-u_p\). Drop primes with \(L_p=0\). The mass \(D=\sum L_p\log p/W\) of the retained layers satisfies \(D>3/4-4/\sqrt m\), since the sum of the three anchor-pair distances equals twice the weight of their ranges. Let \(z_{xp}=\min(a_{xp},v_p)-\min(a_{xp},u_p)\) and consider the vectors
\[
 c_x(p)=(2z_{xp}-L_p)\sqrt{\frac{\log p}{W L_p}} .
\]
These are obtained by orthogonal projection of the threshold signs onto the span of the normalized indicator of each retained block \(u_p<h\le v_p\), in orthonormal block coordinates. Thus the matrix of the \(c_x\) has operator norm at most \(\sqrt2\). For two extensions, nestedness within each prime means that \(\sum_p|z_{xp}-z_{yp}|\log p\) is exactly their distance on the retained layers. Choose two independently and uniformly among the \(s\) extensions. Their total expected distance divided by \(W\) is at least \((1-1/s)/2\); their expected distance on discarded layers is at most \((1-D)/2\), since two independent bits differ with probability at most one-half. Cauchy–Schwarz across the retained primes and Jensen across pairs therefore give
\[
 \mathbb E_{x,y}\|c_x-c_y\|_2^2
 \ \ge\ (D-1/s)_+^2/D\ \ge\ D-2/s.
\]

The two observed valuations, after subtracting constants independent of \(x\), are \(\min(a_x,a_i)-\min(a_x,a_k)+v_p(t_{xi}/t_{xk})\) and the corresponding expression with \(j\) in place of \(i\), prime by prime. Include the zero component for anchor \(k\). Select the components of the largest and smallest anchor at this prime and subtract them. This recovers \(z_{xp}\) up to an additive error \(E_{xp}\) equal to the difference of the two corresponding cofactor valuations. All choices are fixed with the anchors, including a tie convention. Across the retained primes, \(|E_{xp}|\le2B\) since \(p\ge\sqrt m\), and \(\sum_p |E_{xp}|\log p\le3B\log m\). The observed vectors after the same affine normalization as \(c_x\) therefore differ from them by squared Euclidean norm at most
\[
                \frac4W\sum_p E_{xp}^2\log p/L_p
                       \le 24 B^2\log m/W=O_B(1/m).
\]
For the last estimate the \(m\) affinely independent exponent rows require at least \(m-1\) distinct primes, so \(W\ge\log((m-1)!)\).

Valuations take a rank-\(r_0\) coset to an affine space of dimension at most \(r_0\), and the component selections and normalizations just used are affine maps. The \(s\) vectors \(c_x\) are therefore within squared distance \(O_B(1/m)\) of one such affine space. After subtracting the row mean, their total squared norm is at least \(sD/2-1\) by the pair-distance identity. Their total squared distance from the corresponding linear space is at most \(O_B(s/m)\), because centering also centers the errors. At most \(2r_0\) of the total squared norm can lie in that linear space, by the operator norm bound. Rearrangement proves (III.4).

**PROVED — one rational coordinate already costs linear rank.** If only the marks \(i,k\) have small incident cofactors, keep the layers between those two anchors. Now \(D>1/2-o(1)\); the first observed coordinate alone supplies the clipped sum, up to a known prime-dependent sign and the bounded cofactor error. The identical proof gives \(s\le(8+o_B(1))r_0+O(1)\) if the \(F_x\) alone lie in a rank-\(r_0\) coset. Both versions remain valid after multiplying the observed coordinates by individually varying rationals of numerator and denominator at most \(m^{B'}\), for fixed \(B'\). The added valuations have bounded size above \(\sqrt m\) and total absolute log weight \(O(B'\log m)\), exactly the errors already permitted in the proof.

**PROVED — the common-parent triangle equation has an exact rank bound with no asymptotics.** Fix an edge \(ij\). In (II.1) write \(d_{ij}=a_x b_x\) for its factorization induced by an extension \(x\), so
\[
             a_x=\sqrt{d_{ij}d_{ix}/d_{jx}},\qquad b_x=d_{ij}/a_x .
\]
If \(s\) extensions have \(a_x\) in one rank-\(r_0\) multiplicative coset, then **\(s\le8r_0+4\)**. The same holds for \(a_x/b_x\). There is no size condition on the cofactors, or on \(W/(m\log m)\). Indeed \(v_p(a_x)\) is exactly the clipped distance of the exponent of \(x\) from the anchor \(i\) on the interval between \(i,j\). This follows by substituting the three local absolute differences in the square-root formula. When the order of the anchors is reversed, replace the clipped coordinate by its complement; this is a fixed affine change at the prime. Keep all primes, so \(D=w(a_i-a_j)/W>1/2\). Now the block projections in the preceding proof belong *exactly* to an affine space of dimension at most \(r_0\). It follows directly that \(sD/2-1\le2r_0\). The ratio \(a_x/b_x=a_x^2/d_{ij}\) has the same affine valuation rank. Thus even after fixing one eligible divisor and its two endpoints, the factorizations actually used by its many triangle completions have rank proportional to their number. This conclusion also holds for an angle model with the same pair separation.

**PROVED — moving offsets in (II.1).** Dividing that equation by \(a_xu_x^2\) introduces the term \((b_x/a_x)(v_x/u_x)^2=d_{ij}(v_x/(a_xu_x))^2\). Its valuation vector gives the same clipped coordinates after a fixed affine normalization, with errors from \(v_x/u_x\). If both offsets are at most \(m^B\), the one-coordinate proof applies and gives \(s\le(8+o_B(1))r_0+O(1)\) for a coset containing the normalized solutions. Allowing these small offsets to vary does not compress their common rank.

**PROVED — consequence for the suggested use of ESS.** At a growth failure the average \(\log|t_{ij}|\) is \(O_b(\log m)\). Choose \(B=B(b)\) large enough that at most one percent of all edges fail the bound \(m^B\). Averaging over triples of anchors gives a triple with \(\gg m\) extensions having all three incident cofactors this small. By (III.4) their simultaneous rational solutions already require rank \(\gg m\); even a rank-\(o(\log m)\) coset holds only \(o(\log m)+O(1)\) extensions. The uniform ESS bound, exponential in the rank, gives less information than this linear bound in each group. This concerns this unit equation on actual angular cross ratios. It imposes no limitation on an additive system involving other observables.

**COMPUTED.** On all 35, 165 and 969 triples of ordinary Paley anchors at \(m=7,11,19\), respectively, the idealized projected valuation rows before adding the cofactor errors have the full possible affine rank \(m-4\) on the other marks. The check projects to the primes where the anchors are nonconstant and certifies this rank by modular elimination. This checks the information retained by the projection, independently of the existence of genuine small cofactors. [Anchor-rank check](iii_rank_extension.py)

**PROVED — the star identity and the fixed-weight caveat.** Identity (8) in the parallel report is stated for the simple ten-pair family, not all balanced tuples. Dividing by one term and applying \(v\mapsto v/\bar v\) to each multiplicative variable leaves anchored phase ratios up to torsion, even if the rational cofactor coefficients were absorbed into the variables. A rank-\(r_0\) coset gives an affine space of valuation signatures of dimension at most \(r_0\). By real affine independence of the exponent rows, these signatures faithfully record the four formal differences of unit vectors on the five marks. Taking an affine basis of the signatures shows that *every* mark used occurs among at most \(5(r_0+1)\) marks of the basis tuples. The resulting polynomial bound on the number of tuples per coset is already smaller than the quoted ESS bound. Fixing the finite linear forms alone also cannot give an exceptional-line count uniform in the places: for \(j=M,\ldots,2M-1\) use \((2j,1)\), the coordinate forms at infinity, and \(X+iY,X-iY\) at all finite places above divisors of \(\prod_j(4j^2+1)\) in \(\mathbb Q(i)\). With absolute values normalized for the product formula, the product of these form values is \(2j/(4j^2+1)<H(2j,1)^{-1}\). There are \(M\) distinct lines. Individual log heights in this counterexample are small compared with the bound for the parent; the fixed forms by themselves lose this distinction. The fixed-weight version of the parametric Subspace theorem asks for additional common local weights. For the natural two Gaussian forms, the supplementary proof bounds the number of edges on which one fixed ray of local weights can gain \(\delta W\) by \(O(\delta^{-5})\), uniformly in the parent and the varying height parameter. It gives a separate obstruction to batching pairs under that prescription. [Star-rank argument and scope](iii_star_rank.md) · [Twisted-height argument and proof](iii_star_targets.md) · [Parallel report, equation (8)](../lemma/answer.md)

## III.5. What one repeated cofactor can force

**PROVED — density threshold.** Suppose the \(E=\binom m2\) positive integer edge labels have mean logarithm at most \(H\). If the largest multiplicity is \(f\), the \(j\)-th smallest integer in the multiset is at least \(j/f\). Thus \(f\ge E e^{-1-H}\) by the factorial bound. The normalized arc and the mean separation bound give \(H\le W/(4(m-1))+O(1)\). At \(W=(3+o(1))m\log m\), a color has at least \(m^{5/4-o(1)}\) edges. This density does force a ten-cycle by the Bondy–Simonovits bound \(\mathrm{ex}(m,C_{2k})=O_k(m^{1+1/k})\). It does not force a four-cycle by edge count alone; incidence graphs of finite projective planes have order \(m^{3/2}\) edges and no four-cycle. As a calibration, a uniform exclusion of monochromatic ten-cycles would force \(W\ge(16/5)m\log m-O(m)\) for these normalized one-unit clusters. To improve the published constant for arbitrary arcs one must account for any losses in shortening or passing to a unit class. For arbitrary fixed \(b>8\), the log budget alone even permits every edge label to be different.

**PROVED — parity and two equal offsets.** Every fixed-cofactor graph is bipartite. Indeed \(v_{1+i}(z_i-z_j)=2+2v_2(t_{ij})\): the difference is \(2i t_{ij}\) times a Gaussian factor of odd norm, up to sign. For fixed valuation, edges go between the two lifts of one residue class to the next \((1+i)\)-adic level. This excludes every odd cycle, even with only the valuation fixed. If two edges incident with an anchor have the same absolute cofactor \(t\), the edge between the other ends has cofactor divisible by \(2t\). Write each ratio from the anchor as \((s_j+it)/(s_j-it)\) using signed real parts. The norm \(g\) of the Gaussian gcd of the two numerators divides \(t(s_j-s_k)\), by taking determinants after factoring out that gcd. It is odd and coprime to \(t\) by primitivity. Cancellation gives cofactor \(t|s_j-s_k|/g\); the real parts have the same parity. For two fixed anchors with base edge \(d=ab=s^2+w^2\), a common neighbor of color \(t\) consequently has \(h=|w|/t\) integral and, from (II.1),
\[
             c=(a+b\pm2s)/h^2\mid N/d,\qquad a\equiv\mp s\pmod h .
\]
The last congruence follows from \(ac-t^2=(t(a\pm s)/w)^2\) and integrality of a rational square root of an integer. It is the complementary-parent condition that has to be shared by many neighbors. Parity alone admits dense bipartite graphs.

### A balanced four-cycle with the smallest primary-class cofactor

**PROVED using the exact COMPUTED certificate below.** There are one-unit, Gaussian-gcd-primitive four-point clusters on arbitrarily large norms \(N\) with **all seven nonconstant four-mark cuts of weight \((\log N)/7+O(1)\)**, normalized angular width \(\asymp N^{-2/7}\), and the same primitive cofactor **two** on every edge of a \(K_{2,2}\). In the primary prime convention, two is the smallest positive cofactor on a one-unit edge: its primitive factor is a real unit times a product of primary factors, whose imaginary coordinate is even. This gives actual integer common-parent examples at the four-mark exponent obtained by restricting an asymptotically balanced growth-scale profile.

*Proof.* Use the parameters \(x=-5,m=3\) of Part I §6 and clear the root denominators by six. Here is the complete resulting construction. Index the following Gaussian roots starting at zero:
\[
 (0,6),\quad(-135,-36),\quad(-135,-51),\quad(-135,75),\quad
 (-145,81),\quad(-117,-45),\quad(-108,-30),
\]
where \((u,v)=u+iv\). Put \(H_j=n-\rho_j\). The four rows are
\[
 S_0=\varnothing,\quad S_1=\{0,1,2,4\},\quad
 S_2=\{0,1,3,5\},\quad S_3=\{0,2,3,6\}.
\]
Take \(z_i=\prod_{j\in S_i}H_j\prod_{j\notin S_i}\bar H_j\) and divide the four by their common Gaussian gcd. For \(i<j\) their ratio is \(B_{ij}/\bar B_{ij}\), where \(B_{ij}\) uses \(H_k\) on \(S_i\setminus S_j\) and \(\bar H_k\) on \(S_j\setminus S_i\). The two sides have four factors together. Each resulting monic quartic has constant nonzero imaginary part \(b_{ij}\). In lexicographic edge order \(01,02,03,12,13,23\), multiplying the four linear factors gives
\[
 (|b_{ij}|)=(19966500,16110900,14871600,3855600,5094900,1239300).
\]
In particular the angular width is \(\asymp n^{-4}\). On the arithmetic progression
\[
 n\equiv79658575767\pmod{86434978500}
\]
the four gcds \(\gcd(\Re B_{ij}(n),b_{ij})\) for \(ij=01,02,13,23\) are respectively \(9983250,8055450,2547450,619650\). Verification at one residue suffices since the period is the lcm of the four corresponding \(|b_{ij}|\). After dividing out the gcd, each primitive near-real representation has imaginary coordinate of absolute value two and odd real coordinate. Thus the four relevant norms are odd and every cross edge has cofactor two. Ratios of these near-real factors to their conjugates have primary unit one, so connectivity of the four-cycle puts all four points in one primary unit class.

For the cut weights, the seven root signatures run once through all seven nonempty subsets of marks \(1,2,3\), with mark zero as anchor. The seven norms \(\mathrm NH_k\) are \(\asymp n^2\). All fourteen roots including conjugates are distinct. Any common Gaussian factor of two of the corresponding fourteen linear values divides their fixed nonzero root difference. Hence, after discarding a bounded Gaussian factor from each linear value, the fourteen residual ideals are pairwise coprime, conjugation permutes the residual ideals, and every residual prime contributes only the prescribed cut. For example, take a fixed ordinary integer divisible by all the pairwise root differences in \(\mathbb Z[i]\), and remove from each prime valuation the smaller of the valuation itself and that fixed bound. Only bounded factors have been removed. Restoring them affects at most a bounded number of weighted layers in total and changes the common gcd by a bounded factor. This proves \(N\asymp n^{14}\) and gives each cut weight \(2\log n+O(1)\), as claimed.

**COMPUTED.** The stand-alone script multiplies the quartics over \(\mathbb Z[i]\), checks the constant offsets and the progression, then produces the four Gaussian integers at the displayed residue and computes their primitive half-angle tangents directly. Their six absolute offsets are \((2,2,7344,560,2,2)\) and the norm has 135 digits. Rational inequalities certify critical arc constant smaller than \(1/128\). The seven prime blocks can be composite, as in a growth profile; no claim of squarefree values is needed. [Balanced four-cycle certificate](iii_tri_balanced.py)

**PROVED — the numerical calibration.** In a large squarefree-profile model with sixteen unoriented cut weights for five marks (or with fair bits on four marks), a balanced four-subtuple has primitive log norm \((7/8+o(1))W_{\rm parent}\) after deletion of its constant cut. A parent arc of width \(\ll\exp(-W_{\rm parent}/4)\) has the four-mark exponent \(2/7+o(1)\). Our example realizes this scale with exact cofactors on a cycle. Thus a forbidden-*unordered*-four-cycle argument cannot be justified by those four-point consequences of multiplicity, even after imposing all seven cut weights and the integer origin of the star matrix. A restriction on how many such cycles can share vertices of a much larger parent would go beyond this example.

## III.6. An integer matrix retaining half the prescribed divisibility on every pair

**PROVED — partial integral fake.** For each fixed \(L>0\) and arbitrarily large \(m\), there are a fixed squarefree split parent with \(W=O_L(m\log m)\), the Paley profile of III.2, and \(m\) primitive integer columns \(V_i\), with these simultaneous properties:

1. Their doubled arguments are contained in an arc of angular width \(\tfrac12 N^{-1/4}\), and \(\big||V_i|^2/n_i-1\big|<m^{-L}\) at every mark.
2. For all pairs, \([V_i,V_j]/g_{ij}\) is within \(m^{-L}\) of a nonzero signed integer \(b_i-b_j\). All positive absolute values are different. One can require the same error for the scalar angular quantity \(\sqrt{d_{ij}}\sin(\arg V_j-\arg V_i)\).
3. On a single set \(S\) of columns, the Gaussian valuations of *all* the \(V_i\) are exactly the prescribed zeros and ones, on a consistent choice of split factors. In particular
\[
       g_{ij,S}=\prod_{p\in S}p^{\min(a_{ip},a_{jp})}\mid[V_i,V_j],
       \qquad \log g_{ij,S}=\tfrac12\log g_{ij}+O_L((\log m)^2)
\]
for every pair. With a row-dependent choice of omitted first-family constraints, the error in the log weight of the guaranteed subproduct on each pair can be reduced to \(O_L(\log m)\).
4. The prime assignment can be chosen so that **no quartet** admits a matrix satisfying its *full* prescribed divisibilities with pair cofactors of size at most \(3m^3+1\). At this very assignment the original F1 angle model exists. The statements about the columns \(V_i\) use their actual directions. The separate F1 assertion allows its own prime angles.

*Proof.* Write \(K=L+10\), \(A=4L+60\). As in III.2, choose \(2m\) split primes in an interval \([P,P\exp(m^{-K-1})]\subseteq[m^A,3m^A]\). Indeed partition the log interval from \(m^A\) to \(2m^A\) into \(O(m^{K+1})\) bins and use the prime number theorem in the progression \(1\bmod4\). Set \(H_0=P^{m/2}\), \(G=P^{(m-1)/2}\), so \(G/H_0=P^{-1/2}\). Each ordinary row contains \(m\) ones in both families together, and \((m+1)/2\) ones in the first family with diagonal sign positive. Every pair shares \((m-1)/2\) ones across the two families, of which \((m+1)/4\) are in the first family. For the first-family count use the off-diagonal character sum \(-1\) in III.2. Replacing each of its two zero entries by the diagonal sign adds opposite terms since \(\chi(-1)=-1\), so two first-family rows have dot product \(-1\) and each has sign sum one. The two-family count follows from row size and separation. Uniformly,
\[
 N^{1/4}/H_0,\quad \sqrt{n_i}/H_0,\quad g_{ij}/G
       \quad\text{are all }1+O(m^{-K}).
\]

Delete \(O(\log m)\) first-family columns so that every row loses at least three ones. For example sample \(t=\lceil 10\log_2 m\rceil\) columns without replacement. At all \(t\) draws the conditional probability of a zero is at most \(0.6\) for large \(m\), so the chance any of the \(m\) rows has at most two ones among the deletions is at most \(m(1+t+t^2)0.6^{t-2}=o(1)\). Let \(S\) be the remaining first-family columns. At any chosen prime fix a root \(r_p^2\equiv-1\pmod p\). The congruence \(X_i\equiv r_p\pmod p\) prescribes one Gaussian prime above \(p\) in \(V_i=X_i+i\). The two roots are distinct. For each row impose this congruence at its selected ones and impose \(X_i\equiv0\pmod 2\). The CRT period is \(2M_i\), where
\[
 M_i\le (1+o(1))P^{(m+1)/2-3}=(1+o(1))G/P^2 .
\]

Take \(b_i=2m^2i+i^2\) for \(0\le i<m\) and target \(X_i^*=H_0+G b_i\). The positive differences \(b_j-b_i\) lie in disjoint intervals by the value of \(j-i\), since \(0<j^2-i^2<m^2\). With \(j-i\) fixed they determine \(i\). From a CRT representative nearest the target try \(2m+2\) consecutive representatives. At a selected one, removing divisibility by the square of its Gaussian prime forbids one class of the progression parameter modulo \(p\): the derivative of \(X_i^2+1\) is nonzero modulo \(p\), and \(2M_i\) is divisible exactly once by \(p\). At a selected zero, avoiding both Gaussian factors forbids two classes. Each prime exceeds the entire search length, so fewer than \(2m+2\) candidates are forbidden altogether. This produces \(X_i=X_i^*+O(mG/P^2)\). The integers \(V_i=X_i+i\) are primitive.

Their determinants are simply \(X_i-X_j\). The normalized error is \(O(m^{3-K}+m/P^2)\), and the relative error in \(|V_i|^2/n_i\) is \(O(m^3/\sqrt P+m^{-K})\). Multiplying a normalized determinant by \(\sqrt{n_i n_j}/(|V_i||V_j|)\) gives the exact angular quantity in assertion 2, so the additional absolute error is \(O(m^6/\sqrt P+m^{3-K})\). Their doubled angle range is at most \(2(\max X_i-\min X_i)/(\min X_i)^2\ll m^3/(H_0\sqrt P)\). All the claimed inequalities follow with the chosen \(A,K\). The retained common-prime count differs from \((m+1)/4\) by at most the number of deleted columns. Alternatively delete just three marked first-family primes separately in each row and use the constraints shared by the two rows; at most six common primes are lost on any pair. Finally, all choices of prime assignment were allowed throughout this construction. Choose one excluded on every quartet by III.2. Since \((3m^3+1)^2<P/4-1\) eventually, any four columns satisfying the full conditions with these cofactor sizes would pass its test. This proves assertion 4.

**COMPUTED.** The horizontal construction at \(m=19\) gives 19 primitive columns on a 38-prime parent. Eleven columns of the cut array have exactly the prescribed Gaussian valuations. Each pair shares one to five of its nine prescribed primes there. The largest of the 171 normalized determinant errors is at most \(116082/10^9\), all 171 rounded absolute values are distinct, and all 171 full divisibilities fail. Exact rational inequalities verify the doubled arc bound and a 3% upper bound on squared-norm errors. A second construction in a tilted strip gives squared-norm errors below \(10^{-6}\), normalized determinant errors below \(1/990\), and at least one failed full divisibility. [Horizontal integer certificate](iii_horizontal_stars.py) · [Tilted-strip certificate](iii_graph_stars.py)

**PROVED — the full-isotropy endpoint.** If instead all prescribed Gaussian factors divide primitive \(V_i\), their product \(P_i\) divides \(V_i\). If \(|V_i|^2/n_i=1+o(1)\) uniformly, the quotient has positive integer norm \(1+o(1)\) and is a unit. Even a bounded *mean* log slack suffices for a positive fraction: if \(m^{-1}\sum_i\log(|V_i|^2/n_i)\le C\), at least half have quotient norm \(\le e^{2C}\); a positive fraction share a single Gaussian quotient. For these marks \((\prod_p\bar\pi_p^{e_p})P_i/\bar P_i\) are Gaussian integers of norm \(N\), and a common rotation takes their arguments to the doubled arguments of the \(V_i\). Provided the original directions are projectively distinct, these are distinct genuine points in the prescribed arc. Thus a full isotropic integral fake at bounded mean slack would give an actual large cluster. Exact divisibility of minors alone permits arbitrary projective lines at the primes and is a strictly weaker local requirement. For that model the product of all minors and ordinary residue-class collisions still give the leading coefficient three. The supplementary note includes the proof under the aggregate height hypothesis and the precise loss as the discarded isotropic weight increases. [Projective-line estimates and horizontal CRT formulation](iii_star_crt.md)

**PROVED — with exact star norms the local isotropy is automatic.** Suppose primitive integer columns have \( |V_i|^2=n_i\) and satisfy every prescribed minor divisibility. For each split prime \(p\), every column with \(a_{ip}>0\) reduces to one of the two isotropic lines modulo \(p\), since its norm vanishes and its reduction is nonzero. The shared divisibility forces all these columns onto the same line. Primitivity prevents a column from containing both conjugate factors of \(p\), so its entire valuation \(a_{ip}\) lies on the corresponding factor \(\pi_p\). Consequently \(V_i=u_i\prod_p\pi_p^{a_{ip}}\) with Gaussian units \(u_i\), and \((\prod_p\bar\pi_p^{e_p}) V_i/\bar V_i\) are actual points on the common parent, with doubled star angles up to a common rotation. Thus if the proposed integral fake includes the *exact* individual norms, it is equivalent to the original cluster already. With only strip geometry and approximate norms, III.6 leaves the full minor-divisibility case open.

## III.7. Global projective reconstruction for arbitrary large profiles

**PROVED — the growth-specific conclusion.** Let \(m\) be sufficiently large, with an arbitrary prime-power exponent array satisfying \(w(a_i-a_j)>W/2\) for distinct marks. Supply the parent and the exact labeled integers \(d_{ij}\). Up to **one transformation in \(\mathrm{PGL}_2(\mathbb Q)\)** there is at most one projective configuration of \(m\) distinct rational star directions realized by integer columns with
\[
             [V_i,V_j]=g_{ij}t_{ij},\qquad 0<|t_{ij}|\le
                       \sqrt{d_{ij}}/(4 N^{1/4}).
                                                               \tag{III.5}
\]
Here each realization may use its own column lengths, and its prime orientations may be chosen separately provided the exact separation divisors stay fixed. The conclusion even survives deleting a specified collection of constraints in (III.5), with at most \(m/10\) deleted partners at each mark. With \(o(m^2)\) arbitrary deleted pairs it gives uniqueness on \(m-o(m)\) marks. This is a uniqueness assertion whenever a realization exists. It places no height bound on the common projective transformation. The underlying integer divisibilities, together with multiplicity, are essential to the proof.

*Proof of local uniqueness.* For a quartet with a priori bound \(U\) on the absolute values of the three cofactor products and with \(B>2U^2\), there is at most one rational projective triple \([x:y:z]\) satisfying (III.1) with such coefficients. Indeed \(B\mid xz'-x'z\) for any two solutions, since \(A,C\) are coprime to \(B\), and the determinant has absolute value at most \(2U^2\). Thus \(x:z=x':z'\), and the equation fixes the remaining ratio. Consequently the full projective triple of matching products of the six raw minors is fixed, namely \([Ax:By:Cz]\). This gives the rational cross ratio of the four ordered directions. The three matching factors \(A,B,C\) are computable from the six separation divisors by III.1. In particular, the cross ratio is unchanged if another oriented exponent array has the same \(d_{ij}\). This uses their integer factorizations and the small-solution condition.

*Proof that sufficiently many quartets give uniqueness.* On the threshold space with weights \(\log p/W\), let \(f_i\in\{-1,1\}\) be the sign functions of III.4. Put \(G_{ij}=\langle f_i,f_j\rangle=1-2\log(d_{ij})/W\) and fix \(\tau=1/32\). The frame bound proved in III.4 gives
\[
               \sum_i |\langle u,f_i\rangle|^2\le2\|u\|_2^2 .
\]
Call a pair good if \(|G_{ij}|\le\tau\). Each mark has at most \(2/\tau^2\) other partners failing this test. Given three pairwise good marks \(i,j,k\), at most \(8/\tau^2\) extensions \(x\) fail one of the four tests against \(u=f_i,f_j,f_k,f_if_jf_k\). All four functions have norm one. For the other extensions every pair is good and \(|\int f_i f_j f_k f_x|\le\tau\). On these quartets the log weight of each two-two matching is at least \((1-7\tau)W/8\): expand the indicator of a sign pattern together with its negative as in III.2, obtaining six pair moments and one four-sign moment. Also (III.5) bounds every edge cofactor by \(e^{\tau W/4}\), since \(\log(\sqrt{d_{ij}}/N^{1/4})=-WG_{ij}/4\). Taking \(U=e^{\tau W/2}\), we have
\[
              \log B-\log(2U^2)\ge [(1-7\tau)/8-\tau]W-\log2>0
\]
for large \(W\). Large \(m\) forces large \(W\), for example by the distinctness of the exponent rows. Thus the local criterion holds on every one of these quartets for which all six constraints are retained.

*Gluing the directions.* First suppose all constraints are present. Greedily take three pairwise good anchors; every choice rules out only a bounded number of candidates. Identify their three distinct rational directions with \(\infty,0,1\) by a projectivity. For all but at most \(8/\tau^2\) extensions, the single recovered cross ratio with these anchors now fixes the direction in that chart. Given any missed mark \(x\), choose from the recovered marks a \(y\) good with \(x\), then a \(z\) good with both, then a \(w\) passing all four extension tests for \(x,y,z\). This succeeds once \(m\) exceeds an absolute constant. The three already known directions of \(y,z,w\) and this last quartet fix \(x\). Any two realizations therefore agree in the same chart everywhere. The argument can be implemented using just the divisors: test the numerical matching sizes and pair bounds, find the small kernel ray, propagate from three anchors, and check the resulting cross ratios. The sign functions certify the availability of the needed quartets.

If at most \(\eta m\) pairs at each mark were deleted, also exclude deleted partners at each choice. The first chart misses at most \(3\eta m+O(1)\) marks. Choosing \(y,z,w\) for \(x\) inside the recovered marks rules out at most another \(3\eta m+O(1)\) at any step. Thus \(6\eta<1\) suffices, including \(\eta=1/10\). With only a global \(o(m^2)\) deletion bound, discard the \(o(m)\) rows incident with more than \(m/20\) deletions, and apply this argument on the remaining array with the same layer weights and frame inequality. This proves the theorem.

**PROVED — almost every parent in a large pool excludes every assignment.** A separate strengthening of III.2 allows the parent itself to vary. Suppose a narrow bin at scale \(P=m^{O(1)}\) contains a pool of \(K\ge m^{1+\eta}\) split primes, for a fixed \(\eta>0\). Draw \(r=2m\) without replacement. Except for a fraction \(\exp(-\Omega_\eta(m\log m))\) of the resulting parents, **every prime assignment fails on every Paley quartet**, even with independent bounded coefficients in each test. To see this, order the sample and prescribe the three disjoint sets of indices for the matching products. There are at most \(4^r\) choices in all. For any candidate with \(k\ge r/16\) indices in the middle product, expose the primes on its complement. For each nonzero integer coefficient triple of absolute values at most \(U\le P\), the equation forces the complete middle product \((Ax+Cz)/y\). By unique factorization this allows at most one of the \(\binom{K-r+k}{k}\) equally likely remaining subsets. The union bound is at most \(4^r(2U)^3/\binom{K-r+k}{k}\), using the smallest denominator over the allowed \(k\). It has the asserted size. Such bins exist by the same prime-number-theorem averaging used in III.2, on taking the fixed exponent of \(P\) sufficiently large. This still leaves the worst parent open.

**COMPUTED — a planted passing quartet that cannot be extended.** To distinguish an empty local test from an obstruction due to sharing primes, I constructed a second fourteen-prime parent containing three disjoint triple products \(A,B,C\) deliberately chosen to pass. Each prime in a triple is a translate of one of three offsets of sum zero and fixed sum of squares; only the constant terms in the cubic products differ. In the order used in the script the primitive relation is
\[
                   -1145443 A+1927968 B-782525 C=0 .
\]
After scaling by four, each of the three coefficients factors as a product of two even integers at most \(250000\) in absolute value, the common relaxed individual bound. Among all 560,560 unordered triples of three disjoint three-prime products on this **same** parent, exactly **one** admits any signed nonzero solution below the common product bound \(U=62500076862\). Each of the 21 five-mark restrictions of the seven-mark Paley profile requires two three-per-class quartets with different partitions of the primes. The script computes their overlap tables, including the unused primes, allowing free permutations of the three matching names. No pair of passing prime partitions has the required table. Hence some quartet passes, even with individually factorable even cofactors, while **no five marks pass on any assignment** of this parent. Passing here means satisfying the stated relaxed tests; it requires no claim of a genuine four-point cluster. As a further check, at the first fourteen-prime parent of III.2 there are 16, respectively 12, surviving two-per-class partitions when each coefficient must split into two bounded factors, respectively two bounded even factors after common scaling. Local uniqueness justifies testing the primitive coefficient triple: a divisor of a product of two bounded positive integers itself splits into two such factors. [Planted-parent and five-mark certificate](iii_quartet_extension.py) · [First-parent certificate](iii_quartet_exhaust.py)

## III.8. Five-mark compatibility can come for free

**PROVED — cancellation with three independently solved quartets.** Fix distinct \(i,j,k,l,h\). For each of the quartets \(ijkl,ijlh,ijkh\) suppose there is a signed nonzero integer solution of (III.1), with maximum coefficient bounds \(U_1,U_2,U_3\), respectively. These can be independent solutions. For a quartet with anchor pair \(i,j\), denote its inferred raw cross ratio by
\[
 R_{ij;kl}=\frac{g_{ik}g_{jl}}{g_{jk}g_{il}}\frac{a_{kl}}{b_{kl}},
\]
where \(a_{kl},b_{kl}\) are the signed coefficients of the crossed matchings, with signs taken in the order \([ik][jl]/([jk][il])\). Let \(Q\) be the product of the weighted prime layers putting \(i,j\) on one side and \(k,l,h\) on the other, each layer contributing its prime. If \(Q>2U_1U_2U_3\), then
\[
                         R_{ij;kl}R_{ij;lh}=R_{ij;kh}.                  \tag{III.6}
\]

*Proof.* Reduce the displayed quotient of \(g\)-products in each of the three ratios to \(C_s/D_s\). At each prime of \(Q\), the two crossed matchings have the same local exponent in each quartet and the uncrossed matching contains at least the full corresponding additional exponent of \(Q\). Thus \(C_s,D_s\) are units modulo \(Q\), and (III.1) gives \(C_sa_s\equiv D_sb_s\pmod Q\). This uses the identity with signs \([ik][jl]-[jk][il]=[ij][kl]\). The three quotients of \(g\)-products telescope **exactly**: \(C_1C_2D_3=D_1D_2C_3\). Multiply the first two congruences and the reversed third congruence and cancel this common unit. We obtain \(Q\mid a_1a_2b_3-b_1b_2a_3\), whose right side has absolute value at most \(2U_1U_2U_3\). This proves (III.6) without inverting a single cofactor. III.3 compared two quartets to build a modulus for a third. Here the *third small integer solution* is already present and cancels the entire height of the crossed matching factors.

**PROVED — a five-mark consequence.** Call a set regular if all its off-diagonal pair moments and all its four-distinct-sign moments have absolute value at most \(\tau=1/64\), in the notation of III.7. Suppose the five quartets on a regular five-set each pass independently with coefficient bounds \(U=e^{\tau W/2}\); this bound is forced by (III.5). Then the five recovered cross ratios come from one ordered configuration of five **distinct** rational projective points when \(W\) is large enough. Indeed the sign-pattern indicator up to complement on five marks is \(1/16\) times a sum of its constant, ten signed pair moments and five signed four-sign moments. Every two-against-three modulus has log weight at least \((1-15\tau)W/16>3\tau W/2+\log2\) for large \(W\). Thus (III.6) holds with every choice of anchor pair. Fix the first three directions, and place the fourth and fifth by their cross ratios with those three. For each pair of the first three marks, apply (III.6) with the third as an auxiliary. This gives exactly the three remaining quartet cross ratios of the placed directions. They are finite, nonzero, and different from one by the local three-term equations, which also excludes coincident directions. Projective changes of ordering cause no ambiguity because each local solution satisfies the ordinary three-term identity.

**PROVED — gluing on the regular quartets of a large profile.** Under the pair separation and frame hypothesis of III.7, if every regular quartet passes this independent small-integer test, there is a single collection of rational directions realizing all the regular quartet cross ratios. Marks forming a good pair in the sense \(|G_{ij}|\le 1/64\) have distinct directions. This conclusion allows coincidences on untested pairs, and it does not furnish the integral column lengths or the required cofactors.

Here are details of the gluing. Any regular set of a fixed size \(d\) has all but \(2(d+\binom d3)/\tau^2+O(d)\) marks available as regular extensions. Apply the frame inequality to its individual signs and its products of three signs. For each regular triple \(T\) take its rational chart on every regular extension, placing \(T\) at \(\infty,0,1\) in a chosen order. On a finite regular set every five is consistent by the preceding result. Charts on its triples then agree up to the unique transition projectivity: to replace one member of a triple, its union with the next triple has four members, and the consistency on this union together with any additional mark transports each value. Successive single-member replacements connect any two triples of the set. Distinctness in such a finite chart follows from the nonzero three matching products. This reasoning applies also to sets of size three or four directly.

For possibly incompatible triples \(T,T'\), choose an auxiliary triple \(J\) with both \(T\cup J\) and \(T'\cup J\) regular, by the extension count applied to these two sets separately at each draw. Identifying the two charts on the three directions of \(J\) defines their transition. It agrees on all but \(O(1)\) common marks, since all but \(O(1)\) extend both unions regularly and the two finite charts have the same cross ratio on the auxiliary triple with that mark. Its cofinite agreement determines it uniquely: in any cofinite subset of a chart take three marks making a regular set jointly with its anchors by the same greedy extension argument, hence giving three distinct values. Consequently the transitions obey composition, as is seen on a cofinite intersection of any three charts. They also agree at *each* mark \(x\) in the domain of both charts. For this, choose the auxiliary triple by extending the two regular sets \(T\cup\{x\}\) and \(T'\cup\{x\}\) in parallel. The resulting transition is the unique one just determined and it agrees at \(x\). Every mark can be put in a regular triple. The resulting assignment realizes every regular quartet by using three of its marks as a chart, and every good pair extends to a regular quartet. This proves the assertion for sufficiently large \(m\).

**HEURISTIC — consequence for the next experiment.** On balanced five-marks of a growth profile the shared two-against-three block is exponentially larger than the product of three polynomial coefficient bounds. Once each quartet has passed individually, testing the rational cross-ratio identities again on such a five-tuple yields no further rejection. The nontrivial extension test is to fix the *individual edge cofactors*, their divisibility and their heights on a passing base and extend those same integers on every new quartet. Another target is to prove that independent small integer rays fail to exist on some regular quartet of every large parent. The planted example in III.7 illustrates that latter failure with shared primes.

## III.9. A uniform bound for almost all parents in a prime-product family

**PROVED — fixed-size bound allowing every cut profile.** Fix \(D\ge4\). Draw \(r\) distinct primes uniformly from the split primes in \([r^D,2r^D]\), and take \(N\) to be their product. Except for a fraction \(\exp(-c_D r\log r)\) of these squarefree parents for large \(r\), **every normalized cluster has at most 17 marks**. This holds even for arbitrary binary cuts carrying integer star minors as in (III.5). Only pair separation \(\log d_{ij}>W/2\), the height bound there and the quartet relations are needed. For sufficiently large fixed \(D\) the same proof gives 13 in place of 17. This is a density statement in the specified family of parents. No asymptotic counting statement for all integers \(N\) or for their worst case is implicit. The same bound of 17 holds on actual arcs of length \(C\sqrt R\) for each fixed \(C\), by the unit adjustment at the end of the proof.

*Proof.* Consider any \(M\) distinct marks, with \(M\) a fixed prospective forbidden size. Write \(P=r^D\). The available prime pool has size \(K=r^{D+o(1)}\) by the prime number theorem in the split progression, and \(W=(D+o(1))r\log r\). For a column of size \(k\) set \(u=k(M-k)\). Summing pair separation over the \(\binom M2\) pairs gives the prime-weighted mean of \(u\) at least \(M(M-1)/4\). The unweighted mean is at least \(M(M-1)/4-o(1)\), since \(M\) is fixed and all log weights are within \(\log2\). For a random ordered quartet, let \(s\) count the primes in its middle matching class. A given cut contributes to this class with probability \(2u(u-M+1)/(M)_4\), by counting two points on each specified side. Convexity, and monotonicity at the lower bound on the mean for \(M\ge 4\), give
\[
           \mathbb E s/r\ge F_M-o(1),\qquad
           F_M=\frac{(M-1)(M-4)}{8(M-2)(M-3)} .
\]
For this quartet write \(V=|xyz|\), the product of all six cofactor magnitudes. The mean over all edges of \(\log|t_{ij}|\) is at most \(W/(4(M-1))\), since \(u\le M^2/4\) bounds the mean pair distance, and (III.5) applies. Thus \(\mathbb E\log V\le 3W/(2(M-1))\).

Set \(H=\log((K-r)/r)=(D-1+o(1))\log r\). Whenever
\[
                        (D-1)F_M>\frac{3D}{2(M-1)},                       \tag{III.7}
\]
some ordered quartet satisfies \(sH-\log V\ge c r\log r\) for a fixed positive \(c\). To bound the probability that *any* profile can supply it, order the sampled primes and prescribe their four-way allocation to the three matching products and the unused part. Expose all primes outside the middle part, of size \(s\). The number of signed nonzero triples with \(|xyz|\le T\) is at most \(8T(1+\log T)^2\) for \(T\ge1\), by summing \(T/(|x||y|)\). Here \(T=\exp(sH-cr\log r)\); if \(T<1\) there are no candidates. As in III.7, each triple forces at most one subset of the remaining pool to furnish the middle product, with probability \(1/\binom{K-r+s}{s}\le e^{-sH}\). The \(4^r\) allocations cover every choice of the marks, the cuts and the quartets. Summing these bounds gives \(\exp(-\Omega_D(r\log r))\). Inequality (III.7) holds for \(M=18\) and \(D\ge4\), and for \(M=14\) at sufficiently large \(D\). This proves the result.

For an actual arc with several Gaussian units, absorb the unit of each point into its star by multiplying by one of \(1,1+i,i,i(1+i)\). Their quotients by their conjugates are \(1,i,-1,-i\), respectively. The odd shared divisor \(g_{ij}\) still divides the new integer determinant, and the new column lengths increase by factors at most \(\sqrt2\). Distinct points give nonzero determinants. Thus \(1\le |t_{ij}|\ll_C \exp((\log d_{ij})/2-W/4)\), implying pair separation \(\log d_{ij}\ge W/2-O_C(1)\). These changes cost \(O_C(1)\) in the mean of \(\log V\) and \(O_C(1/W)\) in the fixed-\(M\) mean cut estimate. Both are absorbed by the strict margin in (III.7).

**PROVED — extension to lower prime scales and bounded exponents.** Condition (III.7) gives a uniform constant for each fixed \(D>1\), since \(F_M\to1/8\). At \(D=3/2\) it already holds for \(M=38\), giving at most 37 marks for almost every such parent. More generally, on a random support in this dyadic pool the adversary may choose all exponents \(1\le e_p\le E\) after seeing the primes, with \(E\) fixed. The proof works for every such choice simultaneously if \((D-1)F_M/E>3D/(2(M-1))\). Indeed apply the cut calculation across the \(S=\sum e_p\) threshold layers. The middle matching then has expected total exponent at least \((F_M-o(1))S\), and at least its total exponent divided by \(E\) distinct primes. A prime contributes to just one two-two matching of the quartet, by III.1, so specifying that matching and its exponent takes at most \(1+3E\) options per prime. After exposure, the forced middle product still determines its subset of primes uniquely. Finally \(W=(D+o(1))S\log r\), so the same strict entropy margin absorbs the \((1+3E)^r\) possibilities. Increasing exponents without a fixed cap is outside this estimate.

## III.10. What remains after reconstructing the directions

**PROVED — an exact test once the individual cofactors have been chosen.** On a set of distinct recovered rational directions, take any nonzero rational representative \(u_i\) on each. Put \(h_{ij}=[u_i,u_j]/g_{ij}\). The possible families of rational cofactors with these projective directions up to projectivity are exactly
\[
                 t_{ij}=c\,s_i s_j h_{ij},\qquad c,s_i\in\mathbb Q^*.
\]
Indeed a change of the two ambient coordinates contributes its determinant \(c\), and a change of individual representatives contributes the \(s_i\); each rational \(c\ne0\) occurs. Requiring all \(t_{ij}\) to be bounded nonzero integers is therefore a simultaneous condition on one rational scale per mark and one common scale. This gives a concrete meaning to keeping the six base-edge cofactors fixed while choosing the four new cofactors incident to each fifth mark.

If such a family has been found, put \(D_{ij}=g_{ij}t_{ij}\) and \(d=\gcd_{i<j}|D_{ij}|\). Its raw integer minors have an integer two-row realization. Indeed their decomposable rational two-vector determines a rational row plane in \(\mathbb Q^m\). Its intersection with \(\mathbb Z^m\) is a direct summand (the quotient is torsion-free), so a basis matrix \(B\) of this intersection has maximal minors of gcd one. Its wedge is a primitive vector spanning the wedge line. Orient the basis so its minors are \(D_{ij}/d\). Multiplication on the left by an integer matrix of determinant \(d\) realizes the desired minors, and every integral realization is obtained this way since its rows lie in this intersection.

There is also an exact condition for obtaining **primitive columns**. Every column \(B_i\) must be primitive, and for every prime \(p\mid d\) the nonzero reductions of these columns must miss a point of \(\mathbb P^1(\mathbb F_p)\). For necessity reduce the left multiplier modulo \(p\): it is singular, and to preserve primitivity its kernel must miss every column. For sufficiency choose at each \(p\mid d\) a nonzero linear functional vanishing on one missing line. Lift these functionals to a common primitive integer row \((a,b)\): choose a nonzero \(a\) by CRT, then choose \(b\) with its prescribed residues and with \(b\equiv1\) at all additional primes dividing \(a\). Complete this row to a unimodular matrix and multiply its second row by \(d\). This preserves every column's primitivity at the primes of \(d\), and at the other primes the multiplier is invertible. When exact lengths \(n_i\) are imposed, an additional single symmetric form must satisfy \(B_i^t S B_i=n_i\) for every \(i\); necessarily \(S=M^tM\) for the same left multiplier \(M\). Three distinct projective directions already determine that form, since a nonzero binary quadratic has at most two projective zeros.

**HEURISTIC — request for the next check.** I would ask for independent scrutiny of the five-mark cancellation (III.6) and the uniform parent count (III.7). For a new arithmetic test I would search the exceptional parents, starting with shared passing three-per-class partitions in the seven-mark example, then impose the common vertex scales above, primitivity, and the single quadratic form. Counting cross-ratio failures on independent, typical, already-passing five-tuples would miss the remaining obstruction. I have neither an unbounded full integral fake nor a worst-parent saving in the growth rate or its leading constant.

## III.11. A full five-mark integral countermodel at exponent 10/3

**PROVED using the exact coefficient check below.** The exceptional local model can actually be realized on five marks. There are arbitrarily large split parents with ten equal-weight pair cuts up to \(o(1)\) in log weight, valuation index one, and a primitive integer star matrix such that **every** prescribed \(g_{ij}\) divides its minor, every normalized cofactor has bounded size, and
\[
   |V_i|^2=(1+o(1))n_i=(1+o(1))\sqrt N,\qquad
   \operatorname{diam}\{2\arg V_i\}\asymp N^{-3/10}.
\]
All data here belong to the **same** integer matrix on each parent. At least one norm differs from its prescribed integer value, so at least one local line is non-isotropic. Thus retaining integral stars, all divisibilities and asymptotically exact star radii does not improve the five-mark ten-pair exponent. This is a fixed-size model, with possibly composite prime blocks. No F1 realization with these same individual directions is asserted.

*The pair factors.* Take \(s_i=3,7,11,19,23\) and \(u_i=s_i^2\). Write \(J_{ij}=K^2+u_i u_j\). One can take infinitely many positive integers \(K\) such that the ten \(J_{ij}\) are pairwise coprime, odd, have only split prime factors, and each contains a distinguished prime to exact exponent one. Indeed the products \(u_i u_j\) are distinct. Require \(K\) even and divisible by every prime dividing a difference of such products, except the \(s_i\), and require \(K\equiv1\pmod{s_i}\) for all \(i\). Each \(J_{ij}\) is a sum of two coprime squares, and its odd prime factors are split. A common factor of two blocks would divide a difference, hence divide \(K\) or be one of the inert \(s_i\); both alternatives are impossible. Choose ten more distinct split primes larger than all products and assign one to each pair. At that prime choose a square root of \(-u_i u_j\) and a lift modulo its square that is not a root. CRT gives all the conditions on one progression. The same reasoning shows the blocks are coprime to every nonzero difference \(u_i-u_j\), by using the difference of the products with any common third index.

Start with the primitive columns \(U_i=(u_i,u_i^2-K^2)\). Exactly
\[
                  [U_i,U_j]=J_{ij}(u_j-u_i).
\]
Give every prime power in block \(J_{ij}\) the cut marking \(i,j\), and put \(N=\prod J_{ij}\). Each of the ten pair columns occurs with exponent one at some prime. The pair indicator vectors generate the even-coordinate-sum sublattice of \(\mathbb Z^5\): their differences give all coordinate differences and any pair vector has sum two. Adding the all-ones vector gives all of \(\mathbb Z^5\), since five is odd. Modulo the all-ones vector this proves index one for the four anchored differences. Complementing any cuts preserves this index.

*Balancing the actual heights.* Complement the five pair cuts of the cycle \(F=(01,12,23,34,40)\). Every resulting row then has five marked blocks and each pair has two common marked blocks. Put \(D=\prod_{e\in F}J_e\), \(A_i=\prod_{i\in e\in F}J_e\), so the new shared product is \(g_{ij}=D J_{ij}/(A_iA_j)\). Consider the integer row lattice
\[
       L=\{(a,b): a+(u_i+u_j)b\equiv0\pmod{J_{ij}}\ (ij\in F)\}.
\]
It has index \(D\) by CRT. In any basis matrix of this lattice all coordinates of the image of \(U_i\) are divisible by \(A_i\), since modulo \(J_{ij}\) that column is \(u_i(1,u_i+u_j)\). Dividing each image by its \(A_i\) gives primitive integer columns. To check this at a prime power \(p^e\Vert J_{ij}\) of the cycle, the row lattice locally is generated by \((p^e,0),(-(u_i+u_j),1)\) over \(\mathbb Z_{(p)}\). The endpoint columns have first coordinate \(u_i\) or \(u_j\) prime to \(p\), so their images have exact common depth \(e\). Other columns miss this kernel line modulo \(p\), since their determinants with an endpoint are coprime to \(p\) by the preceding block and difference conditions. Away from \(D\) the basis matrix is invertible locally. Its determinant has absolute value \(D\), so this change implements the complementary cuts and retains all ten absolute cofactors \(|u_i-u_j|\).

Here is a basis choice that also puts all five lengths on the required shells to relative error tending to zero. Use polynomials in \(z\) with coefficient triples in increasing order
\[
 f=(3566328745136157,-212247537740,-232553177),\qquad
 h=(-68417652775626,-13737560909,189877).
\]
They have nonzero resultant, and at each of the five cycle edges \(f(-u_i u_j)+(u_i+u_j)h(-u_i u_j)=0\). These six integer checks are in the script. Thus \((f(K^2),h(K^2))\in L\), and its ordinary coordinate gcd is bounded independently of \(K\), by the resultant identity. Dividing it down to a primitive vector \(v\) of the lattice divides by at most that gcd. Its length and the absolute value of its second coordinate are both \(\asymp K^4\). Complete \(v\) to a basis by a row \(w\), reducing its component parallel to \(v\) by a multiple of \(v\). The perpendicular component has length \(D/|v|\), so \(|w|=O(K^6)\). In the basis \((w+jv,v)\) choose integer \(j\) such that the first row's second entry is \(K^7+O(K^4)\). Its first entry is \(O(K^7)\). After applying this matrix to \(U_i\) and dividing by \(A_i=K^4(1+O(K^{-2}))\), each column has first coordinate \(-K^5(1+O(K^{-2}))\) and second coordinate \(O(K^2)\). At least one of three consecutive choices of \(j\) gives the first column a squared length different from \(n_0\): its first coordinate varies with nonzero slope, since the second entry of \(v\) has size \(\asymp K^4\) and \(u_0^2-K^2\sim-K^2\). All three choices have the stated asymptotics. Now \(n_i\sim K^{10}\), \(N\sim K^{20}\) and each minor has magnitude \(g_{ij}|u_i-u_j|\asymp K^4\). Taking sines of the pair angles proves the claimed two-sided angular estimate. Exact isotropy at all the marked primes would imply \(n_i\mid |V_i|^2\) for all \(i\) and hence equality for large \(K\), contrary to the choice of \(j\).

**PROVED — why this construction stops.** With \(m\) fixed and only pair cuts, an edge sees \(2(m-2)\) of the \(\binom m2\) equal-weight blocks. The required strict half-distance with a growing log margin is available for \(m=5,6\) and fails for \(m\ge7\). There is a six-mark example with the full minor divisibilities and the individually prescribed radii: add \(s_5=31\), use the uncomplemented pair cuts, and apply the unimodular shear \((x,y)\mapsto(x+K^{m-3}y,y)\) to the \(U_i\). The lengths are \(\sim K^{m-1}\sim\sqrt{n_i}\), all cofactors remain constant, and the angular width is \(\asymp K^{-2(m-2)}=N^{-4/15}(1+o(1))\) up to fixed multiplicative constants at six marks. The same shear gives the unbalanced five-mark example before the cycle change. In fact the mean weighted edge separation for any profile made entirely of two-against-the-rest cuts is exactly this same fraction of the parent weight, even with unequal block sizes. Each individual cut separates exactly the stated number of pairs. Thus reweighting these blocks cannot produce a seventh normalized mark.

**COMPUTED.** Exact checks at one explicit CRT representative each verify primitivity, all ten or fifteen divisibilities, the shortened arc bound, and squared-norm relative errors below \(10^{-12}\). For the balanced five, the fourth powers of the lengths also agree with \(N\) to relative error below \(10^{-12}\). On the uncomplemented profiles the congruence lattice solver independently recovers every ordered quartet ray; all 120 and 720 ordered five-mark cocycles, respectively, hold exactly. This supplies a passing control for III.8 and exhibits the gap between approximate and exact individual norms. [Full integer pair construction and chart checks](iii_pair_pencil.py)

## Twelve-line summary of Part III

1. **PROVED.** Small integer Plücker relations give a four-mark test using three coprime matching products from the same parent.
2. **PROVED.** In any sufficiently large normalized profile the exact separation divisors determine every realizable rational star direction up to one projectivity.
3. **PROVED.** This reconstruction survives a missing tenth of the pair constraints at every mark.
4. **PROVED.** On regular five-sets, three independently solved quartet congruences telescope to an exact cross-ratio identity.
5. **PROVED.** If all regular quartets pass independently, their ratios already glue; the individual integral cofactor factorizations still matter.
6. **PROVED.** On fixed narrow Paley parents at most one ordering of the two primes in every preassigned twin pair can pass simultaneously.
7. **COMPUTED.** A planted fourteen-prime parent admits a factorizable passing quartet, but no five marks pass under any assignment.
8. **PROVED.** Among parents formed from \(r\) random split primes in \([r^D,2r^D]\), \(D\ge4\), all but an exponentially small fraction have at most 17 points on a fixed critical arc.
9. **PROVED.** A rank-\(s\) coset accounts for at most \(8s+4\) extensions of one edge in the triangle divisor-split equation.
10. **PROVED.** At five marks a primitive integer matrix meets all minor divisibilities at exponent 10/3; at growing sizes I can retain half the log divisibility on every pair.
11. **EVIDENCE.** Exceptional parents still decide the conjecture. The full integral models here have fixed size, and the worst-case growth bound is unchanged.
12. **HEURISTIC.** Please referee the five-mark cancellation, the parent count and the balanced integer countermodel; then test simultaneous individual cofactor scales on a larger exceptional profile.

# Part IV — Can a family of integral matrices produce a growth-scale fake?

**PROVED — outcome of the construction attempt.** I tried to construct (C) by specializing one-parameter families of primitive columns with shared polynomial factors in their minors. Over finite fields of characteristic two this construction attains the growth scale with **every cofactor constant**. Over \(\mathbb C(t)\), **every full construction satisfying the genuine degree bound has at most 24 marks**. This applies to arbitrary cuts, prime powers, and arbitrary projective lines at the marked primes. It also quantifies what has to fail at a single specialization. With 48 marks, at least one raw minor exceeds its critical *degree* budget by at least \(\deg N/136\). A construction with degrees that predict the sizes of its specialized blocks must make that minor unusually small at the chosen integer parameter.

**HEURISTIC — verdict on the question asked.** I do not know whether integral fakes exist over the integers at growth scale. The theorem closes a construction, including a growing degree-controlled version of III.11. It gives no improvement to the coefficient three on arbitrary integer parents. The unproved step for (P) would be to apply a comparable intersection inequality to prime contacts of a *single arithmetic configuration*. The proof below uses an actual surface over \(\mathbb C\). I have no arithmetic substitute for it.

## The construction I tried

**PROVED.** Put \(n_i=\prod p^{a_{ip}}\), \(g_{ij}=\gcd(n_i,n_j)\), and \(D_{ij}=[P_i,P_j]\). The exact arithmetic requirements relevant here are
\[
        g_{ij}\mid D_{ij},\qquad 0<|D_{ij}|\le C_1\sqrt{n_i n_j}N^{-1/4}.                 \tag{IV.1}
\]
The second inequality is the genuine cofactor bound multiplied by \(g_{ij}\). Without a norm condition the cone itself can always be obtained after (IV.1): take an integral unimodular chart in which no second coordinate is zero, choose column signs with positive second coordinate, and apply the common shear \((x,y)\mapsto(x+Hy,y)\) with \(H\to+\infty\). All absolute minors and all column gcds are preserved. Locally, (IV.1) forces the rows marked at a given prime depth onto one projective residue line. The unmarked rows need not occupy a second line.

**PROVED — a seed that attains the requested scale in finite characteristic.** Let \(q=2^d\), \(d\ge3\), and use the \(m=q-1\) nonzero marks \(a\in\mathbb F_q\) and the \(q-1\) linear primes \(t-x\) for \(x\in\mathbb F_q^*\). Mark \(x\) when \(\operatorname{Tr}(ax)=0\). Write \(L_a(t)=a^{-q/2}\sum_{j=0}^{d-1}(at)^{2^j}\), a monic polynomial, and use the primitive columns
\[
                 P_a=(f_a,1),\qquad f_a=(t^{q/2}-L_a(t))/t .
\]
At a marked \(x\) one has \(f_a(x)=x^{q/2-1}\). The two distinct trace hyperplanes for \(a,b\) have exactly \(q/4-1\) common nonzero elements. The difference \(f_a-f_b\) has exactly this degree: its leading coefficient is nonzero by injectivity of \(a\mapsto a^{-q/4}\). Thus its quotient by the product of **all** common marked primes is a nonzero constant. Pair separation is \(q/2>(q-1)/2\), and the log norm of the parent is \((q-1)\log q=m\log(m+1)\). The degree budget for the cofactor is \(1/4\), so constant cofactors meet it. The characteristics of the coefficient field and of the local residue fields are essential here.

**COMPUTED.** At \(q=128\) all 8,001 pairs of the 127 marks satisfy the polynomial factorization exactly. The check multiplies the required linear factors for each pair and compares coefficients. [Trace construction certificate](iv_trace_check.py)

**PROVED — why simply assigning integer primes to this seed fails.** For *any* assignment of split primes to its cut columns, \(d_{ab}\) depends only on \(a+b\), even if the prime weights are unequal. Suppose the resulting array meets pair separation and admits an integer realization at the normalized bound of III.7. For each translation \(u\), restrict to \(a\ne0,u\). The original and translated realizations have identical labeled separation divisors on this restriction. For sufficiently large \(q\), III.7 gives one rational projectivity \(T_u\) taking each of those directions to the direction labeled \(a+u\). Composing on their common domain gives \(T_uT_v=T_{u+v}\), since three distinct directions determine a projectivity. The action is faithful. In \(\mathrm{PGL}_2(\mathbb C)\) an elementary abelian 2-group has order at most four: after putting a nonidentity involution at \(z\mapsto-z\), its centralizer consists of \(z\mapsto cz\) and \(z\mapsto c/z\); there are at most two scalings of order dividing two, and these have index at most two. This excludes that incidence table at *every* integer parent to which III.7 applies. The same rigidity proof works with any fixed multiplier \(C_1\) in (IV.1) under strict pair separation: the local inequality in III.7 has a positive multiple of \(W\) in its log margin and absorbs the fixed multiplier. To try to construct (C), one must at least break this translation symmetry. The obstruction that follows allows **every** incidence table.

## The obstruction for arbitrary polynomial cuts

**PROVED using Langer's orbifold Bogomolov–Miyaoka–Yau inequality.** Take \(N,n_i\in\mathbb C[t]\), with \(L=\deg N>0\), \(n_i\mid N\), and \(\alpha_i=\deg n_i\). Let \(P_i\in\mathbb C[t]^2\) have relatively prime entries and pairwise nonzero minors \(D_{ij}\). For each root of \(N\) and each depth up to its multiplicity, count one layer \(h\); let \(s_h\) be the number of rows marked at that depth. Thus \(\sum_h1=L\), \(\sum_hs_h=\sum_i\alpha_i\). Define the missing degree and the allowed excess over the critical size by
\[
 E_{ij}=\deg\gcd(n_i,n_j)-\deg\gcd(n_i,n_j,D_{ij}),\qquad
 \deg D_{ij}\le \frac{\alpha_i+\alpha_j}{2}-\frac L4+\delta
 \quad(i\ne j),\qquad \delta\ge0 .
\]
For every \(m\ge4\),
\[
 \boxed{\quad
    \frac34(m-24)L+\sum_h(s_h-m/2)^2
       \le m(m+3)\delta+4\sum_{i<j}E_{ij}.
       \quad}                                                        \tag{IV.2}
\]
Full divisibility and the genuine degree bound give \(E_{ij}=\delta=0\), hence \(m\le24\). The proof requires no pair separation, prescribed norms, restriction on the roots, or nonsingularity of the contact pattern.

**PROVED — choosing the surface at infinity.** Regard the columns as sections \(C_i\) of a ruled surface over \(\mathbb P^1\), starting with the product. At finite \(t\), the intersection multiplicity of two sections is exactly the vanishing order of their minor, because their coordinate pairs are primitive. Write \(S=\sum_i C_i^2\) and \(I=\sum_{i<j} C_i\cdot C_j=(m-1)S/2\). Indeed the numerical classes of two sections differ by a multiple of a fiber, so \(2C_i\cdot C_j=C_i^2+C_j^2\).

At infinity, if a point of the fiber contains \(k>m/2\) sections, blow it up and contract the strict transform of the fiber. The \(k\) self-intersections go down by one and the others go up by one. Thus \(S\) decreases by the positive integer \(2k-m\). The identity for \(I\) gives \(S\ge0\), so this terminates. Order the marks in contiguous blocks by their final value on the infinity fiber. Pair every index with the index \(\lfloor m/2\rfloor\) further on in cyclic order. No such pair meets at infinity. For these pairs the entire intersection number is \(\deg D_{ij}\). Summing the degree bounds over this permutation gives
\[
                S\le\sum_h s_h-mL/4+m\delta.                           \tag{IV.3}
\]
This normalization charges cancellations at infinity before any use of the finite divisibilities.

**PROVED — the surface inequality.** Give the sections weight \(a=2/3\). Blow up every collision of \(k\ge3\) sections until they separate into smaller groups, giving each new exceptional curve weight zero. Once only two sections meet, leave a transverse intersection alone. For a remaining tangency of order \(e\ge2\), blow up its full chain of \(e\) centers, including the final separation, assigning all those exceptionals weight \(a\). This also works for a pair separating from a larger group. In a local fiber parameter \(z\) the sections are graphs, and each blowup in their common direction replaces their coordinate \(y\) by \((y-y(0))/z\). A later center meets only its immediate predecessor among the exceptional curves. The final positive-weight boundary has simple normal crossings. Add four generic fibers of weight \(1/2\), avoiding all contacts and blowups. Denote this weighted boundary on the blown-up surface \(Y\) by \(B\).

The divisor \(K_Y+B\) is \(\mathbb Q\)-effective. On the initial ruled surface, \(K_{Y/\mathbb P^1}+C_i+C_j\) is the pullback of the base intersection divisor of the two sections: it has degree zero on fibers, and restriction to \(C_i\) identifies the base line bundle by adjunction. On blowing up a center through \(l\) of those two sections, the discrepancy with their strict transforms is \(1-l\). A vertical multiplicity can decrease only along their shared chain, whose length is at most their initial intersection multiplicity there. This leaves an effective representative after all the blowups. Averaging over pairs proves effectivity with horizontal coefficient \(2/m\). The four half-fibers cancel the base canonical class, and all additional coefficients in \(B\) are nonnegative. Langer's inequality now gives
\[
                   0\le 3e_{\rm orb}(Y,B)-(K_Y+B)^2 .
\]
Here the normal-crossing orbifold Euler number is the ordinary Euler number minus the weighted Euler numbers of the boundary curves, plus the products of coefficients at their crossings. The external result used is A. Langer, *Logarithmic orbifold Euler numbers of surfaces with applications*, Proc. London Math. Soc. 86 (2003), 358–396, Theorem 0.1. The record tested orbifold BMY at seven marks; the estimate here uses growing collision groups and different weights on their exceptional curves. [Earlier seven-mark calculation](/home/user/Jarnik/docs/claude_review/round4/L7.md)

For clarity, the following is the complete accounting. Before blowups the canonical square with section weights included is
\[
          8-4am+(a^2-2a)S+2a^2 I.
\]
A blowup centered on \(k\) sections, with old vertical coefficient \(b_0\) and new coefficient \(b\), subtracts \((1+b-b_0-ak)^2\) from this square. Each such center consumes \(\binom{k}{2}\) pair intersections. In the orbifold Euler number it contributes \(1-2b\) from the surface and new vertical curve, \(bb_0\) from the vertical crossing, and \(ak(b-b_0)\) from moving the horizontal incidences onto the new exceptional. The last contributions telescope down each chain. After including the consumed intersections, its contribution to the difference in the preceding display for BMY is
\[
 F(k,b_0,b)=3(1-2b+bb_0)+3ak(b-b_0)+(1+b-b_0-ak)^2-2a^2\binom{k}{2}.
\]
Every unblown transverse pair contributes \(a^2=4/9\). The remaining terms before adding the generic fibers are \(4-2am+(2a-a^2)S\); each generic half-fiber adds \(-1+am/2\), leaving just \(8S/9\). A center of size \(k\ge3\) has \(b_0=b=0\) and contributes \(4-8k/9\). The first blowup of a two-section tangency contributes \(F(2,0,a)=8/9\), and every later one contributes \(F(2,a,a)=-4/9\). Its total is at most \(4e/9\) for \(e\ge2\).

It follows that one can charge \(4/9\) per unit of actual pair intersection everywhere, and deduct at least \((4/9)\phi(k)\) at an actual contact group of size \(k\), where
\[
             \phi(k)=\max\{0,\tbinom k2+2k-9\}.
\]
For \(k\ge3\) this is exactly the preceding large-center calculation. On a chain of pairs \(\phi(2)=0\), and the bound holds for the whole chain. Intersections at infinity are included in \(I\).

At a prescribed layer of size \(s\), partition its marked rows by their *actual* projective residue class at that depth, with sizes \(s_1,s_2,\ldots\). Write \(E_h=\binom{s}{2}-\sum_j\binom{s_j}{2}\), so \(\sum_h E_h=\sum_{i<j}E_{ij}\). Each group of size at least two belongs to its actual contact group at that depth. Distinct groups at a depth use distinct centers, and distinct depths count separate levels. For positive \(u,v\),
\[
              \phi(u+v)-\phi(u)-\phi(v)\le2uv .
\]
When \(u,v\ge3\) the difference is \(uv+9\le2uv\); when the smaller is one or two, direct substitution gives the inequality (for the other at least three the differences are \(v+2\) and \(2v+5\)). Thus \(\sum_j\phi(s_j)\ge\binom{s}{2}+2s-9-2E_h\), also for an empty layer. Since \(\phi\) is nondecreasing, the surface inequality implies
\[
  0\le \frac89 S+\frac49\!\left(I-\sum_h\binom{s_h}{2}\right)
       +\sum_h(4-8s_h/9)+\frac89\sum_{i<j}E_{ij}.
\]
Substituting \(I=(m-1)S/2\) and (IV.3), the contribution per layer is at most
\(4-m/6-(2/9)(s_h-m/2)^2\), with the additional total \(2m(m+3)\delta/9+8\sum E_{ij}/9\). This proves (IV.2).

**COMPUTED.** The accompanying exact-rational check reproduces the blowup contributions and the completed square. As a check on the infinity terms, it also computes both Chern numbers directly for the graphs \(y=u_i-t^2/u_i\) from the pair construction: each finite pair meets twice, and all sections share two successive centers at infinity. [Intersection check](iv_section_check.py)

## Where this construction stops on the integer problem

**PROVED — fixed families and isolated specializations.** Suppose fixed integer polynomials supply the parent, the row factors and the columns at infinitely many unbounded integer parameters, with a nonconstant parent. Then the integer inequalities in (IV.1) give the degree inequalities with \(\delta=0\), by taking logarithms and dividing by \(\log|t|\). Full divisibility also passes to the polynomial model. Indeed a primitive polynomial gcd of two row factors divides both in \(\mathbb Z[t]\). If \(G(T)\mid F(T)\) on an unbounded integer sequence, divide \(F\) by \(G\) over \(\mathbb Q[t]\) and clear a fixed denominator. The resulting integral remainder has size \(o(|G(T)|)\) and is divisible by \(G(T)\) on that sequence, hence vanishes identically. The same argument applies to \(n_i\mid N\). Common polynomial factors of an input column can first be removed. Subsequent gcds of its two coprime polynomial entries at integers divide a fixed nonzero Bézout constant, so primitive normalization changes no degrees. Thus no full family of this type can have 25 marks. In particular, taking the parameter increasingly large *for each fixed shape* cannot construct (C).

For a single specialization there is a quantitative test. With full polynomial divisibilities, apply (IV.2) to **any 48 marks**, keeping the whole parent. Some edge of that subset satisfies
\[
       \deg D_{ij}\ge (\alpha_i+\alpha_j)/2-L/4+L/136 .
\]
Here \(3(48-24)/(4\cdot48\cdot51)=1/136\). Consequently, for large \(m\), at least a fraction \(1/47-o(1)\) of *all* pairs exceed the critical degree by this amount: the complementary graph is \(K_{48}\)-free, so Turán's bound applies. This is a deterministic statement inside each presentation. Suppose at the proposed integer instance \(W=\lambda L+o(W)\) and \(\log n_i=\lambda\alpha_i+o(W)\) uniformly in the marks, for one positive conversion factor \(\lambda\). On every one of those pairs an actual fake then requires
\[
                \lambda\deg D_{ij}-\log|[P_i,P_j]|
                            \ge (1/136-o(1))W .                         \tag{IV.4}
\]
Thus even predicting all minor logs to \(o(W)\) from their degrees would be too accurate for this construction to succeed. The definition of an integral fake supplies neither that accuracy nor a polynomial presentation.

**PROVED — allowing missing factors does not give a limit construction.** Impose the critical raw-minor degree bound \(\delta=0\), but permit partial divisibility. For \(m\ge36\), let \(\overline E\) be the mean missing degree over all pairs. Apply (IV.2) to 36 uniformly selected marks of this *same* configuration. At a layer of global marked fraction \(x\), the sample size \(s'\) of the marked group has mean \(36x\) and variance \(36x(1-x)(m-36)/(m-1)\), by sampling without replacement. Hence \(\mathbb E(s'-18)^2\ge9(m-36)/(m-1)\), with the minimum at \(x=1/2\). We obtain
\[
        \overline E\ \ge\ \left(\frac1{140}-\frac1{8(m-1)}\right)L .
\]
More generally the first constant from a sample of \(k\) is \((k-18)/(2k(k-1))\), maximized at \(35,36\). If the profile has pair separation at least \(L/2\), summing the separated pairs over layers gives \(\sum_h(s_h-m/2)^2\le mL/4\). It follows by Cauchy–Schwarz that the mean prescribed common degree is \(L/4+O(L/\sqrt m)\). Accordingly these polynomial columns can retain **at most \(34/35+o(1)\)** of the total prescribed common degree at the critical size. The optimal attainable fraction is unknown; \(34/35\) is the upper bound proved here.

**HEURISTIC — the missing arithmetic input.** The integer data give sections of a ruled *arithmetic* surface. Their contact at \(p\) contributes \(\log p\); it supplies no point of a complex base curve at which the corresponding graphs meet. Replacing \(p\) by a polynomial of comparable size preserves products but need not preserve the sums in a two-by-two determinant. This is exactly the information used in (IV.2). I have stopped at this obstruction. The finite-characteristic construction, the symmetry exclusion, and (IV.4) together leave open isolated characteristic-zero configurations with unsymmetric cuts and large specialization discrepancies. There is no implication here that exact Gaussian norms are necessary to prove growth.

**PROVED — correction to Part III.** The exposure in III.9 uses a uniformly random ordering of the sampled primes. Sorting them by size would invalidate the conditional subset probability. No exposure argument is used in Part IV.

## Twelve-line summary of Part IV

1. **HEURISTIC.** I leave (P) and (C) open for integer parents. My concentrated attempt was to construct (C) by specializing a one-parameter family.
2. **PROVED.** The no-norm problem retains the full minor divisibility and size conditions; a common unimodular shear can supply the cone.
3. **PROVED.** Trace hyperplanes over \(\mathbb F_q(t)\) give \(q-1\) marks, weight \((q-1)\log q\), and constant nonzero cofactors on every pair.
4. **COMPUTED.** At 127 marks, all 8,001 full polynomial divisibilities passed exact verification.
5. **PROVED.** The trace incidence table cannot itself lift to a large integer fake at any parent: its translations would act faithfully on a rational projective line.
6. **PROVED.** Allowing arbitrary asymmetric cuts, a full characteristic-zero polynomial model at the genuine degree bound still has at most **24 marks**.
7. **PROVED.** The obstruction comes from orbifold BMY on the surface of the column directions, with prime powers, tangencies and infinity contacts included.
8. **PROVED.** On any 48 marks with full polynomial divisibility, one raw minor overshoots the critical degree by at least \(\deg N/136\).
9. **PROVED.** A single successful specialization with accurately predicted block logs needs a shortfall of at least \((1/136-o(1))W\) in the predicted minor logs on a positive fraction of pairs.
10. **PROVED.** At critical degree, a growing polynomial construction retaining partial factors can keep at most \(34/35+o(1)\) of the prescribed shared degree.
11. **HEURISTIC.** An isolated integer parent need not admit such a lift. The worst-parent growth bound and its leading constant are unchanged.
12. **HEURISTIC.** Please referee the surface inequality (IV.2), particularly the exceptional-curve weights. I would reopen this construction only with a concrete arithmetic source of the discrepancies in (IV.4).

# Part V. An unbounded circle cluster over the fixed field \(\mathbb F_2\)

**PROVED.** For the anisotropic norm from \(\mathbb F_4(t)\) to \(\mathbb F_2(t)\), the answer to (a) is yes. For each \(m=2^r\), \(r\ge2\), there are \(m\) distinct points of degree
\[
             D=m(m-1)(\log_2 m-1)
\]
on one parent of degree \(2D\), with maximum pair-difference degree exactly \(D/2\). The field of constants has four elements upstairs and two downstairs for every \(m\). In coordinates \(z=A+B\theta\), where \(\theta^2+\theta+1=0\), the norm is \(A^2+AB+B^2\). It is anisotropic because it is the norm of a field extension. Part IV's characteristic-zero uniform bound, including the application to genuine circles that I missed, therefore has a different answer over this fixed field.

**PROVED (scope of the degree bound).** The same construction has a quadratic cost even if one chooses the evaluation polynomial, chooses only a subfamily of hyperplanes, reverses selected sides, and divides out their entire common polynomial factor. For \(s\) retained marks the reduced parent degree is \(\Omega(s^2\log s)\) as \(s\) grows. The unrestricted growth-scale question (b) remains open here. The proof below locates the loss at the finite primes and at infinity on the same parent.

## V.1. A common norm and an exact difference

**PROVED.** Put \(d=r+1\), \(Q=2^d=2m\), and \(V=\mathbb F_2[t]_{<d}\). For a binary linear subspace \(H\subset V\), write \(P_H(X)=\prod_{v\in H}(X+v)\). This polynomial is additive in \(X\). Indeed \(P_{\{0\}}(X)=X\) and adjoining one generator \(b\notin H\) gives
\[
      P_{H+\langle b\rangle}(X)
          =P_H(X)P_H(X+b)=P_H(X)^2+P_H(b)P_H(X).
\]
The recurrence proves additivity by induction over \(\mathbb F_2(t)\).

Take the \(m\) linear functionals \(\ell:V\to\mathbb F_2\) with \(\ell(1)=1\), and write \(H=\ker\ell\). Set \(Y=t^L+\theta\), initially with \(L\ge d\). For each \(H\) take the polynomial
\[
       z_H=\frac{P_H(Y)}{Y}
           =\prod_{v\in H\setminus\{0\}}(Y+v).
\]
Conjugation fixes \(v\) and sends \(Y\) to \(Y+1\). Since \(1+H=V\setminus H\), the common norm and its degree are
\[
       N=\prod_{v\in V\setminus\{0,1\}}(Y+v)\in\mathbb F_2[t],
       \qquad \deg N=(Q-2)L,\qquad D=(m-1)L.
\]
This uses the model allowing prime powers.

For two distinct marks let \(J=H\cap H'\), so \(|J|=Q/4\). Write \(H=J+\langle u\rangle\) and \(H'=J+\langle u'\rangle\). In \(V/J\) the class of \(1\) is the nonzero class in neither hyperplane, hence \(u+u'\equiv1\bmod J\). Subtracting the two recurrences gives
\[
       z_H-z_{H'}=P_J(1)\frac{P_J(Y)}{Y},\qquad
       \deg(z_H-z_{H'})
          = (Q/4-1)L+\sum_{v\in J}\deg(1+v).                 \tag{V.1}
\]
The factors on the right are nonzero. In particular the \(m\) marks remain distinct after evaluation.

**PROVED.** The maximum of the sum in (V.1) is \((d-2)2^{d-2}\). Choose a basis of \(J\) with distinct leading degrees \(a_1<\cdots<a_{d-2}\). Then \(a_i\le i+1\). Exactly \(2^{i-1}\) elements have degree \(a_i\), and adding \(1\) preserves all their degrees since \(1\notin J\). The element \(v=0\) contributes zero. Thus the sum is at most \(\sum_{i=1}^{d-2}2^{i-1}(i+1)=(d-2)2^{d-2}\). Equality holds for \(J=t^2\mathbb F_2[t]_{<d-2}\), as seen by taking the two forms extracting the constant coefficient and the sum of the constant and linear coefficients. By (V.1), the half-degree condition is exactly \(2\deg P_J(1)\le L\) for every pair. The least \(L\ge d\) for these marks is therefore
\[
              L=(d-2)2^{d-1}=m(\log_2m-1).                  \tag{V.2}
\]
The extremizing pair attains \(D/2\). Equations (V.1) and (V.2) prove the claimed family.

**COMPUTED.** The accompanying script constructs every mark by multiplication in \(\mathbb F_4[t]\), constructs the parent separately, checks every norm, and subtracts every pair. The scale \(L\) in the script is obtained by enumerating the intersections. The values reproduce (V.2).

| \(m\) | \(L\) | Point degree \(D\) | Parent degree | Maximum difference degree | Pairs checked |
|---:|---:|---:|---:|---:|---:|
| 4 | 4 | 12 | 24 | 6 | 6 |
| 8 | 16 | 112 | 224 | 56 | 28 |
| 16 | 48 | 720 | 1,440 | 360 | 120 |
| 32 | 128 | 3,968 | 7,936 | 1,984 | 496 |

## V.2. Can a better evaluation or more cancellation reach growth scale?

**PROVED.** Fix the same ambient \(V=\mathbb F_2[t]_{<d}\), select \(s\ge2\) distinct forms with \(\ell(1)=1\), and choose either affine side \(S_i\) of each. For any \(T\in\mathbb F_2[t]\) let
\[
       Y=T+\theta,\qquad F_i=\prod_{v\in S_i}(Y+v).
\]
These raw marks all have norm \(P_V(Y)\), hence all have the same degree \(D_0\). If a common polynomial division yields a \(c=0\) cluster, then for \(d\ge48,\ s\ge6\) its **point** degree is at least \(Qs d/16\). In particular the quadratic order in V.1 is sharp for this family. Both the restriction to this \(V\) and the requirement that the divisor be common to the chosen marks are part of the statement.

*Proof.* Let \(G=\gcd_i F_i\) be monic and \(g=\deg G\). Dividing by a common factor of degree \(a\) makes the half-degree condition equivalent to \(2\deg(F_i-F_j)-D_0\le a\) for every pair. Thus the full gcd is the most favorable choice both for contact and for the reduced degree. Put \(B=\max(d,\deg T)\), with \(\deg 0=-\infty\), and \(C=\bigcap_i S_i\), \(h=|C|\). If the selected forms have span of rank \(r'\), then \(s\le2^{r'-1}\), since evaluation at \(1\) is a nonzero linear functional on that span. A consistent system of the selected affine equations has \(Q/2^{r'}\) solutions. Hence \(h\le Q/(2s)\), also when \(C\) is empty.

First, the pair differences force
\[
                        g\ \ge\ Q(d/2-6).                    \tag{V.3}
\]
For any pair choose \(a\in S_i\cap S_j\), which exists because their forms are independent. With \(J\) the common kernel of the two forms, translating by \(a\) gives the raw version of (V.1),
\[
                 F_i-F_j=P_J(1)P_J(Y+a).
\]
The \(Q/4\) elements of \(1+J\) are distinct and nonzero. For \(k=1,\ldots,d-2\) at most \(2^k\) have degree \(<k\). Summing the numbers of elements of degree at least \(k\) gives \(\deg P_J(1)\ge(Q/4)(d-4)\). Also,
\[
                    \sum_{v\in V}\bigl(B-\deg(Y+v)\bigr)\le2Q.
\]
For \(\deg T\ge d\) the sum is zero. Otherwise \(B=d\), and at each \(k=1,\ldots,d\) at most \(2^k\) choices of \(v\) satisfy \(\deg(Y+v)<k\). No \(Y+v\) vanishes. Thus \(\deg P_J(Y+a)\ge QB/4-2Q\), while \(D_0\le QB/2\). Applying the necessary inequality \(g\ge2\deg(F_i-F_j)-D_0\) proves (V.3).

Second, the common factor satisfies
\[
                          g\ \le\ hB+Qd/4.                   \tag{V.4}
\]
Remove the explicit common product \(\prod_{v\in C}(Y+v)\), of degree at most \(hB\). Consider an irreducible \(f\in\mathbb F_4[t]\) in the gcd of the remaining products. Its contraction is generated by a prime \(P\in\mathbb F_2[t]\) of even degree \(e\). Indeed at a root \(\alpha\) of \(f\), the identity \(T(\alpha)+v(\alpha)=\theta\) implies \(\mathbb F_4\subset\mathbb F_2(\alpha)\). Thus \(P=f\bar f\), and for a polynomial over \(\mathbb F_2\), divisibility by \(f^j\) is equivalent to divisibility by \(P^j\).

Among \(v\in V\setminus C\) choose \(v_*\) maximizing \(\operatorname{ord}_f(Y+v)\), and select one side \(S_i\) excluding \(v_*\). For each \(j\) through that maximum, the set of \(v\in V\) with \(f^j\mid Y+v\) is exactly the affine space \(v_*+\{w\in V:P^j\mid w\}\). For \(ej<d\) it has \(Q/2^{ej}\) elements. At most half lie in \(S_i\), since an affine side omitting its center either misses the whole space or meets it in half. For \(ej\ge d\) it is a singleton missed by \(S_i\). Above the maximum no surviving factor contributes. This *one* selected product bounds the valuation of the gcd across all powers of \(f\), by \((Q/2)\sum_{ej<d,\ j\ge1}2^{-ej}\).

For each base prime of even degree \(e\), the two conjugate primes upstairs have degrees summing to \(e\). If \(\pi_2(e)\) counts irreducibles of degree \(e\), then \(e\pi_2(e)\le2^e\) by counting roots. The total residual gcd degree is at most
\[
 \frac Q2\sum_{\substack{2\le e<d\\e\ {\rm even}}}
                   \frac{1}{1-2^{-e}}
 \ \le\ \frac{Qd}{4}.
\]
For the last inequality there are at most \((d-1)/2\) summands and their excess over \(1\) sums to at most \(\sum_{k\ge1}1/(4^k-1)\le4/9\). This proves (V.4).

Combining (V.3) and (V.4) gives \(hB\ge Q(d/4-6)\), forcing \(B\ge s(d/2-12)\ge sd/4\) when \(d\ge48\). For \(s\ge6\), \(B>d\), so every raw factor has degree \(B\), and \(D_0=QB/2\). Even after full division,
\[
       D_0-g\ \ge\ (Q/2-h)B-Qd/4
            \ \ge\ Qd(s-3)/8\ \ge\ Qsd/16.
\]
Finally \(Q\ge2s\) implies \(d\ge 1+\log_2s\), so every sequence in this family with \(s\to\infty\) eventually satisfies the thresholds and has parent degree \(\Omega(s^2\log s)\). \(\square\)

**COMPUTED.** As a separate check of (V.4), at \(d=4\) the script enumerates every selection of at least two distinct forms with either choice of side, for every \(T\) of degree \(<6\), including zero. Translation by \(V\) permutes the sides, so four representatives of \(T\bmod V\) suffice. In 26,176 cases the maximum degree of the gcd after removing the intersection product was **4**. This uses Euclidean polynomial gcds. It also equals the bound obtained before simplifying the prime sum, since the only eligible base prime is \(t^2+t+1\).

**PROVED (odd-characteristic boundary).** The norm identity in V.1 uses a nonzero translation by a base-field element under coefficient conjugation. If \(\bar Y=Y+w\) with \(\bar w=w\), applying conjugation twice gives \(2w=0\). This implementation therefore requires characteristic two. The example rules out uniformity arguments valid unchanged for every fixed finite field and every anisotropic quadratic norm. Its degree scale is compatible with a bound of the form in (c).

## Twelve-line summary of Part V

1. **PROVED.** The answer to (a) is yes for the fixed anisotropic norm over \(\mathbb F_2(t)\).
2. **PROVED.** Every \(m=2^r\), \(r\ge2\), occurs on a parent of degree \(2m(m-1)(r-1)\).
3. **PROVED.** The points have degree \(D=m(m-1)(r-1)\), and some pair attains the exact boundary \(D/2\).
4. **PROVED.** Each mark chooses a hyperplane in \(\mathbb F_2[t]_{<r+1}\), and coefficient conjugation chooses its complementary side.
5. **PROVED.** Additivity of the subspace polynomial gives the full pair difference in (V.1), including its coefficient at infinity.
6. **COMPUTED.** Direct arithmetic verifies every norm and all 650 pair differences for \(m=4,8,16,32\).
7. **PROVED.** For \(s\) chosen sides, the half-degree condition forces common cancellation of degree at least \(Q(d/2-6)\).
8. **PROVED.** Outside their intersection the chosen products have common gcd degree at most \(Qd/4\), for every evaluation \(T\).
9. **PROVED.** These bounds force reduced parent degree \(\Omega(s^2\log s)\) even after optimizing the subfamily and common divisor.
10. **PROVED.** The conjugation-as-translation identity requires characteristic two.
11. **HEURISTIC (research status).** Unboundedness for the ternary circle, growth scale for unrestricted fixed-field parents, and the integer growth rate remain unresolved here.
12. **HEURISTIC (request).** Please have the split-prime valuation bound (V.4) checked independently. For (b) I would test a different additive space \(V\), for which the residue fibers can have different sizes.

# Part VI. Fixed fields, reflection, and growth

**PROVED.** The answer to (b) is **unbounded** over the fixed field \(\mathbb F_3(t)\), for the genuine form \(A^2+B^2\). Each point descends individually to \(\mathbb F_9[t]\). The construction uses a quadratic character on an auxiliary set of roots, multiplicities \(0,1,2\), and a common divisor large enough to reach the critical width. The answer to (c) is also **yes**. Over \(\mathbb F_2(t)\) there are squarefree parents of degree \(O(M^2\log M)\) with \(M\) marks.

**HEURISTIC (research verdict).** I could not decide (a). My strongest general test uses all binary linear combinations of the valuation rows of a single parent. It excludes the record's squarefree skew double core in the binary polynomial model, independently of the prime assignment. A second test gives a quadratic cost for complete Hadamard cores, and a local argument rules out near-linear degree for full hyperplane designs even with nonconstant translation. The gap between \(M\log M\) and \(M^2\log M\) for unrestricted binary parents survives. None of these growth arguments improves the integer bound.

## VI.1. A necessary code condition on every binary parent

**PROVED.** Write \(W=\deg N=2D>0\), \(K=\lceil D/2\rceil\), and \(U_k=1+t^{-k}\mathbb F_4[[1/t]]\). At each split prime \(p=\pi_p\bar\pi_p\) use monic factors and let \(a_{ip}=\operatorname{ord}_{\pi_p}z_i\). Fix an anchor \(i=0\). Of the columns \((a_{ip}-a_{0p})_i\bmod 2\), keep the distinct nonzero patterns. There are \(n\) such patterns. Give a pattern weight equal to the sum of the degrees of its distinct base primes, and form the binary matrix with one column per pattern. If its row span is \(C\) of dimension \(r\), then **every nonzero word of \(C\) has weight at least \(K+1>W/4\)**. In particular, as \(n\to\infty\),
\[
          r\le [1-H_2(1/8)]n+O(\log n),\qquad
          1-H_2(1/8)<0.457,                                      \tag{VI.1}
\]
where \(H_2(x)=-x\log_2 x-(1-x)\log_2(1-x)\). Arbitrary prime powers and unequal weights are allowed.

*Proof.* Set \(R_i=z_i/z_0\). The arc condition says \(R_i\in U_K\). Every product \(R\) of these ratios is still in \(U_K\), so its logarithmic derivative \(R'/R\), taken with respect to \(t\), has order at infinity at least \(K+1\), including the case where it vanishes. In characteristic two,
\[
      \frac{(\pi_p/\bar\pi_p)'}{\pi_p/\bar\pi_p}
          =\frac{\pi_p'}{\pi_p}+\frac{\bar\pi_p'}{\bar\pi_p}
          =\frac{p'}p .
\]
Inert primes have the same valuation on all marks. Thus a combination with binary word \(c\) gives \(R'/R=H'/H\), where \(H\) is the squarefree product of all the base primes in the coordinates of \(c\) equal to \(1\). For \(c\ne0\), \(H'\ne0\), and the order at infinity is \(\deg H-\deg H'\le\deg H\). This proves the weighted distance assertion for the entire span.

The total weight of all coordinates is at most \(W\). Balls of weight radius \(W/8\) about the codewords are disjoint. For \(j=\lfloor n/8\rfloor-1\) and \(n\) sufficiently large, a uniform \(j\)-subset has expected weight at most \(jW/n\). By Markov's inequality at least a fraction \(8/n\) of these subsets have weight at most \(W/8\). Consequently \(2^r(8/n)\binom nj\le2^n\). Taking logarithms and using Stirling's estimate gives (VI.1). \(\square\)

**PROVED (an application to the fake profiles).** No sufficiently large squarefree skew double core from Part III.2 can occur over \(\mathbb F_2[t]\), at any parent degree. On its \(m\) ordinary rows the two binary columns with label \(a\) add modulo two to the coordinate vector \(e_a\), since the flip changes only that row after the border is deleted. The unanchored columns therefore span \(\mathbb F_2^m\). Anchoring leaves rank \(m-1\), unchanged by removing zero or repeated anchored columns. There are at most \(2m\) distinct columns. This contradicts (VI.1) as \(m\to\infty\). The proof also works with disjoint blocks of primes per column provided their odd valuations still realize those patterns. A parent with more distinct patterns or a smaller parity rank can pass (VI.1).

**COMPUTED.** For \(m=167\) ordinary Paley rows the actual anchored binary rank is \(166\) on \(334\) columns. Taking all \(40\)-subsets in the Markov argument already gives
\[
  2^{166}(334-8\cdot40)\binom{334}{40}>334\cdot2^{334}.
\]
Thus this particular profile is impossible on every binary parent, with arbitrary weights. Integer arithmetic gives the same contradiction at \(m=191,251\). The script checks the paired columns and computes the ranks by binary elimination. [Packing certificate](vi_profile_check.py)

**PROVED (a low-rank test of the same derivative condition).** Take \(m=2^d\), \(d\ge3\) odd, and two Walsh cores on the same \(m\) rows. Label their rows by \(x\) and \(x^3\), respectively, for \(x\in\mathbb F_m\). Use all nonzero binary linear forms on each label as the columns, with a distinct split prime of the same even degree \(e\) at each column and with exponent one. The span on the two cores is the direct sum of their \(d\)-dimensional binary row codes. Indeed a linear relation on \((x,x^3)\) would give \(\operatorname{Tr}(A x+B x^3)=0\) for every \(x\). Polarizing in \(x,y\) and using the nondegenerate trace pairing gives \((By)^{2^{d-1}}+By^2=0\) for every \(y\). Squaring yields \(By+B^2y^4=0\) identically on the field of more than four elements, whence \(B=A=0\). Cubing permutes the labels since \(d\) is odd. Each nonzero word on either single core occupies \(m/2\) columns, so the direct sum passes the weighted-distance test at \(W=2(m-1)e\).

Nevertheless the *polynomial* condition in the proof of (VI.1) requires \(e\ge4d-2\). For each nonzero word on the first core alone its prime product \(H_u\) has degree \(me/2\). Consequently \(\deg H'_u\le me/2-K-1=e/2-1\). These \(m-1\) nonzero derivatives are all distinct. Otherwise \(H_u+H_v\) for distinct words would be a nonzero square, of degree strictly below \(me/2\) by monicity. Their \(m/4\) common primes divide this square, so the square of their product divides it, already of degree \(me/2\). Finally \(H'_u\) is itself a square over \(\mathbb F_2[t]\). There are at most \(2^{\lceil e/4\rceil}-1\) nonzero squares of degree at most \(e/2-1\), proving the bound. This forces \(W\ge (2m-2)(4d-2)\) on these profiles. It still permits \(W=O(m\log m)\).

## VI.2. The degree cost of a complete Hadamard core

**PROVED.** Let the \(m\) rows and \(m-1\) nonconstant sign columns \(\chi_v\) form a normalized Hadamard matrix after adjoining the all-ones column, with \(m\ge4\). At label \(v\) let each mark choose \(f_v\) or \(\bar f_v\) according to the corresponding sign. Here \(f_v\) is monic and nonconstant, coprime to its conjugate, and the base-prime supports of distinct labels are disjoint. Each block can have arbitrary degrees and powers and can incorporate independent choices of orientation at its primes. These are conditions on the prime profile, with no additive hypothesis on the roots. Remove any common factor of the polynomial points. If the points form a binary critical cluster, \(m\) must be a power of two, and their point degree satisfies
\[
          D\ \ge\ \frac{m(m-1)}6\bigl(\log_2(m-1)-3\bigr).       \tag{VI.2}
\]
Removing a common factor from a critical cluster preserves the half-degree inequality. The theorem thus also bounds the original degree from below.

*Proof.* Put \(d_v=\deg f_v\), \(D=\sum_v d_v\), \(r_v=f_v/\bar f_v\), and take the positive sign to choose \(f_v\). The marks can be made simultaneously monic since they have the same leading coefficient. Every integral product of the marks whose exponents sum to zero lies in \(U_K\). Orthogonality and balance of the columns give the exact products
\[
       \prod_i z_i^{\chi_v(i)}=r_v^{m/2},\qquad
       \prod_i z_i^{(\chi_v(i)-\chi_w(i))/2}=(r_v/r_w)^{m/4}
       \quad (v\ne w).
\]
Here \(4\mid m\), as follows by counting the four sign pairs in two balanced orthogonal columns. For a principal unit \(r\ne1\), define its depth as \(\operatorname{ord}_\infty(r-1)\). Raising to \(2^a b\), with \(b\) odd, multiplies its depth by \(2^a\). This follows from Frobenius and the nonzero linear term of \((1+x)^b-1\).

If \(u\) is the odd part of \(m\), the first product forces depth \((r_v)\ge uD/m\). The nonzero polynomial \(f_v+\bar f_v\) gives depth \((r_v)\le d_v\). Summing over the \(m-1\) blocks forces \(u=1\). Now both products give
\[
       \operatorname{depth}(r_v)\ge D/m,\qquad
       \operatorname{depth}(r_v/r_w)\ge 2D/m.
\]
Write \(s_v=d_v-D/m\), so \(s_v\ge0\) and \(\sum_v s_v=D/m\). Write \(f_v=A_v+\theta B_v\), with \(\theta^2+\theta+1=0\) and \(A_v,B_v\in\mathbb F_2[t]\). From the first inequality \(0\ne B_v\) and \(\deg B_v\le s_v\). The numerator of \(r_v/r_w-1\) is \(f_v\bar f_w+\bar f_v f_w=A_v B_w+A_w B_v\), before common cancellation. The second inequality gives its degree at most \(s_v+s_w\).

At least \((m-1)/2\) labels have \(s_v\le s=2D/[m(m-1)]\). Among these there are at most \(2^{s+1}\) choices of \(B_v\). At a fixed \(B_v=B\), comparison with one label \(w\) gives \(\deg(A_v-A_w)\le2s-\deg B\), hence at most \(2^{2s+1}\) choices of \(A_v\). Distinct labels have distinct monic factors by their disjoint prime supports. Thus \((m-1)/2\le2^{3s+2}\), which is (VI.2). \(\square\)

## VI.3. Allowing a nonconstant translation

**PROVED.** Let \(V\subset\mathbb F_2[t]\) be an additive subspace of size \(Q\), \(0\ne B\in V\), and \(Y=T+\theta B\) for real \(T\in\mathbb F_2[t]\). Use **all** \(s=Q/2\) forms \(\ell(B)=1\), with either choice of side \(S_\ell\) for each form. Form \(F_\ell=\prod_{v\in S_\ell}(Y+v)\), permit unit twists and common division, and suppose the resulting points are a critical cluster of degree \(D'\). Then \(D'=\Omega(s^2)\). This also covers shared primes among the displayed factors, so VI.2 alone does not give the result.

*Proof.* Dividing \(T,B,V\) by their common real polynomial content changes the raw products by a common factor. We may also pass to the full gcd of the products, since this decreases the degree and preserves criticality. Use the content-free data in what follows. The involution on the labels is \(v\leftrightarrow v+B\). Call such a pair variable if both orientations occur across the chosen sides. At most one pair is fixed. Indeed for two fixed pairs represented by \(v,w\), the value of \(\ell(v+w)\) would be constant on *all* forms taking \(B\) to \(1\), so \(v+w\in\{0,B\}\).

For a label pair define its height
\[
        h_v=\max(\deg(T+v),\deg B)-\deg\gcd(T+v,B).
\]
This is the numerator degree in the reduced unit \((Y+v)/(Y+v+B)\), since its numerator and denominator are conjugate and their common monic divisor is \(\gcd(T+v,B)\). In particular the reduced numerator is supported only above split primes. We claim
\[
          D'\ \ge\ \frac12\sum_{\text{variable pairs}}h_v
             \ \ge\ (Q/8-1/2)\deg B.                           \tag{VI.3}
\]

Fix a split prime \(p=f\bar f\), choose one factor \(f\), set \(b=\operatorname{ord}_p B\), and put \(A_j=\{v\in V:f^j\mid Y+v\}\). The valuation of \(F_\ell/\bar F_\ell\) is the sum over levels \(j\) of \(2|S_\ell\cap A_j|-|A_j|\). For \(j\le b\) this contribution is zero, since \(A_j\) is invariant under \(v\mapsto v+B\). For \(j>b\), the nonempty \(A_j\) are nested affine spaces with direction contained in \(\{v\in V:f^{b+1}\mid v\}=K_f\), and \(B\notin K_f\). If any such level contains a variable label, take \(v_*\) in the deepest such level. A form taking \(B\) to \(1\) and annihilating \(K_f\) exists and occurs among the marks. Its selected side takes all or none of every \(A_j\) from \(b+1\) through that deepest level, with the same sign throughout. Because \(v_*\) is variable, another selected side makes the opposite choice at \(v_*\). It takes at most half of each of these affine spaces if it excludes \(v_*\), or at least half if it includes \(v_*\). At deeper levels only fixed labels remain, with the same contribution for both marks. The range of the valuations across the marks is therefore at least the sum, over variable label pairs, of
\(\bigl|\operatorname{ord}_f((Y+v)/(Y+v+B))\bigr|\): each unit of this last sum is one occurrence of a variable label in an \(A_j\) above level \(b\).

If \(p\) has exponent \(e_p\) in the reduced parent, then \(\operatorname{ord}_f(z_\ell/\bar z_\ell)\) varies over an interval of length at most \(2e_p\). Common division shifts the valuations by the same amount, and units have zero valuation. Multiply the range inequality by \(\deg f\) and sum once per conjugate prime pair. For each primitive block ratio the sum of the absolute valuations times \(\deg f\) is exactly \(h_v\). Since \(D'\ge\sum_{p\ {\rm split}}e_p\deg f\), this proves the first inequality of (VI.3). For \(p\mid B\), at most \(Q/2\) labels satisfy \(p\mid T+v\). Otherwise an affine fiber would be all of \(V\), contradicting common content one. At least \(Q/4-1\) variable pairs therefore contribute the full \((\operatorname{ord}_p B)\deg p\) in the lower bound \(h_v\ge\deg B-\deg\gcd(T+v,B)\). Summing proves (VI.3).

To finish, retain \(s'\ge s/3\) marks with a common unit twist. Take the largest odd \(e\) with \(2^{e+1}\le s'\), giving \(s'/8<2^e\le s'/2\). If \(\deg B\ge2^{e-2}\), (VI.3) gives the conclusion. Otherwise consider odd-degree primes \(p\nmid B\) of this degree. Their degrees sum to at least \(2^{e-1}\) for sufficiently large \(s'\), since among the \(2^e\) elements of \(\mathbb F_{2^e}\), at most \(\sum_{j\le e/2}2^j\) have smaller degree, and deleting primes dividing \(B\) costs at most \(\deg B\).

The following count applies the referees' Theorem E at the good inert primes \(p\). At most \(2^{e-1}\) selected forms vanish on the kernel of \(V\to\mathbb F_2[t]/p\), since these forms factor through the image and send the nonzero class of \(B\) to \(1\). For any two other forms, their restrictions to that kernel attain the pair \((1,1)\): both are nonzero, so they are either equal or independent. If \(J\) is their common kernel on \(V\), at least \(Q/(4\cdot2^e)\) elements of \(J\) therefore equal \(B\bmod p\). Hence \(\operatorname{ord}_p P_J(B)\ge Q/(4\cdot2^e)\).

For the two raw marks the subspace-polynomial identity gives \(F_i-F_j=P_J(B)P_J(Y+a)\), with \(a\in S_i\cap S_j\), exactly as in (V.1). The inert \(p\) divides no \(Y+v\), since \(p\nmid B\). More precisely, after cancelling the displayed common factor \(P_J(Y+a)\), the remaining two factors have the form \(X_0,X_0+P_J(B)\). Thus
\(\delta_{ij}=\deg(F_i-F_j)-\deg\gcd(F_i,F_j)\) is the degree of \(P_J(B)/\gcd(X_0,P_J(B))\). All the counted inert valuations contribute additively to \(\delta_{ij}\). Each prime contributes on \(\binom{s'-2^{e-1}}2=\Omega(s'^2)\) or more pairs. Summation gives \(\sum_{i<j}\delta_{ij}\gg Qs'^2\).

For completeness the common-norm bound in the other direction is \(\sum\delta_{ij}\le s'D'/4\). Common division preserves \(\delta_{ij}\) and the retained units agree. Start each pair gcd count at \(D'\). At each split-prime exponent layer its loss is \(\deg f\) on the pairs separated by the layer, of which there are at most \(s'^2/4\). Inert common factors cause no loss. Hence \(\sum\deg\gcd(z_i,z_j)\ge D'[\binom{s'}2-s'^2/4]\). Using \(\deg(z_i-z_j)\le D'/2\) gives the asserted bound. This proves \(D'\gg Qs'\gg s^2\). \(\square\)

**HEURISTIC (where the growth proof stops).** A sparse selection of forms can omit every form annihilating \(K_f\). At that prime two sides may both take exactly half of each relevant fiber, and the valuation range in (VI.3) can disappear. More generally, (VI.1) permits a low-rate binary code, and the exact integral characters in VI.2 require orthogonality of a complete core. I have no implication forcing any of these profile conditions on an arbitrary growth-scale parent.

## VI.4. Unbounded clusters on the ternary circle

**PROVED.** For every \(m\ge2\) there exist \(m\) distinct monic \(z_j\in\mathbb F_9[t]\) of one degree \(D\), with the same \(z_j\bar z_j=N\in\mathbb F_3[t]\) and \(\deg(z_i-z_j)\le D/2\). One may take \(\deg N=\exp(O(m\log m))\).

*Root labels and descent.* Set \(k=\lceil\log_3(2m+1)\rceil\). Take \(q=3^r>2^{m-1}\), with \(r\) the least positive odd integer satisfying this inequality, so \(q\le9\cdot 2^{m-1}\). Let \(i^2=-1\) and work temporarily over \(\mathbb F_{q^2}\). Define
\[
        V=\left\{\,i\sum_{b=0}^{k-1}x_b t^b:
                            x_b\in\mathbb F_q\,\right\},
              \qquad n=|V|=q^k,\qquad h=n/q .
\]
Choose \(m\) pairwise nonproportional linear forms \(\lambda_j\) in the \(x_b\), with coefficients in \(\mathbb F_3\). There are \((3^k-1)/2\) available, and the chosen forms remain nonproportional over \(\mathbb F_q\). Let \(\chi\) be the quadratic character of \(\mathbb F_q\), extended by \(\chi(0)=0\). For a sign \(\epsilon_j=\pm1\) define
\[
          F_j(X)=\prod_{v\in V}(X-v)^{\,1+\chi(\epsilon_j\lambda_j(v))}.
                                                                    \tag{VI.4}
\]
Here \(\lambda_j(v)\) means evaluation on the coefficient vector \(x\) of \(v\). Each \(F_j\) is monic of \(X\)-degree \(n\), since each value of a form occurs \(h\) times.

The coefficient map \(\alpha\mapsto\alpha^9\) fixes \(i\), sends the vector \(x\) to \(x^9\), and preserves its exponents in (VI.4), since \(\chi(y^9)=\chi(y)\). Its fixed field in \(\mathbb F_{q^2}\) is \(\mathbb F_9\). Thus every \(F_j\) already belongs to \(\mathbb F_9[t,X]\). On that ring coefficient conjugation is induced by the \(q\)-power map, since \(r\) is odd. This sends \(v\) to \(-v\), and \(\chi(-1)=-1\). Paired exponents therefore sum to two, including at the zeros of a form. In particular,
\[
                     F_j\bar F_j=P_V(X)^2,\qquad
                     P_V(X)=\prod_{v\in V}(X-v).
\]
The distinct forms have distinct kernels. At a vector in one kernel outside the other the exponents are \(1\) and \(0\) or \(2\). Thus the formal marks are distinct for every choice of signs.

*Degree of the error before cancellation.* Fix \(j\). Put \(H=\ker\lambda_j\), \(A=P_H(X)\), and \(c=P_H(v_0)\) for any \(v_0\) with \(\epsilon_j\lambda_j(v_0)=1\). The polynomial \(P_H\) is \(\mathbb F_q\)-linear as a map in \(X\). Indeed for a subspace \(J\) and \(u\notin J\),
\[
 P_{J+\mathbb F_q u}(X)=P_J(X)^q-P_J(u)^{q-1}P_J(X),
\]
by \(\prod_{a\in\mathbb F_q}(Z-a w)=Z^q-w^{q-1}Z\) and induction starting with \(P_{\{0\}}(X)=X\). Thus \(P_H(v)=c\,\epsilon_j\lambda_j(v)\), and \(P_V=A^q-c^{q-1}A\). Taking the zero fiber once and the nonzero square fibers twice gives
\[
 \begin{split}
       F_j&=A\bigl(A^{(q-1)/2}-c^{(q-1)/2}\bigr)^2\\
          &=P_V+2c^{q-1}A-2c^{(q-1)/2}A^{(q+1)/2}.
 \end{split}                                                       \tag{VI.5}
\]
Since \(\deg_X A=h\), all pair differences have \(X\)-degree at most \((n+h)/2\). This permits nonconstant cofactor differences under reflection.

*Common division.* Let \(G\) be the product of the same linear factors with the minimum, over all \(j\), of their exponents in (VI.4), and write \(g=\deg_X G\). It too descends to \(\mathbb F_9[t,X]\). Choose the \(\epsilon_j\) independently and fairly, for this existence argument only. If all \(m\) form values at \(v\) are nonzero, the expected minimum exponent equals \(2\cdot2^{-m}\). If exactly \(a<m\) are nonzero, the minimum is \(1\) if all those \(a\) orientations make their character positive, otherwise zero. Its expectation is \(2^{-a}\ge2^{1-m}\). Summing shows that some deterministic choice of orientations has
\[
                      g\ge 2^{1-m}n>h .                           \tag{VI.6}
\]
The formal reduced marks \(Z_j=F_j/G\) have degree \(n-g\), common norm \(P_V^2/(G\bar G)\), and pair differences of degree at most \((n+h)/2-g\). The strict slack below half their degree is \((g-h)/2\). The orientation argument selects a single parent and a single family.

*Specialization.* Each \(Z_j\) is a product of \(n-g\) monic linear factors whose roots have \(t\)-degree at most \(k-1\). Its coefficient of \(X^a\) has \(t\)-degree at most \((k-1)(n-g-a)\). Take \(L=2(k-1)n+1\) and substitute \(X=t^L\). The intervals of possible degrees belonging to different \(X\)-powers are disjoint. All differences remain nonzero, all points have degree \(D=(n-g)L\), and
\[
 \deg_t(z_i-z_j)
       \le \bigl((n+h)/2-g\bigr)L+(k-1)n
       < (n-g)L/2 .
\]
Here \(g-h\ge1\). The substitution commutes with conjugation, so the norm equality is exact after specialization. Since \(k=O(\log m)\) and \(\log q=O(m)\), the parent degree \(2D\le2nL\) is \(\exp(O(m\log m))\). This proves the theorem. \(\square\)

**COMPUTED.** For \(q=27,\ k=2\), the four forms \(x,y,x+y,x-y\) with positive orientations give \(n=729,\ h=27,\ g=92\). At \(X=t^{11}\) the full gcd of the four specialized raw marks has degree \(1{,}038\). After division the four actual polynomials have degree \(6{,}981\), the common parent has degree \(13{,}962\), and the six pair-difference degrees are \(3471,3457,3457,3471,3471,3457\). The script multiplies the polynomials by their coefficient conjugates and checks equality in \(\mathbb F_9[t]\). [Ternary certificate](vi_ternary_check.py)

**COMPUTED.** There is also a five-mark factor certificate, exceeding the exhaustive maximum through \(D=11\). Use \(q=27,\ k=3\) and the forms \(x,y,z,x+y,x+z\) with positive orientations. Enumeration of the \(19{,}683\) root labels gives \(g=1301,\ h=729\). At \(L=67\) the degree of each formally reduced point after substitution is **\(1{,}231{,}594\)**. For a difference, apply the coefficient bound from the proof at \(a\le (n+h)/2-g\). Since \(L>k-1\), its specialized degree is at most
\[
       \bigl((n+h)/2-g\bigr)L+(k-1)(n-h)/2
             =615{,}589 < D/2=615{,}797 .
\]
The script distinguishes all five specialized raw polynomials by direct evaluation of (VI.5) modulo \(t^5-t+1\). Distinctness survives division by \(G(t^L)\). The parent degree is \(2{,}463{,}188\). This five-mark check uses factor exponents and modular arithmetic. The four-mark check above expands the points and norms.

## VI.5. Squarefree binary parents at the same quadratic order

**PROVED.** For \(m=2^{d-1}\), \(d\ge3\), and every integer \(L\ge4m(d-2)\), there is a squarefree \(N\in\mathbb F_2[t]\) of degree \(2(m-1)L\) carrying \(m\) points in a critical cluster. In particular the least degree supplied here is \(8m(m-1)(\log_2m-1)\). Rounding up to a power of two also gives \(O(M^2\log M)\) for arbitrary \(M\).

*Proof.* Let \(U\) be the span of \(t,\ldots,t^{d-1}\) over \(\mathbb F_2\). Define
\[
      V=\{\,a+(f^2+f)^2:a\in\mathbb F_2,\ f\in U\,\}.
\]
It has dimension \(d\). Indeed the map on \(f\) is additive, every nonzero input has positive degree, and the leading degrees of its displayed basis are \(4,8,\ldots,4(d-1)\). All elements of \(V\) are squares. At any even-degree prime \(p\), their residues have absolute trace zero: the trace of \(f^2+f\) is zero, squaring preserves trace, and \(1\) has trace zero in even degree. This leaves room to choose a translate even at the smallest split primes.

Let \(R\) be the product of all distinct even-degree primes dividing at least one nonzero element of \(V\). We have
\[
                           \deg R\le2m(d-2).                       \tag{VI.7}
\]
Indeed \(R\) divides \(H=\prod_{f\in U}(f^2+f+1)\). A prime dividing \(1+(f^2+f)^2=(f^2+f+1)^2\) lies in that product directly. A prime of even degree dividing a nonzero \((f^2+f)^2\) divides \(f\) or \(f+1\), so its degree is at most \(d-1\). The space \(U=t\,\mathbb F_2[t]_{<d-1}\) surjects modulo every such prime, because \(t\) is invertible there. Its residue field contains both roots of \(X^2+X+1\), giving a factor of \(H\) divisible by this prime. In particular \(t^2+t+1\) occurs in at least two factors of \(H\), from its two distinct root residues. Since \(R\) is squarefree and \(\deg H=2\sum_{j=1}^{d-1}j\,2^{j-1}=2+2m(d-2)\), (VI.7) follows.

Use \(\theta^2+\theta+1=0\). Set \(E=t^L+t\) for even \(L\) and \(E=t^L\) for odd \(L\). Seek \(T=E+S^2\). At each \(p\mid R\), prescribe \(\operatorname{Tr}(T)=1+\operatorname{Tr}(\theta)\) in its residue field. Either embedding of \(\theta\) gives the same trace, since its two values differ by \(1\). The absolute trace is surjective because its nonzero trace polynomial has degree smaller than the field size. Squaring is bijective, so an allowed residue of \(S\) exists. The Chinese remainder theorem gives a simultaneous \(S\) of degree \(<\deg R\). By (VI.7), \(2\deg S<L\), hence \(\deg T=L\).

The polynomials \(T+\theta+v\), \(v\in V\), are pairwise coprime over \(\mathbb F_4\). At a common root \(\alpha\) of two of them, \(\theta=T(\alpha)+v(\alpha)\in\mathbb F_2(\alpha)\), so the base prime has even degree. It divides a nonzero difference in \(V\), hence divides \(R\). But the absolute traces of the two sides of \(T(\alpha)+\theta=v(\alpha)\) would be \(1\) and \(0\). Moreover each factor is squarefree. Its derivative is \(1\) for even \(L\), or \(t^{L-1}\) for odd \(L\). In the latter case the factor's value at zero has nonzero \(\theta\)-component.

Now use the Part V marks \(z_H=\prod_{v\in H\setminus\{0\}}(T+\theta+v)\) for the hyperplanes excluding \(1\). Their common norm is the product over \(v\in V\setminus\{0,1\}\), hence is squarefree. Every factor has degree \(L>\max_{v\in V}\deg v\). In (V.1) the intersection \(J\) has dimension \(d-2\) and excludes \(1\). Its echelon degrees are bounded term by term by \(8,12,\ldots,4(d-1)\), the largest \(d-2\) basis degrees in \(V\). Consequently
\[
     \deg P_J(1)\le \sum_{j=1}^{d-2} 2^{j-1}4(j+1)
                  =4(d-2)2^{d-2}\le L/2 .
\]
Equation (V.1) gives distinct points of degree \((m-1)L\) with the required half-degree contact. \(\square\)

**PROVED (sharpness for the profile in VI.2).** Choose a linear complement \(E_0\) of \(\langle1\rangle\) in \(V\). The remaining conjugate pairs are indexed by \(0\ne v\in E_0\). A mark chooses the factor at \(v\) or \(v+1\) according to the value of \(\ell|_{E_0}(v)\). All functionals on \(E_0\) occur, so these are precisely the \(m-1\) nonconstant columns of the Walsh Hadamard matrix. The squarefreeness proof gives disjoint prime supports for its blocks. Together with (VI.2) this determines the least binary parent degree to order \(\Theta(m^2\log m)\) within complete Hadamard cores, along powers of two, even when the primes and orientations within each block can be chosen freely.

**COMPUTED.** At \(m=4,8,16,32\) the CRT moduli have degrees \(6,26,60,176\). The initial squarefree parent degrees are \(96,896,5760,31744\), respectively. The certificate checks \(\gcd(N,N')=1\), every individual norm, and all 650 pair differences by multiplication and subtraction. [Squarefree certificate](vi_squarefree_check.py)

**PROVED (transfer to the integers).** The ternary construction disproves uniformity using a fixed odd characteristic, reflection, finite residue fields and exact anisotropic norms together. Its added input is the finite-field identity (VI.5), with the \(q\)-linearized polynomial and coefficient Frobenius providing the degree drop and descent. Over the binary field, (VI.1) uses the ultrametric multiplicative neighborhood \(U_K\) and a rational derivative that kills squares. In VI.2 the exact amplification of depth under \(2^a\)-th powers and the finite number of bounded-degree traces are essential. VI.3 uses binary additive fibers and even/odd residue degrees. These proofs supply neither an integer growth theorem nor an integer construction. For Part V, the original choice of \(L\) was minimal for its specified undivided family. I missed the unit-twist issue in V.2. The referees' Theorem E handles it, and their degree-\(104\) example improves the eight-mark degree when the sides and common divisor may vary.


## Twelve-line summary of Part VI

1. **HEURISTIC (research status).** I could not determine the unrestricted binary growth rate or improve the integer bound.
2. **PROVED.** At every genuine binary parent, the entire span of the parity rows has weighted minimum distance greater than \(W/4\).
3. **PROVED.** This excludes the record's sufficiently large squarefree skew double-core profile for every assignment of split primes and degrees.
4. **COMPUTED.** Exact weighted packing already excludes the example on 167 ordinary Paley rows and 334 columns.
5. **PROVED.** Complete Hadamard prime profiles require order \(m^2\log m\) in characteristic two, and squarefree Walsh examples attain that order.
6. **PROVED.** Full hyperplane designs with nonconstant translation have degree \(\Omega(m^2)\), including after all common polynomial cancellation.
7. **PROVED.** The ternary circle is unbounded over the fixed field \(\mathbb F_3(t)\), with parent degree at most \(\exp(O(m\log m))\).
8. **PROVED.** Quadratic-character multiplicities give exact reflection norms. A common divisor removes the degree excess, and every mark descends to \(\mathbb F_9[t]\).
9. **COMPUTED.** Five ternary points have degree \(1{,}231{,}594\) with every pair difference of degree at most \(615{,}589\).
10. **PROVED.** For each \(m=2^r\ge4\) a squarefree binary parent of degree \(8m(m-1)(r-1)\) carries \(m\) points at critical width.
11. **COMPUTED.** Direct arithmetic checks squarefreeness, every norm, and 650 pair differences across four binary examples.
12. **HEURISTIC (request).** Please referee the ternary descent and common cancellation first. For growth, I would test the two Walsh cores on \(x,x^3\) at prime degrees \(e\ge4\log_2 m\), where the derivative count still permits growth scale.

# Part VII. Multiplicity two in odd characteristic

**PROVED.** I took (a) first. The multiplicities in VI.4 have two sources. Two copies of a linear factor occur in the character formula itself. Many *different* linear factors also acquire a common zero after specialization. The second source can be removed entirely. For every fixed prime \(p\equiv3\pmod4\), critical clusters on parents with **every exponent at most two** have unbounded size. This holds both for the constant extension and for a fixed ramified quadratic extension with inert infinity. In the latter the involution fixes the constants.

**HEURISTIC (research verdict).** I do not know if the supremum for squarefree ternary parents is finite. The squarefree census first reaches five at point degree 12. Below I give an obstruction to obtaining an unbounded squarefree family by specializing the original VI.4 marks, even with extra common cancellation. For (b) I can exclude a single Hadamard matrix with missing columns and all its rows, but I could neither construct the two-core example nor improve the unrestricted growth bound. I found no improvement over the integers.

## VII.1. Separating the finite zeros

**PROVED.** Fix \(p\equiv3\pmod4\). Take either \(\eta^2=-1\) or \(\eta^2=\Delta(t)\), with \(\Delta\in\mathbb F_p[t]\) of positive even degree and nonsquare leading coefficient, with a finite ramified prime. Write \(\Delta=-1\) in the first case. Let \(M^{(\le2)}(D)\) denote the maximum number of distinct points of pole order \(D\) with common norm in \(\mathbb F_p[t]\), all its prime exponents at most two, and pair differences of pole order at most \(D/2\). Use the order \(\mathbb F_p[t,\eta]\) and the inert place over infinity. For \(a+b\eta\) the pole order is \(\max(\deg a,\deg b+\tfrac12\deg\Delta)\). In either model,
\[
 \liminf_{D\to\infty} M^{(\le2)}(D)\frac{\log\log D}{\log D}
       \ \ge\ \frac12\log_2 p.                                      \tag{VII.1}
\]
In particular every fixed cap \(e_0\ge2\) permits arbitrarily large clusters. The unpadded parents constructed below are coprime to the finite discriminant.

*Proof.* Given \(m\ge2\), choose \(k=\lceil\log_p((p-1)m+1)\rceil\) and choose \(r\) least positive odd with \(q=p^r\ge 2^m\). Write \(n=q^k, h=n/q\). Fix \(m\) linear forms on \(\mathbb F_q^k\) defined over \(\mathbb F_p\), mutually nonproportional, possible by the choice of \(k\). Work temporarily over a constant field containing \(\mathbb F_{q^2}\). Choose a nonzero \(C\in\mathbb F_p[t]\) divisible by
\[
             \Delta(t)\prod_{j=1}^{k-1}(t^{q^j}-t).
\]
One may take this product itself up to a constant. Extra factors for padding will be specified below. The vector space of roots is
\[
 V=\left\{\eta^p C^p\sum_{j=0}^{k-1}x_j t^{pj}\ :\ x_j\in\mathbb F_q\right\}.
                                                                    \tag{VII.2}
\]
It has size \(n\), its elements are \(p\)-th powers, and each has pole order at most \(b=p(\tfrac12\deg\Delta+\deg C+k-1)\) at each geometric infinity. The maps \(\ell_i\) are the chosen forms on its coefficient vector. If two different labels take the same finite value over \(t=\alpha\), then \(C(\alpha)=0\). Indeed, apart from zeros of \(C\eta\), their difference is a unit times a \(p\)-th power of a nonzero polynomial over \(\mathbb F_q\) of degree less than \(k\). Every root of such a polynomial is a zero of one of the displayed \(t^{q^j}-t\). At every zero of \(C\), *all* labels vanish.

For any \(L>p\deg C+1\), put \(s=L-p\deg C\). Take \(T=1+C^p t^s\) if \(p\nmid s\), and \(T=1+C^p(t^s+t)\) otherwise, with \(C\) normalized monic. Then \(T=1\) on the collision set, and \(T'\) is nonzero everywhere off it since \(t\mid C\). The functions \(T-v\) have disjoint finite zero sets. Off \(C=0\) the cover is unramified (even if the defining radicand contains a square factor, this open set omits every zero of the radicand), \(t-\alpha\) is a local parameter after extending constants, and \(d(T-v)=T'\,dt\). Thus every finite zero is simple. At branch points all factors take the value one. Consequently
\[
                      R=\prod_{v\in V}(T-v)
\]
is squarefree in \(\mathbb F_p[t]\). To see the field assertion, the involution sends the set of labels to itself by negation, and the constant-field automorphisms permute the coefficient vectors. Thus the product is an integral polynomial on the base. In the constant case this also follows by combining the \(p^2\)-power map on constants with the \(q\)-power map there. A multiple base zero would remain multiple upstairs, whereas every finite zero just counted has order one.

For completeness, let \(\chi\) be the quadratic character of \(\mathbb F_q\) with \(\chi(0)=0\), so \(\chi(-1)=-1\). For signs \(\epsilon_j=\pm1\), form
\[
 F_j(X)=\prod_{v\in V}(X-v)^{1+\chi(\epsilon_j\ell_j(v))},\qquad
 P_V(X)=\prod_{v\in V}(X-v).
\]
All \(F_j\) lie in the required quadratic order with \(X\) adjoined. In the ramified case constant Frobenius fixing \(\eta\) raises each \(x_j\) to its \(p\)-th power and preserves its character after applying a form. In the constant case use the \(p^2\)-power map fixing \(\eta\). Conjugation can in both cases be extended to fix \(\mathbb F_q\) and negate \(\eta\). In the constant case this is the \(q\)-power map on constants. It follows that \(F_j\bar F_j=P_V^2\).

Let \(A=P_{\ker\ell_j}\), and let \(c=A(v_0)\) where \(\epsilon_j\ell_j(v_0)=1\). For an \(\mathbb F_q\)-space \(H\), its monic root polynomial is \(\mathbb F_q\)-linear as a map. Induct on its dimension using \(P_{H+\mathbb F_q u}=P_H^q-P_H(u)^{q-1}P_H\), starting at \(P_0=X\). In particular \(P_V=A^q-c^{q-1}A\). As \(A(v)=c\epsilon_j\ell_j(v)\), the selected nonzero fibers give
\[
 F_j=A(A^{(q-1)/2}-c^{(q-1)/2})^2
     =P_V+2c^{q-1}A-2c^{(q-1)/2}A^{(q+1)/2}.
\]
Hence \(\deg_X(F_i-F_j)\le(n+h)/2\). In the formal common divisor \(G\) take the minimum multiplicity at each label. At a label with all \(m\) values nonzero, independent fair signs give expected minimum \(2^{1-m}\). With only \(a<m\) nonzero, the expectation is \(2^{-a}\ge2^{1-m}\). Fix one choice with \(g=\deg_XG\ge2^{1-m}n\). In particular \(g-h\ge n/2^m\). The minimum exponents are invariant under the same descent as the marks, so \(Z_j=F_j/G\) are polynomials over the stated order.

Take \(L>b\). Each \(T-v\) has poles of order \(L\) at infinity and no finite poles, hence has finite zeros (on each component in the constant model). Two nonproportional forms have some label with different exponents. The disjoint zeros of the \(T-v\) therefore distinguish all \(z_j=Z_j(T)\), including after common division. These marks have degree \(D_0=(n-g)L\). Their common norm \(N_0\) divides \(R^2\). Some parent prime attains the cap. For any two independent forms, a label can have signed characters \(+1,-1\) respectively, and its negative then has signs \(-1,+1\). Both labels survive in \(R^2/(G(T)\bar G(T))\) with exponent two at their zeros.

The coefficient of \(X^a\) in \(Z_j\) has pole order at most \(b(n-g-a)\), by multiplying its linear factors. Its pair differences have \(X\)-degree at most \((n+h)/2-g\). Since \(L>b\), for every pair
\[
  2\deg_\infty(z_i-z_j)
       \le (n+h-2g)L+b(n-h).
\]
In particular \(D_0-2\max\deg_\infty(z_i-z_j)\ge(g-h)L-b(n-h)\). This estimates the marks after precisely the formal division by \(G\).

To cover every sufficiently large point degree while preserving the cap, include in \(C\) one monic irreducible of each degree \(1,2,4,\ldots,2^{\lceil\log_2 n\rceil}\) over \(\mathbb F_p\). Such irreducibles exist since \(\mathbb F_{p^d}\) has an element of degree \(d\) (for \(d>1\) a generator of its multiplicative group suffices). The added degree is at most \(4n\). The initial product has degree \(O_\Delta(n)\). For target \(D\) put \(L=\lfloor D/(n-g)\rfloor, u=D-(n-g)L<n\). A squarefree product \(e\) of the reserved primes has degree \(u\) by binary expansion. Each reserved prime divides \(R-1\). Therefore \(e^2N_0\) still has exponent cap two, and \(e z_j\) give the required degree. The width holds as soon as the slack above is at least \(u\), in particular if \(L\ge 2^m(b+3)+p\deg C+2\). Since \(b=O_{p,\Delta}(n)\), it suffices that \(D\ge C_{p,\Delta}2^m n^2\) for a fixed constant. Finally \(\log q=m\log 2+O_p(1)\) and \(k=(\log m)/(\log p)+O_p(1)\). The logarithm of this sufficient threshold is at most \((2\log2/\log p+o(1))m\log m\). Taking \(m=\lfloor a\log D/\log\log D\rfloor\) for any \(a<(\log p)/(2\log2)\) proves (VII.1). \(\square\)

**COMPUTED.** With \(p=3,q=27,k=2\), use \(x,y,x+y,x-y\), all signs positive, and take \(C=t^{27}-t, T=1+tC^{36}\). Here \(C\) already vanishes on the branch locus in the ramified example. Both examples have \(L=973,n=729,h=27,g=92\). The script expands \(R\) of degree 709317 over \(\mathbb F_3\) and obtains \(\gcd(R,R')=1\) in both models. The four marks have degree 619801 and common parent degree 1239602. Their certified pair-degree upper bounds are 307762 for \(\eta^2=-1\) and 308815 for \(\eta^2=2t^2+1\). This certificate uses factor exponents for the marks and an expanded gcd for squarefreeness of \(R\). [Cap-two certificate](vii_cap2_check.py)

**PROVED (transfer scope).** Thus ramification, a conjugation fixing constants, and any exponent cap \(e_0\ge2\) can hold simultaneously with unbounded clusters in the function-field model. The decisive identity makes the *same* \(P_V\) a polynomial in \(P_{\ker\ell_j}\) for every \(j\). For comparison, in characteristic zero, if \(P=f(A)=g(B)\) with \(A,B\) of the same positive degree, normalize the inner polynomials to be monic. The outer polynomials have the same degree \(d\) and the same leading coefficient \(c\). At infinity the monic Laurent \(d\)-th root of \(P/c\) has polynomial part \(A+\text{constant}\) and also \(B+\text{constant}\), by expanding in \(1/A\) or \(1/B\). Hence the normalized factors differ by a constant. Our distinct monic root polynomials vanish at zero and have different kernels. The additional device \(dv=0\) in (VII.2) also uses positive characteristic. The referee was right to locate the VI.4 obstruction in the simultaneous right factors. Galois descent alone does not explain it.

## VII.2. What fails at multiplicity one

**PROVED (specialization of the original marks).** Take the original packed space \(V=i\,\mathbb F_q[t]_{<k}\), \(k\ge2\), \(q\equiv3\bmod4,\ q\ge27\), and \(m\ge16\) nonproportional forms with any signs. For **every** nonzero real specialization \(X=T(t)\in\mathbb F_q[t]\), the VI.4 marks cannot have squarefree common norm after division by any common polynomial. This assertion already holds without an infinity condition. It concerns the specified space and exponent formula.

*Proof.* Write \(b=k-1\), \(B=\max(b,\deg T)\) and extend constants to \(\mathbb F_{q^2}\). For an irreducible \(f\) there and positive \(s\) put
\[
    E_s(f)=\{u\in\mathbb F_q[t]_{\le b}: f^s\mid T-iu\}.
\]
The nonempty sets are nested affine \(\mathbb F_q\)-spaces. Call \(u\) mixed if the signed characters of the \(m\) forms at \(u\) include both signs. If \(E_s(f)=\{u\}\) for some \(s\) and \(u\) is mixed, the valuations of the raw marks at \(f\) have range at least two. In fact each nonempty level contains \(u\). The sum of a form's character on that level is zero when the form varies on its direction, and otherwise is the size of the level times the character at \(u\). For a positive and a negative form at \(u\), subtracting their contributions therefore gives a nonnegative difference on every level and two on the singleton level. The common parts of exponents \(1+\chi\) cancel in this comparison. Common division preserves the range. A squarefree norm would give a range of at most one, also after extending constants. Thus **every** occurrence of a mixed label in the divisor of \(T-iu\), counted at every valuation level, has to share that level with another label.

For uniformly chosen \(u\), indicators of a positive character for distinct forms are pairwise independent, each with mean \(a=(q-1)/(2q)\). Chebyshev bounds the probability of no positive sign by \((1-a)/(ma)\), and the same holds for a negative sign. Thus the fraction unmixed is at most \(2(q+1)/(m(q-1))\le4/m\). The fraction of labels of degree less than \(b\) is \(1/q\). For degree exactly \(b\), \(T-iu\) has degree \(B\), including when \(\deg T=b\), since \(i\notin\mathbb F_q\). Dividing all incidence counts by \(q^k\), the shared levels must consequently contribute degree at least
\[
                      (1-4/m-1/q) B.                               \tag{VII.3}
\]

For the upper bound group the \(f\)'s over an irreducible \(P\) of \(\mathbb F_q[t]\), of degree \(e\). The direction of \(E_s(f)\) consists of the polynomials of degree at most \(b\) divisible by \(P^s\), since \(f^s\mid u-u'\) for a real difference forces \(P^s\mid u-u'\). Shared levels therefore require \(es\le b\) and occupy at most the fraction \(q^{-es}\) of labels. For even \(e\), the two upstairs degrees sum to \(e\). Since \(e\,\pi_q(e)\le q^e\) by counting roots in \(\mathbb F_{q^e}\), all even-degree contributions together are at most
\[
       \lfloor b/2\rfloor+
          \sum_{\substack{e\ge2\\e\ {\rm even}}}
                 \frac{q^{-e}}{1-q^{-e}}
       \ \le\ b/2+\frac{q^{-2}}{(1-q^{-2})^2}.
\]
For odd \(e\), \(f=P\) stays irreducible. The divisibility \(P^s\mid T-iu\) implies \(P^s\mid T\) and \(P^s\mid u\) by taking real and imaginary parts. The total shared incidence here is at most \((\deg T)/q\), using \(q^{-es}\le1/q\) at each such level of \(T\). Together the bounds would imply
\[
     (1-4/m-1/q) B \le (1/2+1/q) B +q^{-2}/(1-q^{-2})^2 .
\]
This fails for \(m\ge16,\ q\ge27,\ B\ge1\). At \(T=0\) the raw marks vanish. If one first divides away the formal zero factor, every remaining linear factor evaluates to \(i\) times a real polynomial. All reduced marks lie on the same real line up to constant units, so common norm and a common leading coefficient give at most one distinct mark. \(\square\)

**PROVED (how far the cap-two examples are from squarefree).** In the unpadded VII.1 examples, the primes of exponent two account for at least the fraction \(1-2(q+1)/(m(q-1))\) of the parent degree, counting multiplicity. Indeed the same pairwise-independence calculation just used applies for each fixed orientation in VII.1. At any mixed label \(v\), the minimum exponent is zero both at \(v\) and at \(-v\). The product of \(T-v\) over all mixed labels is invariant under descent and reflection, is squarefree by VII.1 and has degree \(L\) times the number of mixed labels. Its square divides \(N_0\), whose degree is \(2(n-g)L\le2nL\). This proves the fraction bound. \(\square\)

**PROVED (cost of editing the signs on one copy).** In the formal variable \(X\), take two independent signed forms on any anti-invariant \(\mathbb F_q\)-space \(V\) as in VI.4. On any common symmetric subset of \(V\setminus\{0\}\), choose two squarefree selectors with exponents \(a_i(v),a_j(v)\in\{0,1\}\), one root from each reflected pair. Set both exponents to zero outside that subset. Define
\[
 E=\{v:2(a_i(v)-a_j(v))\ne\chi_i(v)-\chi_j(v)\},
 \qquad d_*=\frac{n(q-1)^2}{4q^2}.
\]
Here \(\chi_i\) is the signed character, including its zeros. Critical contact after any common division requires \(|E|\ge d_*/5\). In particular, editing only the union of all \(m\) kernels fails when \(mh<d_*/5\), as happens for the large \(q\) in VII.1 as \(m\) grows.

*Proof.* Let \(U_i,U_j\) be the selected products. The monic rational function \(Q=(U_i/U_j)^2/(F_i/F_j)\) is supported on \(E\). It is nonconstant since its exponent is odd on a label in exactly one of the two kernels. Its exponents sum to zero and each has absolute value at most four. Its reduced denominator degree is therefore at most \(2|E|\), which bounds the order of \(Q-1\) at infinity. There are \(d_*\) labels with signs \(+,-\). On every such label outside \(E\), \(U_i\) takes the root and \(U_j\) does not. Thus their degree after common division is at least \(d_*-|E|\). If they are critical, their ratio and its square have depth at least \((d_*-|E|)/2\). By the displayed formula for \(F_i\) in VII.1, \(F_i/F_j\) has depth at least \((n-h)/2\ge d_*/2\). Dividing these units gives \(2|E|\ge(d_*-|E|)/2\), as required. \(\square\)

**PROVED (arbitrary choices inside a single kernel).** On \(V\setminus\{0\}\), use the positive signed character off \(H_j=\ker\ell_j\), then fill the reflected pairs in \(H_j\setminus\{0\}\) arbitrarily. Two nonproportional forms cannot yield a critical pair for \(q\ge11\), including after common division. In fact \(U_j^2/(F_j/X)\) is the ratio of the selected kernel product to its conjugate. Its numerator and denominator have degree \((h-1)/2\), and its signed support is exactly \(H_j\setminus\{0\}\). Thus \(Q=(U_i/U_j)^2/(F_i/F_j)\) is nonconstant of denominator degree at most \(h-1\). On the \(d_*\) labels where the original signs are \(+,-\), the two selections differ, so critical contact would give \(Q\) depth at least \(\min(d_*/2,(n-h)/2)>h-1\) for \(q\ge11\). This contradicts the denominator bound. \(\square\)

**PROVED (a squarefree replacement on two conjugate translates).** There is also a formal obstruction which does not require the packed space. Let \(V\) be any \(k\)-dimensional anti-invariant \(\mathbb F_q\)-space in a field of odd characteristic with an involution fixing \(\mathbb F_q\), \(q\equiv3\bmod4,\ q\ge11\). Replace the two occurrences of each original root by two translates \(V+\tau,V-\tau\) with \(\bar\tau=-\tau\) and \(\tau\notin V\). Their union has squarefree monic root polynomial. On each translate and off \(\ker\ell_j\), select exactly the roots whose original label has positive signed character. On the two kernel translates fill the remaining conjugate pairs in **any** way, selecting one per pair. Then even two nonproportional forms cannot give critical contact in \(X\), after any common division.

*Proof.* Write \(h=|V|/q,\ n=qh,\ s=(q-1)/2\), \(A_j=P_{\ker\ell_j},\ b_j=A_j(\tau)\ne0\), and choose \(c_j\) as in VII.1. The kernel selections contribute a monic polynomial \(K_j\) of degree \(h\) with \(K_j\bar K_j=A_j^2-b_j^2\). The whole selected polynomial is
\[
 F_j=K_j\big((A_j-b_j)^s-c_j^s\big)\big((A_j+b_j)^s-c_j^s\big).
\]
At \(X=\infty\) the two linear terms in \(b_j/A_j\) cancel. All terms involving \(c_j^s/A_j^s\) are \(O(X^{-2h})\), and \(P_V/A_j^q=1+O(X^{-(q-1)h})\). Consequently
\[
                       F_j/P_V=K_j/A_j+O(X^{-2h}).
\]
For two independent forms there are \(n(q-1)^2/(4q^2)\) labels with prescribed signs \(+,-\). Both copies of each are selected on the first mark alone. Thus after common division the remaining degree is at least \(d_*=n(q-1)^2/(2q^2)\). Critical contact of the divided marks, which have equal leading coefficient, would give \(\operatorname{ord}_\infty((F_i-F_j)/P_V)\ge d_*/2\ge2h\). Hence \(K_iA_j-K_jA_i\) has degree at most zero. It vanishes at \(X=0\), since both \(A_i,A_j\) do, and is therefore identically zero. Taking norms of \(K_i/A_i=K_j/A_j\) gives \(b_i^2/A_i^2=b_j^2/A_j^2\). The \(b\)'s are nonzero and the \(A\)'s monic. This would give identical kernel polynomials. \(\square\)

**HEURISTIC (scope).** The specialization obstruction permits arbitrary extra common zeros at finite places and prescribes the original factors and their exponents. The editing bound measures the number of departures from those signs in the formal variable. The one- and two-translate bounds permit arbitrary kernel choices and prescribe the signs off the kernels. General squarefree ternary clusters need not carry these signs.

## VII.3. Exact squarefree and ramified tables

**COMPUTED.** For the constant ternary extension, with \(D\) the degree of a mark, I obtain the following exact maxima.

| \(D\) | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| \(M_3^{\rm sf}(D,0)\) | 2 | 4 | 3 | 4 | 3 | 4 | 4 | 4 | 4 | 4 | 4 | **5** | 4 | **5** | 4 | **5** |

**COMPUTED.** A squarefree parent at the first degree admitting five marks is
\[
                  N=t^{24}+t^{20}+2t^{18}+t^{16}+2t^{10}+2t^2+2 .
\]
The ten pair differences have degrees \(5,6,6,6,6,6,5,5,6,6\). A separate FLINT script checks all sixteen witnesses, including equality of norms, squarefreeness, distinctness and critical contact. To reproduce the upper bounds, compile the C++ source with C++17 and run it with argument 16. [Census source](vii_squarefree_table.cpp) · [Witnesses and enumeration counts](vii_squarefree_table.out) · [Independent witness check](vii_squarefree_witness.py)

**PROVED (why the census is exhaustive).** Normalize the common leading coefficient to one. In this constant extension every irreducible of a squarefree norm splits since an inert factor has even norm valuation. A mark takes exactly one member of each conjugate pair. In five marks, each such binary cut is constant on at least one of the ten triples. Weight the cuts by the degree of one upstairs factor. Some triple therefore has a common block \(G\) of degree at least \(\lceil D/10\rceil\). On the remaining cuts write the triple as
\[
          z_i=G U_i\prod_{j\ne i}\bar U_j,\qquad i=1,2,3.
\]
Here \(U_i\) chooses the factor of the singleton on each cut where \(i\) is the exceptional row. The four blocks have squarefree, mutually coprime norms. With \(K=\lceil D/2\rceil\), the critical width condition is \(\operatorname{ord}_\infty(z_i/z_j-1)\ge K\). In the displayed factorization it says that the three units \(U_i/\bar U_i\) agree to depth \(K\). Writing \(U=A+iB\) monic, this is equivalent to equality of the slopes \(B/A\) modulo \(t^{-K}\), by the invertible Cayley transform in odd characteristic. The nonzero numerator \(U_i\bar U_j-U_j\bar U_i\) also shows that \(\deg U_i+\deg U_j\ge K\). Thus
\[
 \max_i\deg U_i\le D-\lceil D/10\rceil-K,\qquad
 \deg G\le D-\lceil3K/2\rceil.
\]
These bounds are at most six for the range where the search uses them.

The program enumerates the blocks by degree, truncated slope and disjoint norm support. For each triple and common block it counts every choice of the individual conjugate factors on that parent matching the triple's phase. Slopes add under multiplication by \((a+b)/(1-ab)\). Switching a factor changes its phase by twice its original phase in this group, with sign depending on the original choice. Doubling is invertible because the truncated group has odd order \(3^{K-1}\). Thus subset sums of *single* phases suffice to count switches. Matching phases at equal monic norm gives matching prefixes. Indeed the ratio \(u\) of the two marks has norm one and \(u/\bar u=u^2\), so \(u^2-1\) and \(u-1\) have the same order. Every cluster of at least five is counted in its own phase class by one of the triples above. A four-mark witness then suffices to certify a maximum of four if the search finds no larger class. The code uses all parents at \(D\le5\) and an unrestricted triple fallback if needed to establish the smaller lower bounds. Ordinary polynomial multiplication independently groups the marks by their actual leading coefficients each time the recorded maximum increases. \(\square\)

**PROVED.** The five-mark lower bound recurs at every \(D=12L,\ L\ge1\). Substitute \(t^L\) if \(3\nmid L\) and \(t^L+t\) if \(3\mid L\). In the first case the only possible finite critical value is zero and \(N(0)=2\); in the second case the derivative is one. The parent stays squarefree and all nonzero polynomial degrees scale by \(L\).

**COMPUTED.** For the ramified order \(\mathbb F_3[t,\eta]\), \(\eta^2=2t^2+1\), with all nonzero leading coefficients of \(N\) allowed, the exact maxima are as follows. Here \(D=\max(\deg A,1+\deg B)\) for \(A+B\eta\).

| \(D\) | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| all parents | 2 | 3 | 2 | 4 | 3 | 4 | 3 | 4 | 4 | 4 |
| squarefree parents | 1 | 2 | 2 | 3 | 3 | 4 | 3 | 4 | 3 | 4 |

**PROVED (coverage of the ramified computation).** A cluster has one common pair of coefficient strings above pole order \(\lfloor D/2\rfloor\). Within each pair of strings the program enumerates all low coefficients and sorts the exact norms \(A^2-(2t^2+1)B^2\). All distinct pairs in a norm group give distinct marks with the required contact. Independent sign changes of the entire \(A\)- or \(B\)-coordinate preserve every test, so the leading coefficient pairs \((1,0),(0,1),(1,1)\) cover all cases. The last type has parent leading coefficient two and is included. The optional squarefree filter is Euclidean gcd with the derivative. Run the script with argument 10. [Ramified census](vii_ramified_check.py)

## VII.4. A single incomplete binary Hadamard matrix

**PROVED.** Take a Hadamard matrix of order \(m\ge8\) with an all-ones column. Retain all rows and a proper subset of the other \(m-1\) columns. At each retained column let a mark choose the monic polynomial \(f_v\in\mathbb F_4[t]\) or its conjugate by the sign in that column. Assume \(f_v\) is nonconstant and nonreal, and that different columns have disjoint base-prime supports, as in VI.2. Then the \(m\) products cannot be distinct marks of a critical cluster, at any degrees. One may first remove common factors, including common parts of a block and its conjugate. The statement applies to Walsh matrices without a linearity assumption on the polynomials.

*Proof.* Work with the reduced blocks and write \(n\) for their number, \(d_v=\deg f_v>0,\ D=\sum_v d_v,\ K=\lceil D/2\rceil,\ r_v=f_v/\bar f_v\). Deleting any now constant blocks still leaves a missing column. The cluster condition is \(z_i/z_j\in U_K=1+t^{-K}\mathbb F_4[[1/t]]\). Thus any product of the marks with integral exponents of sum zero lies in \(U_K\). Let \(w\) be a missing nonconstant column. For each retained \(v\), orthogonality and balance give
\[
               \prod_i z_i^{(\chi_v(i)-\chi_w(i))/2}=r_v^{m/4}.
\]
In fact the signed exponents are integral, their sum vanishes, and summing them on the positive rows in any retained column gives \(m/4\) at \(v\) and zero elsewhere. Here \(4\mid m\) by orthogonality of the original sign columns. Write \(u\) for the odd part of \(m\). The depth \(\operatorname{ord}_\infty(r-1)\) of a nonidentity principal unit is multiplied by \(2^a\) on raising to \(2^a b\) with \(b\) odd. This follows by Frobenius and the linear term of the odd power. Thus
\[
                    d_v\ \ge\ \operatorname{ord}_\infty(r_v-1)
                      \ \ge\ (4u/m)K\ \ge\ 2uD/m.
\]
The first inequality uses the nonzero polynomial \(f_v-\bar f_v\). Summing gives \(n\le m/(2u)\).

Form the real \(m\times n\) matrix \(S_{iv}=\chi_v(i)\sqrt{d_v}\). The Gram matrix \(SS^\mathsf T\) has diagonal \(D\) and zero row sums. Its off-diagonal entries are nonpositive. Indeed the product of common choices on a pair divides \(z_i-z_j\), hence their total degree is at most \(D/2\), and the entry is twice that degree minus \(D\). The Gram matrix is a weighted graph Laplacian. Its nullity is the number of components, since its quadratic form is the sum over edges of the positive edge weights times squared coordinate differences. Each vertex has positive diagonal and thus each component has size at least two. Its rank is at least \(m/2\). Its rank is also exactly \(n\), by orthogonality of the retained columns. We get \(u=1,\ n=m/2\), all \(d_v=d=D/n\), and every component is a pair. On each such pair the off-diagonal entry is \(-D\), so the two rows have opposite signs in every retained column.

Fix a retained column \(w_0\) and choose from each pair the row positive at \(w_0\). The restricted \(n\times n\) sign matrix is Hadamard because any product of two column signs is identical on the two rows of each old pair, and hence its sum on the chosen rows is half the old sum. On these rows \(w_0\) is constant. The first depth inequality was sharp, so \(f_v+\bar f_v=1\) at every label. Write \(f_v=A_v+\theta\) with \(A_v\in\mathbb F_2[t]\) and \(\theta^2+\theta+1=0\). For \(v,v'\ne w_0\), distinct, the exponents \((\chi_v-\chi_{v'})/2\) on the restricted matrix now give \((r_v/r_{v'})^{n/4}\in U_K\). Importantly the original depth \(K\) is still valid on the subset of marks. As \(n\) is a power of two and \(n\ge4\), this gives depth at least \(4K/n\ge2d\) for \(r_v/r_{v'}\) (infinite if equal). Its numerator after subtracting one, before cancelling, is \(A_v+A_{v'}\) and its displayed denominator has degree \(2d\). Thus \(\deg(A_v+A_{v'})\le0\). At most two such \(A_v\) can be distinct, whereas the \(n-1\ge3\) blocks are distinct by disjointness. This proves the contradiction. \(\square\)

**PROVED (scope for the binary growth question).** Combined with VI.2, this leaves the full \(m-1\) columns as the only possible support for a single complete set of \(m\ge8\) Hadamard rows, and costs \(\Omega(m^2\log m)\) there. The Gram-matrix step uses pairwise divisibility and real linear algebra. The other step uses the ultrametric depth of a product, exact doubling of depth under squaring in characteristic two, and the two available constants in \(\mathbb F_2\). It supplies no corresponding integer bound. A sparse set of rows inside a larger Walsh matrix, or two Walsh matrices on differently permuted rows, need not satisfy the hypotheses of the proof.

**COMPUTED (the first two-core test).** For eight rows labelled by \(x,x^3\) as in VI.1, the two-core profile with one distinct split prime of common degree \(e\) per column is impossible at \(e=10,12,14,16\). Ten was the first degree the previous test permitted. The new calculation tests all assignments, even before choosing any factor above a base prime. It does not settle \(e=4d-2\) for general odd \(d\). [Exact binary search](vii_double_core.cpp) · [Counts](vii_double_core.out)

| \(e\) | available primes | admissible four-prime sets | two-dimensional subspaces | three-dimensional subspaces |
|---:|---:|---:|---:|---:|
| 10 | 99 | 67 | 10 | 0 |
| 12 | 335 | 319 | 34 | 0 |
| 14 | 1161 | 4382 | 226 | 0 |
| 16 | 4080 | 40546 | 1453 | 0 |

**PROVED (criterion for the computation).** By VI.1 the isolated words of the first Walsh core belong to the span of the two-core rows. Their seven nonzero supports are four-prime sets forming a three-dimensional binary space under symmetric difference. For each product \(H\) on such a support, \(H'/H\) must vanish to order at least \(7e/2+1\), by logarithmically differentiating a product of ratios in \(U_{7e/2}\). Thus \(\deg H'\le e/2-1\). To generate all eligible sets, write two disjoint monic pair products as \(P=A^2+tB^2,\ Q=C^2+tE^2\). Then \(\sqrt{(PQ)'}=AE+BC\) and \(\deg A=\deg C=e\). With \(b=\lfloor(e/2-1)/2\rfloor\), the derivative bound is equivalent to
\[
 \left\lfloor t^{2e-b-1}B/A\right\rfloor
       =\left\lfloor t^{2e-b-1}E/C\right\rfloor ,
\]
where the floor means the polynomial part at infinity. Indeed the difference has order strictly positive exactly when \(\deg(AE+BC)\le b\). The program sorts every pair product of distinct degree-\(e\) primes by this key. Each eligible four-set is obtained from precisely its three pairings. Finally the program lists all lines \(\{a,b,a+b\}\) in the eligible set system and tests each outside \(u\) for the three further sums \(u+a,u+b,u+a+b\). This exhausts the three-dimensional subspaces. At \(e=10\) the script also tests all \(\binom{99}{4}=3{,}764{,}376\) products by directly forming the polynomial derivative. \(\square\)

## Twelve-line summary of Part VII

1. **PROVED.** Over \(\mathbb F_3[t]\), unbounded critical clusters exist for every fixed parent multiplicity cap \(e_0\ge2\).
2. **PROVED.** This extends to every fixed \(p\equiv3\bmod4\) and every ramified quadratic model with inert infinity specified in VII.1.
3. **PROVED.** The involution in the ramified example fixes the constants, and the unpadded parent avoids every branch prime.
4. **PROVED.** With cap two the maximum at every sufficiently large degree \(D\) is at least \((\tfrac12\log_2p-o(1))\log D/\log\log D\).
5. **COMPUTED.** In both models a degree-619801 example has four marks, with a degree-709317 squarefree polynomial whose square contains their norm.
6. **COMPUTED.** For squarefree ternary parents the maximum first reaches five at \(D=12\), is five again at 14 and 16, and is four at 13 and 15.
7. **PROVED.** No specialization of the original packed VI.4 family with \(m\ge16,\ q\ge27\) gives a squarefree cluster after common division.
8. **PROVED.** Two reflected translates with arbitrary choices on the kernel fibers cannot split the formal double factors at critical contact for \(q\ge11\).
9. **PROVED.** A binary Hadamard matrix with all \(m\ge8\) rows and any missing nonconstant column cannot give a critical cluster.
10. **COMPUTED.** The eight-row \(x,x^3\) binary profile is impossible for each even prime degree from 10 through 16.
11. **HEURISTIC (status).** I could not decide squarefree odd-characteristic uniformity, unrestricted binary growth, or the integer growth rate.
12. **HEURISTIC (request).** Please referee the finite-place separation in VII.1 first. I would next ask whether cap-two clusters remain unbounded when a fixed positive fraction of parent degree comes from primes occurring once.

# Part VIII. Low-height primes on integer circles

**PROVED (choice and result).** I chose (c). The low-height window contains nonarithmetic fakes with \(P_+(N')\le m^{23/12+o(1)}\) and \(W\le(23/3+o(1))m\log m\). One can simultaneously impose all monomial gaps, point residues through \(2m\), and separation \(\gg1/m\) of the individual prime angles. The construction and the distinction between point residues and prime-factor residues are in VIII.3–4.

On genuine integer circles I get two restrictions, uniform in the parent. The first is
\[
                    M\le 1+\omega(N')+K_C.                         \tag{VIII.1}
\]
Here \(N'\) is the norm after common Gaussian gcd division and \(K_C\) depends only on the arc constant. In particular, on circles with four-full norm the coefficient \(2/3\) can be replaced by \(1/2\). The second restriction, (VIII.6) below, rules out a sequence of arbitrarily large clusters throughout the low-height window whenever even one coordinate has finite order \(\exp(o(M))\) in the quotient by the row lattice. The rank estimate is the conclusion that round 20 obtained conditionally from a moving-line height comparison. The proof here needs only finitely many inert residue fields. [Round 20, §§3 and 6](</home/user/Jarnik/docs/claude_review/round20/wild.md>) · [Round 21, §§0–1](</home/user/Jarnik/docs/claude_review/round21/README.md>)

## VIII.1. Recovering affine rank at fixed inert primes

**PROVED.** Fix \(C>0\). Choose any finite set \(S\) of rational primes \(\ell\equiv3\pmod4\) with
\[
            \sum_{\ell\in S}\frac{\log\ell}{\ell+1}>\log C,
               \qquad K_C=\sum_{\ell\in S}\ell .                  \tag{VIII.2}
\]
Such sets exist by the prime number theorem in the progression \(3\bmod4\) and partial summation. For \(C<1\) one may use the empty set. For an arc of length at most \(C\sqrt R\) on \(|z|=R\), divide the marks by their common Gaussian gcd. Write \(N'=\prod p^{e_p}\) for the resulting norm, \(W=\log N'\), \(r=\omega(N')\). Only split rational primes remain. Fix a choice \(\pi_p\mid p\), and let \(a_{ip}=v_{\pi_p}(z_i)\) on the normalized marks. The \(M\) rows \(a_i\) have affine rank at least \(M-1-K_C\). More precisely, appending the residue-class indicator columns of the normalized marks at the primes of \(S\) makes them affinely independent.

*Proof.* Common division is a similarity, and the new arc constant is at most \(C\). At a common mark norm \(N'\) the other valuation is \(e_p-a_{ip}\). Thus for \(i\ne j\), putting \(g_{ij}=\gcd(z_i,z_j)\), \(b_{ij}=(z_i-z_j)/g_{ij}\), and \(d_{ij}=\sum_p|a_{ip}-a_{jp}|\log p\), we have
\[
 \log|g_{ij}|=(W-d_{ij})/2,\qquad
 d_{ij}\ge W/2+2\log|b_{ij}|-2\log C.                            \tag{VIII.3}
\]
These statements include the individual Gaussian units of the marks. For each \(\ell\in S\) the residues lie in one nonzero norm fibre in \(\mathbb F_{\ell^2}\). It has \(\ell+1\) elements by cyclicity. Append the coordinate vector of the residue in this fibre to each row.

An affine dependence of the enlarged rows would give probability measures \(P,Q\) on disjoint nonempty sets of marks with the same barycentre. Take \(i\sim P, j\sim Q\) independently. For each \(p\), the variables \(X=a_{ip}/e_p,Y=a_{jp}/e_p\) have the same mean \(\mu_p\). The inequality \(|X-Y|\le X+Y-2XY\) on the unit square gives \(\mathbb E|X-Y|\le 2\mu_p(1-\mu_p)\le1/2\). Therefore \(\mathbb E d_{ij}\le W/2\). At every \(\ell\in S\) the two residue distributions are also equal. Their collision probability is the sum of their squared masses, at least \(1/(\ell+1)\). On a collision \(\ell\mid b_{ij}\), since \(\ell\nmid N'\). For each pair the colliding rational primes divide \(b_{ij}\) simultaneously, and the remaining Gaussian integer has modulus at least one. Consequently \(\mathbb E\log|b_{ij}|\ge\sum_{\ell\in S}\log\ell/(\ell+1)>\log C\), contradicting (VIII.3). The affine dimension contributed by each block of indicators is at most \(\ell\). This proves the rank assertion and (VIII.1). \(\square\)

**PROVED (consequences for parents).** Write \(Y=P_+(N')\) when \(N'>1\). Then
\[
 M\le\pi(Y;4,1)+1+K_C,\qquad Y\ge(2-o_C(1))M\log M\quad(M\to\infty).
                                                                    \tag{VIII.4}
\]
In particular for each \(\eta>0\) the size is uniformly bounded in the subwindow \(Y\le(2-\eta)M\log M\). Suppose each prime active after normalization occurs to exponent at least \(k\) in the **original** norm \(N=R^2\). For fixed \(k\) and \(C\), along any such sequence with \(M\to\infty\),
\[
 \log N\ge(k-o_C(1))M\log M,\qquad
 M\le (2/k+o_{C,k}(1))\frac{\log R}{\log\log R}.                  \tag{VIII.5}
\]
Indeed the \(r\) active primes are distinct, so their log sum is at least \(\log((r+1)!)\). Use (VIII.1), Stirling's estimate, and inversion of \(x\log x\). Formula (VIII.4) follows by inverting the prime number theorem \(\pi(y;4,1)\sim y/(2\log y)\). In particular (VIII.5) holds whenever every prime divisor of \(N\) has exponent at least \(k\). Its improvement on \(2/3\) starts at \(k=4\). No claim about a four-full *normalized* norm is needed.

**PROVED (field test).** The rank argument works in both fixed-field constant quadratic extensions. Replace the real absolute value by \(|f|=p^{\deg f}\), and \(\log\ell\) by \((\deg\lambda)\log p\) at an inert irreducible \(\lambda\). The odd-degree base irreducibles are inert, with norm fibre of size \(p^{\deg\lambda}+1\). Their contributions in (VIII.2) diverge because at each odd degree \(d\) the number is \((p^d+O_p(p^{d/2}))/d\), by the formula for irreducible polynomials. All remaining inequalities in the proof are unchanged for the corresponding chord bound. On the squarefree binary family VI.5, the \(M-1\) nonconstant Walsh columns are mutually orthogonal and are present on disjoint nonempty factor blocks. The affine rank is exactly \(M-1\). On the separated ternary family VII.1 the column on a root label, after subtracting its constant part, has entries proportional to \(\chi(\ell_i(v))\) up to the fixed mark signs. For distinct nonproportional forms the pair of values is uniform on \(\mathbb F_q^2\) as the label varies. The nonzero character rows are therefore mutually orthogonal, since the character extended by zero has mean zero. Common division translates columns, preserving affine rank. Both constructions satisfy the rank bound.

## VIII.2. A uniform height and order inequality

**PROVED.** Use the round-21 normalization \(B_{ip}=2(a_{ip}-a_{1p})\), \(2\le i\le M\), and put \(L=\operatorname{row}_{\mathbb Z}B\subseteq\mathbb Z^r\). For \(v\ne0\) put \(w(v)=\sum_p |v_p|\log p\). If \(v\) has finite order \(D_L(v)\) modulo \(L\), then on every genuine integer cluster with \(M\ge M_0(C)\),
\[
             (1+w(v))\log(2D_L(v))\ \ge\ c_* W                    \tag{VIII.6}
\]
for an absolute effective \(c_*>0\). Finite order includes order one. In particular each finite-order prime coordinate in \(Y\le M^{2-\delta}\), for fixed \(0<\delta<1\), has order at least \(\exp(c_\delta M)\) for sufficiently large \(M\). This last assertion follows from \(W\ge3M\log M-O_C(M)\) and \(w(e_p)=\log p\). It applies even if all other coordinates are free. Thus for a specified subexponential upper envelope on one finite order the cluster size is uniformly bounded in the window. [Growth bound, Theorem 7.5](</home/user/Jarnik/docs/claude_review/two_thirds/proof.md>)

*Proof.* Write \(t_p=2\arg(\pi_p)/\pi\). The arc gives an integer vector \(k\) and a real error vector \(u\) with \(Bt=k+u\) and \(\|u\|_\infty\le (2C/\pi)e^{-W/4}\). This uses angular differences from the first mark, so its common phase may drift arbitrarily. Let \(s=\operatorname{rank} B\) and let \(H\) be the maximum absolute value of a square minor of \(B\), including the empty minor. Column \(p\) of a \(j\)-minor has Euclidean length at most \(2\sqrt j\,e_p\). Hadamard and \(\sum e_p\le W/\log5\) give
\[
 |\det|\le j^{j/2}\left(\frac{2W}{j\log5}\right)^j.
\]
Uniformly for \(j\le M-1\), the logarithm of the first factor is at most \(M\log M/2\). For \(M\) large the logarithm of the second is increasing in \(j\le M\), since \(W/M\ge 3\log M-O_C(1)\). Its maximum is \(o_C(W)\), as \(\log x/x\to0\) at \(x=W/M\). Hence
\[
                    \log H\le(1/6+o_C(1))W.                       \tag{VIII.7}
\]
Choose a nonsingular \(s\)-minor \(T=B_{I,J}\), put \(d=|\det T|\), and correct only the coordinates \(J\) of \(t\), by subtracting \(T^{-1}u_I\), to obtain \(t^*\). If \(s=0\) there is no nonzero finite-order vector. Otherwise \(B_I t^*=k_I\). For any row \(b\) of \(B\), its coordinates in the basis of rows \(I\) are \(b_J T^{-1}\). Multiplying these coordinates by \(d\) gives integers of magnitude at most \(H\) by the row-replacement minor formula. Thus \(d(bt^*-k_b)\) is an integer of magnitude at most \((s+1)H\|u\|_\infty<1\) once \(M\) is sufficiently large. This proves \(Bt^*=k\) on *all* rows, including rows outside the minor.

Finite order is equivalent to membership in the rational row space. Applying the same row-basis formula to \(v\) shows \(dv\in L\), so \(D=D_L(v)\) divides \(d\), and \(D v\cdot t^*\in\mathbb Z\). Every adjugate entry of \(T\) is bounded by \(H\). It follows that
\[
  \|Dv\cdot t\|_{\mathbb R/\mathbb Z}\le (D/d)s\|v\|_1H\|u\|_\infty .
\]
For \(w(v)\le W\) this is at most \(\exp(-(1/12-o_C(1))W)\), since \(\|v\|_1\le W/\log 5\). For \(w(v)>W\), (VIII.6) is immediate upon decreasing \(c_*\).

Set \(\alpha=\prod_p(\pi_p/\bar\pi_p)^{v_p}\). Its numerator and denominator in lowest terms have norm \(\exp(w(v))\), its absolute Weil height is \(w(v)/2\) in logarithmic normalization, and it is nontorsion by unique factorization. Moreover \(\alpha^{2D}=\exp(2\pi iD v\cdot t)\). For any nontorsion \(\alpha\in\mathbb Q(i)\) on the unit circle and any integer \(n\ge1\), the two-logarithm Baker–Wüstholz bound gives
\[
 \log|\alpha^n-1|\ge -K(1+h(\alpha))\log(2n)
\]
with absolute effective \(K\). To see the uniformity, use the principal logarithm of \(\alpha\) and \(\log(-1)=i\pi\) and subtract \(2j\log(-1)\) from \(n\log\alpha\) for the nearest winding \(j\). The form is nonzero, \(|j|\le n\), the two logarithms have bounded imaginary parts, and the field degree is two. Distance to one bounds the reduced form from below up to the absolute factor \(2/\pi\). This is a fixed two-logarithm application even when the number of active primes grows. Combining it with the preceding upper estimate proves (VIII.6). \(\square\) [Baker and Wüstholz, *Logarithmic forms and group varieties*](https://doi.org/10.1515/crll.1993.442.19)

**PROVED (field failure, with an explicit contradiction).** The parameterization by a real angle modulo a fixed period, the use of its integer windings, and the complex logarithm bound are the steps without the required local-field analogues. Indeed (VIII.6) itself fails on VI.5. Take one nonconstant Walsh factor block and let \(v\) be its signed indicator, orienting each of its prime factors consistently with the column. It has \(w(v)=W/(M-1)\) by equal block degrees and squarefreeness. In \(B\), the entries on any block are its signed Walsh column minus the value of that column at the first row. Summing the row differences against the chosen nonconstant Walsh character gives \(Mv\), by orthogonality and zero mean. Thus \(D_L(v)\mid M\), and the left side of (VIII.6) divided by \(W\) tends to zero. In either positive characteristic, raising a nontrivial principal unit to \(p^a\) multiplies its depth by \(p^a\), which precludes a two-logarithm estimate of the displayed form with only logarithmic cost in the exponent. This test does not use an assertion about the squarefree odd-characteristic problem.

## VIII.3. A fake below the quadratic prime cutoff

**PROVED.** For arbitrarily large \(m\) there is a normalized angle model with \(m\) distinct marks on a critical arc of constant one, with
\[
 P_+(N')\le m^{1.995},\qquad
 W=\log N'\le(7.96+o(1))m\log m.                                  \tag{VIII.8}
\]
It satisfies both sharp monomial gaps for every integral exponent vector, including the rational characters of the profile. It also admits compatible point residues satisfying round 5's axiom (R) at every odd prime power up to \(2m\). The successful angles have positive Haar measure on a fixed component of a realization coset. In particular one can choose nonarithmetic angles. This answers round 21's question (Q) for a fixed power saving in the prime cutoff. The construction uses a sparse interval of primes and an exponent cap of two. It gives an upper bound on the weight of a **fake**. Equation (VIII.1) remains the result above that applies to every genuine circle. [Round 21, §6.3](</home/user/Jarnik/docs/claude_review/round21/fakes.md>) · [Point residue axiom, §1.2](</home/user/Jarnik/docs/claude_review/round5/construct.md>)

*Proof.* A profile has point phases \(\langle 2a_i-e,\phi\rangle\) up to an individual Gaussian unit. Its two gaps at a nonzero integer \(v\) require distances at least \(\arcsin(e^{-w(v)/2})\) from \((\pi/2)\mathbb Z\), and \(\arcsin(e^{-w(v)/2}/\sqrt2)\) from \(\pi/4+(\pi/2)\mathbb Z\). We will use a Legendre matrix of order \(q\) with slightly fewer than half its rows. Its integer kernel supplies vectors congruent modulo \(q\) to the power sums of the column labels. Newton's identities will give the first usable power sum for each short monomial. For longer monomials it is essential to count how many can have a long initial string of zero power sums, and then to bound the probability that their integer lifts also vanish.

Use \(q\to\infty\) through primes \(3\bmod4\), \(d=1.99,\ P=q^d,\ H=\lfloor q^{.08}\rfloor,\ A=q^{.12},\ B=\lceil q^{1.708}\rceil,\ J=\lfloor q/4\rfloor\). Turn angles \(\theta_j=\phi_j/(\pi/2)\) take values modulo one. In this coordinate it suffices for each nonzero \(v\), \(h=\|v\|_1\), to avoid intervals about \(0,1/2\) of radii \(2\exp[-w(v)/2]\), which is stronger than both required arcsine bounds for large \(q\). Except for pairs of rows we upper bound each radius by \(2P^{-h/2}\).

### The matrix and primes

Let \(K_{ij}=\chi(j-i)\), \(i,j\in\mathbb F_q\), with the quadratic character extended by zero. In this proof \(K_S\) denotes the submatrix on rows \(S\), and products restricted to columns \(I\) are interpreted by extending the other vector by zeros. Its row and column sums vanish, \(K^T=-K\), and \(KK^T=K^TK=qI-\mathbf1\mathbf1^T\). For the inner products off the diagonal, the map \(j\mapsto(j-i)/(j-i')\) takes all values but one as \(j\ne i'\) varies, giving \(-1\). Remove the column indexed by zero, writing \(D=\{0\}\) and \(I=\mathbb F_q\setminus D\). The class \(\mathcal A\) of rows taking the value \(+1\) on that column has size \((q-1)/2\). Fix any real \(0.99<\rho<1\). The marks will be a subset \(S\subset\mathcal A\) of size \(m=(\rho/2+o(1))q\). Set \(a_{ij}=1+K_{ij}, e_j=2\) on \(I\), so \(2a-e=2K\). Removing this column preserves every difference of selected rows. For any nonidentical \(i,i'\) the \(\ell^1\) distance of their full \(K\)-rows is \(q+1\). Indeed, there are \((q-1)/2\) opposite signs where both are nonzero and two singleton zeros.

Take \(q-1\) distinct split primes in some interval \([P_0,P_0 e^\sigma]\subset[P,2P]\) with \(\sigma\ll q\log(q)/P\), by the prime number theorem modulo four and pigeonholing that dyadic interval. Randomly permute them on \(I\). Simultaneously for every pair of rows of \(\mathcal A\), one can achieve
\[
 w(K_i-K_{i'})-\tfrac12 W \ \ge\ 2\log P-o(1).
\]
Indeed, the coefficients \(|K_{ij}-K_{i'j}|-1\) on \(I\) lie in \([-1,1]\) and sum to two. The weighted sum has mean twice the average log prime. Its deviation exceeds \(C\sigma\sqrt{q\log q}\) with probability at most \(2q^{-4}\) for an absolute \(C\). Under sequential exposure, swapping two remaining values changes the total by at most \(2\sigma\), giving the bounded-difference martingale inequality. The deviation bound tends to zero. This works also when \(S\) is chosen later. In particular \(W=(2d+o(1))q\log q\) before any further normalization.

### An uncertainty estimate for congruences

For primitive integral \(v\) on the \(q\) indices with \(\sum v_j=0\), if \(h=\|v\|_1\le q^{.10}\), every odd prime \(\ell\le h\) has at least \(\sqrt q-h-1\) nonzero coordinates of \(Kv\bmod\ell\), and \(Kv\bmod4\) has at least this many. Here "primitive" means the gcd of all coordinates is one.

Take the Paley skew conference matrix \(C\) of order \(q+1\) whose finite part is \(K\), \(C_{i\infty}=1\), \(C_{\infty i}=-1\), and diagonal zero at infinity. Augment \(v\) by zero there. We have \(C^2=-qI\), \(Cv\equiv v\bmod2\) since \(\sum v=0\). Put \(g=\sqrt{-q}\). As \((g-1)/2\) is an algebraic integer, the two eigenvectors \(y_\pm=(Cv\pm g v)/2\) are integral over \(\mathbb Z\), and their difference is \(gv\). Modulo any prime of \(\mathbb Q(g)\) above \(\ell\ne q\), at least one reduces nontrivially. For odd \(\ell\) its reduced support lies in the union of those of \(v\) and \(Cv\bmod\ell\), and for \(\ell=2\) use \(Cv\bmod4\).

Every such nonzero reduction of an integral eigenvector has at least \(\sqrt q\) nonzero coordinates. In fact use \(s_i=(1,i),s_\infty=(0,1)\), so \(C_{ij}=\chi(\det(s_i,s_j))\). Each element of \(\operatorname{SL}_2(\mathbb F_q)\) permutes the projective indices and rescales their representatives. Its permutation with these scalar quadratic characters as signs commutes with \(C\). Transitivity lets us move a unit coordinate of the eigenvector to infinity. If \(f\) is its finite part and \(\lambda=\pm g\) its eigenvalue, then \(-\sum f=\lambda y_\infty\) is a unit. At nonzero additive Fourier frequencies on \(\mathbb F_q\), the equation \(Kf+\mathbf1 y_\infty=\lambda f\) puts the Fourier support of \(f\) on just one of the two quadratic classes because the convolution eigenvalues of \(K\) are proportional to \(\chi(\xi)\) by a change of variables, and their squares equal \(-q\) by \(K^2=J_{\rm ones}-qI\). These assertions are over characteristic zero. Convolving \(f\) with its argument dilation by a nonsquare gives the constant \((\sum f)^2/q\). This identity in the original number field reduces modulo the chosen prime, and the constant stays nonzero. A cyclic convolution of two vectors each with support \(b\) can be nonzero at at most \(b^2\) places. Thus \(b^2\ge q\), proving the claim also in residue characteristic two.

Independently select each member of \(\mathcal A\) into \(S\) with probability \(\rho\), and no other indices. With probability tending to one we have \(\big||S|-\rho(q-1)/2\big|\le q^{2/3}\) and, simultaneously for all primitive balanced \(v\) of \(h\le q^{.10}\),
\[
  \gcd((Kv)_{S^c})\le2,\qquad
  \|(Kv)_{S^c}\|_2^2\ge (1-\rho)q\|v\|_2^2/2 .
\]
For the gcd test, failure means all coordinates detecting some odd \(\ell\le h\) or detecting \(4\) were sampled. For each obstruction this has probability at most \(\rho^{\sqrt q-h-1}\). The number of \(v\)'s is at most \((2q+1)^{q^{.10}+1}\). For the norm test, the total squared norm is \(q\|v\|_2^2\), and each summand at most \(h^2\). The expected retained squared norm on \(S^c\) is at least \((1-\rho)q\|v\|_2^2\). The multiplicative lower-tail bound for independent variables in \([0,1]\), applied after dividing by \(h^2\), gives failure probability at most \(\exp(-c_\rho q\|v\|_2^2/h^2)\le\exp(-c_\rho q/h)\), which also admits this union bound. Fix such an \(S\). Moreover for every nonzero \(v\bmod q\) supported on fewer than \(q/2\) indices, \(Kv\bmod q\) has at least \((q+1)/2\) nonzero coordinates, hence at least \(c_\rho q\) on \(S^c\) for a fixed \(c_\rho>0\). Indeed \((Kv)_i=\sum v_j(j-i)^{(q-1)/2}\) modulo \(q\), a polynomial of degree \((q-1)/2\). It is nonzero by the Vandermonde matrix in the distinct labels, since the binomial coefficients in that power are all units. A zero polynomial would force all moments through degree \((q-1)/2\) to vanish.

### Forced characters, offsets and point residues

In calculating characters extend any input \(v\) on \(I\) by zero on \(D\). Every zero-sum rational row coefficient \(\lambda\) on \(S\) with \(v=K_S^T\lambda\) integral on \(I\) is integral. Indeed its image on \(D\) vanishes by constancy, \(q\lambda=Kv\) on all indices (extend \(\lambda\) by zero), and \(K^T(q\lambda)=0\bmod q\) forces \(q\lambda=0\bmod q\) by the preceding polynomial observation. Thus this covers the whole rational span of row differences. Put \(s=\|\lambda\|_1\). Apart from \(s=2\), the pairs, all nonzero characters have \(\|v\|_1\ge1.2q\). For \(\lambda\) primitive modulo \(q\) in the sense of having a nonzero entry there, \(v\) has at least \(q/2\) nonzero entries. On each use \(y^2\le(s+1)y-s\) for integer \(1\le y=|v_j|\le s\). Then
\[
 (s+1)\|v\|_1\ge q\|\lambda\|_2^2+s q/2\ge 3qs/2 .
\]
For \(s\ge4\) this suffices. If all entries of \(\lambda\) are multiples of \(q\), divide out their common \(q\)-power first.

We next choose offsets that preserve the point residue condition as well. This is where the order of *auxiliary* genuine Gaussian directions matters. Choose \(m\) distinct auxiliary split rational primes \(2m<\ell_i\le q^{1.01}\), possible by the prime number theorem in progressions, and put \(L_0=\max\ell_i\). Choose Gaussian factors \(\omega_i\) with arguments \(\beta_i\in[0,\pi/4]\). Distinct prime norms imply nonparallel directions, hence
\[
             |\beta_i-\beta_j|\ge|\sin(\beta_i-\beta_j)|\ge1/L_0.
\]
The second inequality uses the nonzero **integer** determinant of the two factors. Set \(\Delta=e^{-W/4}/10\) in quarter-turn units and choose the \(\epsilon_i\) independently uniform on the intervals of radius \(\Delta/(100L_0)\) about \(\Delta(1/4+\beta_i/2)\). They lie in \([0,\Delta]\), and every difference satisfies \(|\epsilon_i-\epsilon_j|\ge .48\Delta|\beta_i-\beta_j|\). Put
\[
            t_j=(2q)^{-1}\sum_{i\in S}K_{ij}\epsilon_i\qquad(j\in I).
\]
The frame identity and deletion of a constant column give \(2K_S t=\epsilon+\mathrm{constant}\) across the marks, with \(|t_j|\le\Delta\). All subsequent corrections have constant \(K_S\)-image on these rows, as real vectors. Thus the windings can all be chosen zero, and \(v\cdot\theta=\lambda\cdot\epsilon/2\) for every character in the quarter-turn normalization. The physical arc width is \((\pi/2)\Delta<e^{-W/4}\).

Here are the residue data and their full cost. For each odd prime power \(b\le2m\), choose compatible norm-\(N\) residues \(s_0\), for instance by taking the reduction of a fixed Gaussian integer of norm \(N\). Define the point residues by \(s_0\omega_i/\bar\omega_i\pmod b\). All denominators are units at these primes. For \(\lambda\ne0\) integral of sum zero and any Gaussian unit \(u\), let \(D_u\) be the product, over the underlying rational primes, of the maximal powers \(b\le2m\) where the product of these residues to powers \(\lambda_i\) equals \(u\). The product of the chosen moduli divides
\[
 A-u\bar A,\qquad
 A=\prod_{\lambda_i>0}\omega_i^{\lambda_i}\prod_{\lambda_i<0}\bar\omega_i^{-\lambda_i}.
\]
It is nonzero by unique factorization, and divisible by \(1+i\) since \(A\) has odd norm. Consequently \(D_u\le\sqrt2 L_0^{s/2}\), where \(s=\|\lambda\|_1\). For a pair and \(u=1\) we have the stronger \(D_1\le|\operatorname{Im}(\omega_i\bar\omega_j)|\le L_0|\beta_i-\beta_j|\), since \(D_1\) is odd. Axiom (R) requires
\[
 |\sin((\pi/2)v\cdot\theta-\arg(u)/2)|\ge D_u e^{-w(v)/2}/\sqrt2 .
\]
For a pair, the margin \(w(v)-W/2\ge2\log P-o(1)\) and \(L_0=o(P)\) prove it for \(u=1\) by the lower bound for \(|\epsilon_i-\epsilon_j|\). For each other unit the sine is bounded away from zero, while the required right-hand side tends to zero. The pair's ordinary monomial gaps also hold since \(|\epsilon_i-\epsilon_j|\gg\Delta/L_0\).

For all other characters we can impose the stronger condition of avoiding both \(0\) and \(1/2\) modulo one by \(2 L_0^{s/2}P^{-h/2}\), \(h=\|v\|_1\). This suffices for (R) for all four units and the ordinary gaps, by the elementary lower bound on sine by distance to its zero set. Put \(b_0=\log L_0/\log P<.51\). We claim that
\[
           h-b_0s\ge1.06q,\qquad h-b_0s\ge.74h
             \quad(s\ge4,\ q\text{ sufficiently large}).           \tag{VIII.9}
\]
The support bound above also applies if all entries of \(\lambda\) are multiples of \(q\), by first dividing them by their largest common power of \(q\), which preserves support. Hence \((s+1)h\ge q\|\lambda\|_2^2+s q/2\) in all cases. For \(s=4,6\) use \(h\ge1.2q\). For \(8\le s\le q/2\) use \(h\ge (8/9)(3q/2)\). For \(s\ge q/2\) use \(\|\lambda\|_2^2\ge s^2/m\ge2s^2/q\) to get \(h\ge(8/9)(2s+q/2)\). These three bounds give the first inequality. Also \(sh\ge\|v\|_2^2=q\|\lambda\|_2^2\ge q s^2/m\), giving the second.

Exposing one nonzero coefficient of \(\lambda\), the conditional density of \(\lambda\cdot\epsilon/2\) modulo one is \(O(1+L_0/\Delta)\). There are at most \(2^q\binom{h+q}{q}=\exp(O(h))\) integer vectors of length \(h\ge1.2q\), and each determines its row coefficient uniquely. For \(1.2q\le h\le2q\) use the first bound of (VIII.9), for larger \(h\) the second. Their combined failure probability tends to zero since \(1/\Delta\le10(2P)^{q/2}\). Fix offsets obeying all the character inequalities.

### Exact integer moment lifts

For each \(1\le s\le J\) there is an integral \(q\)-vector \(d^{(s)}\) congruent to \((j^s)_j\bmod q\), of sup norm \(O(q^{3/2}\sqrt{\log q})\), with \(K_S d^{(s)}\) constant. To see this, choose independently mean-zero integer representatives \(d_j\equiv j^s\bmod q\) of magnitude at most \(q\) (round between adjacent representatives). The vector \(\tau=Kd/q\) is integral since modulo \(q\) summing \(j^s(j-i)^{(q-1)/2}\) gives zero since its degree is below \(q-1\). Choose the representatives with \(\|\tau\|_\infty\ll\sqrt{q\log q}\) by Hoeffding. Keep \(\tau'_i=\tau_i\) for \(i\in S\). Modulo \(q\) complete on \(S^c\) so that \(K^T\tau'=0\bmod q\). This is possible since any at least \((q+1)/2\) distinct columns of \(K^T\) span its image, as one sees by expanding \((j-i)^{(q-1)/2}\) and applying Vandermonde as before. Independently round the chosen residue representatives on \(S^c\) with zero mean and magnitudes at most \(q\). The deterministic contribution of \(S\) to \(K^T\tau'\) is at most \(O(q^{3/2}\sqrt{\log q})\), and with positive probability the random contribution is of the same order, uniformly in the coordinates. Then \(d^{(s)}=d-K^T\tau'\) has the required residues and height. The identity \(KK^T=qI-\mathbf1\mathbf1^T\) gives constancy on \(S\). Fix such lifts, then set
\[
       D_s=d^{(s)}+q\,K_{S^c}^T\xi_s
\]
with independent integral coordinates of each \(\xi_s\) distributed proportionally to \(\exp(-n^2/(2A^2))\), independently across \(s\). The added vector still has constant \(K_S\)-image and vanishes modulo \(q\). With probability \(1-o(1)\) all digits have sup norm at most \(E=q^{1.625}\). Indeed a sum of up to \(q\) independent signed lattice Gaussians of width \(A\) has tail \(O(\exp[-c x^2/(qA^2)])\), since each summand has moment generating function at most \(\exp(A^2 u^2/2)\), and independence bounds the product by \(\exp(qA^2u^2/2)\). The one-variable bound follows by completing the square, because a shifted periodized Gaussian is maximized at an integer by positivity of its Fourier coefficients. Use a union bound on at most \(q^2\) sums. The fixed initial lifts also fit this bound. In particular \(B/(HE)\to\infty\).

Two concentration estimates used below follow from the same periodized Gaussian Fourier formula. For a nonzero integral vector \(w\), with \(n_0\) nonzero entries,
\[
 \sup_b\Pr(w\cdot\xi_s=b)\ll (A\sqrt{n_0})^{-1}.
\]
Indeed its atom is at most the integral over a period of \(\prod|\psi(w_i x)|\), for the characteristic function \(\psi\) of one coordinate. Hölder and invariance under multiplication by each nonzero integer bound this by \(\int|\psi(x)|^{n_0}\ll1/(A\sqrt{n_0})\). To justify this last estimate uniformly here, Poisson summation gives \(\psi(x)=\sum_k e^{-2\pi^2 A^2(x-k)^2}/\sum_k e^{-2\pi^2 A^2k^2}\). On \(|x|\le1/4\) its nearest term gives \(\psi(x)\le e^{-2\pi^2 A^2 x^2}(1+O(e^{-cA^2}))\). Outside this interval it is \(O(e^{-cA^2})\). We use \(n_0\le q\).

If \(T_0=\max|w_i|\le q^{.10}\) we also have
\[
 \sup_b\Pr(w\cdot\xi_s=b)\ll \frac{\gcd_i(w_i)}{A\|w\|_2}+O(e^{-c(A/T_0)^2}).
\]
Where some \(\|w_i x\|\ge1/(3T_0)\), one factor in the integral has the indicated exponential bound. Else nearest-integer ratios \(k_i/w_i\) all coincide, since \(|w_i k_j-w_j k_i|<1\). Their common value is an integer multiple of \(1/\gcd_i(w_i)\). Integration near these values gives the first term by the explicit Gaussian bound, with a harmless \(1+o(1)\) product correction. In particular, for a balanced primitive \(v\) of length \(h\le q^{.10}\), the probability \(v\cdot D_s=0\) is at most \(O(q^{-.12-.5}h^{-1/2})\) by our choice of \(S\). For all balanced nonzero \(h<q/2\), it is \(O(q^{-.12-.5})\) by the polynomial support estimate. Both bounds use the unrestricted distribution of the digits, with no conditioning on their heights.

### Completion and the tail estimate

We have \(1.708>1.625+.08\). Independently of the integer digits, choose \(X_s\) uniform in \([1,2]\), \(c_0\) uniform in \([1/(20H),1/(10H)]\), and \(\zeta\) with independent uniform \([0,1]\) coordinates indexed by \(S^c\). On \(I\) set
\[
 \theta_j=t_j+c_0
       +\frac1{100H^2}\sum_{s=1}^{J}D_{s,j}X_s B^{-s}
       +B^{-J-2}(K_{S^c}^T\zeta)_j .
\]
All terms except \(t\) have constant row evaluation by \(K_S\), also after deleting \(D\) since its columns were constant on \(S\). Hereafter extend \(v\) by zero on \(D\). Call \(v\) balanced if \(\sum v_j=0\). For a balanced nonzero \(v\) of length \(h<q/2\), put \(t=h/2\). The two multisets of \(t\) labels for its positive and negative entries are different over \(\mathbb F_q\). Their power sums cannot agree for \(s=1,\ldots,t\), since Newton identities recover the coefficients of a monic degree-\(t\) polynomial from those sums, with all denominators units. Hence some \(v\cdot D_s\ne0\) by the congruences, with \(s\le t\le J\). The same bound uses \(h/(2g_0)\) in place of \(h/2\) if \(v\) is \(g_0\) times a primitive vector.

On the event \(\max|D_{s,j}|\le E\), **all** \(h\le H\) are safe. The full phase of each has magnitude at most \(1/10+o(1)\), so the half-turn causes no problem. An unbalanced \(v\) has its constant term of magnitude at least \(1/(20H)\), and the total variation from the digits is at most \(2hE/(100H^2(B-1))=o(1/H)\), with negligible contributions from \(t\) and the final term. For balanced \(v\) let \(k\) be its first nonzero digit and use integrality. Relative to \(B^{-k}\) the modulus of the later digit sum is at most \(2hE/(B-1)=o(1)\), while the leading digit times \(X_k\) has magnitude at least one. Thus the gap is at least \(B^{-h/2}/(200H^2)\), after absorbing \(h\Delta+h q B^{-J-2}\). This exceeds \(2P^{-h/2}\) since \(h\ge2\), \(P/B\gg H^2\).

It remains to show, **without** conditioning on the bounded-digit event, that the probability of any bad noncharacter of length \(h>H\) tends to zero. The number of integer vectors of length \(h\le q\) is at most \((Cq/h)^h\) for an absolute \(C\), by distributing positive and negative tokens among \(2q\) bins. For \(h\ge q/2\) one may instead use \(\exp(O(h))\) as in the character estimate.

* Unbalanced vectors, any \(h>H\). Condition on everything except \(c_0\). The nonzero integer multiplier \(\sum v_j\) gives a density modulo one at most \(1+20H\). For \(H<h<q/2\) the union bound is \(O(H)(Cq/h)^hP^{-h/2}\) at length \(h\), summable to \(o(1)\) since \(h>q^{.08}/2\). For larger \(h\), use \(\exp(O(h))P^{-h/2}\).
* Balanced noncharacters with \(h\ge q/2\). At least one \((Kv)_i\ne0\) for \(i\in S^c\). Otherwise \(qv=K^TKv\) gives a zero-sum rational combination of the rows in \(S\), since \(Kv\) sums to zero. Use the corresponding free \(\zeta_i\) to bound density modulo one by \(1+B^{J+2}\). The tail is \(o(1)\) since \(J\le q/4,\ h\ge q/2,\ B<P/q^{.28}\).
* Balanced nonprimitive \(H<h<q/2\). There are at most \(h(Cq)^{h/2}\) such vectors, since the primitive length is at most \(h/2\). The first nonzero digit \(k\) is at most \(h/4\). Conditional on all digits and all other continuous parameters, \(X_k\) bounds the risk by \(O(H^2 B^{h/4}P^{-h/2})\). The sum tends to zero since \(1/2+1.708/4<1.99/2\).

For the remaining primitive balanced vectors, \(H<h<q/2\), partition by the first nonzero digit \(k\le h/2\). For \(k-1\) prior digits to vanish, \(k-1\) modular power sums must vanish. There are at most
\[
 \binom{q+t-1}{t}\binom{q+t-k}{t-k+1}
       \ \le\ \exp(O(h))\, (q/h)^{h-k+1}
\]
balanced vectors eligible at \(h,k\). To see this, specify all \(t\) positive labels as a multiset, and any predetermined \(t-k+1\) of the negative labels, for example the first in increasing order. The remaining \(k-1\) labels are recovered by Newton identities. The binomial bound follows from \(n!\ge(n/e)^n\) and \((h/u)^u\le\exp(O(h))\) for \(0\le u\le h\). Among the eligible vectors, the probability for the prior \(k-1\) integer digits to vanish is at most \(\mathfrak a^{\,k-1}\), by independence across digits and the atom estimates above. Once these digits are fixed and \(v\cdot D_k\ne0\), continuous averaging over \(X_k\) gives risk \(O(H^2 B^kP^{-h/2})\). This last bound is independent of all integer digits, and requires no bound on \(|v\cdot D_k|\) beyond integrality and nonvanishing.

Write \(\gamma=\log h/\log q,\ \nu=k/h\in[0,1/2]\). The sum of risks at \(h,k\) is at most \(q^{h F}\) where the \(o(1)\)'s below tend to zero uniformly since \(h\gg q^{.08}\).
For \(h\le q^{.10}\) use \(\mathfrak a\ll q^{-.62}h^{-1/2}\). This gives
\[
  F\le (1-\gamma)(1-\nu)-.995+(1.708-.62-\gamma/2)\nu+o(1)
    \le .049-.75\gamma+o(1)<-.009 .
\]
For \(q^{.10}<h<q/2\), use \(\mathfrak a\ll q^{-.62}\). Then
\[
  F\le (1-\gamma)(1-\nu)-.995+(1.708-.62)\nu+o(1)
    \le .049-.5\gamma+o(1)<-.0005.
\]
Thus summing over \(h,k\le q\) gives \(o(1)\). Combining this with the bounded-digit event and the already fixed good splitting proves simultaneous avoidance of the two gaps for all \(v\).

If a column of \(I\) lacks an endpoint, subtract its common minimum, replace \(e_j\) by the attained range, and discard constant columns. This decreases \(W\) and adds a common phase to the marks. Every remaining gap follows by zero extension. Character images are unchanged after zero extension, since the coefficients have sum zero. For the new parent one changes only the common norm representative \(s_0\) in the residue construction. Its cancellation in character products preserves all \(D_u\). The deleted primes had no effect on pair distances. Each pair still has \(\ell^1\) distance \(q+1\) with entries of magnitude at most two, so at least \((q+1)/2\) active primes remain. In particular \(\sum p^{-1/2}\gg q^{.005}\) even after normalization.

The positive-probability construction yields positive Haar measure on the fixed realization component for the chosen offsets. Indeed the continuous perturbations in the all-ones direction and the directions \(K_{S^c}^T\) span its real tangent space after restriction to \(I\). To check the annihilator, extend any real vector by zero on \(D\). If its sum is zero and \((Kv)_{S^c}=0\), then \(qv=K^TKv\) with \(Kv\) supported on \(S\) and of sum zero. Conversely every zero-sum combination of the selected rows annihilates these directions. Surjectivity of this linear map sends the absolutely continuous box law of \((c_0,\zeta)\) to an absolutely continuous law on the component, even conditional on the digits and their multipliers. Projection onto the retained coordinates has the same property after normalization, by the same annihilator argument or by projecting this torus. In particular the good set survives deletion of the countable set of arithmetic angle choices.

Finally \(P_+\le2q^{1.99}\le m^{1.995}\) eventually, and \(W/(m\log m)\le 4d/\rho+o_\rho(1)\). Taking fixed \(\rho\) tending to one and a diagonal sequence of sufficiently large \(q\) gives (VIII.8). No uniform starting threshold in \(\rho\) is being asserted. The real-angle construction and its Lebesgue density estimates take place at the archimedean place. The finite field \(\mathbb F_q\) here indexes the sign matrix and the digits. The parent is an integer. The model imposes neither genuine Gaussian prime arguments nor a factor-by-factor realization of its point residues. \(\square\)

**PROVED (range of the construction).** In fact there is a sequence satisfying the same point residue axiom for which
\[
 \limsup\frac{\log P_+(N')}{\log m}\le\frac{23}{12},
 \qquad \limsup\frac{W}{m\log m}\le\frac{23}{3}.                  \tag{VIII.10}
\]
In particular this gives fakes in every fixed window \(P_+\le m^{2-\delta}\) with \(0<\delta<1/12\). Here are the parameter changes and the estimates at each affected range.

*Proof.* Fix \(0<\xi<1/12\), put \(d=2-\xi\), and choose numbers \(\kappa>2\xi,\ \tau>\kappa+\xi,\ a>\tau\) with \(1.5+a+\kappa<d\) and \(\tau<1/2\). This is possible since at the three lower endpoints the required strict inequality reads \(1.5+6\xi<2\). Choose \(\varepsilon>0\) small enough that
\[
 b=1.5+a+\kappa+2\varepsilon<d,\quad
 \kappa>2\xi+4\varepsilon,\quad \tau>\kappa+\xi+2\varepsilon .
\]
Use \(P=q^d, A=q^a, E=q^{1.5+a+\varepsilon}, B=\lceil q^b\rceil, H=\lfloor q^\kappa\rfloor\), and the stronger atom estimate up to \(q^\tau\). The support union bounds hold since \(\tau<1/2\) and the Gaussian estimate holds since \(a>\tau\). The permutation error is \(O(q^{1.5-d}(\log q)^{3/2})=o(1)\). In (VIII.9) replace \(.51,.74\) by \(.53,.73\), since \(1.01/d<.53\). The same three cases still give its first bound of \(1.06q\), with middle-case lower bound \((4/3-.53/2)q>1.06q\).

In the angle formula drop the factor \(1/(100H^2)\) from the series. Choose \(\xi/2<\lambda<\kappa/2\), and sample \(c_0\) using \(H_u=\lfloor q^\lambda\rfloor\) in place of \(H\). For balanced \(h\le H\) the tail relative to the leading digit is \(O(hE/B)=o(1)\), the whole series is \(o(1)\), and the lower bound \(B^{-h/2}/2\) beats \(2P^{-h/2}\). For unbalanced \(h\le H_u\) the digit series is \(O(H_uE/B)=o(1/H_u)\), so its common-offset test still works. For every unbalanced \(h>H_u\) the density cost is \(O(H_u)\), and the shell bound starts with exponent per token \(\xi/2-\lambda+o(1)<0\). Balanced noncharacters at \(h\ge q/2\) are handled by the final continuous perturbation since \(b<d\). For nonprimitive balanced words of \(H<h<q/2\) there are at most \(h(Cq/h)^{h/2}\) choices, by applying the token count to length at most \(h/2\). The digit cost is at most \(O(B^{h/4})\), so the base-\(q\) exponent per token is at most \((1-\kappa)/2+b/4-d/2+o(1)<0\).

For primitive balanced words write \(\gamma=\log_q h\). The same moment count gives an exponent per token bounded by
\[
 \frac{\xi+\kappa+2\varepsilon}{2}-
 \begin{cases}3\gamma/4,& h\le q^\tau,\\ \gamma/2,&h>q^\tau\end{cases}
 +o(1).
\]
Indeed the atom costs \(a+1/2+\gamma/2\) and \(a+1/2\) respectively per erased integer digit, in base-\(q\) logarithms, and the worst case has \(k/h=1/2\). Both displayed bounds are strictly negative uniformly in their ranges by the choices of \(\kappa,\tau\). This accounts for all changed estimates. First fix \(\xi,\rho\) and take \(q\) large. Then diagonalize with \(\xi\uparrow1/12,\ \rho\uparrow1\). The prime-height exponent tends to \(2-1/12\) and the weight coefficient is at most four times it, proving (VIII.10). \(\square\)

## VIII.4. What prime-angle information is still available

**PROVED (a test for genuine integer clusters).** Suppose all active primes lie in a single interval \([P,Pe^\sigma]\) and their chosen Gaussian factors have arguments in an interval of length \(\eta\), with \(0\le\sigma,\eta\le1\). Then
\[
 M\ \le\ 1+K_C+O\big((1+\sqrt P\,\sigma)(1+\sqrt P\,\eta)\big). \tag{VIII.11}
\]
Indeed put unit squares around the distinct Gaussian factors. They have disjoint interiors and their union lies in the annular sector enlarged by an absolute distance. At radius comparable to \(\sqrt P\) the enlargement adds \(O(1)\) to the radial thickness and \(O(P^{-1/2})\) to the opening angle. Its area is bounded by the product in (VIII.11). Bounded \(P\) is absorbed into the absolute constant. Apply (VIII.1). A bounded number of associates or conjugate sectors only changes the absolute constant. This step counts planar integer points using both the radial and the angular interval, a count with no verbatim counterpart at ultrametric infinity.

For the first construction the angles themselves all lie within \(O(q^{-.243})\) of a common direction, since the digit series is \(O(E/(H^2B))\) and the two perturbations are smaller. Its radial bin has \(\sigma\ll q\log(q)/P\). Formula (VIII.11) allows only \(O(q^{.757}\log q)\) marks there for fixed \(C\), whereas the fake has order \(q\). In the modified parameters of (VIII.10) its angular diameter is \(O(q^{-\kappa-\varepsilon})\) and the same area test again gives \(o(q)\). Thus the near-axis feature of the most economical construction would be detectable even by counting all Gaussian integers in its joint height-angle range.

**PROVED (spread at the same weight).** The models in (VIII.10) can simultaneously be required to have distinct prime angles separated by \(\gg1/m\), even up to conjugation. Thus the angular concentration just exhibited is avoidable with the same asymptotic prime-height and weight bounds.

*Proof.* Use any fixed parameter choice for (VIII.10), including fixed sampling density \(\rho\). We first build a rational background in the realization torus, with denominator \(12q\). Write \(J_1=(0,1,\ldots,q-1)\) in the column order. We seek
\[
    b_j=1/12+(j/q+z_j)/3,\qquad z_j\in\{0,1,2\},\qquad
    (2K b)_i-(2K b)_{i'}\in\mathbb Z\quad(i,i'\in S).
\]
Define \(T\) to be the labels for which some choice of \(z_j\) places \(b_j\) within \(10P^{-1/2}\) of \(\tfrac12\mathbb Z\). Each of three equally spaced grids meets those neighborhoods in \(O(1+qP^{-1/2})\) positions, so \(|T|=O(q^{\xi/2}+1)\). For \(j\in T\) prescribe \(z_j=0\), giving \(1/12\le b_j\le5/12\).

With probability tending to one over the original row sampling, \(K_{S^c,T}\) has full column rank over \(\mathbb F_3\). For a nonzero vector on \(T\) modulo three take an integer lift \(v\) supported on \(T\) and augment it by zero at infinity. Use the two integral eigenvectors \(Cv\pm gv\) of the Paley conference matrix from VIII.3. At a prime above three their difference \(2gv\) is nonzero, since \(q\ne3\). The convolution argument above gives support at least \(\sqrt q\) for one reduction. It follows that \(Kv\bmod3\) has at least \(\sqrt q-|T|-1\) nonzero finite entries. Here the one additional index accounts for infinity, so this argument also covers unbalanced \(v\). The probability all those entries are sampled is at most \(\rho^{\sqrt q-|T|-1}\). A union bound over \(3^{|T|}\) vectors proves the claim. Choose \(S\) satisfying this and all previous erasure estimates.

The power-sum calculation gives \(\tau_1=KJ_1/q\in\mathbb Z^q\). Work modulo three, put \(u_S=(\tau_1)_S\), and choose \(u_{S^c}\) so that \((K^Tu)_T=0\), possible by the full column rank just proved. Set \(z=-q^{-1}K^Tu\pmod3\) and take its stated representatives. On \(S\), \(\tau_1+Kz\) is constant modulo three by \(KK^T=qI-\mathbf1\mathbf1^T\). This gives the required windings, including after removing column zero. Distinct grid positions differ from an integer by at least \(1/(3q)\) when subtracted, since their column labels differ modulo \(q\). Their sums have odd numerator over \(6q\) and are at distance at least \(1/(6q)\) from an integer. This proves separation also after changing a direction to its conjugate.

Use \(b+t\) as background plus offset correction. Divide the remaining perturbation in the angle formula of (VIII.10), including \(c_0\) and the final continuous term, by \(24q\). Choose a fixed integer \(n_0\ge3\) with \((d-b_*)n_0>2\), where \(b_*=1.5+a+\kappa+2\varepsilon\) is the constant base exponent. Reorder the first digits as \(2,3,\ldots,n_0,1,n_0+1,\ldots,J\). All marks have the same offsets as before. A character acquires a shift of \(\lambda\cdot k/2\) in quarter-turn units from the integer mark windings \(k\). This permutes the two forbidden centers. In the round-5 residue axiom the prescribed phase \(\arg(u)/2+(\pi/4)\lambda\cdot k\) cancels this extra term in radians. All character and residue estimates survive.

For short balanced words \(h\le H\), or short unbalanced words \(h\le H_u\), a background value outside \(\tfrac12\mathbb Z\) misses it by at least \(1/(12q)\). The perturbation contributes at most \((.1+o(1))/(24q)\) on the bounded-digit event. This suffices for \(h\ge2\). At \(h=1\) the construction of \(T\) already gives the required larger gaps. At a center a short unbalanced word is protected by \(c_0\), giving \(\gg1/(qH_u)\), enough for \(h\ge2\) since \(\lambda<\kappa/2<d-1\). For a balanced word at a center, multiplying the background by \(12q\) and reducing modulo \(q\) shows that the first power sum vanishes. Write \(n=h/2\). If \(n<n_0\) then \(n\ge2\), and Newton's identities give a nonzero digit by position \(n-1\). If \(n\ge n_0\), there is one by position \(n\). The leading-digit argument now gives distances \(\gg B^{1-n}/q\) and \(\gg B^{-n}/q\) respectively. Since \(P\gg q\) and \((P/B)^{n_0}\gg q\), both suffice.

At larger lengths the continuous density bounds lose a factor at most \(O(q)\). For nonprimitive words the first usable index increases by at most \(n_0\), costing at most \(B^{n_0}\). For primitive balanced words whose first nonzero digit lies at \(k>n_0\), the first \(k-1\) moments are the same set as before, so the erasure count and the atom estimates are unchanged. For \(k\le n_0\) the unrestricted bound \((Cq/h)^h P^{-h/2}\) times a fixed power of \(q\) already suffices when \(H<h<q/2\). The large tail still uses \(\zeta\). The strict negative exponents per token therefore survive all these polynomial factors for each fixed parameter choice. The successful bounded-digit perturbation shifts each coordinate of the rational grid by \(o(1/q)\). At least \((q+1)/2\) primes survive normalization. Their angular diameter is bounded below by an absolute positive constant by the spacing. Positive conditional Haar measure follows from the same tangent-space calculation. Lastly take the diagonal sequence sufficiently far out for each fixed parameter choice. This proves the spread assertion at (VIII.10). It gives no estimate comparing an assigned angle with the true Gaussian angle at that particular prime norm. \(\square\)

**COMPUTED (arithmetic checks).** On the ten-mark integer example \(N=1176852625\), at arc constant \(C<33/4\), the exponent affine rank is five. Appending residue indicators modulo 3 raises it to six and then modulo 7 to nine. The Smith invariants of \(B\) are \(2,2,2,2,24\) and the coordinate orders are \(24,8,6,4,24\). The corrected exact phase in quarter-turns is \((7/24,3/8,1/6,1/4,11/24)\). At \(q=383\) the script uses 189 of the 191 positive rows of column zero. After deleting that column, all 17,766 pair differences have length 384 on 382 labels, the exact two-unit surplus used in the proof. All 95 moment lifts have the required congruences and constant row evaluations, and the largest coefficient in absolute value is 7,558. In the spread construction, exactly 124 labels require prescription at this order. The restricted matrix has rank 124 modulo three, and the resulting grid meets all the one-prime clearance, separation and row-winding tests. These are finite checks of the inputs to the construction. [Integer marks and moment lifts](viii_window_check.py)

**COMPUTED (the erasure step).** At \(q=127\) the minimum support of \(Kv\) over primitive balanced length-four words is 73 modulo 3, 56 modulo 4, 74 modulo 5 and 74 modulo 7. Replacing the mod-four check by a mod-two check already fails for length two, where the minimum is 2 and the argument needs 10. For two disjoint three-element multisets with repetition over \(\mathbb F_{11}\), the exact ordered-pair counts with the first zero, one, two and three moments equal are 40,370, 3,520, 440 and 0. The script explains its affine symmetry reduction. [Erasure and Newton checks](viii_erasure_check.py)

**HEURISTIC (verdict).** I would close the function-field programme. For the integer problem the low-height restriction alone no longer singles out an obstruction to building an angle model. A genuinely different input is the association of each individual prime norm with its actual Gaussian argument or its prime-level residues. The height-order inequality (VIII.6) gives one provable consequence of that association. I did not turn it into a growth bound on the whole window, or improve the coefficient \(2/3\) for unrestricted integer parents. The spread example also leaves open the sharp union sieve with the unchanged prime angles and any factor-by-factor residue model.

## Twelve-line summary of Part VIII

1. **PROVED.** The low-height window contains nonarithmetic angle models of growth-scale weight.
2. **PROVED.** One sequence has \(P_+\le m^{23/12+o(1)}\) and \(W\le(23/3+o(1))m\log m\).
3. **PROVED.** Every monomial satisfies both sharp gaps, including all rational characters of the profile.
4. **PROVED.** Every odd prime power through \(2m\) satisfies the point residue axiom (R), simultaneously for every character and unit.
5. **PROVED.** The gap-admissible angles have positive Haar measure on a fixed realization component.
6. **PROVED.** At the same weight and height, the assigned prime angles can be spaced by \(\gg1/m\) even after conjugation.
7. **PROVED.** Genuine integer clusters satisfy \(M\le\omega(N')+1+K_C\) by affine independence after adding finitely many inert residue columns.
8. **PROVED.** This gives coefficient \(1/2\) on four-full integer norms and bounds \(M\) in \(P_+\le(2-\eta)M\log M\).
9. **PROVED.** On a genuine low-height cluster every finite-order coordinate in the exponent quotient has order at least \(\exp(c_\delta M)\) for large \(M\).
10. **COMPUTED.** At \(q=383\), all 95 integer moment lifts and all 17,766 pair margins pass exact checks.
11. **HEURISTIC (scope).** The integer growth rate and the unrestricted coefficient \(2/3\) remain open to improvement.
12. **HEURISTIC (request).** Please referee the mod-four support argument, the intermediate-length erasure count, and the ordered residue witnesses first. These are the three steps on which the point-residue conclusion depends.
