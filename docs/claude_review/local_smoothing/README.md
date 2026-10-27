# OpenAI "local smoothing in all dimensions" (Oct 26 2026): relevance to Cilleruelo-Granville

Judge's summary of facts.md, chains.md, lifts.md, steelman.md and skeptic.md (same directory).
**Verdict: NOT USEFUL.** Nothing survives that justifies a round.

## 1. What was claimed, and verification status

**Primary source (read):** github.com/openai/math at commit 7ef5c0cc (2026-10-26), CONTENTS.md entry 079, and the preprint
`preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/` (paper.pdf plus LaTeX sources).
The Lean side is in `lean/ComparatorChallenges/LocalSmoothing.lean` and `lean/OAI/Analysis/Fourier/LocalSmoothing/`.
We could not fetch openai.com, X, arXiv or any news host (DNS failures). Everything from those hosts comes from search snippets only.

**Exact statement (Thm 1.1).** For n >= 2, 2 <= p < infinity and eps > 0,
||e^{-iT|D|} f||_{L^p(R^n x [1,2])} <= C ||J^{s_+(p)+eps} f||_p, where s_+(p) = max(0, (n-1)/2 - n/p).
eps = 0 is allowed only at p = 2.
- This is the Euclidean, constant-coefficient Sogge conjecture. It is sharp up to eps, not eps-free. No FIO or manifold version is claimed.
- The core result is a critical annular estimate at p_c = 2n/(n-1) with a (log N)^{2/p_c} loss.
- The paper also proves that some log loss at p_c is necessary (lower bound (log N)^{1/p_c}), so an eps-free endpoint at p_c is false.

**Claimed consequences.** All go through known implication chains, with explicit log powers:
- Bochner-Riesz, and Stein's square function for p >= 2.
- Restriction for the sphere, paraboloid and cone in the open range only.
- Kakeya/Nikodym maximal with a (log 1/delta)^{4/p_c} loss, and Kakeya sets of full dimension.

No decoupling theorem is claimed, and the paper contains nothing arithmetic.

**Method.**
1. Imaginary translation along light rays gives a positive semidefinite, divergence-free tensor.
2. A fibrewise quaternionic Monge-Ampere Dirichlet problem is solved, with an Aubin-Mabuchi-type energy.
3. Two circle-mean comparisons follow, each losing one log.

**Verification.**
- Not peer reviewed. It is released as model output.
- The Lean comparator states six eps-level targets with standard axioms. The comparator file has the usual `sorry` placeholders. A text scan of the proof library (about 2,261 modules) found no sorry, axiom or admit.
- **We did not compile it.** The formalized log power, (1+log 1/h)^{2n+5}, is much weaker than the paper's. The sharp log refinements are unformalized.
- We found no expert assessment of the proof and no identified gap. The reactions we saw (snippets only) are general comments on OpenAI's mass releases.

**Already known in dimension 2 (our dimension):**
- Local smoothing for all p with an eps loss: Guth-Wang-Zhang, Ann. of Math. 192 (2020).
- Circle restriction and Bochner-Riesz: Fefferman, Zygmund, Carleson-Sjolin.
- Kakeya: Davies 1971, and Cordoba 1977 with a better log power than this paper's 2-D result.
- Cone restriction in R^3: Barcelo 1985.
- Decoupling: Bourgain-Demeter 2015.

Every eps-level planar statement was already a theorem before our program began. The only new things in 2-D are an explicit (log N)^{1/2} at p = 4 and a new proof.

## 2. Verdict for our problem: not useful

These are the decisive reasons. Each is proved in the companion notes, and the skeptic re-checked them.

1. **One circle makes the estimates vacuous (lifts F1).** All lattice points on x^2+y^2 = N lie in one Laplace eigenspace. On such data the wave propagator and J^s act as scalars, so every local-smoothing and Bochner-Riesz inequality becomes an identity. This includes the paper's periodic version on R^n/(100Z)^n.
   - This stays true for every time T, so a "long-time torus version" would not help either.
   - The paper's imaginary-time filter does act at our arc scale R^{-1/2} (steelman S1; skeptic confirmed it in barrier.tex and energy.tex). But the resulting energy is L^4 averaged over directions, and gives at best M <~ R^{1/4+o(1)}.
2. **L^p and energy estimates cannot see a cluster at our scale (chains Prop 2, lifts F2).**
   - An M-point arc cluster has |F| >= M/2 only on a set of area about C^{-3} R^{-1/2}.
   - So any bound ||F||_p <= K_p ||F||_2 gives only M <~ R^{1/p}.
   - Any L^p, Lambda(p), energy, decoupling or discrete-restriction chain is therefore capped at Theta(log R). That is worse than our proved (2/3+eps) log R/loglog R.
   - The obstruction is a volume factor, not a loss. An eps-free or log-free estimate changes nothing.
3. **Fakes and generic real configurations satisfy the estimates (chains Prop 1).**
   - The inequalities are Euclidean and stable under perturbation.
   - Our fake circles satisfy them, and so does the transcendental set Z = {R e^{ik/R}}, which has about C sqrt R points on the arc.
   - So they fall under the fake-circle cap and the round-10 Thm B cap.
   - At our scale the arc is a single canonical sector, where these estimates carry no information.
4. **The open routes are not reached.**
   - **Route 5 (MIN^pair)** needs a pointwise lower bound about discrete polynomial forms on a slab-cut ternary cube, in growing dimension, at precision e^{-W/4}. The paper gives averaged upper bounds in fixed dimension with curvature. Every one of these features mismatches. Its "determinant" (Monge-Ampere) has nothing to do with the record's interpolation determinant.
   - **P3-all*** is a finite combinatorial/Kronecker (L^infinity) question with no Euclidean structure.
   - **H^gen/L7** concerns Gaussian primes. The paper has no primes. Repository entry 003, a quasi-RH claim, is about Q(sqrt -3) and finite-order characters, which is the wrong field and the wrong characters. It is also a robust input capped by round-10 Thm B.

**Rejected claims from the companion notes:**
- The facts note proposed revisiting if the method is extended to a "periodic/long-time T^2" setting. The skeptic refuted this: by item 1, single-shell data stays scalar for all T.
- The lifts note called "single-circle Lambda(p) at p ~ log R/loglog R" its least remote link. The skeptic showed it is a cap: even an optimal Lambda(p) bound gives M <= 4 K_p^2 log R.
- The reverse direction (our arc bound feeding toral eigenfunction or nodal estimates à la Bourgain-Rudnick) helps other problems, not ours.

## 3. Surviving mechanisms worth a round

**None.** The strongest candidate (steelman S3/S4) is recorded here as a precise question, but only to log it. It does not justify a round triggered by this news.

> **Q (discrete anti-concentration, AC(u)).** Fix u in the route-5 gap region. Is it true that no polynomial of degree m' ~ sqrt(d log r) with leading form <k',u>^{m'} has size <= d^{m'} eps_0 on half of the slab-cut ternary cube? The only exceptional directions should be the collision strata.

Why it is not worth a round:
- The continuous version is classical and holds with room to spare (Carbery-Wright 2001). Even that room is only asymptotic: by skeptic_checks/cw_room.py, with Carbery-Wright constant 1 it fails for n below about 6,000-12,000.
- The discrete version is the record's existing J(A) problem, and fake Thm 2.1 forces any proof of it to decide the barrier.
- AC is also stronger than MIN^pair needs: it quantifies over every half of the cube, not just A. It controls kappa_V, while Thm 4.1 needs kappa_H <= kappa_V.
- The paper contributes nothing to either the continuous or the discrete half. Its Monge-Ampere and Aubin-Mabuchi tools are the classical pluripotential tools (Zaharjuta, Berman-Boucksom). The finite set A is pluripolar, and k_x^3 - k_x kills every norming inequality at degree >= 3.

**Cost:** if Q is ever pursued, it belongs to route 5 on its own merits, not to this announcement. It would be a multi-round effort with a low chance of success (the steelman estimated under 5% overall).

## 4. What would change the verdict

These would change it:
- A **discrete** anti-concentration or inverse Littlewood-Offord theorem for polynomials of degree ~ sqrt(d) on slab-cut ternary cubes, with collision strata as the only exceptional directions. Equivalently, a discrete analogue of the paper's Monge-Ampere/energy method that controls J(A). This would reopen route 5.
- A technique that turns the positivity method into **pointwise lower bounds at a single point** that use integrality. Even then, it would have to avoid the fake-circle cap, so it must use something fakes lack.
- A harmonic-analysis statement about the **supremum** of a single-shell count, not an average, that genuinely uses arithmetic (Gaussian primes or integrality across many shells). No such statement exists or is claimed.

These would not change it:
- An eps-free or log-free form (the obstruction is volume, not loss).
- A periodic or long-time T^2 version (still scalar on one shell).
- Confirmation of the n >= 3 results, or independent verification of the proof.

If peer review confirms the paper, record it as background (log-loss Kakeya and restriction) only.
