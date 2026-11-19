# Referee report, block "bmy": the proof of inequality (IV.2) in Astra's Part IV

Scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4_bmy/`:
* `bmy_symbolic.py`: sympy derivation of every local and global constant from intersection numbers.
* `bmy_global.py`: a lattice-level engine. For explicit primitive columns in Q[t]^2 it builds the whole surface Y and recomputes everything.
* `stress.py`: pencils, nested trees and long tangency chains.
* `genuine.py`, `cube_test.py`: genuine function-field clusters from the record.

All of this code is mine. I did not use Astra's `iv_section_check.py`, except to compare its output with mine afterwards.

## Verdict

**(IV.2) is CORRECT.** I checked the proof line by line, and every constant comes out the same in my own symbolic and lattice-level computations. The corrections are minor and none of them weakens the statement:

* m >= 3 suffices, not only m >= 4.
* The hypothesis delta >= 0 is never used.
* For the "some minor overshoots" consequence, 49 marks give the slightly better constant 75/10192 instead of 1/136.
* I could check Langer's Theorem 0.1 only through secondary sources (details in (a)).

**Significance, which Astra missed.** (IV.2) applies verbatim to *genuine* clusters of the record's function-field model over C(t) or Q(i)(t). It gives:

* M(D,0) <= 24 for every D;
* sup_D M(D,c) < infinity for every c.

The record lists both as OPEN (round 29 ff.md §1, round 30). See §6. This is the first exceptional-set-free uniform bound in any model of the problem. It has no arithmetic transfer.

## (a) Langer's theorem and the pair (Y,B)

**The statement.** arxiv.org is not reachable from this session (fetch failed: host not resolvable). I could not read Langer's paper (arXiv math/0012180; PLMS 86 (2003)). Secondary sources found by WebSearch quote [Langer, Thm 0.1] as follows:

> X a normal projective surface, D a boundary Q-divisor (coefficients in [0,1]), (X,D) log canonical, K_X+D **Q-effective**. Then (K_X+D)^2 <= 3 e_orb(X,D).

The equality case forces K_X+D to be nef. The sources are Pokora, arXiv 1612.05141, and the abstract of the paper itself ("log canonical case ... BMY type inequality"). This is the form Astra uses: Q-effective, not nef. The Q-effective hypothesis is the one I verify in (b). **Caveat:** this is unverified against the original text.

**e_orb.** On a smooth surface with an SNC boundary with coefficients a_i in [0,1), the local orbifold Euler numbers are:
* 1 at an interior point;
* 1-a_i on D_i;
* (1-a_i)(1-a_j) at a crossing.

Summing these gives exactly Astra's description:

e_orb = e(Y) - sum a_i e(D_i) + sum_{crossings} a_i a_j.

I derived this identity myself. That Langer's analytic local e_orb equals these values in the SNC case is standard, but I did not check it in the original.

**SNC.** I checked the procedure in local coordinates (z, y), where the sections are graphs.
* A blow-up in the common direction gives y' = (y - y(0))/z. The strict transforms stay graphs, so they are transverse to E.
* The strict transform of E_{nu-1} meets E_nu at the "vertical" point, and that point lies on no section.
* Hence the positive-weight boundary consists of:
  * the sections (weight 2/3), crossing only transversally in pairs, away from every positive-weight exceptional;
  * chains E_1 - ... - E_e (weight 2/3), where E_e meets the two separated sections at distinct points;
  * four generic fibres (weight 1/2).

  This is SNC. Weight-0 curves (the exceptionals at k >= 3 centres, and the special fibres) are not in B.
* With coefficients < 1 the pair (Y,B) is klt, hence log canonical.

The engine checks this independently: it asserts that the lattice value C_i'.C_j' equals the combinatorial count of the remaining transverse points, on every example.

## (b) Q-effectivity of K_Y+B

**On X = F_e.**

K_{X/P^1} + C_i + C_j ~ (b_i + b_j - e) F = (C_i.C_j) F.

I take the representative V_ij = sum_p (C_i.C_j)_p F_{pi(p)}.

**Blow-ups.** A blow-up through l of the two sections adds (1-l)E. The invariant to carry is: "at every common point q, mult_q V_ij >= (C_i.C_j)_q". It is preserved, because when l = 2:

(C_i.C_j)_c = 1 + sum_{q in E} (C_i'.C_j')_q.

So V_ij stays effective.

**Averaging.** Average over all pairs, then add the remaining terms:
* (a - 2/m) sum C_i, which is >= 0 iff m >= 3;
* the exceptional weights b >= 0;
* pi^*K_{P^1} + 2F ~ 0, which is where the four half-fibres come in.

This gives K_Y+B ~_Q effective.

**Machine check.** `bmy_global.py` verifies, on every example, two things:
* the class identity K_Y - pi^*K_{P^1} + C_i' + C_j' = (C_i.C_j)F + sum_c (1-l_c) E_c^tot;
* the nonnegativity of every coefficient of V_ij on the strict-transform components.

**Correction:** m >= 3 is enough.

## (c) Normalisation at infinity, (IV.3), and I = (m-1)S/2

* **The identity for I.** On any F_e, 2 C_i.C_j = C_i^2 + C_j^2, so I = (m-1)S/2.
* **Elementary transformation at a point with k > m/2 sections.** The self-intersections are -1 on the k sections and +1 on the others, so S drops by 2k - m.
* **Termination.** I >= 0 gives S >= 0, so the process terminates.
* **The cyclic shift sigma: i -> i + floor(m/2).** Every final block has size <= floor(m/2), so i and sigma(i) lie in different blocks. Hence C_i.C_sigma(i) = deg D_{i,sigma(i)}: the finite contacts are unchanged by the transformations at infinity.
* **Summing.** sum_i C_i.C_sigma(i) = S, because sigma is a permutation. Summing the degree bounds over sigma gives

  **(IV.3) S <= sum_h s_h - mL/4 + m delta.** Correct.

* **Machine check.** The engine implements the transformation as an ultrametric rule on the contact orders at infinity: c-1 for pairs inside the group, c+1 for pairs outside, 0 for mixed pairs. It verifies sum over fibres of the contacts = C_i.C_j, 2I = (m-1)S, and (IV.3) on every example. Elementary transformations occurred in many of the tests, up to 5 of them.

## (d) The accounting

**One blow-up.** I derived the update from classes: C' = C - E, V' = V - E, K' = K + E.
* The change in (K+B)^2 is -(1+b-b0-ak)^2.
* The change in virtual e_orb is 1 - 2b + b b0 + ak(b-b0) - a^2 C(k,2).
* Re-adding a^2 per consumed pair, which is Astra's convention of charging a^2 per final transverse pair, gives exactly Astra's F(k,b0,b).

**Values at a = 2/3:**
* F(k,0,0) = 4 - 8k/9;
* F(2,0,a) = 8/9;
* F(2,a,a) = -4/9;
* a chain of length e totals 4/3 - 4e/9 <= 4e/9 for e >= 2;
* each half-fibre contributes -1 + am/2;
* after the four fibres, the constant and S terms leave 8S/9.

Total:

3e_orb - (K+B)^2 = 8S/9 + 4T/9 + sum_c F_c.

**The deduction phi.** At a centre with k >= 3,

F - (4/9)C(k,2) = -(4/9)(k+6)(k-3)/2 = -(4/9) psi(k), with psi(k) = C(k,2) + 2k - 9.

So phi is exact for k >= 3. For a chain, F - (4/9)e = (12-8e)/9 < 0, consistent with phi(2) = 0.

**Superadditivity.** psi(u+v) - psi(u) - psi(v) = uv + 9 <= 2uv for u, v >= 3. The bound phi(u+v) - phi(u) - phi(v) <= 2uv was checked exhaustively for u, v < 300.

**Layers.** A class at layer (t0, r) is the r-th infinitely near point on the fibre t0. Distinct classes, depths or roots give distinct centres, so there is no double counting. The actual group at a centre may contain unmarked rows; phi is monotone, so this only helps. Centres at infinity and other unmatched centres have -(4/9)phi <= 0, so they may be dropped. Finally:
* sum_h E_h = sum_{i<j} E_ij, from the gcd definition of E_ij;
* sum_j phi(s_j) >= psi(s) - 2E_h, which also holds for an empty layer.

**Global check, independent of Astra's formulas.** `bmy_global.py` builds Pic(Y) = <C0, F, E_c>, the classes of every section, exceptional curve and fibre, K_Y, B, e(Y) = 4 + n, e_orb by the SNC formula and (K_Y+B)^2. It then checks four things:
* the ledger identity, exactly;
* slack >= 0, consistent with Langer;
* (IV.3);
* (IV.2), with delta and E computed from their definitions.

It ran on these examples:
* about 40 random Hermite-interpolated configurations: m = 4..8, prime powers up to 3, groups up to size 8, partial divisibility;
* pencils y = 1 + u_i t^r, where the slack is exactly 4r, so BMY is nearly tight and the per-layer constant 4 is attained;
* nested trees and long tangency chains (50 blow-ups);
* the record's genuine clusters: the cube at D = 8, c = 0, and the vertical-line clusters M(7,1) >= 8 and M(8,1) >= 14.

All checks passed.

## (e) The final substitution

* The coefficient of S is 8/9 + (4/9)(m-1)/2 = 2(m+3)/9 > 0, so substituting the upper bound (IV.3) is legitimate.
* The per-layer identity is exact:

  (2(m+3)/9)(s - m/4) - (2/9)s(s-1) + 4 - 8s/9 = 4 - m/6 - (2/9)(s - m/2)^2.

* The delta term is 2m(m+3)delta/9 and the E term is 8E/9.
* Multiplying by 9/2 gives exactly **(3/4)(m-24)L + sum_h (s_h - m/2)^2 <= m(m+3)delta + 4 sum E_ij**. The constants are confirmed.
* At delta = E = 0 this gives m <= 24.
* The value 3(m-24)/(4m(m+3)) at m = 48 is 1/136, as stated. Its integer maximum is at m = 49 (75/10192), but the difference is negligible.

## 5. Novelty

The record's only BMY computation is round4/L7.md §7. It is a "no obstruction" check for seven balanced curves in M_{0,7}, with Kummer weights on a fixed cubic model. Part IV uses a different surface and weights, and a growing collision tree whose layer and prime-power bookkeeping is new. The orbifold-BMY bound for sections of a ruled surface with prescribed contacts is not in the record. It is in the spirit of Miyaoka's bounds for curve configurations, which are known in general, but this application is new here.

## 6. Significance: a consequence for the record's function-field model

Take a cluster z_a = prod_k (t-rho_k)^{a_k} (t-conj rho_k)^{e_k-a_k} of the round-29/30 model over C (or over Q(i), embedded in C). Put:
* P_a = prod_k (t-rho_k)^{a_k}, with columns (Re P_a, Im P_a) in R[t]^2;
* n_a = P_a conj(P_a), dividing N;
* deg N = 2D.

These data have the following properties:
* The columns are primitive: gcd(P_a, conj P_a) = 1.
* E = 0: n_{min(a,b)} divides P_a conj(P_b) and is real.
* z_a - z_b = 2i conj(Q_full) Im(P_a conj P_b) / (conj P_a conj P_b). Hence the cluster condition deg(z_a - z_b) <= K is exactly deg D_ab <= (alpha_a + alpha_b)/2 - L/4 + delta, with delta = K - D/2 <= c.

This works for non-squarefree N too. It is confirmed by machine on the cube cluster (delta = 0, E = 0) and on the vertical-line clusters (delta = 1/2 and 1).

Applying (IV.2), also to subsets of 25 or 48 points, gives:
* **M(D,0) <= 24 for all D**. The record has M(8,0) = 8 and sup >= 8.
* M(D,c) <= 47 for D > 68c, and M(D,c) <= 24 for D > 467c.
* Hence **sup_D M(D,c) <= max(47, 2^{2c+2}(68c+2)) < infinity for every c**. This is the record's open question, answered over C. It is not answered over F_{q^2}: BMY fails in characteristic p, consistent with Astra's trace construction.

This is a theorem of exactly the "truncated function-field Vojta" shape the record predicted, and (R_4) becomes moot for finiteness.

For the integer problem it changes nothing. The proof uses an actual complex surface. The arithmetic analogue is an arithmetic BMY / Vojta-type statement on a blown-up ruled surface, so it sits behind F2 and Question A. Astra says the same.

Sources: [Langer, arXiv math/0012180](https://www.arxiv.org/abs/math/0012180) (not fetched); [Pokora, arXiv 1612.05141](https://arxiv.org/pdf/1612.05141) (statement quoted via search results).
