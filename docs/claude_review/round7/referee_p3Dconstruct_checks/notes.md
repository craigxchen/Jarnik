# referee notes p3Dconstruct (scratch)
## Findings so far (step 1, after reading sections 0-4)
- F1 (major, overstatement): summary/table "every completion has W>=..." but Cor 1 is probabilistic (completion prob >= exp(-M^{eta/4})). Every label vector is in the support of the uniform design, so literal "every completion" would be a full barrier. 
- F2 (major, scope): Thm1 hypothesis is lambda_j >= (log M)^-C (depends on Legendre symbols (p_j/q) at k(q)>=3), summary/Consequence 2 replace it with reciprocal mass of Q*. For Q* = {q = +-1 mod 2^K} (all k>=3), columns with (p_j/q)=-1 never collide singly. Fixable via 2-column monomials e_j - e_k with equal parity (check).
- F3 (major, gap): "Answer to suggestion" dichotomy random->Thm1 / coherent->Thm2 is not exhaustive. Suggestion = coherent non-twins + designed (solved) twins => twins deterministic, neither (Ind) nor G-coherent. Thm 2 needs coherence on twin columns (phi_t, step 5 includes all physical columns). Without twin coherence |T*| unbounded.
- F4 (major, false exact claim): §3.2 "Matching parities costs 2^K >~ X/(log X log p_j)" labelled exact; zeta_j = pi_j * rho^2 (Norm = p_j * square) matches ALL parities at zero cost (sharp §6.6 lists this). So "parity is genuine obstruction on q=+-1 mod 16" holds only for zeta_j not of that form.
- F5 (moderate): Consequence 4: (Ind) as defined (conditional on data of columns outside J) FAILS for twin columns: given non-twin labels and kappa, t_x is deterministic. Needs unconditional variant (Ind').
- Check: Thm 2 proof steps 1-8 rederived OK. Cor 2 rates OK (M^{1/2-o(1)} correct but weak; with C small get M^{2/3}).
- Minor: "Gaussian primes of odd norm" includes inert primes (zeta/conj zeta=1): should be split primes.
- Minor: "non-root repair ... (exact)" uses GRH; quadratic bound is an upper bound on the certificate, not a lower bound on cost.
## step 2
- F6 (moderate): (Ind) conditional on data of columns outside J is FALSE for sharp's random fibre with J = non-twin cols (twin data reveal A^{nontwin} k^{nontwin} up to shift -> M linear constraints). So Consequence 1 (Prop 6.5 as case) and Cons. 4 mis-verify (Ind). Proof of Thm 1 works with UNCONDITIONAL independence (Ind'); fix statement.
- F7 (minor): Thm1 via m^gen + aggregated (1.2) loses factor 4: C=0 gives (1/2)log^2/loglog not Prop 6.5's 2 (per-unit). "Prop 6.5 is the case C=0" true in rate only.
- F8 (minor-moderate): Open pt 3 "Any construction at O(M log M) must control smooth parts of binom(M,2) integers of size exp(M^{1/2-o(1)})" -- Cor 2 only for G-coherent on Paley, and only AVERAGE height.
- Todo: SNF of Paley 0/1 rows subsets (Cons. 4 surjectivity "fails mod 2"); rerun root checks independently; check T2 numerically small; check matching claim with zeta=pi.
## step 3
- paley_snf.py / paley_snf2.py: P=(1+H')/2 (M x q0). rank_F2 = M-1 (q0=3 mod 8), M/2 (q0=7 mod 8). Greedy finds Z-surjective row subsets |R| = M/2 for ALL q0 in {43,47,59,67,71,79,83,103,107}. => Cons. 4 caveat "fails mod 2 for Paley" is wrong for subsets; Cons 4 holds with eps=1/2 (given (Ind') fix) for these q0. BUT Cons 4 redundant: if non-twin labels random, Thm 1 applies to them directly.
- F9 (major, overstatement): "robust" only w.r.t. prime set/weights, not label distribution: independent labels uniform on parity class MINUS the 4-torsion -> zero single-column hits at every q. So "independent randomness on any prime set is fatal" (summary) and Cons. 3 overstated. Robust replacement: two-column monomials e_j - e_k; for i.i.d. labels w/ common distribution per (q,parity), Pr[ell_j-ell_k in T4] >= sum_C mu(C)^2 >= p(q) (Cauchy-Schwarz over T4-cosets).
- F10 (moderate): weights fixed before labels. Adaptive design (choose p_j after seeing hits; at k=2 classes hits independent of parity via coupling u -> u + (m/4)[par != nu]) pays only sum_j H_j log X = O(M log M) for single columns. Single-column mechanism not robust to adaptivity; pair-of-columns argument (each j has a partner k with ~log M/loglog M hits; vertex cover >= r/3) plausibly restores rate. Not in write-up.
## step 4 (refines F4)
- Matching on ALL q forces p_j N(zeta_j) square (non-square is QNR mod some prime) -> N(zeta_j) in p_j Q*^2, |zeta_j|^2 >= p_j; then certificate m_(v,s) <= nu^-1 |gamma_v| ~ e^{w/2}: slab measure O(1), completion sum of Prop 6.6 diverges (needs m_v <= M^{beta||v||}, beta log M < tau/2; here beta log M = 2 tau). So full matching is incompatible with certified completion (zeta=pi_j is the actual circle, which pins phi). Author's density "2^{-|S|}" + "2^K >~ X/(log X log p_j)" is pseudosquare heuristic, not exact. Under GRH least q in a 2-adic class with (n/q)=-1 is (log n)^{2+o(1)} 2^{O(K)} -> matching on S impossible for n = M^{O(1)}; unconditionally only Linnik-type polynomial.
- Thm 2 re-derived: all steps OK. Cor 2: with C <= K' log M actually Hbar >= M^{2/3-o(1)} (stated 1/2 is weaker, true).
- Thm 2 vacuous for G = I (Hbar ~ bM/2) — only says pair certification by size fails; fine.
## step 5
- root_ref.py: 89 prime powers, 168562 collisions, all Q | Im(gamma^(2^k)); f odd label, f=+-i iff k=2, classes k=2,3 OK, odd-label order divisible by 2^k OK. (Im=0 cases = repeated zeta in random choice -> Z_v=1, consistent with 3.3(4)).
- tail_ref.py: M=64: RAND 8.96/3.73 COH 6.28/3.22; M=256: 10.72/3.82, 6.83/3.12; M=1024 (2 seeds): 11.14/3.86, 7.33/3.00 and 10.83/3.87, 7.24/3.00. Matches author's table within noise.
- snf_cross.py: sympy SNF confirms half-size row subsets Z-surjective (q0=43,47).
- F2 repair: pairs e_j - e_k with equal parity; greedy matching >= |J|/12 disjoint pairs with lambda_jk >= Lambda*/4 (density >= 1/3 graph). Unconditional in Legendre symbols, constant loss.
- Minor: Thm 2 "uncertainty" = Donoho-Stark/Elad-Bruckstein support uncertainty; cite.
- Minor: Cons 2 "stays at Theta(...)": upper bound for the hybrid not shown.
## step 6
- pairs_ref: M=44 5 seeds D1 max 4.22-6.13, mean .911-.977; M=132 max 6.01,5.56 mean .952,.961 -> author's maxima single-seed noisy.
- M4 sharpened: positive n = p_j N(zeta_j); quadratic chars mod 2-power: chi_{-1}, chi_2, chi_{-2}; n>0 so exceptional n in {1,2}Q*^2. Under GRH, n not in {1,2}Q*^2 has (n/q)=-1 for ~half of primes of each 2-adic class in [M^{1-eta},X] once X >= (log n)^{2+eps}. => matching any class k>=3 forces N(zeta_j) in {1,2} p_j Q*^2 => |zeta_j|^2 >= p_j => certificate >= trivial e^{w/2}: no certification. n=2*square matches ALL k>=3 classes (zeta=(1+i)pi_j) with repair only at k=2.
- M2 pair repair: k' threshold per pair with w_j+w_k <= T_M; matching n/12; W >= (eps/12 - o(1)) M T_M.
