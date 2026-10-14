# notes for fatpoint.md (working; not a report)

## Sources read
- round12 README (all), fake.md §1-2, §5; refJ.md §0, §2 (fat-point inclusion; Fano spans in Bs(L_8));
  loctheory.md (all: alpha_F criterion, parity locality Thm 3.4, string thm, Sylvester images Lemma 5.1,
  complete blocks Y_t non-local, rho(Y_t) bounds); loccompute.md §0-2 (Paley-12 Loc fails b<=5 over Q).
- H = U^2, U = |A| (t1rigour l.79).

## Key ideas (parked)
1. TC(A,m) = {u : l_u^m in J(A)_m}. Bs of local system = {u : l_u^m in sum_F (I_F)_m} subset TC unconditionally.
2. Single family: Bs(V_F) = {u in span F : coords in TC(B_F,m)} (tensor/bigrading argument). Linear => degenerate.
3. Joins: u = sum u_F, sum nu_F(u_F) < m => u in Bs (multinomial).
4. Rational-span: lin. forms l vanishing on W = sum span F: l^m in every V_F. So Bs(subsystem) subset W.
   Nondegenerate pt => Q-irreducible component Z spanning everything.
5. PRODUCT CERTIFICATES (generalised zonotopes): prod_i <z_i,u> in J(A)^perp_m iff Pol(f)(z_1..z_m)=0 for all f in J_m.
   Then TC subset union of hyperplanes {<z_i,u>=0}, and kappa >= |prod|/(|c|_1 d^m).
   - b0=1 (width) condition: count(v) := #{i : <z_i,v> != 0} <= tau_A(v) for every v (necessary always, since J_w subset J).
   - family F: restricted a_i = pi_F(z_i); sufficient: N_F := #{i: a_i != 0} <= alpha_F - 1 (+ width counts).
     necessary: N_F <= rho(B_F)  (J(B)_e = S_e for e > rho).
   - braid products impossible in Sylvester-type: each pair split by ~half the columns => column average m/2 = 4j > 2j.
   - coset differences 1_{x+H}-1_{y+H}, |H|=2^s: column load m/2^{s+1}; but subgroup-indicator directions 1_{U^perp}
     (tau ~ 4j(1-2^-t)) are overloaded at t = s+1 (load ~ m/2 = 4j). Mixed s_i might fix on average.
   - complete t-blocks: restricted product lives on class-sum box C_t (self-similar). Zonotope-sufficient condition there
     is X-type (sum |z(S_a)|) for generators with H_i cap U^perp = 0. If alpha(Y_t) ~ 4j for all t (heuristic), and
     generators touch most large blocks, the N_F <= alpha-1 criterion fails; need finer J(Y_t)^perp analysis.
   - Arithmetic: X_z - W/2 ~ abar for coset differences in Sylvester (Plotkin-tight) => product bound is cheap.
     CAVEAT: non-pair characters may vanish exactly for actual clusters => degenerate case; need family of certificates.
6. BLOCK-PRODUCT (tower) reduction, Sylvester: A contains prod_b Syl_{c_b+U}(N0, j_b) with sum j_b <= j (disjoint cosets).
   nu additive on products, kappa transfers (sub-box monotone via duality; product c = tensor). Needs blocks with
   m_b/j_b > 8 (LP-supercritical blocks: N0 >= 32 per zetc). Reduces Sylvester MIN^pair to ONE fixed block statement.
   For generic profiles, small blocks have full traces => block box = small ball, rho/j_b ~ 2 => useless.
   Full cube on X needs |X| <= 4j (pair-split count) but rho = 2(|X|-1) >= 8j needs |X| > 4j: LP-borderline failure.

## Computations so far
- BUG fixed in fplib.tangent_space_mod (int64 overflow in V @ (Em.T % p)); first sdc7 runs (m=6,7) invalid, rerun.
  Validation: zcheck.py (explicit zonotopes are found in T_m), zdebug.py.
- sdc7 L=2 m=8: no product in pool passes width counts. m=7: braid products (7 pairs) in T_7 EXACTLY (= braid
  zonotopes of length 7 found by referee_fake). LP-subcritical, consistent with Thm 3.1 (needs m >= 4w).
- yblock t=2 (single line, parity): count-type up to N = jc-1 (Prop 4.1), at e = jc only centred monomials
  (lambda_j sign: j=3 kills (3,3,3)).
- yblock t=3 j=2: e=5: all; e=6: all parity-feasible (lines<=5) heavy products in J^perp regardless of weighted cost
  (count-type); e=7: heavy2 products in J^perp only with max axis count 2 (weighted 4 = 2j): 4/5; count>=3: 0/~3000.
  light: 24/394 at e=7.  Exhaustive heavy-only runs at e=7 launched (yblock_exh).
## Dichotomy reasoning (parked)
- If blocks impose weighted (zonotope) cost on generators heavy in the block, average over blocks containing column a:
  fraction of blocks U (dim t) with H cap U^perp = 0 is ~ prod_{i}(1-2^{i-t}) -> 1 for t >> s. Then
  sum_i |z_i(S_a)| <~ 2j/c per column => sum X_z <= 2jW/c vs >= m'W/2 = 4jW: contradiction if c > 1/2.
- With FIXED b0 (blocks t <= log2(b0+1)), very heavy generators (s > t_max+1) rarely touch small blocks (Fourier
  support n/2^{s+1}, U cap H^perp = 0 generically) => escape. So: Loc(b0) fixed => no obstruction seen;
  blocks of growing t needed (Loc fails) => weighted cost => certificates impossible (if weighted-cost holds for all t).
## Tower theorem (parked derivation) -- see fatpoint.md §6.2
- locsys: M8 m=8: b0=1 deficit 21, b0=3 exact (validates against judge J6c). sdc7 m=8: Loc_8(1) HOLDS EXACTLY
  (J(A)_8 = (J_w)_8, H=1) => TC(sdc7,8) = Bs_R(L_8) = {Phi=0}, which has certified collision-free regular zero
  (fake §3.2, sep >= 0.2496) and nearby nondegenerate points (R11 Prop 6.2). So a width fat-point system CAN have
  nondegenerate real base points (underdetermined, LP-subcritical, ratio 2).
- yblock_exh e=7 lines<=5: heavy count>=3: 0/88200 (rigorous), count 2: 1536/3072. mixed: 0/4705 above zonotope cost.
- Theorem 6.4 (tower) + Lemma 6.5 (BLOCK <=> ZE-TC^braid of one block, Lojasiewicz) + Cor 6.6 written.
- Final: small2 classification: sdc7 m7 24/24, pal7 m7 57/57, syl8 L2 m7 696/696 certificates, all zonotopes.
  pal11 pairs-only m=8: 0 multisets pass width counts. M8' m=10,11: Loc(1) exact.
- yblock_exh e=7 lines<=6 run did not finish within the session (not used).
