# kill_lit working notes
## Normalisation check (from gauss.md)
- h_{-K} = W' = log R''^2 (Lemma 1.1, local heights at pi, pibar normalised 1/2 each); Galois-invariant toric divisors = E_A + E_{A^c}; on balanced points all toric heights ~ (sum a_A) w.
- A-point: lambda_E >= (1/4+delta/2) h_B. dist^gamma H bounded for gamma = 4/(1+2 delta). delta = 1/30 -> 15/4; delta=1/10 -> 10/3. OK.
- (L5circ) counterexample (eta0,c): alpha <= 15/4 + O(eta0+c) [exponent of diam <= 7/15 + O(eta+c)]. So target "alpha_ess,bal > 15/4" must be read: exists eta0 with alpha_ess on eta0-balanced > 15/4 + kappa(eta0). Imprecise but ok direction.
- Vojta D=E: lambda + h_K >= 4 lambda - h_B -> alpha_ess >= 4 - eps. OK.
- diagonal conic points are NOT distinct tuples -> not A-points at all (lit says "Z excludes anyway": minor imprecision).
## Findings so far (reading)
- F1 (MAJOR, math error): lit Lemma 1.1 / summary: "-K_X ample, X toric Fano (permutohedral variety is Fano)". For k>=4 the
  linear functions m_sigma = e_{sigma1} - e_{sigmak} on max cones are NOT strictly convex: ray 1_{13} lies outside chain
  {1}<{12}<{123} yet <m,1_{13}> = 1. So -K nef & big, not ample (k=4,5). Sibling toric.md says the same (correct). Fix:
  replace by "big and nef, positive on every curve meeting T"; harmless for negative verdict (toric heights comparable on T(Q)).
  VK attribution therefore wrong as stated. -> verify by script.
- Prop 1.3(b) logic: forced-divisibility section bound = width-based toric height h_Q = (1/2) sum_layers w_A log p; on
  balanced, (w/4) sum_A w_A / ord >= 2^{k-3} w = 4w. OK. Covers non-symmetric supports too (widths translation/sign invariant)
  although (a),(b) are stated only for Galois-invariant L. Minor. It is row 20 restated (honestly flagged).
- Only "forced" (profile) divisibility is covered; per-point auxiliary sections (Siegel lemma per point) and multi-point
  (Thue-Roth) not covered by Prop 1.3 per se; Roth-type covered separately by beta<=tau (for fixed S={inf}).
- "diagonal conics ... which any admissible Z contains": they are not distinct tuples => never A-points. minor.
- "conditional ones ... rest on Vojta, i.e. on V_E(5) itself": imprecise (Vojta for THOSE varieties). minor.
- §8 "gives (L_5^circ) for that delta" — (L_5circ) has no delta; bookkeeping OK (choose eta', c' small). minor.
- split dP6 "alpha_ess = tau = 3" asserted in §8 without proof in lit.md; I checked: x=(P1P2,P2P3,P1P3), P_i = n + c_i
  gives dist~1/n, H=lcm~n^3 -> gamma 3; union over c dense -> alpha_ess <= 3 = tau. True (my derivation).
## Computations
- fan.py/fan.out: k=3 ample; k=4,5,6 nef, NOT ample (strict-convexity violations 24/480/7920; e.g. sigma=id, A={0,2}, <m,1_A>=1).
  Pic^Gal rank 15 (k=5) OK; widths of P(-K) all 2; Thm-B per-m bound 15/4 OK; -K.C_v = 2(maxv-minv) OK; vertices = e_i-e_j.
- reg2.py (flint, binomial basis, last coords, primes 1000003, 998244353): k=3 m<=5: 3m; k=4 m<=4: 3m. Reproduces lit table.
## Searches
- Q1 "Voskresenskii Klyachko toroidal Fano varieties root systems permutohedral..." -> VK 1985 Math USSR Izv 24, 221-244:
  finiteness + classification of centrally symmetric toric Fano; no "permutohedral is Fano" statement found.
- Q2 '"permutohedral variety" "is Fano" if and only if dimension' -> snippet (Abe-Horiguchi-Masuda?, "Fano and weak Fano
  Hessenberg varieties", Michigan Math J 73(3), arXiv 2003.12286): "permutohedral variety ... not Fano unless n<=3, always
  weak Fano". Confirms F1.
- Q3 MR EJM Liouville: snippets confirm EJM 2016, arXiv 1306.2977; "generalize Liouville to arbitrary varieties, extension
  involving the asymptotic (also: stable) base locus; cubic surfaces all nef L, x not on a line". Exact hypotheses NOT seen.
  On P^1: alpha_x(O(1)) >= 1/d (Liouville) and >= 1/2 (Roth): normalisation dist^gamma H bounded consistent with lit.
- Q4 MR Invent: Invent. Math. 200 (2015) 513-583 confirmed. "alpha >= eps/2 (Thm 6.3)" NOT confirmed by my searches.
- Fetch attempts: link.springer.com, www.math.uwaterloo.ca: DNS fail (network-blocked).
- reg2_k5.out: k=5 m=1..5 max ord mod p (two primes) 3,6,10,13,16 -> reproduces lit exactly (|S| 21,131,471,1251,2751).
- vander.out: exact: Vandermonde shifted support (120 pts) in 3P not 2P; moments deg<10 vanish, some deg-10 nonzero: ord 10 exact.
- Q5 Huang ANT: arXiv 1809.02001 abstract confirmed (split toric, Salberger torsor, best approx on min-degree rational
  curves, McKinnon conj); "Def 4.7 essential constant first introduced by Pagelot [Pag08]" CONFIRMED by snippet. ANT 15
  461-512 not confirmed by my search (lit18 had same; [unverified pages]).
- Q6 McKinnon JAG 16 (2007) 257-303, DOI 10.1090/S1056-3911-06-00458-9 confirmed; Conj 2.7 content (paraphrase) confirmed.
  Snippet: Castaneda thesis extended to split smooth rational surfaces (lit's "San19" item probably this; lit unsure).
- Q7 McKinnon-Satriano TAMS 374 no.5 (2021) 3557-3577, DOI 10.1090/tran/8318 CONFIRMED; erratum 03-29-23: "Vojta does not
  imply canonical boundedness in general; does for toric; corrected Prop 1.3".
- Q8 Batyrev-Tschinkel IMRN 1995 no.12 591-635 CONFIRMED (arXiv alg-geom/9411009, not 9510014 as lit's primary pointer).
  Abstract: "all K-rational points of X are contained in T(K)"; Poisson formula; prove Batyrev-Manin conjecture.
- Q9 offline mirror (scratchpad/lit18/tools/db.json, 92926 items): 2403.02480, 2504.05213, 2509.03998, 2308.07911,
  2111.01509 (v3 2026-07-19), 2510.09954, 2509.05300 (math.GM) abstracts all match lit's descriptions. Keyword scan
  (approximation constant/essential constant/anisotropic tor/best approximation/McKinnon): nothing on anisotropic tori.
- Q10 Monahan-Satriano 2308.11847 confirmed (Vojta => curve through P better than any Zariski dense sequence, split
  horospherical; Vojta => canonical boundedness).
- Q11 McKinnon 2310.01640 confirmed (cubic hypersurfaces with a rational line; Coba).
