# refJ working notes

## J5 re-derivation (done in head, to be checked by script)
- u-coords: u_base = 0; <k',y> for y = 1_{S} - 1_S(base) equals k(S) since sum k = 0. Sylvester column a: <k,beta_a> = 2k(S_a), |k(S_a)| <= j.
- Lemma 5.0 (nearstrata): rho = any interpolation degree of pi(A); image of F_b in {-j..j}^b so rho <= 2jb; o_b = m'-2jb.
- 8j < m' automatic (w irrational => 8 floor(w/2) <= 4w < ceil(4w)).
- Case A: any L in TC (o>=1) with P(L) < N2/Q0: generic y in L lies in no Delta_pi with P(pi) > P(L) (would need pi refine pi(L)).
  dist(y,Delta_pi) >= gap/2. Fails once sigma <= gap/2.
- Case B: dist(u, Delta_pi) >= eta (n-1)/(Q0 n) for refinements of pi_1; mixing >= (1-eta)/2; normalisation |u| <= 1+eta.
  failure iff abar >= g*(Q0) = 2[m' ln4 + o ln(cQ0) - ln(d/(m'+1)) + ln(1/C)]/(o Q0 - d), c = n/(n-1).
  d g*/dQ sign = o(1 - ln cQ) - d/Q - a < 0 for cQ >= e. (a includes ln(1/C) >= 0 when C<=1.)
  eta at threshold <= exp(-a/o) automatically; actually failure ineq itself forces eta^o <= eps0 4^{-m'} <= 4^{-m'} if eps0<=1.
- C-dependence: README says "C <= 1" but g* formula is C = 1. For C < 1 thresholds shift by 2 ln(1/C)/(oQ-d). Minor correction.
- MIN* => MIN needs C <= 1 (C' = 1): MIN* rhs (d/(m'+1))e^{-d g/2} at g = abar = W/(2d) + 2 ln C' equals d e^{-W/4} C'^{-2d}... check.
  Actually d*abar/2 = W/4 + d ln C'. so rhs = d e^{-W/4} C'^{-d}/(m'+1). For C<=1: = d e^{-W/4}/(m'+1) >= eps0. OK.

## J7 hierarchy: checked in head, all inclusions OK (G subset {Gamma<=abar} subset R_per subset R_unif).

## J5 DONE (section 1 written). j5_recheck_out.txt all confirmed.

## J6(a) notes
- Fat point at [v] lives in u-space; D_v^{tau+1} Phi = 0 <=> mult >= m - tau at [v]. tau_A(v) = #values - 1 <= h_A(v) (strengthening).
- linear span of v_1..v_s in Bs iff (sufficient) sum h < m: Phi|_Lambda has deg_{t_i} <= h_i. Converse false in general
  (BDP linear cycles are BY DEFINITION spans of subsets; a line with >=3 fat points of total excess > m is in Bs too).
- "exactly" = definitional match with BDP linear cycles only.

## J6(b) notes
- sl2 proof: V = k[x]/(x^3), E = x, F(1)=0,F(x)=2,F(x^2)=2x, H = diag(-2,0,2). B = tensor, E_tot = mult by sum x_i.
- telescoping: h_i >= h_{i-N} for i <= s+L, <= for i > s+L. total = sum_{|i-s|<=L} h_i.
- signs: automorphism x_i -> eps_i x_i. disjoint cuts: product, dims multiply.

## J6(c) plan
- M8: base 0; columns c_a = -1_{T_a}, T_a = {x != 0 : a.x = 1}; cut |k'(T_a)| <= 2.
- Fano line triple {a,b,a+b}: c_a+c_b+c_{a+b} = 2 * (-1_{nonbase \ {a x b-common pt}}), half is integer direction z; width of z <= 6.
- tests: exact J(A), exact J_w, exact J_w+col pairs (deg 6), triples sandwich, Fano vs non-Fano triples, extended pairs.
- tau vs h: primitive v with max|v_x|>=4 has tau >= 6 (values 0,+-M,+-b,+-(M+-b)); so |v_i|<=3 enumeration complete up to degree 6.

## M8 results so far (exact)
- |A|=1037, H_A=[1,7,28,77,161,259,287,175,42], dimJ(A)=[0,0,0,7,49,203,637,1541,2961] exact Q.
- width set: 7 (h=2), 49 (h=4), 302 (h=6); no tau<h; max|v_i|=2. J_w exact Q all deg<=8: [..,616,1485,2940].
- column pairs: image = full grid {-2..2}^2, J(Y)=(y1^5,y2^5), rho=8 (NOT < m): pullbacks are width powers -> trivially nothing.
- non-line triples: full grid 5^3, J=(y^5), add nothing. Fano-line triples (7): |Y|=57, rho=6<8, gens deg5 (6) deg6 (3).
  J_w + 7 Fano triples = J(A) in all deg<=8 (exact Q deg<=6, mod p 7,8 -> sandwich needed).
- planes inside Fano 3-space U: contributing planes = pencil through z=(c_a+c_b+c_{a+b})/2; gain 2 per line; total 630 (R<=5). 7 short.
- U_l in Bs(L_8) via 6 rational collinear triples (plane-curve peeling) though no cheap spanning family (all dirs in U have h>=4).

## step7: 51933 planes from pairs of 358 low-width dirs, 3213 contribute, total 630 < 637 (exact). b0=2 fails on that family.
## NEW (big): parity box P_j = {y in [-j,j]^3, sum even}: rho = 3j (exact j<=4). PROOF of rho<=3j: slice by z=(y1+y2+y3)/2,
  order z=0,1,-1,2,...; slice H_m = hexagon, rho(H_m) <= 3j-2|m| by diagonal (y1+y2=c) slicing ordered by |c|. Peeling lemma:
  rho(Y1 u Y2) <= max(rho Y1, deg g + rho Y2), g=0 on Y1, nonzero on Y2.
- Sylvester Fano triple: pi_U(A) subset P_j (parity: sum = 2k(S_a cap S_b)), so o_U >= m' - 3j, P(U)=P(L_2) (4 cells).
- per-subspace: k=12: oQ-d = 1413*4.0029-4095 = +1561; k=13: 2066*4.0015-8191 = +76 >0 -> per-subspace J5 extends to n=8192!
- next: 8-cell XOR family rho (M8-type box with unbounded integer k), c_3 constant.
## persub_fano: k=12 Fano span U fails per-subspace from abar>=5.26 (W>=0.96 W_min) -- all near-extremal W! k=13 from abar>=158.
## Combos: s Fano blocks + c columns (independent a's): rho <= (3s+2c)j (product of P_j's and intervals), cells 2^{2s+c}.
 - per-subspace best (2,0): o=m'-6j, Q~16: k=13,14 fail (proved). k=15 no.
 - uniform: (2,1) has cost 8j<m' and 32 cells -> Case A for Q0<Q5~32; Case B then needs o_1 Q5 > d: extends to k~20.
 - t=3 XOR block rho: j=1:3, j=2:<=8 (mod p upper bound rigorous: rho_Q <= rho_p). conjecture 4j. not used.
