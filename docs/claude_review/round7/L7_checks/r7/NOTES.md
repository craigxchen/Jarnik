# Round-7 L7 working notes (scratch, not the write-up)

## Plan
1. cub13 in Kapranov psi_7 model, P^4, p1..p5 = e0..e4, p6 = (1,...,1).
   Curve c(t) = (h, X2, X3, X4, X5) cubics. Project from p1 -> C' in P^3 (X2..X5).
   Planes with 1 -> lines in P^3: 123,124,125,134,135,145 = 6 coordinate edges => X_d = lam_d prod_{e!=d}(t - r_de).
   126: X3=X4=X5 at u ; 136: X2=X4=X5 at v.  v = other root of psi(v)=psi(u), psi=(s-r24)(s-r34)/((s-r25)(s-r35)).
   Planes without 1 (all contain 6): 246:{h=X3=X5}, 256:{h=X3=X4}, 346:{h=X2=X5}, 356:{h=X2=X4}, 456:{h=X2=X3}.
   Contact params: third roots (rational) except 456 where root w is free param and r23 solved (Moebius).
   Final: 5 interpolation conditions on cubic h (4 coeffs) -> divided-difference DD = 0.
2. Check if DD linear in some param -> W Q-rational.

## Results so far (session 1)
- kap1/kap2: param exact (sympy). Gauge u=t126=inf, r45=0, r24=1. Free (r25,r34,r35,w). r23 Moebius-solved.
- kap3/kap5: DD numerator on coordinate lines: genuine factor G has degrees (7,5,5,5) in (r25,r34,r35,w);
  ALL roots of G genuine (19 contacts, 13 planes), all other factors degenerate. (NUMERICAL genuineness)
- So W = {G=0} birational to cub13 family; projection to any 3 coords has degree 5 (or 7).
- cub12 family (drop plane 456): Q-RATIONAL 4-dim, explicit (r,u) param, h by Lagrange.
- fibre126.py: contact map W -> E_126 has degree >= 4 (numerical, 4 points found).
- Torus (Hadamard) structure: Kapranov psi_7 with p1..p5 coord pts: T = torus; 30 toric divisors,
  10 pair hypersurfaces x_i=x_j, 16 hidden (p6-centres). Toric-enrichment lemma: C (.) f^{e_A}, f=1 at hidden
  contact params, preserves hidden contacts, adds deg(f)*gamma_A. Sequential-linear count: balanced needs
  39q hidden conditions vs 30q capacity: deficit 9q (same as -9q lifting).

## Session-1 conceptual conclusions (to write up)
- Kapranov psi_7 model with p1..p5=coord pts, p6=1: point (c1:..:c5) <-> configuration (c1,...,c5, 0, inf) (labels 6->0, 7->inf).
  So a cubic curve c(t) <-> y(t) = (h, X2, X3, X4, X5, 0, inf).
- Enrichment lemma (torus T_ab): C -> C (.) f^{e_A}, f = 1 mod (hidden contact poly), class + d*gamma_{a,b,A}.
  gamma classes: f_i (B+8K=-10) or |A|=2,3 type (pair 6, triple 2, B+8K=-8). Cone of enrichments excludes beta.
  qbeta = gamma0 + sum d gamma_enr forces gamma0 ~ (1-small) qbeta: base family must already be large & triple-rich.
- Hadamard products: hidden contacts only intersect (monotone). Dilution argument for P1xP1 product charts.
- Mover deficit: free mover in anchor frame: 25q join conditions vs 16q+1 capacity (k=7); k=6: 10 vs 8; k=5: 3 vs 4.
  Boolean (pencil) movers free but r Boolean movers -> (Z/2)^r cover, genus 5 for r=4.
- Toric sequential interpolation capacity: sum over distinct last splits (2d+1) + 4 >= 39q needed; <= 20q+14. Impossible.
- FIBRATION CRITERION (to prove): Q-conic bundles phi_k: X -> Z_k with fibre classes g_k, beta in cone(g_k)
  => iterated sections (ruled surfaces phi^{-1}(phi(C)) with section C) give Q-dense classes g0 + sum m_k g_k => (L7) false.
  Fibre classes have -K.g = 2 (pair 6). Rich ones need triple>=10, all psi-degrees>=3. cub13 is NOT a fibration (N>=8).

## Session-2 (after exact cub13 verification)
- verify_cub13.py: EXACT member over quintic field (r25,r34,r35)=(3,-2,5/3), w root of
  329715w^5-8401893w^4+44833850w^3-89897718w^2+73951827w-20297381: class cub13, 19 distinct transversal
  contacts, W smooth there, P^1xW -> X differential rank 4. => Theorem 1 done.
- projdeg.py: projections (r25,r34,phi) of W have degree >= 4 for all natural phi (numerical).
- CORRECTION: fibration fibres need NOT have -K=2 (lines through a general point: -K=5).
- KEY NEW IDEA: "addability": family V addable iff ev: U_V -> X has a rational section over Q.
  Fibrations: ev birational. cub13: ev degree N>=8, monodromy transitive -> not addable.
  Rich 4-dim family Q4: quartics with all 10 coordinate planes (c_m = lam_m prod_{e!=m}(t-r_me), K5 edges)
  + 7 hidden planes (clusters = triples of [5] not containing {1,2}); pair 9, triple 17, B+8K=+2.
  If Q4 Q-rational and ev_{s0}: V -> X birational (fixed intrinsic parameter s0) => addable => (L7) FALSE
  (beta in cone(sym(Q4), sym(f)) since (21,35) in cone((9,17),(6,0))).
