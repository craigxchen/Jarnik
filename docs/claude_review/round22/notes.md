# notes (round 22, toric)

## N1. Basic toric facts (to verify by script)
- N = Z^5/Z1, M = {sum m = 0}. Rays 1_A (30). Cones = chains. Galois acts by -1 on N: 1_A -> 1_{A^c}.
- P_{-K} = {sum x = 0, <x,1_A> <= 1 all A} = {sum x=0, ||x||_1 <= 2} = conv(e_i - e_j) = root polytope A_4.
  Vertex of max cone sigma (perm): e_{sigma5} - e_{sigma1}; non-simple => -K nef, big, not ample (k>=4).
- -K = H + H' with H = pi_1^* O_{P^4}(1) (simplex Delta), H' = pi_2^* O(1) (-Delta).
- width_v(R) = 2(max v - min v); B.C_v = 2(max v - min v) for 1-PS (Q-form: norm-one subtorus).
- eps_e(-K) = 2: upper by C_{1_A} (in diagonals), lower by H + H' superadditivity.
- tau_e(-K) = lim reg(mR cap M)/m, in [10/3 (Vandermonde, m=3, ord 10), 15/4 (construct Thm B)].
- If tau > 10/3: toric Liouville gives A_5(delta) for delta > 2/tau - 1/2 (all profiles), off base locus.
- Balanced Liouville: alpha_bal >= 15/4 (pair product, mean width equality).
- Class bookkeeping on Bl_e: gamma (E_A.gamma = 1/2 all A), -K.gamma = 15; l_E line in E (E.l_E = -1).
  A-points class gamma - mu l_E per w, mu = (15/4)(1+2delta). mu=15/4 counting/Vojta-critical, mu=4 Liouville.
  (K_Y+E).(gamma - mu l_E) = -15 + 4 mu.
- covering rational families: (k-1) mu_0 <= B.C => 4 mu <= 15 at k=5 (mu <= 15/4); at k=4: 3 mu <= 7, Liouville mu=2.
- expected dim of balanced rational curves of class a(gamma - mu l_E) with one-branch contact: a(15-4mu)+2.
  mu = 4, a = 2d: 2 - 2d. d=1 square system 28x28.

## N2. Results so far
- reg(mR_5) = 3,6,10,13,16 (m=1..5, reg.py); m=6: full rank at 20 (reg_cert_m6.out); m=7 full at 23 (reg_cert_m7.out).
  => t_max(m) = floor(10m/3), m<=7. k=4: t_max = 3m, m<=9.
- k=4 Thm 4.1 curves: E_A.C=1 all 14, mult 4, ratio 7/2, V-null (k4class.py).
- k=5 d=1 V-null balanced curves: system = 4 coupled ideal complex PTE (size 8, deg 7): Im p_n(R_j)=0.
  LM (k5curves, kpte 300 starts), Rabinowitsch (k5pte 300), symmetric (12)(34) ansatz (k5sym 400),
  complex Newton (400) & Newton homotopy (300 paths): no nondegenerate solution; all converged endpoints have c_j=0
  (R_j conj-closed: all points coincide). Inconclusive. Calibration: kpte at k=4 finds 55/200 nondegenerate.
- Argument principle: c_j != 0 => P_j has exactly 4 roots in each half plane (arg P_j in (0,pi)).
- alpha_bal(k) = 4(1-2^{1-k}); covering threshold k-1.
