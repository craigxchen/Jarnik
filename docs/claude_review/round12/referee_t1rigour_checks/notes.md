# referee notes (t1rigour)
Checked by hand: Thm 2.1(a)-(f) proofs, Lemma 3.1, Thm 3.1, Cor 3.2 algebra, Remark 3.3, Lemma 4.1, Thm 4.2,
Props 5.1-5.3, Lemmas 6.1-6.2, Thm 6.3, Cor 6.4, Prop 6.5(a,b), Lemma 6.6, Thm 6.7, Thm 6.8, Lemma 6.9 proof,
Lemma 6.10, Thm 6.11 (a),(b) logic, product TC formula, J2 Chebyshev.
FINDINGS so far:
F1 (major): Prop 5.5 uses o_F <= m' - 2js - 2e, i.e. rho_A(F) >= 2(js+e). Unproved; false for near-twin columns:
   S1,S2 sign columns, |S1\S2|=|S2\S1|=delta: image {(t1,t2): |ti|<=j, |t1-t2|<=2delta}, Newton => rho <= 2j+2delta,
   while (s,e) of F=(1_S1,1_S2) is (2,0) [claims >=4j] and of F'=(1_S1,1_S1-1_S2) is (1,2delta) [claims >= 2j+4delta].
   Also rho(Lambda) can be < rho(projection) (example {(0,0),(1,1),(2,0)}), so even rho_A(F) >= 2j for s>=1 not immediate.
   => 'n >= 5185 scale compatibility' unsupported. Recompute with only proved bound (o_F <= m'-1).
F2 (minor): Thm 6.11(a) toy: zeta(A2(4,1) x cube^2) >= 8 >= m (not zonotope-free), w not from formula (w=sqrt(48 ln4)=8.2),
   m < 4w. Only (b) is regime evidence. (b) zonotope-free by Hadamard L1 argument: sum L_i <= 2w r/n = 2w < m'.
F3 (minor): §7.1 T=1 J1 cost: single direction 1_{V_c} = (1/4)(1+e1 chi_a + e2 chi_b + e1e2 chi_{a+b}) has width
   <= (lambda m' + 4j)/2 on T=1 A_{Q,lambda}; family (b,a+b,1_{V_c}) cost <= 6j + lambda m'/2 (+O(1)) -> cheap for m' > ~12j
   (lambda=1), in window from n=2^13, not 2^15. 'iff m' > n/2+4j' false.
F4 (minor): Def 5.4 needs finite set; (Q2) uses all rational subspaces of TC_H (possibly infinite). Fix: maximal ones / closure.
Verified exactly: m_max table [3,12000], thresholds 2899/3622/5808 last failures, 2966/5930, interval constants,
 m_max(2^12..2^16)=2631,5260,10518,21034,42065; Lemma 6.9 s<=6; toy dims 193/96/33/6 and zeros; zeta(toy)>=8.
Consistency: author's script reproduces 5184/3661 but uses exponent m'-d_F; with proved rho_A>=2min(j,n-largest class)
 (zonotope cert, matching) threshold ~255100 (coef), ~250600 at W_min (r=n,C=1, double).
F5 (structural): ZE-TC^rat trivially true at rational zeros (R u is a rational subspace of the cone). (Q1)+(Q2) imply
 ZE-TC_eps relative to finite family L_coll = {Delta_pi : P >= binom(n,2)/Q0}; that single hypothesis suffices for Thm 6.8.
F6 (minor): Prop 6.5(c) 'iff k>=16' verified only k<=22 by author; tail is easy (k w < 0.6417d-2.2137 for k>=17).
