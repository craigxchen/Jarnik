# Referee notes (ze.md)

## Lit
- Philippon Bull SMF 114 (1986) 355-383 confirmed by search (eudml 87515, numdam bsmf.2060). arxiv/numdam fetch: arxiv DNS blocked.

## Lemma 3.1 check
- Step 1 re-derived: Norm F >= N^{jD} prod p^{-j h_S(a_p)} (min over supp; inert/1+i full power); |F| <= N^{jD/2}|c|_1 (DjTheta/2)^mu/mu!.
  Gives (mu - jLambda)W/2 > 2ln|c|_1 + 2mu ln(DjC'/2) - 2ln mu!. MATCHES (3.1) at s=0. OK.
- Step 2: theta^alpha P_c = P_{c k'^alpha}; ord drops by |alpha|; |c|_1 grows T_A^|alpha|. OK.

## MAJOR GAP 1: "route 5 certifies ONLY powers Z^j"
- Same argument applies to y_x = prod_y z_y^{M_xy} (conj for negatives), any integer n x n M with equal |row sums|
  (or pad by N powers): point Y = Z^M := image of Z under endomorphism M of G_m^d (reduced coords),
  cost Lambda_M(c) = sum_p omega_p h_S(M beta_p), angle bound |<k,M delta>| <= |M| D Theta/2.
- Examples: (a) sigma in Aut(profile) (perm preserving weighted multiset {beta_p} up to sign): cost = Lambda(c);
  Sylvester: AGL(k,2) (size ~ n^{log2 n}); and A, V_H sigma-invariant -> P_c(Z^sigma)=0 for all c in V_H already.
  (b) rank-deficient profiles: v in left kernel of [a_q (q != p), 1]; M = 1 m^T + u v^T: Y = (Z^v)^{u'} lies on 1-dim
  subtorus t -> t^{u'}, cost = omega_p |<v,beta_p>| h_S(u) (one light prime) -> long progressions (Z^v)^{j u'}.
- Prop 3.5 as stated assumes Sigma subset <Z>; Thm 3.6 application counts only Y^j. Interpretation "all route-5 data"
  unproved. Need: is barrier robust to Aut-orbit? compute.

## Checked OK
- Prop 3.3 proof (Omega(S)=E_eps h_S(eps), construct.md Cor B1) correct. Cor 3.4 correct (W>=40).
- Cor 3.2 J_lemma table: hand-checked n=2966 E=1.69, J values all consistent with m-2jw>E.
- Prop 3.5(1): Sigma(d) contains d-multiples -> d<=2J. At T=1, T_0=0 forced (d T_0 < m_max < d). margin ~ ln(d! 2^d).
- Prop 4.4 numbers: 22371, 19114 hand-check OK.
- Prop 4.5 hyperplane ratio 2d/(m-1+d)>1 ok; uses "#distinct values of <k,a_p> on A = L+1" (interval property unproved; minor).
- Lemma 4.2 ok; Prop 4.3 ok (degree-level, geometric estimates).

## Minor
- Prop 5.1 proof: "|<v,delta>|<2pi" false for large v; conclusion true by unique factorisation + affine independence.
- Prop 5.2 constant: from 2e^{-X/2} and |delta'|_inf <= Theta (differences) get 2/(C'|u|), not 4/(C'|u|).
- Thm 3.6 numbers only at n in {2966,3000,4096,8192}, not "n=2966..8192"; monotone in mu OK.
- Thm 3.6 witness c is complex, not integral: no conclusion about small-value estimates (need integer P, heights).
- lambda=n/(n-1) = max possible lambda (trace bound) = best case; no Hadamard at n=2966 (2 mod 4); summary omits profile restriction.

## MAJOR GAP 2: dual picture uses V(A,m), route 5 knows only V_H
- Record: V_H = span of integer c in V(A,m') with |c|_inf <= H=U^2; dim V_H >= |A|/2 only ([D] l.142, 287; R12 README l.96).
  dim V(A,m) >= |A| - binom(m+d-1,d) ~ (3/4)|A| at m_max. So V_H may be a proper subspace.
- Rigidity gives Z in Bs(V_H) only, NOT Z in Bs(V(A,m)). Hypothesis "Q vanishes on phi_Z(A)" is not known at Z.
- ZE (route 5) = Z notin Bs(V_H): STRONGER than Z notin Bs(V(A,m)). So "ZE^all(A,m) => ZE" (l.548 "implies ZE^all, hence ZE")
  is unproved; needs V_H (x) Q = V(A,m) (x) Q, not on record (and Siegel basis heights astronomically above Liouville budget).
- Correct target: ZE^all_H: Bs(V_H) has no nondegenerate point. Not group-structured.
- Effect on verdict: the dual picture is even less applicable; negative conclusions survive; map statement must be corrected.

## Lit: eudml, arxiv, numdam? DNS blocked (arxiv, eudml). Shape matches my memory of Waldschmidt 2000 Thm 5.1:
  binom(T+l*,l*) card((Sigma+G*)/G*) H(G*;D) <= H(G;D), vanishing order > nT on Sigma(n), H(G_a^d0 x G_m^d1)=(d0+d1)!/d0! D0^d0 D1^d1.
- report sections 1.1-1.9 written; next: 2 recomputations, 3 lit, 4 issues, 5 corrected, 0 verdict
