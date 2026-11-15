notes round29 spectral
## Li-Liu pp1-20
- Lemma 3.1: G=D-P exact finite orthogonality (odd chars mod N); 0<=P<=D<<log N unconditional => trivially |G-cS|<<log N.
- GRH used only for holomorphy of L'/L(1/2+s,chi) on Re s>0 to shift contours (Lemma 3.5) -- no cancellation.
- Beta model = continuous HL main term over families m=n+kN, n=m+kN, n+m=kN with density S(kN); sum S(kN)k^-w = S(N)zeta(w)P_N(w).
- T_N kernel E_N(s,zbar)/((1+s)(1+zbar)) on two lines alpha_N=1/logN, beta=1/8; HS.
- r(N)=<T v,v> where v = Riesz vector of l_N(h)=(1/2pi i)int N^s Phi~(s)(1+s)h(s)ds; ||v||^2<<1. This is the definition, unpacked.
- X_N(q)=prod(1-4mu_j^2/delta^2 q(1-q)) = det(I-4delta^-2 q(1-q)T^2); zeros in Re q>1-eta <=> ||T||>delta/(2sqrt(eta(1-eta))) (exact dictionary p11). S <=> (2.28), R <=> (2.30) ||T||<=delta_N.
- GRH alone: ||T|| <= ||T||_HS << (log N)^2 (Thm 6.2) -- WORSE than trivial unconditional |G-cS| << log N.
- Hyp operator norm bound sup over all h is strictly stronger than the one test vector v_N => hypothesis implies conclusion and more.
## Li-Liu pp21-63
- Lemma 6.4 averaged L'/L bound: uses only diagonal of a Dirichlet polynomial (orthogonality) + GRH zero-count comparison. size = diagonal size. No GRH cancellation.
- Prop 6.5: P,D,K each HS << (log N)^2, triangle inequality. No cancellation between components measured.
- Prop A.4 (paper's own): fixed N, c->0: ||T_{N,c}|| >= z/c -> infinity; R(N,c) FALSE for small c. Norm hypothesis depends on auxiliary strip parameter.
- App B: N=30030 numerics: exact count 905 recovered (unconditional, = finite orthogonality); finite model ||C||/delta_N = 0.83 (delta_N=0.43). Authors: does not certify R.
- Sec 8: "central advance is the exact identification c S - G = <T v, v>".
- Model K at real s=sigma: (pi/N)^{2 sigma} ~ e^{-2a}, a=sigma log N; K ~ -e^{-2a}/(4 sigma^2): matches Goldston-Gonek-Montgomery diag/offdiag structure (memory; verify).
PLAN numerics: E_N(sigma,sigma) for N=2p, sigma=1/9 and 2/log N; lower bound ||T||>=8 d(s)|E(s,sbar)|/|1+s|^2 (5.8). Scaling in N.
## numerics E_N(sigma,sigma) (N=2p): sigma=1/9: E=0.52(202),0.24(2018),0.21(6002),0.23(20014),-0.06(60022); sigma=.05: 4.2,0.84,1.7,-0.52
  lb from (5.8) << delta_N. R not contradicted; no clear N^{-1/2} decay; GRH bound (logN)^2 ~100.
## Transplant ideas
- Hecke: lambda_k(N)=prod_p U_{e_p}(cos4k phi_p): angle-only -> fake-blind (identity holds for fakes).
- Aux modulus q>2W: (1/phi(q)) sum_chi |S_chi|^2 = sum_n r_I(n)^2 exact Parseval; = #strip + E_off; E_off>=M(M-1).
- E_off = sum_{P primitive, |arg P|<=Theta/2} #(Pbar Z[i] cap strip) (equal-norm pairs: w/z = P/Pbar). model ~ (W/R)(6/pi^2) C^2 log R-ish (to verify numerically).
- SCALE OBSTRUCTION: trivial M <= r2(N) <= exp((log2+o(1)) log N/loglog N) = R^{o(1)}; any input with Q^eps or exp(c logQ/loglogQ) loss is vacuous (Lindelof under GRH: exp(c logq/loglogq)). Li-Liu has N^{1/2} room; arcs have only subpolynomial room.
- ratios conj error O(main * family^{-1/2+eps}): W=R: X~C sqrt R, q~R -> error ~ C R^eps. vacuous unless eps=0.
