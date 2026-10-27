# notes (round 16 antitheory)
## Literature located (WebSearch snippets only; arxiv.org DNS fails for WebFetch)
- MNV: Meka, Nguyen, Vu, Anti-concentration for polynomials of independent random variables, Theory Comput 12 (2016) art 11, arXiv 1507.00829.
  rank r = max # disjoint size-d sets S_j with |a_S|>=1 (multilinear, deg d). arXiv Thm 1.6: bound ~ d^{4/3} sqrt(log r)/r^{1/(4d+1)} + exp(B d^2 (loglog r)^2)/sqrt r ;
  published: polylog with exp(B d log d loglog r)/sqrt r. "non-trivial for degrees up to c log r/loglog r"; "no non-trivial bound possible at degree ~ log_2 r".
- Kane, structure theorem for poorly anticoncentrated polynomials of Gaussians, Ann Probab 45(3) (2017) 1612-1679, arXiv 1204.0543: diffuse decomposition; degree fixed; size bound effective but large.
- FKS: Fox, Kwan, Sauermann, Combinatorial anti-concentration inequalities, with applications, Math Proc Cambridge 171 (2021) 227-248, arXiv 1905.12142. Bernoulli, nonneg coefficients etc.
- KS20: Kwan, Sauermann, algebraic inverse theorem quadratic LO, Discrete Analysis 2020:12, arXiv 1909.02089: point prob >> n^{-1+o(1)} => close to low-rank quadratic form.
- KS25: Kwan Sauermann, Resolution of the quadratic Littlewood-Offord problem, Compositio Math 2025: O(1/sqrt n).
- JKSW: Jin, Kwan, Sauermann, Wang, Algebraic aspects of the polynomial LO problem, arXiv 2505.23335 (2025) (appendix Kwan-Sah-Sawhney): MNV bound improvable unless f close to special algebraic structure; robustly irreducible.
- GK: Grebennikov, Kwan, Geometric LO problems via lattice point counting, arXiv 2505.24699: bounded Chow rank (polynomial of boundedly many linear forms).
- Grebennikov, Optimal bound for the polynomial LO problem, arXiv 2610.08708 (6 Oct 2026) -- only via MathWorld: C r(F)^{-1/2}, absolute C; NOT verified; attribution to GPT-6 Pro session per MathWorld.
- MOO: Mossel, O'Donnell, Oleszkiewicz, Ann Math 171 (2010) 295-341, arXiv math/0503503: invariance for low-influence bounded degree.
- DSTW: Diakonikolas, Servedio, Tan, Wan, regularity lemma for low-degree PTFs, CCC 2010; Theory Comput 10 (2014) 27-53.
- Filmus-Mossel: Harmonicity and invariance on slices of the Boolean cube, PTRF 175 (2019) 721-782, arXiv 1507.02713: harmonic multilinear deg o(n): slice ~ cube distribution.
- Razborov-Viola, Real advantage (ECCC 2012/134): deg ~ (1/2)log2 log2 n polys uncorrelated with parity; uses Costello-Tao-Vu style anticoncentration with r disjoint max-degree monomials.
- Nguyen-Vu, Optimal inverse LO theorems, Adv Math 226 (2011) 5298-5319, arXiv 1004.3967.
- Carbery-Wright 2001 (Math Res Lett 8): |{|f|<=t}| <= C min(d,n) t^{1/d} (||f||_2)^{-1/d} on convex body vol 1 (as quoted in later papers).
- Becker, Klein, Slote, Volberg, Zhang, Dimension-free Remez inequalities and norm designs, arXiv 2310.07926; Klein-Slote-Volberg-Zhang ITCS 2024 arXiv 2301.01438 (dimension-free Remez on hypergrids; constants depend on degree).
- Glazer, Mikulincer, Anti-concentration of polynomials: dimension-free covariance bounds and decay of Fourier coefficients, arXiv 2108.04268.

## computations so far
- params_out.txt: CW proxy (optimistic C=1, W=W_min) nontrivial only from n ~ 2e4 (log10 bound +3.6 at 2900, +2.0 at 1e4 r=n, -6.4 at 2^16). MOO error >= m'(1 - ...) ~ m' always. MNV range log r/loglog r ~ 4-8.
- barrier_product_out.txt: Q = prod_{i<=m'} k(B_i)/sd, b = floor(d/m'), Pr(Q != 0) = 10^{-328} (n=2900), 10^{-2088} (n=1e6). p_b >= 0.4/sqrt b for 2<=b<=3000 exact; p_1=1/3.
- act_strata_out.txt: partition strata q_max >= 49 (n=2900) ... 86 (n=1e8) exceptional for AC_T(1/2).
- act_bits_theta_out.txt: bit strata (|X_i| = n/2, n=2^k, r=n-1):
  theta=1/2: q_max 11..21 (k=12..60); K*o>d for some q up to k=36; Gamma<=abar failures at W=3n ln n for k=12..22 (all checked).
  theta0 (profile-free mu(A) bound 1-(1/4)(4r)^{-1/2}): q_max 7..9; K*o>d up to k=19; Gamma-fail at 3n ln n/4n ln n for k=12..17.
  thetaS (Sylvester mu(A) bound 1-2r exp(-3w^2/8n)): q_max = 4 for all k (= A's cheap range); no failures.
  At W = W_min: no Gamma-failures in any case (k<=22).
## record facts
- determinant.md Lemma 3.2 proof: sum v_x delta_x = <2u,phi> - (pi/2)<v,kappa>; so for v in ker profile, sum v delta = 0 -> hat u in col span. Prop 5.2 low-rank invisible. Known.
- round13: Fano block cost 3j; 32-class span L_{2,1} in TC(A,m') every k>=5; coll_sigma false k=12..14.
