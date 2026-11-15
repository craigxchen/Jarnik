# notes, lane prove round 30
## verified so far
- S1: z_i - z_j = A(B - Bbar) = 2i A Im B (sign: +2i with B = factors of z_i on S). deg Im B <= w-L-1 = w - D/2 + c (even D). OK.
- S2: arg(t-rho) increasing -pi -> 0 (deriv y/|t-rho|^2). theta = sum delta_k phi_k. Im B = |B| sin theta.
  zeros >= |p-n|-1 by IVT => |p-n| <= w-L. Generalises to ANY integer kernel vector a: U_a = prod (t-rho)^{a+} (t-conj rho)^{a-};
  |sum a| <= ||a||_1 - L, i.e. min(a^+,a^-) >= L/2.
- S3: Psi_a = sum a_k phi_k; Psi_a' = sum a_k y_k/|t-rho_k|^2 = W/N, deg W <= 2D-2-L. eta_a = 2i Psi_a' dt.
- S5 horizontal line: rows span x^0..x^{L-1}; PTE => |T_i \ T_j| >= L (not (L+1)/2), all |T_i| equal; Bonferroni => M<=2 for D>6c+6.
- Cube locus: zero-sum sets F <-> (a,b), a=c+-c-, b=c++c-. m=2: sum b eps=0 combinatorial; m=1,3: sum a y=0, sum a(3/y-y^3)=0.
  Extra F (a not 0, not +-1) => 2 eqns on 2-dim surface: finitely many special members. COMPUTE.
- cube_special2.py: for base Q, extra zero-sum (a,b) with valid b: NO nondegenerate (|y| distinct) complex solution. So cluster through Q is exactly 8 on every nondeg member.
- cube_fibres.py running: all base points (d in {-2..2}^4).
- Symplectic form: m = a+ib <-> (a,b) in R[t]^2, omega(m_i,m_j)=Im(conj(m_i) m_j)=a_i b_j - a_j b_i.
  Pair condition (w = distance): Im(B) deg <= w-L-1. For a triple: blocks A,O_1,O_2,O_3 (odd one out), m_i = majority poly on O_i,
  omega(m_i,m_j) has deg <= beta_i+beta_j-L-1 and Pluecker: om23 m1 + om31 m2 + om12 m3 = 0.
  Equidistant (minimal) triple at c=0: omegas constant; kappa23 m1+kappa31 m2+kappa12 m3=0 (real linear relation, coeffs sum 0).
- Min-weight zero-sum F (w=L+1): q_F = P_F + i kappa_F; roots = P_F^{-1}(-i kappa): UHP count floor/ceil(w/2).
- S3 correction: near-identities hold for Blaschke beta_v = n_v/conj(n_v), NOT for N_v (s3_check.py: N+/N- = 1+337/36 s+...).
- cube special members: real ones all share |y| set {0.1238,0.9112,2.6097,3.3971}; max fibre 8 (cube_fibres2/3). generic 8.
- equidistant D=8: 5-point (rank 4) found, local dim 2, residual 1e-51 (evidence); 6-point: none nondegenerate (150 starts x3 classes).
- equidistant phase ordering lemma: D=0 mod 4, c=0: theta_ij in (0,pi) or (-pi,0) for all real t (Im B=kappa no zeros, p=n).
- doubled fakes (s,-s) with rows (x, x+e), e in annihilator: all differences balanced p=n, M ~ D/2 -> balance does not separate.
- Literature (abstracts only; arxiv blocked): LRT 2304.05438 Osaka J Math 62 (2025): smooth Fano X, Mor(B,X), closed V containing all non-dominant comps. LT 1702.08508, LST 1805.10580.
- Deformation on Y: d(L+1)+1 <= Bd_Y.C if K_Y+Bd psef (Iitaka, subvar of torus) -> dim Y <= 3 at c=0 (sketch); rigid configs have Y=C, invisible.
