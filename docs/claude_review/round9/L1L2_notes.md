# L1/L2 working notes (not the write-up)

## N1. Reduced chord is a rational integer times a fixed small Gaussian integer
Normalised cluster: z_x = i^{k_x} prod_j pi_j^{a_xj} conj(pi_j)^{e_j-a_xj}, N' split-only.
v = a_x - a_y, A_v = prod_{v_j>0} pi_j^{v_j} prod_{v_j<0} conj(pi_j)^{-v_j} = X + iY.
g_xy = gcd = prod pi^{min} conj pi^{min(e-a)}; z_x/g = i^{k_x} A_v, z_y/g = i^{k_y} conj(A_v).
beta = z_x/g - z_y/g = i^{k_x}(A - i^s conj A), s = k_y - k_x mod 4:
 s=0: 2i Y ; s=1: (1-i)(X-Y) ; s=2: 2X ; s=3: (1+i)(X+Y).
So beta = i^{k_x} gamma_s B_xy, B_xy = D_s(A_v) in Z\{0}, |gamma_s|^2 in {4,2,4,2}.
Odd q not dividing N: q^a | z_x - z_y  <=> q^a | B_xy.  (q does not divide g.)

## N2. Margin identity
Norm g_xy = N e^{-d_xy}. t_xy := |z_x - z_y| N^{-1/4} (exact chord, <= C).
|beta|^2 = t^2 e^{mu}. So mu_xy = log|gamma_s|^2 + 2 log|B_xy| - 2 log t_xy.

## N3. Logical structure
* L2 (k window primes divide B_xy) => |B_xy| >= M^{(1-eta)k} => mu_xy >= 2(1-eta) k log M + log 2 - 2 log C.
* near-extremal => mean mu <= W/(2(M-1)) <= c0 M log M/(2(M-1)).
* L1 (K-uniform) + near-extremal => every |B_xy| <= C M^{K c0/4 (1+o(1))} => each B_xy has
  <= K c0/(4(1-eta)) + o(1) window primes.  So under near-extremality L1 => not L2 (quantitatively);
  i.e. Thm C is a proof by contradiction; neither lemma can be "true and informative" alone.
* Route 1 really needs (GM): mean over pairs of log|B_xy| >= f(M) log M. Exact: W >= 2(M-1) mean mu
  >= 2(M-1)(2 mean log|B| + log 2 - 2 log C).  Uniformity is a device to transfer max -> mean.

## N4. Code view of L2
For q in window, q not | N: z_x mod q in torus T_q (order q -+ 1); log: L_q(x) = sum_j a_xj ell_j(q) + c_q,
ell_j(q) = dlog(pi_j/conj pi_j) mod (q -+ 1)  (careful with units / i^{k}).
L2 = the linear code x -> (L_q(x))_q has two words agreeing in >= k coordinates.

## N5. Bisector lemma (lattice lines, outside.md Prop C.2)
z_x z_y = eps g^2 Norm(alpha) => midpoint angle = arg g + arg(eps)/2 mod pi => line spanned by
v = g or g(1+i), |v| <= sqrt2 |g|, lies inside the arc => (C.2, eta = 0) M <= sqrt2 C^2 |g_xy| + 2.
Normalised: mu_xy <= W'/2 - 2 log((M-2)/(sqrt2 C'^2)). Complementary pair only if M <= sqrt2 C'^2 + 2.
Data: records & families all have max mu = W/2 (complementary pair; they sit on symmetry lines).

## N6. Point totals
D_x <= sum_y (mu_xy + c_C) = W/2 - (1/2) sum lambda_s s_x sigma(s) + (M-1) c_C;
sum lambda |sigma| <= sqrt(W * 4Q), Q < MW/4 - (M^2/2) log M + O(M^2) (stability 4.1(c))
=> near-extremal: lambda_max(D) <= max D_x <= (1/2) sqrt(c0(c0-2)) M^{3/2} log M (1+o(1)).

## N7. Exact relation lattice
phi(c) = sum c_x z_x on {sum c = 0} (rank M-1) -> Z[i]; kernel rank >= M-3, covol <= sqrt(M) * M C^2 N'^{1/2}.
Minkowski: relations with |c| <= sqrt(M) (M^{3/2} C^2 e^{W/2})^{1/(M-3)} = M^{c0/2 + 1/2 + o(1)} near-extremal.
Content = lattice-ness of points (Pluecker); not used further.

## N8. Fraction barrier
Markov: #{mu > K mubar_+} < binom/K automatically. Thm C needs a heavy pair in the K-good set.
Poisson: heavy-pair frequency at L2 threshold is M^{-2+o(1)}; so uniformity is needed for all but M^{o(1)} pairs.
