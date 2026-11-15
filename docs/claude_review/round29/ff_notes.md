# ff lane notes (working)
- record check: round18/theory.md 342 (covering families genus<=1, function-field), L8.md 690 (ff analogue of L8, TSMT no contradiction),
  families.md §4 Lemma P (polynomial templates Z[i][u]; power sums; exhaustive table D<=8: max 6 bounded C at D=4, 5 at D=6, 4 at D=8).
- Setup: z monic WLOG (leading coeff units must agree). Squarefree normalised N: roots rho_k (upper half plane), points <-> s in {0,1}^n (box prod[0,e_k] in general).
  Arc: deg(z_i-z_j) <= K:=floor(D/2+c); L := D-1-K equations p_m, m=1..L (Newton). Condition sum_k delta_k Im(rho_k^m)=0.
- Coordinator (i) off by one: equations m<=L=D-1-K (= D/2-c-1 for D even).
- Newton: two multisets of size w with equal p_1..p_w are equal => every nonzero difference has ||delta||_1 >= L+1 (stronger than Vandermonde 2s>=L+2).
- L<=2: moment columns arbitrary (y>0, 2xy): problem = Littlewood-Offord.
- cyclotomic t^{2D}+1: M=4 (D even), 2 (odd), D<=20 (brute.py).
- zero-sum reformulation: base z_1=Q; cluster <-> subsets F of roots of Q with sum_F v(r_k)=0, v(r)=(Im r^m)_{m<=L}; = V' cap {0,1}^D, V'=ker W linear.
- generic count r*L <= 2D-2 => rank<=4 for D>=7; LP (rank4codes.py): rank<=4 & dist>=D/2 => M<=8 (rank3 => 6).
- cube (3 coordinate splits of F_2^3) realised numerically at D=8,c=0 (zerosum.py): M=8 numerically, 2-dim family mod gauge.
- c=0 D=2d minimal structure: z=P_i G_i, P,G in R[t]+iR (Im constant), s factorizations, M=2s+2; equivalently N=(X-2mu_i)^2+S_i^2 (real parts differ by constants).
  count s(d-1)<=3d-1 -> s<=3 (d>=4).
- L2 search D=6: max 4 for integer vectors box<=3 (families.md had 5 at D=6 for Gaussian-int roots? check).
- EXACT: symmetric cube family: roots +-r_k (k=0,1,4,5), x_k y_k = eps_k (+,-,-,+), sum y=0, prod y=1 (cube_rational.py: numerator -3(a+b)(a+c)(b+c)(abc(a+b+c)+1)).
  Q = Qtilde(t^2) with Im w_k = +-eta (Sperner middle layer, 6 pts) + extra face {r0,r1,r4,r5} (2 pts) = 8. Within the family M=8 exactly.
- gauss8.py B=7: Gaussian-integer roots: max 2 splits (M=6) — consistent with families.md polynomial record 6.
- cube_exact.py: EXACT over Q(i): y=(-9/4,-2/3,1/4,8/3) etc: cluster size exactly 8 among all 256 transversals, max deg diff 4, N squarefree deg 16.
- cube_Fq.py: over F_q (q=19,23,31,43) exact cluster size 8 (q=7,11: no nondegenerate y in F_q).
- Theory: image curve in T_Gamma (rank r), -K.C = 2D, contact m=L+1; dominating => r m+1 <= 2D (generic freeness, char 0) => r<=3 at c=0.
  ff-Vojta (truncated) would give (r-1)m <= |S|-2 => r<=4/5. Exceptional set = rigid families (cube). Wronskian gives only m <= (r+1)/2 (|S|-1): useless.
- Mason on pairs vacuous for squarefree (= record's abc item 6). Newton (ultrametric, no derivation) gives distance >= L+1.
- fakeexact.py (exact for L<=2: free columns = true problem): D=5,L=2 (c=0): M=3. D=6,L=2: >=5 (running) — matches families.md "5 at D=6"; my integer-box searches (B<=3, random B<=8) only found 4: box too small.
- cube_template.py: Gaussian-int roots scale 72: 8-pt polynomial template, chord/sqrt(R) -> 3.18e9 (before gcd); mechanism only.
- Thm D proof (generic freeness): F=f^*T_X(-m inf); generically generated subsheaf G rank>=r-1, df gives section of F(-1) => deg F >= 1 => -K.C >= r m + 1. -K.C = 2#{k: column nonzero} (cochar content 1 for 0/+-1 rows).
- fakeexact D=6,L=2: M=5 exact (245105 closures). D=5: 3. D=7,c=0 (L=3): rank-3 4-pt clusters realised numerically (d7b.py), no rank-4 realisation found (20 starts/code): evidence M(7,0)=4; rigorous 4?<=M<=A(7,4)=8.
- rank5.py: D=8 rank-5 spans in H8 (175 with >=9 words; 40 tested) none realised numerically.
- LP thr 4/7: rank4 ->5, rank3->4.
- fe_7_2 stopped at best 6 (lower bound only) for D=7,c=1
