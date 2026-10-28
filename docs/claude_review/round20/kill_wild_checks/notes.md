# kill_wild notes
- KEY: (S_min) at n=2 is refuted by (1:g:g^2), g=A+i, A->inf: min-anchor sum = L=log(A^2+1) (anchor g), h(ell)=L. ratio 1 > 2/3.
  Proposal only took K->inf in Prop 6.6 and claimed F1,F2 do not refute S_min. Wrong for n=2 with A->inf.
- Also genuine short-arc rank-deficient 3-cluster family: z=((b-i)^2, b^2+1, (b+i)^2), b even, N=(b^2+1)^2, C~4,
  z1 z3 = z2^2 (rank 1), Plucker gcd = 2b (!), h(ell)=log(2(b^2-1)), min anchor sum = log(b^2+1). ratio->1 > 2/3.
  => also refutes "Plucker gcd = 1 on actual clusters" as a general pattern (only evidence claim, labelled).
- n=3: (1:g:g^2:g^3) min anchor 2L, h(ell)~2L+log3 -> ratio 1 = n/3 > (1/3-delta)*3.
- subsets family u_S = prod_{j in S} g_j, g_j=A+j+i: r=2 (n=3) min-anchor 1.5L vs h(ell)~L ?  check by code.
- Proposal's own data: N=2576450045 k=4 min-anchor ratio 1.0263 > n/3=1.
- Checked: Lemma 6.1 OK, Thm 6.2 OK (X_M matches overview 4.1), Lemma 6.3 OK, Prop 6.4 OK, Prop 6.6(d) algebra OK,
  Lemma 6.7 OK (N odd after normalisation), Prop 6.8 algebra OK.
- Cor 6.9 "ceiling of lens" only for subgroup-membership; cones heuristic -> overclaim in summary bullet (4) "even an ideal UI input".
- BIG: (S_min) false for EVERY n: u=(1, X_1+i,...,X_{n-1}+i, (X_1+i)(X_2+i)), X_x=T+x: dependent (u_n=u_1u_2), 
  h(ell)~2logT (Plucker: 1, X_x-X_y, ~T, ~T^2), min-anchor ~ n logT => ratio ~ n/2 > n/3. Rank n-1 (deficiency 1).
  Pair-products version: k generic + all k(k-1)/2 pair products: rank k ~ sqrt(2n), still ratio ~ n/2.
