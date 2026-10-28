# Working notes (round 21, fakes in the window)

## Key facts collected
- Gap set G(p) subset T^r: dist(<v,phi>, pi/2 Z) >= arcsin(e^{-w/2}), dist(<v,phi>, pi/4+pi/2 Z) >= arcsin(e^{-w/2}/sqrt2), all v != 0.
- Known fakes: phi Haar on realization torus; non-character v handled by union bound Pi-1, Pi = prod(1+2/(sqrt p_j -1)) (outside.md A.1(ii), needs p_j >= 21 r^2; lw.md 6.1(b)).
- Capped input row 3 (stability Thm 3.3(3)): omega(N') > M/4 (C<=1). So in window Lambda = sum p^{-1/2} >= M^{delta/2}/4.
- Equal-weight window [P, P e^{1/(4r)}] is NOT the essential obstruction (groups of primes per cut fix weights); the union bound is.

## Heuristic picture
- Supercritical range: #{v: w(v)<=X} > e^{X/2} for X < X_0 ~ sum 2 log p/sqrt p.
- Transplants: phi_j = arg beta_j, beta_j Gaussian integers norm n_j <= p_j, pairwise coprime norms => gap-admissible; realization <=> actual cluster on N' = prod n_j^{e_j} with constant C(N'/N)^{1/4}.
- Hybrid genuine background G + Haar free F: free primes need log p_f >~ 2 log Pi_G and sum_F p^{-1/2} <~ 1 with |F| >~ M-1.
- m=0 (rank L = r) rational-angle fakes: discrete Haar over characters of Z^r/L; same union bound.

## Plan
- Thm: Haar scarcity on realization torus (product over coordinates independent mod L_Q).
- MC: decay law of mu(G(p)) vs Lambda, Pi for small prime sets. exp(-c Lambda) vs exp(-c Pi)?

## MC results (mc1_out.jsonl, mc2_out.jsonl) -- upper estimates of mu(G(p)) (truncated at T)
- {13}: .302; {13,17}: .0361; {13,17,29}: .00315; {13,..,37}: .00024; {13..41}: 0/1e5
- 6 primes 97..149 (Pi 3.06): .0805; 9 primes (Pi 4.86): .0102; 8 primes 193..277 (Pi 2.85): .1226; 12 primes (Pi 4.52): .0213
- -ln mu/(Pi-1): 1.1-1.2 for larger primes, 1.55-1.8 for smallest primes. Poisson law mu ~ exp(-c(Pi-1)).
## next_angle (p=101, genuine bg): measure .773,.629,.442,.305,.075,.0072 for Pi_bg=1,1.77,2.90,4.22,8.06,17.4
  -> ln(meas/.773) ~ -0.29(Pi_bg-1). Haar bg at Pi=17.4: .025 and transplant dirs not special.
## cascade: G_1(5) = {arctan 1/2, arctan 2} (n=1 windows tile: arctan(1/2)+arctan(1/3)=pi/4).
  bg{5}->13 pinned; bg{5,13}->17 pinned; bg{5,13,17}->29 pinned (float). bg{5}->17 not pinned (contains arg(2+3i)).

## Theorem C (initial-segment rigidity) -- PROVED (Minkowski + tight lines)
S contains all split primes < p, phi_s genuine (any conj) for s in S; gaps for v=(v',+-1) => phi_p = +-arg pi_p mod pi/2.
Step1: gaps => |sin(theta - arg z)| >= (p Nz)^{-1/2} for all primitive z, Nz < p (odd: axis gap; (1+i)z0: diag gap).
Step2: Minkowski box |Re(z e^{-i theta})|<=(1+eps)sqrt p, |Im|<=(1-eps/2)/sqrt p, area>4 => theta within arcsin(1/p) of +-arg pi_p mod pi/2.
Step3: z_+- on lines Im(z conj(pi~)) = +-1 with R = Re(z conj pi~) in (0,p): kills eta in (0, 2 arctan(1/R)) resp (-2arctan(1/R),0); 2arctan(1/p) > arcsin(1/p).
Cert (cascade_cert.py, mpmath 50 digits, T=20): 13,17,29,37,41,53 all branches certified (1+2+4+8+16+32 runs).
## seq builds: 8 primes: 4/4 verified to w<=34 (360k monomials). 10 primes T=27: 1/3 verified at T2=31. 12 primes T=26: 0/3 complete (died at steps 9-11); T=22: complete but fail at 25.
