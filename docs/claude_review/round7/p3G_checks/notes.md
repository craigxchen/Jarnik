# p3G scratch notes (not the write-up)

Key idea: genuine angles theta_j = arg pi_j satisfy every (1.2) (planted point).
Tight relation (v,s): m_(v,s) = |D_s(A_v)| (D_s odd, divides odd lcm(1..X), coprime to N).
Then (1.2) <=> dist(<v,phi>, s pi/4 + pi Z) >= dist(<v,theta>, same)
<=> x = <v, phi - theta> in sigma [0, pi - 2d] + pi Z.
Circuit lemma: sum n_e sigma_e v_e in Pi, sum n_e a_e < pi => all x_e in pi Z.
Exact tests: sigma = sign(XY) (s=0), sign(Y^2-X^2) (s=1), -sign(XY) (s=2), sign(X^2-Y^2) (s=3).
sin^2 d = nu_s^2 D_s^2/|A|^2; a+a' < pi <=> S+S' > 1; 2a < pi <=> S > 1/2.
Collapse: all columns pinned mod pi/2 => fake = actual cluster with unit-adjusted primes.
Dichotomy: c with v(c)|_F = 0 => |D_s'(A'_v(c))| <= (nC/(4 nu)) e^{s_c/2}.
Rank: any member has rank V >= 2M/3 - 1 (from cut identity + Gram).
Heuristic: free columns need log p >~ c H, so W_fake >~ (2/3) c M H.
