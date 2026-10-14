# referee notes (nearstrata)
- Lemma 1.1: checked by hand (Hoeffding 2exp(-w^2/8d)=1/2r; rigidity inequality; Step 0 exact via theory.md L2.4 sharp remainder). OK.
- Prop 2.1, L2.2, L2.3: OK.
- Thm 2.4: OK. Thm 2.5: gadget hexpair 8 points checked by hand; case analysis checked (tau=2sigma, 2L+2t, s+L even/odd); (b) checked. OK.
- Lemma 2.5 (III) Newton interpolation degree s+L: OK.
- Prop 2.6, 2.7: OK as math; "Reading" is interpretive.
- L3.1, L3.2, L3.3: OK. Thm 3.4 OK. Cor 3.5 OK. Cor 3.6 OK.
- Cor 3.7: rho_min uses exponent m'-d_F <= o_F, so rho_min is a LOWER bound for true min radius; sigma>rho_min is necessary not sufficient. Wording "beats the smallest radius any near-strata statement can have" inaccurate unless o_F = m'-d_F.
- Prop 5.1: necessary condition o N2/P <= d is for MIN* at ALL large g (g->inf); Gamma_L depends on delta_y (uncontrolled). At fixed g=abar no such conclusion. Summary says "(MIN) forces" -> overstatement.
- MAJOR: §4.3 minimality. d_xy (hence eps_xy) determined by profile+primes; individual gaps |u_x-u_y|>=e^{-eps_xy/2} are pair data. MIN* discards them. MIN^pair weaker. Equal-weight Sylvester: eps_xy = abar all pairs => feasible set empty if e^{-abar/2}>2/(n-1).
