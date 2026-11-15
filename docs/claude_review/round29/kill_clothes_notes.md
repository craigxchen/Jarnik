# kill_clothes notes
- Re-derived Thm 4.1 (i)-(iv), lattice/form statement (correct; form primitive since values N and q^2 coprime; converse needs gcd(n,q)=1 and h(Z[i])=1, argument: L_f Z[i]=gamma Z[i], |gamma|=1).
- Prop 4.2 correct. Thm 4.3(i) correct (re-derived; exact check 1396 classes).
- KILL: Prop 4.4 false. Inside a class mod q, l^t||q, residues mod l^(t+b) lie in exactly l^b classes (Lemma 3.2 lifting). So v_l(c/q) pigeonhole gives E(M,l^b). Sec.6 of proof.md needs only depth one: inert E(M,l)>=E(M,l+1); split E(M,p)>=M^2/(2p)-M/2>=c_pM^2/4-M/2 since c_p<=1/(p-1)<=2/p. => D(M) bound, Thm 7.5 with C=2/lam for EVERY odd q coprime to N. Hypothesis "no prime factor <=M" unnecessary. Lead (smooth conductor regime) collapses.
- exact: (A) lifting count 4528 cases 0 fail; (B) 1396 classes, confinement/pigeonhole/Thm4.3(i) 0 fail.
- twin numbers reproduced: lam=1 max 5 (N=2576450045,q=269), arc C=2 max 3; lam=1/2 max 8 (N=1105,q=3 and N=32045,q=7), arc C=4 max 5. The 8 at q=3 is trivial pigeonhole 32/4. 12 of 24 norms skipped (r>1100). smooth-q <=4 only at lam=1; at lam=1/2 smooth q=3 gives 8.
- Minor: Sec 2(b) "exponent 1/2-1/8=3/8 at k=4" wrong: CC k=4 exponent is 1/2-1/10=2/5 (conclusion <=4 still ok since 3/8<2/5). Also with p chosen larger the box/residue identity holds for any box; H^{o(1)} is then just divisor-bound strength.
- 2^{omega(q)+1} could be 2^{omega(q)} (c^2=1 mod q has 2^omega(q) roots, q odd).
