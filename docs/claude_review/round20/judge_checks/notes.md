# Judge notes, round 20
## Read so far
- gauss: killed as route. survives: Thm 3.4/3.4' (A_k => uniform), Thm 4.1 (k=4 refuted; + killer corollary k=3), Thm 4.2 cond V_E(5), Thm 4.3(ii); 4.3(i) narrowed. Lemma K1 (killer): gcd rescaling invisible. (L_5^circ) open, capped, Vojta-strength, same standing as (L_8).
- arakelov: killed; nothing survives as lead. Thm A dictionary exact; Thm B = Fekete w/ crude (1+i); Prop B.2 needs lower.md witness 48 M log M at X>=2M; Prop 3.3(b)=row 8; Thm C = row 2. Totals 29/25.67.
- composite: cap for robust composite inputs survives (Thm M, P; A>max(4,1/g) lw / A>max(2,1/g) P3char). Killed: "no separating power", "only one leftover". Unclassified: (a) sharp prim LS w/o mult = max-scale cross-gap; (b) with mult = product-scale gap (heur. collapsing); (c) F3 smoothed sector counts w/ BA edges; (d) effective smooth Hecke sums w/o complete additivity.
## Judge idea J1 (to check): F3 WITHOUT multiplicativity is fake-satisfiable by a "swap extension":
 for each S-smooth d (moved to psi_d), pick distinct protected-free ideal alpha_d, N alpha_d within H=X^{1/2+e} of Nd and genuine angle within eta=X^{-1/2+e} of psi_d on same side of edges; give alpha_d the angle theta(d). Multiset of (norm,angle) changes only by norm shifts -> F3 error change <= ||f'|| * sum |Nd-Nalpha_d|/X <= Psi_S(2X) X^{-1/2+e} = o(1) for A>2.
 Edge-side issue for a whole class of BA edges: F_K0 nowhere dense, dim<1; Borel-Cantelli over psi_d for Haar phi with A large (K0 fixed before A). 
 Not for F4 at super-poly precision (swap error ~ X^{1/A-1/2}).
 Plan: exact/evidence check reusing killer c3b setup: genuine vs minimal graft vs swap extension.
## J1 result (j1_out.txt): swap extension restores F3 exactly to genuine level on killer's 30 sectors:
 X=1.5e6: genuine mean .141, minimal 2.298, swap .141, max|swap-genuine| .0014 (bound .144). X=2.5e6 similar (.0015, bound .089). 0 partner failures.
## J1 statement (to write): finite sector family, f in C^1_c, A>4, primes >= M^A: change <= C_f X^{2/A-1/2+o(1)}; H=X^{1/2+1/A}, eta=X^{-1/2+2/A}.
 Corollary sketch: Lipschitz-in-angle test fns (smooth Hecke sums, fixed k): eta=X^{-1/4}, H=X^{3/4} -> change <= C(1+|k|)||f||_{C^1} X^{1/A-1/4+o(1)}. So F4 capped at precision X^{-B}, B<1/4-1/A; super-poly F4 still open.
 Uniform-in-type class version: heuristic (Borel-Cantelli via dim of bounded-type set), not proved.
