# Referee notes (search.md)
R1. Lemma 1 split analogue: in split form, sigma(p) has vector v(p) for ALL sigma, so ray-generic places are Q-rational
    (split Thm 2(i): 16 rational places). So "rationality of places" is NOT twist-specific. Attribution of finiteness to
    anisotropy suspect. Check: (u,v) invariant under b -> i b (e2->-e2, e3->-i e3, e4->e4): u=e4/e2^2 inv, v=e3^2/e2^3 inv.
R2. Task Result line says "only finitely many of its [class] members over Q (proved)" -- search.md proves it for M only.
R3. Need: degree of m -> (u,v) param (birational?) for Thm B step 3.
R4. Need: other factors of disc of L/Y in xsym (other x-structured components?).
R5. Y = b^2/(-v): W = Y^2 - Y + u, so b rational in Y: splitting fields equal. G_m derivation checked by hand.
R6. Parity argument (odd ord disc => odd monodromy) is rigorous; RH g >= |N|+1 >= 3 ok.
R7. c1_sym.py: Im P_j const (j=1,2,3), Im n in {X^3,X}, L even, L(0)=0, no t: (a),(b) sufficiency OK.
R8. c1b_disc.py: disc(L/Y) = -[4 cubics = |p_r|^2 (sigma=0)] * F24 * [3 quartics]^2 * b^4 * (sum)^20. Quartics divide A0 and
    appear squared => L = A3 Y^2 (Y+e)^2, D=0: reducible conic (two lines, single rational branch each) -> degenerate. F24 only relevant. OK.
R9. c2_phi.py: F24 == F(e2,e3,e4) exactly; Phi quintic irreducible over Q; Phi(u(m),v(m))=0; pencil disc (k-4)^5(9k-4)^3;
    gcd test: m->(u,v) has degree 1 (birational over Q). Step 3 of Thm B OK.
R10. c3_genus.py: (my normalisation) disc = -2^47 (m-3)^12 (m-1)^16 (m+1)^2 (m+3) c^3 quart; lc = 256 (m-3)^4 (m-1)^4.
    Exponents differ from search.md (42,34) by even amounts (normalisation) -- parity same. 8 odd points distinct, lc != 0 there,
    m=inf order 4 (even). G derivation: Res_W = G. Genus >= 3 argument correct.
R11. c4_construct.py (own code, 60 dps): m=-1.01,-1.1,-1.333,-1.5,-2,-2.5,-2.95: 30 places, all 30 rays once, contact slopes 4.000 at w->inf and w->0 for all u_j.
R12. c5_jac.py (80 dps): residual of (3.1) 5e-47 (limited by double-root precision), J rank 28 (smallest sv 0.089/0.011 vs max 1e12/3e9),
     kernel = span(family tangent, gauge), independent. M is locally a full component (numerical, wide margin).
R13. Thm B step 2 re-proved: sum of the 8 base x-values is real and in Q(i) => centre rational; scale real => Y_r in Q. OK.
R14. c6_qsearch.py H=60: only m=-3 has a rational root (b=0: degenerate endpoint; author excluded m=+-1,+-3 silently). Galois group of
     G_m over Q is S4 at m=-2,-3/2,5,7/3 => generic group S4; geometric N normal in S4 with odd element => N=S4.
R15. SPLIT TWIN: complexified places +-(t/b_r) + i b_r. Split Q-structure: each place Q-rational => after rational centring,
     both t/b_r and i b_r in lambda*Q => b in (lambda/i)Q, Y_r = b_r^2/(-e2) in Q, (u,v) scale-invariant in Q => m in Q, G_m splits.
     Same Faltings finiteness. So Thm B is form-independent; NOT an anisotropy effect. Only anisotropy input remains parity.
R16. c7_scan.py: 399 grid points m in (-3,-1) (40 dps): all 30 places, 30 distinct rays; e>0, D>0, e<D, sigma<0 throughout.
     Reduction to proof: X-values +-t/b_r +- i b_r distinct iff b_r^2 distinct & nonzero (disc_Y != 0, u != 0 on (-3,-1));
     Re X0 != 0 => unramified, H vanishes at exactly one point (beta(X0)!=0); kappa_j != 0 automatic. So nondegeneracy on the
     whole interval follows (labels locally constant). c8_signpat.py: 278 patterns, 15 classes: confirmed.
R17. Thm B step 2 for complex points of M_C: partner places conj in X' => X'_0 in Q, b_r in lambda*Q. OK. Same in split form.
