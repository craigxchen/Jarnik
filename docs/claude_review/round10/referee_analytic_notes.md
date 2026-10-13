# Referee notes (analytic.md)

## First pass observations
- Thm B uses lw.md Thm 6.1 with P = M^A, "A >= 4"; lw P_1 = max((40M/C_0)^4, (64M)^(2/eps)) -> at A = 4, M^4 < (40M/C_0)^4. Need A > 4. CHECK lw.md.
- Sec 6 "open regime P_+ < M^(2-delta)": but Thm B only established for A>4 (lw) / sharp.md Z >= 3M^A (log M)^2 with A ? -> regime M^2..M^4 also not covered unless another fake exists. CHECK construct.md / exact.md prime sizes.
- Def 3.3: robust inputs are PRIME-supported only. Composite-supported (Hecke large sieve over all ideals, sums over all ideals) need multiplicative graft; R2-type survives via Euler factor ratio (|prod_{p in S}(1 -+ p^{-1/2})^{-+2}| = 1 + O(M^{1-A/2})) but not proved in the write-up. Scope limitation.
- Lemma 5.3 only needs LCG_n(phi; n) (Norm alpha <= n), not n+1.
- LCG_n(phi; n) <=> phi = arg(beta), beta primitive, Norm beta <= n (converse trivial). So LCG is literally "phi_j is a lattice direction of norm <= p_j": the input IS the problem (author agrees).
- GPA: |sin(phi - arg w)| < (h Nw)^(-1/2) <=> dist(w, line R e^{i phi}) < h^(-1/2): a strip, not a sector. Equivalent to "every line through 0 not of low-height lattice direction passes within h^{-1/2} of a Gaussian prime".
- Prop C correct (trivial). Extends to distinct-tuple inclusion-exclusion.
- Lemma D.2 correct.
- (2.1) and Remark A.1 Parseval correct (4 r(N')).

## To test
- explicit fakes (construct.md, exact.md skew double core): are their phi_j lattice directions of norm <= p_j? (expect no -> LCG violated). Compute.
- Does any fake class satisfy LCG at the "monomial" level only (outside Lemma 0.1)?

## Second pass (after reading lw.md 6.1, sharp.md 5.1, outside.md Lemma 0.1, stability Cor 3.6, round9 README)
- CONFIRMED: lw.md 6.1 needs P >= (40M/C_0)^4 (and sqrt P' >= 1600 M^2) => A > 4 strictly; Thm B says "A >= 4" and summary "any A > max(2,1/gamma)" with W <= (2A+o(1)) M log M: wrong for 2 < A <= 4. For 2<A<=4 only sharp.md (W <= (bA+o(1))M log M, b>=33). Constants only.
- sharp.md: Z = max(3M^A, 21 r^2)(log M)^2, r = b(M-1): fakes exist down to primes ~ M^2 log^2 M (P3-char). 
- Round-9 R1 caveat: caps with p_j > X = 3M^A only bar arguments fixing A before c_0. Thm B's P3-char clause inherits this; not stated. Also gamma must be fixed before target constant (fake constant 2 max(4,1/gamma)). Published inputs have gamma >= ~0.3 absolute, so OK for those; mention.
- Prop A slogan "input that separates theta from phi must involve primes outside N" FALSE: fibre ARITHMETIC inputs (N^{2k} S_k(N) in Z; sum_{|z|^2=N} z^m in Z) separate; e.g. e=1: S needs prod_j 2 p_j^{2k} cos(4k phi_j) in Z. These are P3-all-type ("residues of all monomials at all moduli" = actual primes, outside.md 2.4), NOT capped. Prop A precise statement (identities valid for all real theta) is fine.
- Remark 3.4.2 + Sec 6 table: "sharp (X+K) large sieve" listed under pointwise prime separation "yes under GPA (Prop D.1)". Unjustified: sharp LS (const c(X+K), K -> inf) only yields spacing |sin(phi_p - theta_w)| >~ 1/max(p, N w) (WCG), weaker than PCG's (p N w)^{-1/2}. Lattice directions beta with p < N beta <= p/c^2 satisfy WCG but violate LCG. GPA/Prop D.1 do not apply. Status of sharp LS: open (neither capped nor shown fake-killing).
- GPA region is a STRIP of half-width h^{-1/2} around the line (not a sector): |sin(phi-arg w)| < (h N w)^{-1/2} <=> dist(w, R e^{i phi}) < h^{-1/2}.
- "GPA is not implied by GRH" should be "not known to follow from GRH".

## Computation results
- lcg_fakes_out: construct q=11 b=8 (r=100), skew dc q=3,7,11 (r=6,14,22, P=(40M)^4): LCG violated at ALL primes with Norm alpha <= p_j; ratios min .0004-.096, max .71-.94; genuine angles min |m|=1 exact; gaps evidence fine; fibre integrality dist to Z .095-.417.
- sharp_ls_check: 12546 composite beta, p<N beta<=4p (p=101,1009,10009): LCG violated all, WCG_1/2 holds vs all primes N<=1e5, PCG violated (witness <=1e5) 9360. Optimal LS constant: lattice graft == genuine; near-coincidence graft ~2x at K<=1/t.
- gpa_lemmaD2: h=29,53,197 all lattice directions h<N<=20h have witness (max ratio 77.1); Lemma D.2 lines: all 8977 split p<=2e5, max 113 steps.
- composite_graft: X=3000 all Gaussian integers; graft at 101,197 changes Lambda by up to 7.5% (K=30); Lemma-3.2-type bound useless (exceeds Lambda).
