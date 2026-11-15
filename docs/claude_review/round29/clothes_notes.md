# notes (clothes lane)
- started; skeleton written.
- KEY OBSERVATION (to verify): residue twin z = z1 mod q  <=>  z/z1 = prod u_p^{a_p-a1_p} = 1 in (Z[i]/q)^*, u_p = pi_p/conj(pi_p).
  So twin = box points in a COSET OF A KERNEL (exact subgroup) instead of a Bohr set (approximate subgroup).
  Twin Liouville: q | Im(P_v) => w(v) >= 2 log q; same pair/cut inequalities as arc.
- Form (a): m=N-Re(z conj z1)=|z-z1|^2/2, m(2N-m)=v^2, twin: q^2|m.
- dictionary.py: (a) m(2N-m)=v^2, 2N=g(a^2+b^2) verified 3792 pairs; (d) z==z1 mod q <=> z conj z1 == N mod q <=> q | a, 99900 pairs, 0 failures.
- twin.py / twin_smooth.py: twin counts on 19 N's: <=5 at lam=1 (arc(C=2) <=3), <=8 at lam=1/2 (arc(C=4)<=5). Smooth q never larger. EVIDENCE only.
- RECORD ALREADY HAS: overview Thm 3.10 (round8/stability Thm 3.1): M < 2W/(log P_+ - 4log^+C) [= the peeling bound (f)], Thm 3.11 near-extremal on M^{O(1)}-smooth norms; fakes share it.
  round4/hardness.md Prop 6.1: arc clusters are pairwise incongruent mod sigma with |sigma|>C N^{1/4}; "Lenstra at alpha=1/4 is the non-archimedean twin", "archimedean analogue of Lenstra at critical exponent is U itself".
  => my twin (d) = Gaussian Lenstra problem at EXACTLY alpha=1/4 with rational modulus, vertex divisors. New to record: the twin as a problem in its own right; ring class group form.
- Twin Thm1.1 transfer: identical (2.5) with C=2/lambda, but pigeonhole at primes l | q is lost. For q divisible by all l<=M (possible iff M <~ W/4) the method gives nothing.
- Peeling+Plotkin (elementary, no CC): m_a <= 1 + W_1/(e log p - 4 log C), M <= (e+1)(...). Rediscovers Thm 3.10 up to constant.
- Literature: CC (search snippet via CG): arc sqrt2 R^{1/2-1/(4[k/2]+2)} has <= k points. Lenstra / CHN O((alpha-1/4)^{-3/2}) for alpha>1/4 (search snippet). Blomer-Granville Duke 2006 averaged r_f(n) uniform D=o(x).
