# referee notes (rigtheory)
- Thm 3.4 proof re-derived: OK. Lemma 3.3 OK. Cor 3.5(3) OK. Prop 3.7 OK. Thm 2.5 OK (q_core not monotone, only cores).
- ISSUE A (major): associate twist pi_j -> i pi_j: labels L -> L + (m/2) t (all levels), phi -> phi + (pi/2) t,
  theta0 -> theta0 - (pi/2)<e,t>, k_x -> k_x + 2<a_x,t>. Member -> member, same W. Strictness (u = 1 mod 4, eps only) not invariant.
  => RigS as stated implies class empty in regime (since twist of any member is non-strict) => Thm 4.1 hypothesis contains conclusion.
  Fix: strict up to (eps,t) group (4^r label systems).  Need computational demo.
- minor: log i = m/4 at every level not compatible with reduction for inert q, even a (log i = 3m/4 then). Author code uses independent generator per level; harmless.
- minor: |sin((delta - arg eta)/2)| >= 1/sqrt2 slightly false (>= 1/sqrt2 - O(Delta)).
- RigS's K dangling.
- ILP infeasibility = HiGHS float, not exact certificate (task standard).
- Independent: genuine identity 198,560 checks OK; forced sets K=2..5 = 7/17/40/65 (certified to 4A_K+5, B=8 Sidon {1,2,4,8});
  one-modulus K=2..5 cores 7/41/79/157 match; 167/199/211 at K=5 96/40/832 F1 match; witnesses (22 files) verified, 0 viol.
- twist.py: twisted genuine system non-strict at 14/17, 37/40, 61/65 forced moduli (K=3,4,5). Lemma: n odd & (p_1/q)=1 => non-strict.
- cover.py: covered moduli (balls/double balls) K=5 up to 131 with gaps (73,83,97,103 not covered); 137-157 NOT covered, yet rigid.
  K=3 core reaches n=10>|S_1|=9. So 'anticode threshold explains' is coincidence/evidence, not theorem (onset at 163 is F1-only).
- r>=M^gamma not part of regime (1.4); if r=M^{o(1)}, bounded-K core <= O(r^{K/2}) < M^{1-eta}: Rig(K) can't give RigS at bounded K.
- joint K=4/K=5 data conditional on fixing genuine labels at 3,7 (and top >163 at K=5): 'joint core grows' not unconditional.
  Unconditional: q_core^{J,eps}(3)=37 (if ILP right), 37 <= q_core^{J,eps}(4) <= 61.
