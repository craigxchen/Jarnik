# Notes (round 14, C1 construction side) -- working file

## Plan of results
1. Break analysis (regime 1.4): random fibre (norm 1), torsion-avoiding (norm 2), L4K Thm 2 / Prop H systems
   (norm K+1, demand ~ theta(X)/r^h), coherent+root repair (deep classes uncertified), flexible tops (genuine core pins).
2. Pinning Prop P: per-unit genuine collisions of e_j at q | D_s(pi_j), s=1,2,3 (all q <= sqrt(2p_j) < X) pin
   phi_j = theta_j mod pi; all columns pinned => actual cluster. Consequence: C1 fakes must be non-genuine at small moduli.
   Check derivation: x^g in (0, pi/4): s=1,2,3 one-sided => x = x^g.
3. Linnik deep-class matching (unconditional).
4. Rank-uncertainty: Paley twin profile, G-coherent at q=3 mod 8 <= X, X >= c b M^2 C W: C + 2Hbar >= (M-1)/8.
5. Shallow certified labels, zeta norms in deep classes => X-uniform.

## Derivation scraps
- per-unit budget B(v,s) = w(v)/2 + [s odd] log2/2.
- e_j collision prob under uniform parity class: q=+-3 mod 8: 4/m_q; q=+-1 mod 8: 8/m_q if nu_j=0 else 0.
- Mertens in APs: sum_{q<=x, q=a mod 8} 1/q = (1/4) loglog x + B_a + O(1/log x).

## Status at end of round (summary of what was done)
- Write-up sections 0-6 complete. Theorems: B1, B2, B3 (break analysis), Lemma P0, Prop P, Cor P1, P2 (pinning),
  Lemma L, Prop G, Thm C (coherence), Thm R (rank-uncertainty, linear), Thm S + Cor S1 (shallow certification).
- Not done: Rig_pin finite test (recipe in open point 2); DeepBal; PairCompat.
- Key logic point: case (b) (no fakes) makes the class exact; growth then = L7 for actual data.
- Lemma P0 derivation: s=2 gives x in [-theta, theta]; theta<pi/4: s=3 forces x=theta; theta>pi/4: s=1 gives
  x=theta or x<=pi/2-theta, s=3 then forces x=pi/2-theta.
