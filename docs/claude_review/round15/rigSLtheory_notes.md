# notes (rigSLtheory, round 15)

## Derivations parked
- J0 normalised: a_s = m_s^2/(nu_s^-2 N). Box: cos y in [2a_2-1, 1-2a_0], sin y in [2a_3-1, 1-2a_1].
  J0 <=> (OPP) a0+a2<=1, a1+a3<=1; (ADJ) max(0,2a0-1,2a2-1)^2+max(0,2a1-1,2a3-1)^2<=1; (MAX) max(|2a2-1|,|1-2a0|)^2+max(|2a3-1|,|1-2a1|)^2>=1.
- Genuine identities: D0^2+D2^2 = N, D1^2+D3^2 = 2N; D0 D1 D2 D3 = Im(Z^4)/4; SL: prod m_s <= (N^2/4)|sin 4x|.
  m0 m2 <= (N/2)|sin 2x|, m1 m3 <= N |cos 2x|.
- One-modulus extra (v,s'): PU: q0^2 <= U_{s'}^2/sin^2 theta; OPP: q0^2 <= U_{s'}^2 (1 + cot^2 theta (1 - U_{s'+2}^-2)).
  U_{s'+2}=1 => OPP <=> q0 <= U_{s'} <=> fits at genuine phi.
- Transport (quadratic): m'_s^2 + m'_{s+2}^2 <= D_s^2 + D_{s+2}^2.

## Scale
- stability Thm 3.3(3): omega > W'/(4W'/M + 8 log^+ C) for every cluster M>=4 (no near-extremality). Members: via Prop 7.1
  inputs + same-level clause; sub-models of members are members (minus (G)), residues shifted by a common factor.
- K=2: every q with n_q < 2r+1 forced (pigeonhole on S_1). r > M/4 => all q <= 2M (n<=M/2) forced.

## Computations (exact)
- classify.py: 83(K4), 163(K5) first4: all non-class PU survivors have OPP-failing extras, all with U_{s'+2}=1; no A1.
  107: 108 all-C (SL at genuine), 36 with OPP (U=1). rand4@31: 8 ANG-ok each with 4 A1; bal4@41: 8 with 2 A1; last4@31 all C.
- slab_iv.py (mpmath iv, rigorous): rand4@31 all 8 INFEASIBLE (21 nodes), bal4@41 all 8 INFEASIBLE (49 nodes).
- core_iv (fixed bnb): rand4@31 core = {e2 (U=(2,1,1,1)), e4 (U=(1,1,2,1)), v=(0,2,0,1) A1 at unit 0 (U2=802)}.
  bal4@41 core = {e3, e1, e1+e4 (all odd-saturated), v=(1,0,-1,2) A1 at unit 2}. Mechanism: odd-saturated v pinned.
- Lemma B1: OPP equality => sin 2x_v = sin 2x_0 (odd pair) / cos 2x_v = cos 2x_0 (even pair). Label-independent set
  P_0 = {v: U_1=U_3=1 wrt F'}: c' colliding at odd unit => A0 => infeasible; else pinned. => Phi finite.
- all-moduli-charged one-modulus scan first4 K4: SL-rigid all q0<=113 (109 PU-flex SL-rigid), 127 SL-flex (512 at phi0).
- affordable.py: first4 K5 F_5: PU/ANG/C affordable extras per q0: e.g. 47: 65/41/32; ANG-C small. all-moduli: 127+: <=3.
- Theorem B (index [Z^r:<T>]=1) at every forced q0>=7 for first4/rand4/bal4/last4, K=4 and K=5 (tg_*.txt).
- one-mod SL exact (angle-free): rand4 K4 rigid <=47, flex 59 (lifted witnesses); bal4 rigid <=41, flex 43; last4 flex 31.
- strictness: K4 one-mod class survivors at 79,107,113 include local twists/unit shifts surviving SL at phi0 (e.g. 113:
  96 class SL survivors, 24 genuine-strict, 48 strict only with local t, 24.. non-strict). K5: 107,113 all genuine-strict.
- Theorem P (P3-G collapse): member strict-genuine (eps,t) on prime factors of Re(Z'_v^2), v in generating triple family
  => phi = phi0' mod pi => points are Gaussian integers => actual cluster. TC all-moduli holds first r<=64 (K=3), loose sets.
