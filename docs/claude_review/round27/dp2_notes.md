# dp2 working notes (parked derivations)

## Geometry at P (affine z=1, P=(0,0), w(P)=1, pi etale at P)
- w = 1 + w1 + w2 + w3 + ... local branch. D: w - T2(x,y) = 0, T2 = 1+w1+w2 (2-jet). Unique member of |-2K| with
  mult>=3 at P (6 conditions on h^0=7). Tangent cone: w3(x,y)=0 (binary cubic). pi(D) = {f - T2hat^2 = 0}: quartic with
  triple point at p, rational, param by slope of lines through p (over Q).
- Branch fields at P = splitting of w3. MR Thm 2.16: alpha_D = 4/(r m): real irrational root -> 2; rational root -> 4.
- eps = 4/3 (L-4/3E = L/9 + 4/9 Dt), tau=3/2 (Dt^2=-1), Vol: 2-t^2 on [0,4/3], 2(3-2t)^2 on [4/3,3/2]; beta = 17/18.
- Curves with alpha_C<2 through general P: need k=2 conjugate real branches mult m with 3m<=d<4m, or rational branch
  mult m, 1.5m<=d<2m. d=3,m=1: C in |-K+E| (E a line) nodal at P = tangent plane section of cubic surface Bl-down of E;
  Q-rational iff E Q-rational. Generic S: none over Q.

## Integer form
Q=(a:b:c:s), f(a,b,c)=s^2, H=|c|, A=max(|a|,|b|)= cX. u := s - T2hat(a,b,c) in (1/N)Z, |u| ~ c^2 |w3(dir)| X^3.
Liouville: |u|>=1/N => X >~ c^{-2/3}. Window gamma in (3/2,2) <-> c in (A^2, A^3).
m = gcd(a,b), (a,b)=m(a',b'), A' = H(l) = A/m. t = m/c on fiber E_l, X = A' t.
m=1 points: s^2 = G(c) monic quartic in c, coeffs f_j(a',b') ~ A'^j: Runge bound c <~ A'^3 = Liouville(D).

## Fibration Bl_{P,P'} S -> P^1 (lines through p), fibre F = -K-E_P-E_P'
N_k = -kK-(2k-1)E_P, fibre degree 1, M_k = pi_* N_k line bundle of degree mu_k = 2k-k^2 = -k(k-2).
k=2: mu=0, section = D (consistent: D|E_l = 3P + R_2(l), R_2(l) traces D).
F_k(l): unique section of L(k S_inf) on E_l with ord_P = 2k-1; extra zero R_k(l); locus Gamma_k.
Heuristic: primitive integral F_k has real-bounded norm lambda ~ A'^{k(k-2)}.
Bound: X^{2k-1} c^k A'^{k(k-2)} >~ 1 for Q not in Gamma_k. With A' = c^theta: X >~ c^{-e}, e=(k+theta k(k-2))/(2k-1).
e < 2/3 iff theta < 1/(3k). k=3: alpha >= 5/(3+3 theta) > 3/2 iff theta<1/9.  (TO VERIFY NUMERICALLY)

## Bounded-u points <-> integral points on twisted K3 double covers
V_u: zeta^2 = (s - T2hat)/u, double cover of S branched on D (D in |-2K| => K3 with D_4 point at the triple point).
Points with u fixed = integral points of V_u w.r.t. ramification curve, accumulating at the D_4 point.
Resolution boundary: Rt + E' + A_1+A_2+A_3 (A_i over the cubic branch field). Lang-Vojta predicts degeneracy.
