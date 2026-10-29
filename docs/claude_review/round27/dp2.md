# Round 27, lane dp2: alpha^ess > tau at a rational point with tau < dim, on a degree-2 del Pezzo surface

Status: complete for this round (author-checked only; **unrefereed**). Outcome type: obstruction + restricted theorem +
conditional reformulation. No unconditional `alpha^ess > tau` for generic sequences was found.

## 0. Summary of results

**Setting.** `S : w^2 = f(x,y,z)` (explicit, §1), `P = (0:0:1:1)`. Certified: `S` smooth; `P` off the 56 lines and
off the ramification curve; `Pic(S) = ZK` (Frobenius at 3 has eigenvalue 1 of multiplicity one on `NS`).

**Part (1), verified exactly.** The member `D in |-2K|` with a triple point at `P` exists, is unique, irreducible and
rational (Lemma 2.1); `eps_P(-K) = 4/3`, `tau_P(-K) = 3/2 < 2`, `beta_P(-K) = 17/18`; the tangent field of `D` at `P` is
a cubic field with one real place, and `alpha_{P,D} = 4/2 = 2` (McKinnon-Roth Thm 2.16). Over `Q` no curve of degree
`<= 6` through `P` has `alpha < 2` (census §2.4: the cuspidal fibre `E_0` and `D` have `alpha = 2`; the four members of
`|-3K|` with square tangent cone at `P` are conjugate over a quartic field). Odd degrees are impossible.

**Part (2).** No unconditional `alpha^ess > tau` result at a rational point of a dP2/dP1 (or any variety with
`tau < dim`) was found; Broustet confirms `eps = 4/3`; Castaneda's thesis covers split surfaces of Picard rank `<= 5`.

**Part (3), main findings.**
1. *Elliptic pencil = two-point Liouville (Prop A1).* The pencil through `P` is the rational elliptic surface
   `Bl_{P,P'} S` with `MW(Q(t)) = Z sigma'`; the multiples `k sigma'` trace rational curves `Gamma_k` of class
   `(k^2-k)L - (k^2-1)E_P - (k^2-2k)E_{P'}` (`Gamma_2 = D`). Fibre-wise Liouville with the natural integral
   normalisation (height cost `A'^{-mu}`, `mu = (2k^2+1-a^2-b^2)/2`) is *identical* to global Liouville with forced
   divisibility at `P` and `P'` for `Gamma_k`.
2. *Theorem A (complete two-point section cap).* Every Liouville inequality from any `Q`-section of any multiple of
   `-K`, with forced divisibility at `P` and `P'`, follows from those of the `Gamma_k`; the exact exponent function is
   `x_sec(nu_+,nu_-)` (§4.3). It beats `3/2` **only** for approximants that are `p`-adically close to `P`
   (`nu_+ > 0`, trivial) or to `P'` on a scale `m_- > H^{2/9}` (via `Gamma_3`), and never for `gcd(a,b) = 1`. This is a
   rigorous restricted theorem (`alpha > 3/2` for sequences in `{m_- >= H^{2/9+eps}}` off `Gamma_3`), but it is
   a section/Liouville statement (T4 class).
3. *Barrier for fibre-wise methods (Prop A3).* On the smooth dP2 `S_0 : w^2 = z^4 + z f_3 + x^4` at a point off the
   lines, with `tau = 3/2`, Liouville is attained: a cuspidal cubic carries points with `X^{3/2}H -> 1`, and they are
   independent of `P' - P` on their fibres. Hence no inequality depending only on fibre data
   `(A', c, m_+-)` (elliptic logarithms, uniform Siegel, Mumford gaps, Runge/Masser exponents) can beat `3/2`; the
   needed statement (UE) is a Runge-exponent improvement `log c <= (3 - eta) log A'` specific to the pencil of `S`.
4. *Square condition.* `D` contributes finitely many window points (Roth, Prop B1). Sieves cannot produce an empty
   window; descent via bitangents does not see proximity; the only non-averaged content is the integral-point problem
   (5.3).
5. *Theorem C (conditional, new reformulation).* Near-tight Liouville points (`|s - T(a,b,c)|` bounded) are
   `Sigma`-integral points of the log general type pair `(Bl_P S, Dt + E)` (`K + Dt + E = pi^*(-K) - E` big), equivalently
   integral points on quadratic twists of a K3 double plane branched along `D`. Lang-Vojta for this one pair implies
   that `X^{3/2}H <= K` holds only on a proper Zariski-closed set, for every `K`. Known Subspace-type criteria fail:
   for the relevant component `beta(L, E) <= 17/18 < 1` (it *is* McKinnon-Roth's `beta`), and on the K3 model
   `beta(H, Lt_+) <= 2/3`.
6. *Transfer to `X` at `e`.* The pencil/section part transfers and gives only the record's cap. The integral-point
   reformulation does **not** transfer: on `X` the near-tight points are integral for the ten diagonals `Z_ab`, and
   `K + sum Zt_ab + E = pi^*(2(-K) - B_{1,4}) - 6E` with `tau_e(2(-K) - B_{1,4}) = 5` exactly for `k <= 3`
   (`toric_order.py`), so the pair is not of log general type (numerically); balance must enter on `X`.

**What is new mathematics here, and what is not.** New: the identification of the dP2 pencil with a rational elliptic
surface whose `Q`-sections `Gamma_k` determine the complete two-point section cap (Theorem A); the calibration surface
`S_0` showing that Liouville is attained at some dP2 points off the lines and that fibre-local bounds are impossible
(Prop A3); the reduction of the top of the window to Lang-Vojta for one log general type pair (Theorem C), with the
explicit failure of Subspace criteria; the computation that this reduction does not transfer to `X`. Not new / mapped:
everything that is a Liouville inequality (Theorem A is one); square sieves (averaged); Roth at a rational target.

## 1. The concrete surface S and point P; genericity

**The surface.** `S : w^2 = f(x,y,z)` in `P(1,1,1,2)` over `Q`, with

    f = x^3 y + 2x^2 y z - 2x y^3 - 2x y z^2 + 2x z^3 + 2y^4 + 2y^3 z + y^2 z^2 - y z^3 + z^4,

`P = (0:0:1:1)`, `p = pi(P) = (0:0:1)`, `P' = sigma(P) = (0:0:1:-1)` (`sigma` = Geiser involution `w -> -w`).
`-K = pi^*O(1)`, so `H_{-K}(a:b:c:s) = max(|a|,|b|,|c|)` for a primitive integral representative.
(Found by `dp2_checks/search_f.py`, seed 1; a first trial surface was singular at `(-1:0:1)` and was discarded.)

**Certified properties** (all exact; scripts in §9):

| property | certificate | script |
|---|---|---|
| `f` smooth over `Q`, good reduction at 3 and 5 | Groebner bases of `(f_x,f_y,f_z)` over `Q`, `F_3`, `F_5` contain pure powers of `x`, `y`, `z` among their leading monomials (`x^3,y^3,z^7` over `Q` and `F_3`; `x^3,y^5,z^7` over `F_5`), so the singular locus is empty | `smooth_bitan.py` |
| `P` not on the ramification curve | `f(0,0,1) = 1 != 0`, so `pi` is etale at `P` and `(x,y)` are local coordinates | -- |
| `P` on none of the 56 lines | no bitangent of `B={f=0}` passes through `p`: the two bitangency conditions along the pencil of lines through `p` have resultant `-80353/1024 != 0`, and the vertical line fails them (`19/8, 119/64`) | `smooth_bitan.py` |
| `rank Pic(S)^{Gal} = 1`, i.e. `Pic(S) = Z K` | point counts `#S(F_{3^n})`, `n=1..7`, give traces `0,2,3,2,0,-1,0` of Frobenius on `NS(S_{\bar F_3})`; the characteristic polynomial is `(X-1)(X+1)(X^6-X^3+1)` (`Phi_1 Phi_2 Phi_18`), so eigenvalue 1 has multiplicity one. The Galois-invariant part of `Pic(S_{\bar Q}) = Pic(S_{\bar F_3})` lies in the Frobenius-invariant part, hence has rank 1, and `K` is primitive. Consistency check at `n=8` in §9. | `ffcount.py`, `frob_pic.py` |
| tangent cone of `D` at `P` | `F_3 = 2w_3 = (8x^3 + 20x^2y - 10xy^2 + 19y^3)/8`, irreducible over `Q`, discriminant (in the slope) `-1667/4 < 0`: one real root, two complex | `geom.py` |
| `D` irreducible | `Res(F_3, F_4') = -80353/1024 != 0` (Lemma 2.1) | `geom.py` |
| lines through `p` tangent to `B` | discriminant of the pencil = `t^2 * (irreducible degree-10 polynomial)`: the line `y = 0` is a flex tangent (`f(x,0,z) = z^3(2x+z)`), the other 10 tangent lines are conjugate over a degree-10 field | `smooth_bitan.py` |

**Why this is general enough.** What the arguments below use is: `S` smooth; `P` off the 56 lines and off the
ramification curve; `D` irreducible with an ordinary triple point whose tangent field is a cubic field with a real
embedding; every fibre of the elliptic pencil through `P` irreducible (no bitangent through `p`); and
`Pic(S) = ZK`, which rules out `Q`-lines and `Q`-conics, so the nodal tangent-plane cubics (the only curves through a
general point with `alpha = 3/2`, see §2.4) are not defined over `Q`. One special feature is visible and harmless: the
fibre over the flex line `y=0` is a cuspidal (Kodaira type II) rational curve through `P`; it is irreducible, so the
Mordell-Weil lattice of the pencil is still `E_8` over `\bar Q`.

**Notation used throughout.** For a rational point `Q = (a:b:c:s) != P` near `P` on the `P`-sheet, with `(a,b,c)`
primitive and `c > 0`: `H = H(Q) = c`, `X = X(Q) = max(|a|,|b|)/c` (comparable to `dist(P,Q)`),
`m = gcd(a,b)`, `(a,b) = m(a',b')`, `A' = max(|a'|,|b'|)` (the height of the line `l(Q) = pQ`), so `X = A' m / c`.
For odd `p | m` one has `s^2 = f(a,b,c) == c^4 (mod p)`, so `s == +c^2` (`Q` is `p`-adically near `P`) or `s == -c^2`
(near `P'`); `m_+`, `m_-` are the corresponding parts of `m` (the prime 2 is absorbed in constants).
Exponents are written `X = c^{-x}`, `m_+- = c^{nu_+-}`, `A' = c^theta`, so `theta = 1 - x - nu_+ - nu_-`, and a
sequence with exponent `x` has approximation constant `alpha = 1/x` for `-K`.

## 2. Part (1): the triple-point member D of |-2K|, eps_P, tau_P, alpha along D; other low-degree rational curves

**2.1 The curve `D`.** `H^0(S,-2K)` is spanned by the six quadrics in `x,y,z` and `s` (`h^0 = 7`). Near `P` put
`z = 1`; since `pi` is etale at `P`, `s = w(x,y) = sqrt(f(x,y,1)) = 1 + w_1 + w_2 + w_3 + ...` (`w_i` homogeneous of
degree `i`). A section `lambda s + q(x,y,z)` has order `>= 3` at `P` iff `q = -lambda(1 + w_1 + w_2)` (six linear
conditions), so there is exactly one member

    D : s = T(x,y,z),   T = z^2 + z(x - y/2) + (-x^2/2 - xy/2 + 3y^2/8)   (the homogenised 2-jet of w),

and its local equation is `w - T = w_3 + w_4 + ...` with
`w_3 = F_3/2`, `F_3 = x^3 + (5/2)x^2y - (5/4)xy^2 + (19/8)y^3`. So `mult_P D = 3` exactly and the tangent cone
`F_3 = 0` consists of three distinct lines (`disc != 0`). `pi(D) = {f - T^2 = 0} = {z F_3(x,y) + F_4'(x,y) = 0}` is a
plane quartic with an ordinary triple point at `p`, and `pi|_D` is birational (`s` is a polynomial in `x,y,z` on `D`).

*Lemma 2.1.* `D` is irreducible and rational.
*Proof.* A plane quartic with a point of multiplicity 3 is either irreducible or contains a line through that point:
the components through `p` have multiplicities summing to 3, an irreducible conic is smooth and an irreducible cubic
has multiplicity `<= 2`, so without a line through `p` the only splittings (conic+conic, cubic+line not through `p`)
have multiplicity `<= 2`. A line `(x,y) = m(a',b')` is a component iff `zF_3(a',b') + mF_4'(a',b')` vanishes
identically, iff `F_3(a',b') = F_4'(a',b') = 0`; `Res(F_3,F_4') = -80353/1024 != 0` and the point at infinity fails
too. An irreducible plane quartic with a triple point has geometric genus `3 - 3 = 0`. `QED`

The same resultant shows that no line through `p` is a bitangent of `B`: on a bitangent `l` through `p`, `f|_l` is a
square and its square root through `P` is `T|_l`, so `l` would be a component of `pi(D)`.

**2.2 Seshadri constant, `tau`, volume, `beta`.** On `Y = Bl_P S` let `L = pi^*(-K)`, `E` the exceptional curve, and
`Dt = 2L - 3E` the strict transform of `D` (`L^2 = 2`, `E^2 = -1`, `Dt^2 = 8 - 9 = -1`, `Dt` irreducible by 2.1).
* `eps_P(-K) = 4/3`: `L - (4/3)E = (1/9)L + (4/9)Dt` is nef (it is `>= 0` on `E`, on `Dt` it is `4/9 - 4/9 = 0`, and on
  any other irreducible curve both terms are `>= 0`), and `Dt` shows it is not nef beyond `4/3`. This agrees with
  Broustet's value `4/3` for `-K` at a general point of a degree-2 del Pezzo surface (Broustet 2006, see §3).
* `tau_P(-K) = 3/2 < 2`: `L - (3/2)E = Dt/2` is effective; since `Dt^2 < 0` and `Dt` is irreducible,
  `h^0(nDt) = 1` for all `n`, so `Vol(L - (3/2)E) = 0`; if `L - tE` were effective for some `t > 3/2`, then
  `L - (3/2)E` would be a positive combination of the big class `L` and an effective class, hence big. Contradiction.
* Zariski decomposition: `Vol(L - tE) = 2 - t^2` on `[0,4/3]`; on `[4/3,3/2]` the positive part is
  `(9-6t)L - (12-8t)E` (negative part `(3t-4)Dt`), of square `2(3-2t)^2`.
* `beta_P(-K) = (1/2) Int_0^{3/2} Vol = (1/2)(152/81 + 1/81) = 17/18`. So the Faltings-Wuestholz/McKinnon-Roth
  Theorem 5.1 bound `alpha^ess >= beta` gives only `17/18`, below Liouville's `3/2` (the usual pattern at a rational
  target, papers_check README §1.2).

**2.3 `alpha` along `D`.** Normalise `D` by the slope `t` of lines through `p` (each line through `p` meets `pi(D)` in
`p` (multiplicity 3) and one further point, so `t` is a rational parameter over `Q`). The three points of `P^1` over
`P` are the roots of `F_3(1,t)`, i.e. of `19t^3 - 10t^2 + 20t + 8` (irreducible; discriminant `-1667 * 2^10`, `1667`
prime, so the Galois closure is an `S_3`-extension with quadratic resolvent `Q(sqrt(-1667))`): one real cubic
irrationality and two complex ones. McKinnon-Roth Thm 2.16
(`alpha_{x,C} = min_q d/(r_q m_q)`, `d = -K.D = 4`, all branches smooth so `m_q = 1`; `r = 2` for the real branch,
`r = 0` for the complex ones) gives **`alpha_{P,D}(-K) = 4/2 = 2`**: the points of `D` approximating `P` are the
Dirichlet approximants of the real cubic slope, `H(Q) ~ H(t)^4`, `X ~ H(t)^{-2} ~ H^{-1/2}`. In the notation of §1:
these points have `A' = H(t) = c^{1/4}`, `u := s - T(a,b,c) = 0`.

**2.4 Other rational curves through `P` of low degree.** Because `Pic(S) = ZK`, every curve on `S` defined over `Q`
has class `-nK`, hence even degree `2n`.

*Which curves can have `alpha < 2` at all.* Let `C != D` be an irreducible `Q`-curve through `P` of degree `d`, with
a Galois orbit of `k` branches at `P` of multiplicity `m` realising the minimum in Thm 2.16. `Dt.Ct >= 0` gives
`2d >= 3 mult_P C >= 3km`. If the branch is rational (`r = 1`): `alpha = d/m < 2` forces `d < 2m <= 2 mult_P C <= 4d/3`,
so `m >= 2` and `3m/2 <= d < 2m`. If it is real irrational (`r = 2`, `k >= 2`): `alpha = d/(2m) < 2` forces
`3km/2 <= d < 4m`, so `k = 2`, `3m <= d < 4m`. (`r = 0` gives `alpha = infinity`.) So `alpha < 2` needs a very
singular point at `P`; the threshold `alpha = 3/2` needs `Dt.Ct = 0`.

*Expected dimensions (over `\bar Q`; HEURISTIC: I did not prove that the conditions are independent, and a fixed
rational `P` is never "very general").* Rational curves of degree `d` move in a `(d-1)`-dimensional family, and a branch of
multiplicity `m` through a fixed point costs `2m - 1` conditions. With two conjugate branches (`r = 2`) the cost is
`4m - 2 <= d - 1 < 4m - 1`, so `d = 4m - 1` and the curves are rigid: `d = 3` (`alpha = 3/2`: the member of
`|-K + E|` nodal at `P`, i.e. the tangent-plane section of the cubic surface obtained by contracting the line `E`),
`d = 7` (`alpha = 7/4`), ..., `alpha = 2 - 1/(2m)`. With a rational branch (`r = 1`) the cost `2m - 1 > d - 1`
always: none. So (heuristically) for general `P` the only curves with `alpha < 7/4` are the 56 nodal cubics (one per line), and
these are defined over `Q` only when the line is. On our `S` there is no `Q`-line.

*Exact census for our `P`, degrees 2, 4, 6:*

| class | curves through `P` defined over `Q` | `alpha` at `P` |
|---|---|---|
| `-K` (deg 2) | the fibres `E_l = pi^{-1}(l)`, `l` a line through `p`. They are smooth of genus 1 except over the 12 tangent lines from `p`; over `Q` only the flex line `y = 0`, whose fibre `E_0` is a cuspidal rational curve, smooth at `P` (the cusp lies over `(1:0:0)`) | `E_0`: `2/1 = 2`. Smooth `E_l`: `infinity` (Siegel; a rational point of an elliptic curve) |
| `-2K` (deg 4) | `D` (triple point); the rational members with a node at `P` (expected: a 1-parameter family over `\bar Q`); rational members smooth at `P` (expected: 2-parameter family) | `D`: 2; nodal: `4/2 = 2` (real-quadratic tangents) or 4; smooth at `P`: 4 |
| `-3K` (deg 6) | `alpha < 2` would need a 4-fold point with tangent cone `kappa q^2` (or `l^4`). The members of the net `|-3K - 4P|` (spanned by `xD`, `yD` and one more section) with square tangent cone are exactly 4, conjugate over a quartic field (the eliminant is an irreducible octic, even in `q_2`, in every chart); none is defined over `Q` | all `>= 2` |

So on our `S` **no curve of degree `<= 6` defined over `Q` passes through `P` with `alpha < 2`**, and `alpha = 2` is
attained on `D` and on `E_0`. The Liouville exponent `3/2` is not attained by any curve of degree `<= 6`; by the
(heuristic) dimension count it is attained by no curve at all for general `P` (nodal cubics excluded by
`Pic(S) = ZK`).
The window is therefore exactly the one the coordinator described:

    3/2 = tau_P(-K) <= alpha^ess_P(-K)   (Liouville with base locus, MR-L Thm 3.3, base locus D),
    alpha^ess_P(-K) >= 2 = dim S         (predicted by Vojta: erratum Prop 0.1, LMS Prop 1.4),

and `alpha_{P}(-K) >= 3/2` with no curve known below 2 (`deg6.py`, `geom.py`, `smooth_bitan.py`).

## 3. Part (2): literature

Searches (WebSearch, October 2026; arxiv.org and uwspace hosts are not reachable from this session, so only search
snippets were read; nothing below is a full-text verification):

* **No source computes `alpha` or `alpha^ess` at a rational point of a degree-2 or degree-1 del Pezzo surface.**
  Queries combining "approximation constant", "essential", "degree 2 / degree two del Pezzo", "McKinnon",
  "curve of best approximation" returned only: McKinnon-Roth (arXiv 1306.2976, 1306.2977), McKinnon-Satriano
  (2004.05212), Lehmann-McKinnon-Satriano (2403.02480), Huang (1809.02001 and the dP6 paper 1610.02542), McKinnon's
  cubic hypersurfaces paper (2310.01640: `alpha` computed when a rational line exists; Liouville-type mechanism, the
  infimum `alpha`, not `alpha^ess`), and weak-weak-approximation papers on dP2 (0809.1251, 2303.09299), which are about
  density, not approximation constants. A search snippet of Huang 1809.02001 says earlier computations "mostly
  consider rational surfaces, often (weak) del Pezzo surfaces of degree >= 3".
* **Castaneda, PhD thesis, Waterloo 2019** ("Rational approximations on smooth rational surfaces", advisor McKinnon;
  UWSpace record 10012/14859, from search snippets). Per LMS 2024 (as quoted in a snippet) it proves McKinnon's
  conjecture for many split rational surfaces, "including all split smooth rational surfaces of Picard rank at most
  five". A dP2 has geometric Picard rank 8 and our `S` is maximally non-split (`Pic(S) = ZK`), so it is outside. The
  method (fibres of Hirzebruch-surface morphisms, nef cone) is curve/section based.
* **Broustet 2006** (arXiv math/0601225, "Constantes de Seshadri du diviseur anticanonique des surfaces de del
  Pezzo"; snippet): `eps(-K_{X_7}; x) = 4/3` at a general point of the blow-up of `P^2` in 7 general points, 1 at other
  points. This confirms §2.2. It concerns Seshadri constants only.
* **Pagelot** (essential constant), **Huang** (toric surfaces; dP6): `alpha^ess = tau >= dim` in every computed case
  (papers_check README §1.2, read_huang_local.md).
* **LMS 2024**: essential boundedness is a hypothesis (their Thm 1.5); unconditional results are for split surfaces and
  concern curves realising `alpha <= alpha^ess`, not lower bounds above `tau`.

**Conclusion of Part (2).** I found no unconditional statement `alpha^ess_P > tau_P` for any rational point of any
dP2 or dP1, nor for any variety with `tau < dim` (this repeats, for surfaces, the papers_check finding). Memory-quotes are
labelled where used (§4.5).

Sources (search results only): [Broustet, arXiv math/0601225](https://arxiv.org/abs/math/0601225);
[Castaneda thesis record](https://uwspace.uwaterloo.ca/handle/10012/14859);
[LMS 2403.02480](https://arxiv.org/pdf/2403.02480); [McKinnon 2310.01640](https://arxiv.org/pdf/2310.01640);
[Huang 1809.02001](https://arxiv.org/pdf/1809.02001); [McKinnon-Satriano 2004.05212](https://arxiv.org/pdf/2004.05212);
[weak weak approximation dP2, 2303.09299](https://arxiv.org/pdf/2303.09299).

## 4. Part (3a): the elliptic pencil mechanism

**4.1 The pencil as a rational elliptic surface.** Lines `l` through `p` pull back to `E_l = pi^{-1}(l) in |-K|`,
which all pass through `P` and `P'` (`(-K)^2 = 2`). So `St = Bl_{P,P'} S -> P^1` is a rational elliptic surface with
fibre class `F = L - E_P - E_{P'}` (here `L = pi^*(-K)`), zero section `O = E_P` and a second `Q`-section
`sigma' = E_{P'}`. All fibres are irreducible (an `E_l` splits only over a bitangent through `p`, excluded in §1), so
`MW(St/\bar Q(t)) = E_8` (Shioda) and the height is `<s,s> = 2 + 2 s.O`. Since `Pic(S) = ZK`, the Galois-invariant
part of `Pic(St)` is `<L, E_P, E_{P'}>`, so `MW(St/Q(t))` has rank 1; `sigma'` has height 2, the minimum of `E_8`, so
`MW(St/Q(t)) = Z sigma'`.

*Lemma 4.1 (the curves `Gamma_k`).* For `k in Z` let `Gamma_k` be the image in `S` of the section `k sigma'`. Its class
on `St` is

    [Gamma_k] = O + k^2 F + k phi',   phi' = sigma' - O - F = 2E_{P'} - L,
              = (k^2 - k) L - (k^2 - 1) E_P - (k^2 - 2k) E_{P'}.

So `Gamma_0 = E_P`, `Gamma_1 = E_{P'}`, **`Gamma_2 = D`** (`2L - 3E_P`), `Gamma_{-1} = D' = sigma(D)` (`2L - 3E_{P'}`), and
`Gamma_k^2 = -1`, `K.Gamma_k = -1`, `F.Gamma_k = 1` for all `k` (`fibersec.py`). On each fibre, with `P` as origin,
`Gamma_k cap E_l = [k]P'`. `Gamma_k` is a rational curve of degree `2(k^2-k)` with `k^2 - 1` branches through `P` (one for each fibre on which
`k sigma'` meets `O`, i.e. `P'` is `k`-torsion), so `alpha_{P,Gamma_k} >= k^2 - k` (MR Thm 2.16).
*Proof.* Shioda's formula `[s] = O + (1 + s.O)F + phi_s` with `phi_{k sigma'} = k phi'` and
`2 + 2(k sigma').O = <k sigma', k sigma'> = 2k^2`. `QED`

**4.2 Fibre sections, and why they are global sections in disguise.** On `E_l` (line `(x,y) = m(a',b')`) use
`t = m/c`; the curve is `s^2 = g_l(m,c) = sum_j f_j(a',b') c^{4-j} m^j`, `f_0 = 1`, `f_j` a binary form of degree `j`,
and near `P` (resp. `P'`) `s = +c^2 w(t)` (resp. `-c^2 w(t)`), `w(t) = sum_n w_n(a',b') t^n`, `w_n` a form of degree
`n`. The space `L_k = H^0(E_l, -kK|_{E_l})` has the integral basis `m^i c^{k-i}` (`0 <= i <= k`),
`s m^j c^{k-2-j}` (`0 <= j <= k-2`), dimension `2k`. For `a + b = 2k - 1`, `a >= k+1`, the sections with
`ord_P >= a`, `ord_{P'} >= b` form a line, spanned by `F_{k,a,b}(l) = sum alpha_i m^i c^{k-i} + s sum beta_j m^j c^{k-2-j}`
whose coefficients are maximal minors, i.e. binary forms in `(a',b')`.

*Lemma 4.2 (height cost).* `beta_j` and `alpha_i` are forms of degree `-mu + j`, resp. `-mu + i`, with

    mu = mu(k,a,b) = (2k^2 + 1 - a^2 - b^2)/2   ( = deg pi_* O_St(-kK - aE_P - bE_{P'}) ),

without common factor (checked for `(k,a,b) = (3,5,0),(3,4,1),(4,7,0),(4,6,1),(5,9,0)` in `fibersec.py`; the degree
count is general: row `i` of the minor matrix has entries `w_{i-j}` of degree `i - j`). Hence the primitive integral
`F^{int}` has real-normalised size `||F^{int}|| := max_i |alpha_i| A'^{-i}, |beta_j| A'^{-j} ~ A'^{-mu}` (the content
is bounded because coprime forms take values with bounded gcd at coprime arguments).

*Fibre Liouville with forced divisibility.* If `Q` is on `E_l` near `P` and `F^{int}(Q) != 0`, then
`|F^{int}(Q)| <= C ||F^{int}|| c^k X^a` (Taylor at `t = 0` in the rescaled coordinate `A't`, where the fibre is
uniformly bounded), and `m_+^a m_-^b` divides `F^{int}(Q)` up to a bounded factor (at odd `p | m_+`, `Q` is `p`-adically
on the branch `s = +c^2 w(t)` with `|t|_p = |m|_p`; similarly at `P'`). So

    (4.1)    X^a c^k A'^{-mu} >= kappa m_+^a m_-^b.

**Proposition A1 (fibre Liouville = two-point global Liouville).** Substituting `A' = Xc/m`, (4.1) is *identical* to

    (4.2)    X^{a-mu} c^{k-mu} >= kappa m_+^{a-mu} m_-^{b-mu},

which is Liouville with forced divisibility at `P` and `P'` for the global section of the class
`(k-mu)L - (a-mu)E_P - (b-mu)E_{P'}` (unique, since `pi_* O(-kK-aE_P-bE_{P'}) = O(mu)`). For `(k,a,b) = (k,2k-1,0)`
(`mu = -k(k-2)`) this class is exactly `Gamma_k`, and the exceptional point `R_k(l)` of `F_{k,2k-1,0}(l)` is
`Gamma_k cap E_l = [k]P'`. *Proof:* (4.1) times `(m/(Xc))^{-mu}`; and the class computation of 4.1. `QED`

So the elliptic pencil, used through Liouville on the fibres with the natural integral normalisation, produces
nothing that a global section does not: the height cost `A'^{-mu}` is exactly the twist that turns the fibre section
into a global one. (For non-extremal fibre sections, `a + b < 2k - 1`, the analogous statement goes through the
Harder-Narasimhan slopes of the vector bundle `pi_* O(-kK - aE_P - bE_{P'})`; I only sketch this and do not use it.)

**4.3 The complete two-point section cap.**

**Theorem A (two-point Liouville cap on `S`).** Let `G` be any nonzero section over `Q` of `-dK`, with
`ord_P G >= a`, `ord_{P'} G >= b`, and `Q` a rational point near `P` with `G(Q) != 0`. Write the resulting
inequality `c^d X^a >= kappa m_+^a m_-^b` in exponents: `a x <= d - a nu_+ - b nu_-`. Then this inequality is a
consequence of the inequalities attached to the curves `Gamma_k`, `|k| >= 2`, and to the fibres, at every `Q` lying on
none of the (at most three) `Gamma_k` used. Consequently the best exponent obtainable from Liouville with forced
divisibility at `P` and `P'`, with any family of sections, is

    x <= x_sec(nu_+, nu_-) := min( 1 - nu_+ - nu_- ,  min_{|k|>=2} [ (k^2-k) - (k^2-1) nu_+ - (k^2-2k) nu_- ] / (k^2-1) )

(for `k >= 2` the bracket is `k/(k+1) - nu_+ - (k(k-2)/(k^2-1)) nu_-`; the `k <= -2` terms exceed `2/3` for
`nu_- <= 1/3`). In particular `x_sec(0, nu_-) = 2/3` for `nu_- <= 2/9`, and
`x_sec(nu_+, nu_-) < 2/3` iff `nu_+ > 0` or `nu_- > 2/9` (the `k = 3` curve).

*Proof.* (i) *Cone.* `dL - aE_P - bE_{P'}` is a `Q`-effective class on `St`. Over `\bar Q`, by the cone theorem
(`K_St = -F`, `F` nef, fibres irreducible) `NE(St)` is the closed cone spanned by `F` and the `(-1)`-curves, which are
sections. Averaging over Galois is the orthogonal projection `pi_V` onto `V = <L,E_P,E_{P'}> = <O, F, phi'>`. For a
section `s` with `s.O = n`, `pi_V[s] = O + (n+1)F + tau phi'` with `tau = -<phi_s,phi'>/2 in (1/2)Z` and
`tau^2 <= n+1` (Cauchy-Schwarz in the negative definite `E_8`). If `tau = k in Z`:
`pi_V[s] = Gamma_k + (n + 1 - k^2)F`. If `tau = k + 1/2`: `n + 1 >= k^2 + k + 1` (an integer `> tau^2`), and
`pi_V[s] = (Gamma_k + Gamma_{k+1})/2 + (n + 1 - k^2 - k - 1/2)F`. So the `Q`-effective cone in `V` is contained in the
cone spanned by `F` and all `Gamma_k`, and it contains it. (ii) *Inequalities.* For fixed `(x, nu_+, nu_-)` the
inequality `a x + a nu_+ + b nu_- <= d` is linear in `(d,a,b)`; it holds for `F = (1,1,1)` (because `A' >= 1`), for
`E_P, E_{P'}` trivially, and for `Gamma_k` at points off `Gamma_k` (Liouville with forced divisibility for the section
defining `Gamma_k`, a fixed polynomial with integer coefficients). By Caratheodory, `(d,a,b)` is a nonnegative
combination of at most three generators. (iii) The formula: `Gamma_k = (k^2-k, k^2-1, k^2-2k)`. `QED`

*Reading.* The section class, even with the second point `P'` and its forced divisibility, beats `3/2` only for
approximants that are also `p`-adically close to `P` (`nu_+ > 0`, a trivial gain: then `X >~ m_+ c^{-2/3}`) or
`p`-adically close to `P'` on a scale `m_- > c^{2/9}` (through `Gamma_3`, of class `6L - 8E_P - 3E_{P'}`). This is a
restricted theorem that narrows the window on the set `{Q : m_-(Q) >= H(Q)^{2/9 + eps}}` (I did not determine whether this
set contains Zariski-dense sequences tending to `P` with `x` near `2/3`; heuristically such points are rare, since a
window point with `m_- >= c^{2/9}` needs a line of height `A' <= c^{1/9}`):

    alpha({Q_i}) >= 1 / (3/4 - (3/8) nu_-)  > 3/2   along sequences off Gamma_3 with m_-(Q_i) >= H(Q_i)^{nu_-}, nu_- > 2/9.

It is honest to say what it is: a Liouville bound with forced divisibility at two rational points (T4 class); the
pencil only helped to find the curves `Gamma_k`. It says nothing about the generic approximant, which has
`m = gcd(a,b) = 1` (§5.1 count: the window points, if any, have `m` bounded on average). For `m = 1` the cap is
exactly `x <= 2/3`, i.e. `alpha >= 3/2`.

**4.4 Calibration: a dP2 point where Liouville is attained by a curve (the fibrewise barrier).**

*Proposition A3.* Let `S_0 : w^2 = z^4 + z f_3(x,y) + x^4` with `f_3 = 3x^3 + 2x^2y + xy^2 + 2y^3`, `P_0 = (0:0:1:1)`.
Then `S_0` is smooth (Groebner), `P_0` is off the ramification curve and off the 56 lines
(`Res(F_3, F_4') = 1`), `D_0` is irreducible, `tau_{P_0} = 3/2`, and the cuspidal cubic

    C_0 = { f_3(x,y) = 2x^2 z,  s = z^2 + x^2 }   (cusp at P_0 with rational tangent x = 0)

carries the points `Q_b = (1 : b : N : N^2+1)`, `N = f_3(1,b)/2` (`b` odd), with `u(Q_b) = s - T = 1` and
`X^{3/2} H -> 1` (`universal.py`: `0.99999...` at `b = 100001`). So **`alpha_{P_0}(S_0) = 3/2 = tau`**, attained on a
rational curve (MR Thm 2.16: `d/(r m) = 3/(1*2)`). `S_0` has the `Q`-lines `E_0 = {z = 0, s = x^2}` and `sigma(E_0)`;
`C_0 + E_0` is the member `{s = z^2 + x^2}` of `|-2K|`, so `C_0 in |-K + sigma(E_0)|` is the tangent-plane section of the
cubic surface obtained by contracting `sigma(E_0)`, cuspidal because `P_0` is special.

*Consequence (barrier for fibrewise methods).* The fibre of `S_0`'s pencil through `Q_b` is
`E_N : s^2 = c^4 + 2N c m^3 + m^4`, and `Q_b = (m,c) = (1,N)`. `C_0` is a section of `St_0` of odd degree, so it is
not a multiple of `sigma'` (all `Gamma_k` have even degree); its height is `2 + 2 C_0.O = 6` and `sigma'` has height 2,
so `C_0` and `sigma'` are independent in `MW(St_0/Q(t))` (a dependence would force `3 = lambda^2`, `lambda` rational).
By Silverman's specialization theorem, for all but finitely many `b`, `Q_b - P` is independent of `P' - P` in
`E_N(Q)`. Hence: *no inequality `X(Q) >= kappa c^{-2/3+eta}` (`eta > 0`) can hold for all rational points `Q` with
`m = 1` of all fibres of all such pencils, even after discarding on each fibre the subgroup generated by `P' - P` and
torsion.* Any argument that treats the fibres one at a time, with constants depending only on the fibre data
`(A', c, m_+-)` and exceptional points defined from `(E_l, P, P')`, is therefore capped at `alpha = 3/2`. The
exceptional set must be global and surface-specific (on `S_0` it must contain `C_0`; on our `S` no such curve exists,
since `Pic(S) = ZK` forbids odd degrees).

**4.5 What a uniform elliptic-logarithm bound would have to say.** For `m = 1`, `Q` is an integral point of the
quartic model `s^2 = G_l(c) = c^4 + f_1 c^3 + f_2 c^2 + f_3 c + f_4` (`f_j = f_j(a',b')`, `|f_j| <~ A'^j`), whose two
points at infinity are `P, P'`; and `X = A'/c`. Runge's method (equivalently Liouville for `D`) gives `c <~ A'^3`, and
by `E_N` the exponent 3 is sharp in the weighted family (the weighted analogue of Masser's bound `|x| << H^3` for `y^2 =` monic quartic, known to me only from a search
snippet; the sharpness of the weighted exponent is our own example `E_N`). The window
`3/2 < alpha < 2` is `A'^2 <~ c <~ A'^3`. So the needed statement is

    (UE)   for l outside a thin set and Q in E_l(Q) with m = 1, off a fixed curve:   log c <= (3 - eta) log A'.

In elliptic-logarithm language: with `P` as origin, `|z(Q - P)| ~ X = A'/c`, the curve has `h(E_l) ~ log A'`
(discriminant a form of degree 12), and `h(Q) ~ log c` lies between `2h` and `3h`. A lower bound for `|z|` beats
`3/2` iff it has the form `log(1/|z(R)|) <= (2/3 - eta) log H(R)` *in the regime `h(R) ~ (2..3) h(E)`*: linear in
the heights with a sharp constant. Known lower bounds for elliptic logarithms (David; Hirata-Kohno; memory-quote, not
re-checked) are of the form `log(1/|z|) <= C(r) h(E)^{O(1)} (log h(R))^{r+2}` with `C(r)` exponential in the rank `r`
of `E_l(Q)`, which is unbounded along the pencil; at `h(R) ~ h(E)` they are worse than trivial. Mumford's gap
principle and uniform Siegel bounds count points per fibre and cannot exclude one point on each of infinitely many
fibres. The Siegel descent route (`Q - P = 2R + T`, then Roth at the algebraic target `-T/2`) needs approximants of
height much larger than the target and curve heights; here (scales only, heuristic) `h(R) ~ h(Q)/4 ~ (1/2..3/4) h(E_l)`, below both, so every
Thue-Siegel-Roth gap principle is vacuous. And Prop A3 shows that any bound that depends on the fibre data only
cannot exist. **Verdict on (a): the pencil reduces Question B(dP2) to (UE), a uniform Runge-exponent improvement
in a one-parameter family of genus-one quartic models; no known bound has the needed shape, and fibre-local bounds
provably cannot.**

## 5. Part (3b): the square condition

**5.1 Elementary form of Question B(dP2).** Write `T = c^2 + c L(a,b) + Q(a,b)` (`L = x - y/2`,
`Q = -x^2/2 - xy/2 + 3y^2/8`) and `u := s - T(a,b,c) in (1/8)Z` (the value of the section defining `D`). Then

    (5.1)   u (2T + u) = f - T^2 = c F_3(a,b) + F_4'(a,b) = m^3 ( c F_3(a',b') + m F_4'(a',b') ),
    (5.2)   |u| <= C c^2 X^3,   and   |u| ~ c^2 X^3 |w_3(direction)|.

`u != 0` off `D`, and `|u| >= 1/8` is Liouville (`x <= 2/3`). The window `3/2 < alpha < 2` is
`1 <~ |u| <~ c^{1/2}`, equivalently `A^2 <~ c <~ A^3` with `A = max(|a|,|b|)`. Solving (5.1) for `c`: with `U = 8u`,
`16U c^2 + (2U L_8 - 64F_3) c + 2U Q_8 + U^2 - 64F_4' = 0` (`L_8 = 8L`, `Q_8 = 8Q`), so `c` is rational iff

    (5.3)   Delta_U(a,b) := (64F_3 - 2U L_8)^2 + 64U (64F_4' - 2U Q_8 - U^2)   is a perfect square.

`Delta_U` is an inhomogeneous sextic in `(a,b)` whose top form `(64F_3)^2` is a square. Heuristic count
(**heuristic**): `A^2` pairs, `U` values, square probability `A^{-3}`: `U/A` points per dyadic box. At exponent `gamma`,
`U ~ A^3/c` and `c ~ A^{gamma/(gamma-1)}`, so the total over dyadic scales is finite for every `gamma < 2` and is
dominated by `m = gcd(a,b) = O(1)`. This matches "H^3 r^2 H^{-2}" with critical `r = H^{-1/2}`.

**5.2 `D` is harmless in the window (a theorem).**
*Proposition B1.* For every `gamma < 2` only finitely many `Q in D(Q)` satisfy `X(Q) <= H(Q)^{-1/gamma}`.
*Proof.* The line of slope `tau = b'/a'` meets `pi(D)` again at `(F_3(a',b')a' : F_3(a',b')b' : -F_4'(a',b'))`, divided by
a gcd that is bounded (`F_3, F_4'` coprime). Near the real root `theta` of `F_3(1,t)`, `F_4'(1,theta) != 0`
(`Res != 0`), so `H ~ H(tau)^4` and `X ~ |F_3(a',b')|/|a'|^3 ~ |tau - theta|`. `X <= H^{-1/gamma}` becomes
`|tau - theta| <~ H(tau)^{-4/gamma}` with `4/gamma > 2`, finitely many by Roth (`theta` a real cubic irrationality). `QED`
(These points have `m = |F_3(a',b')|/g ~ A' ~ c^{1/4}`, all on the `P`-side: `nu_+ = 1/4`, `x = 1/2`; see T5.)

**5.3 Sieves and descent.**
* *Square sieve / quadratic large sieve in a shrinking box.* They bound
  `#{|a|,|b| <= A, c ~ C : f(a,b,c) = square}` from above by (number of candidates) times a saving coming from
  character sums of `f` (Weil). The best conceivable output is `O(A^2 C^{1/2+eps})` (sieve level `<= C^{1/2}`) or a
  main term of the expected size `A^2/C` plus an error `>= C^{-o(1)}`. To prove `alpha^ess > 3/2` one needs the window
  to be empty off a fixed curve for all large scales, i.e. a count `< 1`; an upper-bound sieve cannot give this
  (record barrier §5 item 2). The calibration surface `S_0` (Prop A3) has about one Liouville point per dyadic scale,
  on a curve, so "at most `C^{-o(1)}` points per scale" is genuinely the wrong currency: what must be shown is that the
  few points that exist lie on finitely many curves.
* *2-descent over the splitting field of `f`.* For each of the 28 bitangents `l_i` (over its field `K_i`),
  `f = q_i^2 + l_i g_i` and `(s - q_i)(s + q_i) = l_i(a,b,c) g_i(a,b,c)`, with `gcd | 2q_i`. These factorizations
  describe the descent data of `Q` (which torsor under the Neron-Severi torus it lifts to), but near `P` all the
  factors have their generic size (`l_i(Q) ~ c`, `g_i(Q) ~ c^3`, `s +- q_i ~ c^2`) because `P` is on no line; they do
  not see the proximity to `P`. The only factorization that sees it is (5.1), coming from `D` itself, and `D` is twice
  `-K`: its descent is the double cover branched along `D`, which is §6.1.
* *Geiser involution.* `sigma` maps window approximants of `P` to window approximants of `P'` with the same height and
  exchanges `m_+` and `m_-`; with the Bertini-type involution of §6.2 it generates a group acting losslessly on the
  section inequalities. It does not create new inequalities.

**Verdict on (b): every arithmetic input of square-sieve or descent type either is averaged (cannot produce an empty
window) or is the factorization (5.1) of the section `D`, i.e. Liouville. The non-averaged content of the square
condition is the integral-point problem (5.3), treated in §6.1.**

## 6. Part (3c): other ideas

**6.1 The top of the window is an integral-point problem on twisted K3 surfaces (new reformulation).**
Fix `U = 8u in Z \ {0}`. By (5.3) the rational points `Q` of `S` near `P` with `u(Q) = U/8` inject into the integral
points `(a,b,v) in Z^3` of the affine surface `Y_U^o : v^2 = Delta_U(a,b)`. Homogenising with `zeta`
(`Delta_U` has only terms of degree 6, 4, 2, 0 in `(a,b)`), `Y_U` is the double plane

    Y_U : v^2 = (64F_3 - 2U L_8 zeta^2)^2 + 64U zeta^2 (64F_4' - 2U Q_8 zeta^2 - U^2 zeta^4)   in P(1,1,1,3),

branched along a sextic `B_U`, so `Y_U` is a K3 surface with ADE singularities, and `Y_U ->` `S` (via `c` and `s`) is
generically finite. The integral points are integral with respect to the line at infinity `zeta = 0`. There
`v^2 = (64F_3)^2`, so the boundary splits into two rational curves `L_+-: v = +-64F_3(a,b)`; they meet over the three
roots of `F_3`, where `B_U` has nodes (`64^2 U F_4'(theta) != 0` because `Res(F_3,F_4') != 0`), i.e. `A_1` points of
`Y_U`. On the minimal resolution `Yt_U` the total transform of the boundary is the reduced divisor

    Dt_U = Lt_+ + Lt_- + N_1 + N_2 + N_3   (five (-2)-curves; Lt_+ and Lt_- each meet every N_i once; N_i conjugate over
                                             the cubic field of F_3),   Dt_U ~ H = pi^*O(1),  H^2 = 2,

so `K_{Yt} + Dt_U = H` is nef and big: **`(Yt_U, Dt_U)` is of log general type**. Up to the twist `(a,b,v) -> (-a,-b,-v)`
all `Y_U` are quadratic twists of one surface (`Delta_{lambda^2 U}(lambda a, lambda b) = lambda^6 Delta_U(a,b)`).

*Two equivalent models of the near-tight points.* (i) On `Y = Bl_P S` the divisor `Dt + E` is simple normal crossings
(`Dt = P^1`, meeting `E` transversally in the three tangent directions) and `K_Y + Dt + E = pi^*(-K) - E` is nef and big
(`(L-E)^2 = 1`): **`(Bl_P S, Dt + E)` is of log general type.** A point with `u(Q) = U/8` is `Sigma_U`-integral with
respect to `Dt + E` for `Sigma_U = {infinity, 2} cup {p | U} cup {bad primes of the model}` (it reduces into `pi^* D = Dt + 3E` exactly at the primes
dividing `u`). (ii) By Chevalley-Weil (`D in |2(-K)|`), such points lift to integral points of finitely many quadratic
twists of the double cover branched along `D`, i.e. of the K3 surfaces `Y_U` above.

*Theorem C (conditional; new).* Assume the Lang-Vojta conjecture for the log general type pair `(Bl_P S, Dt + E)`:
for every finite set of places `Sigma`, the `Sigma`-integral points are not Zariski dense (equivalently, by (ii), for the
pairs `(Yt_U, Dt_U)`). Then for every `K > 0` the set
`{Q in S(Q) near P on the P-sheet : X(Q)^{3/2} H(Q) <= K}` is contained in a proper Zariski-closed subset `Z_K` of `S`.
Equivalently: Liouville's inequality `X^{3/2} H >= kappa` is never tight on a Zariski-dense set.
*Proof.* By (5.2), `|u| <= C X^3 c^2 <= C K^2` on that set, so it lies in `D` union the finitely many sets
`{Q : u(Q) = U/8}`, `0 < |U| <= 8CK^2`; each is a set of `Sigma_U`-integral points of `(Bl_P S, Dt + E)`, not Zariski
dense by hypothesis, and a finite union of non-dense sets is not dense. `QED`

This is the dimension-2 analogue of the coordinator's "Question B at the top of the window is a stability statement
for near-tight Liouville": *near-tight Liouville points are integral points on twists of a K3 double cover branched
along `D`*. It is strictly weaker than `alpha^ess > 3/2` (it does not exclude `|u| -> infinity` slowly), and it is
consistent with `S_0`, where the near-tight points (`u = 1`) do exist but fill a curve (`C_0`), which is exactly what
Lang-Vojta allows.

*Why known integral-point theorems do not prove the hypothesis.* Subspace-theorem methods (Corvaja-Zannier, Levin,
Autissier, Ru-Vojta) need the boundary components to be "large". Ru-Vojta's inequality
`sum_i beta(L, D_i) m_{D_i}(x) <= (1+eps) h_L(x)` (off a proper closed set), with
`beta(L, D_i) = Int_0^inf Vol(L - tD_i) dt / Vol(L)`, gives degeneracy of integral points only when
`sum_i beta(L,D_i) h_{D_i} > (1+eps) h_L` along them; when one component carries the height this needs
`beta(L, D_i) > 1` for that component.
* On `(Bl_P S, Dt + E)`: near-tight points approach `E` (at a free direction), so `h_E` carries the height. For
  `L = pi^*(-K) - sE`, `beta(L, E) = Int_s^{3/2} Vol(pi^*(-K) - rE) dr / Vol(pi^*(-K) - sE)`, maximal at `s = 0` where
  it is **`beta_P(-K) = 17/18 < 1`** (§2.2). So the Subspace route on this pair is exactly McKinnon-Roth Thm 5.1, which
  gives `alpha^ess >= 17/18`; it does not even reach Liouville. (`L` outside `<L, E>` was not tried.)
* On the K3 model: the points approach `L_+` (from (5.3), `v = 32Uc + 2UL_8 - 64F_3 ~ +64F_3`). With `L = H`, by the
  Hodge index theorem `Vol(H - tLt_+) <= ((H - tLt_+).H)^2/H^2 = (2-t)^2/2`, and `H - tLt_+` is not big for `t > 2`,
  so `beta(H, Lt_+) <= (1/2) Int_0^2 (2-t)^2/2 dt = 2/3 < 1` (inside the boundary configuration the Zariski
  decomposition gives `Vol = 2 - 2t - t^2/2` on `[0,2/3]`, `4(1-t)^2` on `[2/3,1]`, `beta = 4/9`).

The hypothesis of Theorem C is therefore open; the reduction moves the top of the window from "rational points near a
point" to "integral points on a log general type surface", a standard Lang-Vojta problem for which the available
criteria fail by an explicit margin.

*Evidence (exact search, `usearch.py`).* All points of `S` with `max(|a|,|b|) <= 400`, `c > max(|a|,|b|)`,
`0 < |8u| <= 24` were found by solving (5.3) exactly: there are four, all with `max(|a|,|b|) <= 11`
(`(-1,0,2)`, `(-11,4,246)`, `(9,-3,43)`, `(-11,4,18)`; the best Liouville ratio `X^{3/2}H` among those with
`x > 1/2` is `2.33` at `(-11,4,246)`, `x = 0.564`). None in `11 < A <= 400`. The heuristic predicts
`sum U/A` over dyadic `A`, i.e. `O(1)` points in total; the data are consistent with this (and say nothing more,
by the record's lesson that failed searches are not evidence of emptiness).

**6.2 A lossless symmetry: the Bertini-type involution fixing `P`.** Fibrewise inversion `iota` (origin `P`) is an
automorphism of `St` (relatively minimal elliptic surface), hence of `Bl_{P'} S` (a degree-1 del Pezzo surface when
`P'` is general, and `iota` is its Bertini involution, the base point of `|-K_{Bl_{P'}S}|` being `P`). It fixes `P`
with `d iota_P = -id` and acts on `V` by `O -> O`, `F -> F`, `phi' -> -phi'`, so

    iota^* L = 3L - 4E_{P'},     iota(Gamma_k) = Gamma_{-k}.

For an approximant `Q` with data `(x, nu_-)` (in exponents w.r.t. its own height), `iota Q` is an approximant of `P`
with the same `X` and height `c^{3 - 4 nu_-}` (up to constants): `alpha(iota Q) = (3 - 4kappa) alpha(Q)`,
`kappa(iota Q) = (2 - 3kappa)/(3 - 4kappa)` for `kappa = nu_-`, an involution with fixed points `kappa = 1/2, 1`.
Since `iota` permutes the generators of the effective cone (`Gamma_k <-> Gamma_{-k}`), the family of two-point section
inequalities of Theorem A is `iota`-invariant: the symmetry is lossless, like Galois averaging in the record's
anisotropy analysis. It maps generic approximants (`kappa = 0`) to bad ones (`alpha` tripled) and cannot be used to
improve exponents.

**6.3 Where an exceptional approximant must live (structural).** If `Q in E_l(Q)` and `Q - P` lies in
`Z sigma'(l) + E_l(Q)_tors`, then `Q` lies on some `Gamma_k` (plus torsion translates, finitely many curves), where
`alpha >= k^2 - k >= 2` for `|k| >= 2`. Otherwise either `rank E_l(Q) >= 2` (a rank jump over the generic rank 1) or
`k(Q - P) = n sigma'(l)` with `k >= 2`, i.e. `Q` is on one of countably many division curves; for `k = 2` these are
4-sections of genus `<= 9` by Riemann-Hurwitz over the twelve singular fibres (finitely many rational points each,
Faltings). So window approximants off the `Gamma_k` live on rank-jump fibres or on sporadic division curves. This is
consistent with all heuristics but proves nothing: the record's barrier list (§5 item 9) already says that dense
supplies with finitely many points per curve are not excluded by curve-by-curve finiteness.

**6.4 Small points (Lang's conjecture).** If window approximants had `h(Q - P) = o(h(E_l))`, Lang's conjecture
(proved by Hindry-Silverman for bounded Szpiro ratio) would exclude them. They do not: `h(Q - P)` is comparable to
`log c`, between `2` and `3` times `log A' ~ h(E_l)` (scales only). Lang's conjecture does not bite.

**Verdict on (c): the one genuinely new structural statement is the reduction of the top of the window to integral
points on the log general type pair `(Bl_P S, Dt + E)`, equivalently on twisted K3 double covers branched along `D`
(Theorem C, conditional on Lang-Vojta for that pair), with an explicit reason why Subspace-type integral-point theorems
do not apply (`beta(L, E) <= 17/18` on the pair, `beta(H, Lt_+) <= 2/3` on the K3 model).**

## 7. Pre-registered tests T1-T5 and the dimension-2 controls

The claims of this note are: Lemma 2.1 and the census of §2 (exact); Lemma 4.1, Lemma 4.2, Prop A1, **Theorem A**
(two-point section cap, with the restricted gain for `nu_- > 2/9`); **Prop A3** (calibration surface `S_0`, barrier for
fibrewise methods); Prop B1 (`D` finite in the window); **Theorem C** (conditional on Lang-Vojta for
`(Bl_P S, Dt + E)`). No growth or uniformity claim for the circle problem is made.

**T1 (fake test).** No growth claim, so not applicable as stated. Its dimension-2 analogue is the calibration surface
`S_0` of Prop A3, which agrees with `S` in every fibre-local and averaged statistic one could feed a sieve (smooth dP2,
`P` off lines, `tau = 3/2`, same weighted shape of fibres) but has a curve of exact Liouville points. Theorem A passes:
on `C_0`, `x = 2/3`, `nu_+- = 0`, and Theorem A allows exactly `x <= 2/3` (tight, as it must be). Theorem C passes:
on `S_0` the near-tight points (`u = 1`) fill the curve `C_0`, which a Lang-Vojta exceptional set may contain.

**T2 (`k = 4` negative control) and its dimension-2 analogues.**
* *`P^2`, `O(1)`, rational `p`* (`alpha^ess = tau = 1`, attained densely by `|a|,|b| <= c^eps`). Sections of `O(d)`
  vanish to order `<= d` at `p`; there is no second point; the pencil of lines has genus-0 fibres with
  `mu_k = deg pi_* O(kH - kE) = k > 0` (no cost, no gain). Theorem A's analogue gives `alpha >= 1` only. Correct.
  Theorem C's analogue ("`X H <= K` is not dense") is *true* on `P^2` (such points lie on the finitely many lines
  through `p` of height `<= K`) although `alpha^ess = tau`; this shows that Theorem C's conclusion is strictly weaker
  than `alpha^ess > tau`, and the note does not claim more.
* *Huang's `Y_4` at `[1:1]x[1:1]`* (`alpha^ess = tau = 2`, attained densely by the Pell curves `C_{a,b}`). The Pell
  points with `d^2 H <= K` lie on the finitely many `C_{a,b}` with `ab` bounded (Huang I Thm 5.5), so the Theorem-C
  type conclusion again holds while `alpha^ess = tau`: consistent. A Theorem-A type section cap is a family of Liouville
  inequalities for sections; for points with no `p`-adic closeness to other fixed points (`nu = 0`) it is capped by
  `tau = 2`, and Huang's Pell points, which attain `2` on fixed curves, are on the boundary of every such
  inequality but not excluded (each `C_{a,b}` is a curve of `alpha = tau`).
* *A mechanism that would prove the `k=4` analogue.* None of the statements here is an unconditional
  `alpha^ess > tau` statement for generic sequences, so none can prove a false analogue.

**T3 (split test).** Not applicable to Question B as stated. The analogue for dP2: on a surface with `Q`-lines the
tangent-plane cubics through `P` are `Q`-rational with `alpha = 3/2` (or a cuspidal one, `S_0`), so any proof of
`alpha^ess > 3/2` must place them in `Z`; Theorem C's hypothesis allows exactly that. Theorem A as proved uses
`Pic(S) = ZK`; on a split dP2 the `Q`-effective cone of `Bl_{P,P'}S` is larger (all 240 height-2 sections are rational)
and the cap must be recomputed. No Erdos-type consequence is involved.

**T4 (criticality).**
* *Lemma K1 analogue.* Forced divisibility at `P` exactly compensates the gcd: the `k = 2` inequality is
  `x <= 2/3 - nu_+` (in the reduced coordinates it is the same Liouville). This is the change-of-units phenomenon.
* *Section cap.* Theorem A **is** a section/Liouville statement with forced divisibility at two points; its content
  is that the complete cap is computed exactly (all `Q`-sections, all multiples, via the cone of `Bl_{P,P'}S`), and
  that a second point `P'` (`p`-adic closeness to `P'` while really close to `P`) is the only way sections can beat
  `3/2`. It does not beat `3/2` for the generic approximant (`m = 1`).
* *Pencil.* Prop A1 proves that fibre-wise Liouville on the elliptic pencil is the same as global two-point
  Liouville (the height cost `A'^{-mu}` is the twist making the fibre section global). So the "fibration by genus-one
  curves" mechanism, used through Liouville, is in the record's mapped class; used through Mordell-Weil it needs (UE),
  and Prop A3 shows that no fibre-local version of (UE) exists.
* *Class-L / multi-point.* Not used. Products of sections at several points `Q_1, Q_2` of `S` at the rational target
  `P` have index filtrations that are tensor products of the one-point filtrations; the extremal sections are pure
  tensors `D^{d_1/2} x D^{d_2/2}`, so Dyson-type arguments return Liouville (as in round 25).
* *Theorem C* is not a section bound: its input is a Lang-Vojta statement for integral points on a different
  (log general type) pair. Its Subspace-type attack collapses to McKinnon-Roth Thm 5.1 (`beta = 17/18`), which is the
  record's "Roth at a rational target gives `beta < tau`".

**T5 (points that exist must not be excluded).**
* *`D`* (`alpha = 2`): points have `x -> 1/2`, `nu_+ -> 1/4`, `nu_- = 0` (`dpoints.py`: e.g. `c = 5.2e28`,
  `x = 0.476`, `nu_+ = 0.254`). They lie on `Gamma_2 = D` (excluded from the `k=2` inequality by construction), and the
  `Gamma_3` inequality `x <= 3/4 - nu_+` holds at all 13 computed points and is asymptotically **tight**
  (`x + nu_+ = 0.664, 0.681, ..., 0.730 -> 3/4`). A nice consistency check of Theorem A.
* *`E_0`* (cuspidal fibre, `alpha = 2`): its points are `((f^2-e^2)/2 : 0 : e^2 : e^3 f)`; with `f = e + 2j`,
  `m = 2j(e+j)`, and the primes of `e + j` are on the `P'`-side (`s = e^3 f == -c^2`), those of `j` on the `P`-side. So
  `x = (1-delta)/2`, `nu_+ = delta/2`, `nu_- = 1/2` for `j = e^delta`; Theorem A requires
  `x <= min(1 - nu_+ - nu_-, 2/3 - nu_+, 3/4 - nu_+ - (3/8)nu_-) = (1 - delta)/2`: satisfied, with equality from the
  fibre term (these points lie on a fixed fibre). Passed.
* *`S_0`, `C_0`* (`alpha = 3/2`): tight, see T1.
* *Round-24 real twisted family, split Pell curves, 8- and 10-point families*: they concern `X`, not `S`; nothing here
  constrains them.

## 8. Transfer to the 4-fold X at e

**8.1 Covering families of higher-genus curves through `e`.** They exist: any 4-dimensional subspace
`V subset H^0(X, -K (x) m_e)` gives `X --> P^3` whose fibres are curves through `e` (complete intersections of three
anticanonical members through `e`, of large genus); one-parameter subgroups give rational ones. The mechanism of
Prop A1 is formal: a fibre section of `-kK|_{C_b}` with order `a` at `e`, normalised integrally, costs `H(b)^{-mu}`
with `mu` the degree of the corresponding direct image, and the resulting inequality is the Liouville inequality of a
global section of a twisted class on the total space, i.e. of a section of some `n(-K)` on `X` vanishing at `e` and
at the base points of the family. So fibre-wise Liouville on any such family is inside the record's section cap
(`tau_e = 10/3`, `15/4` on balanced profiles), up to forced divisibility at the other base points; balanced
near-tight A-points have no reason to be `p`-adically close to those (the analogue of `nu_- = 0`). Used through the
group structure of the fibres (genus 1) or Faltings (genus `>= 2`), the family needs a uniform statement of type (UE),
and the analogue of Prop A3 is the record's real twisted family and split Pell curves at exactly `15/4`. **No gain
transfers.**

**8.2 The analogue of Theorem C fails to transfer (computed).** On the dP2 the near-tight points are integral points of
the log general type pair `(Bl_P S, Dt + E)`. On `X` the near-tight points are the tuples with all ten pair cofactors
`|t_ab| <= K`; they are integral with respect to the ten "diagonals" `Z_ab = closure{x_a = x_b}` (`x_j = zeta_j/zeta_5`,
`x_5 = 1`), whose sum is the divisor of the Vandermonde minus a boundary part: with
`V = prod_{a<b}(zeta_a - zeta_b)/prod_j zeta_j^2`, `ord_{E_A} V = binom(|A|,2) - 2|A| in {-2,-3}`, so

    div_{3(-K)}(V) = sum_{a<b} Z_ab + B_{1,4},    B_{1,4} = sum_{|A| in {1,4}} E_A,

and A-points are not integral with respect to `B_{1,4}` (balanced profiles have every `n_A` large). Each `Z_ab` is smooth
through `e`, so on `Bl_e X`

    K + sum Zt_ab + E = pi^*(K + sum Z_ab) - (10 - 3 - 1) E = pi^*L' - 6E,    L' = 2(-K) - B_{1,4}.

Bigness needs `tau_e(L') > 6`. Exact computation of the maximal order at `e` of Laurent polynomials with Newton
polytope in `kP_{L'}` (`toric_order.py`, integer ranks): **orders 5, 10, 15 for `k = 1, 2, 3`** (calibration: the same
code gives 3, 6, 10 for `-K`, i.e. `floor(10k/3)`, matching round 23). Each value is exact (the next order is excluded by
a full-rank integer matrix), so for `k <= 3` the order per unit is exactly 5; I did not prove `tau_e(L') < 6` for all
`k` (for `-K` the per-unit order jumps from 3 to `10/3` at `k = 3`; a jump from 5 to above 6 would be needed). If, as the
data indicate, `tau_e(L') = 5 < 6`, the log canonical class of `(Bl_e X, sum Zt_ab + E)` is not even effective (moreover the boundary is the braid arrangement
through `e`, not normal crossings, and resolving it only subtracts more). **Lang-Vojta predicts nothing for this
integral-point problem.** On `X`, unlike on `S`, near-tight Liouville points are not controlled by integrality alone:
the top-of-window stability must use balance (the 15 cut heights being equal), which is exactly what Vojta's `D = E`
term `V_E(5)` and the heuristic count see. Likewise the four characters: with `Z_a = closure{x_a = 1}`,
`Z_a ~ sum_{A ni 5, A not ni a} E_A` and `K + sum_a Z_a` has empty polytope (singletons force `m_i >= 1`, `i <= 4`, so
`m_5 <= -4`, while `A = {5}` forces `m_5 >= -3`): not effective either.

This is a clean structural difference between the dimension-2 calibration and `X`. In general, for a divisor `D`
with multiplicity `mu` at the target on an `n`-fold, `K + Dt + E = pi^*(K + D) - (mu - n)E`. On `S`: `n = 2`, `mu = 3`,
`K + D = -K`, so the pair is of log general type iff `1 < tau_P(-K) = 3/2`: yes. On `X`: `n = 4`, `mu = 10`,
`K + sum Z_ab = L'`, so it would need `6 < tau_e(L')`, and `tau_e(L') = 5` (exact for `k <= 3`, not proved for all `k`): no.

**Answer to "how would the mechanism transfer":** the pencil/section part transfers and gives nothing beyond the record's
cap; the integral-point reformulation (the only new structural input here) does not transfer, by an explicit
computation.

## 9. Scripts and certificates

All in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/dp2_checks/`
(python3, sympy 1.14, numpy, python-flint 0.9; exact arithmetic unless stated).

| script | what it certifies | output |
|---|---|---|
| `geom.py` | `f`, 2-jet `T`, `F_3 = 2w_3`, `F_4'`, triple point, `F_3` irreducible, disc `-1667/4`, one real root, `Res(F_3,F_4') = -80353/1024` | stdout |
| `smooth_bitan.py` | Groebner smoothness over `Q`, `F_3`, `F_5` (pure powers of `x,y,z` among leading monomials); no bitangent through `p`; tangent-line discriminant `t^2 * (deg 10 irreducible)` | stdout |
| `search_f.py` | the search that produced `f` (seed 1) | stdout |
| `ffcount.py`, `frob_pic.py` | `#S(F_{3^n})`, `n <= 8`; traces `0,2,3,2,0,-1,0,2`; char. poly `(X-1)(X+1)(X^6-X^3+1)`; `n = 8` trace predicted 2, counted 2 | `frob_seed1_n8.out` |
| `deg6.py` | net `|-3K - 4P|` (dimension 3), tangent cones, square-tangent-cone members: eliminant irreducible octic in every chart (no `Q`-member) | stdout |
| `fibersec.py` | fibre-section minors: degrees `-mu + j`, no common factor; `Gamma_k^2 = -1`, `K.Gamma_k = -1`, `F.Gamma_k = 1` | stdout |
| `universal.py` | `S_0` smooth, `Res = 1`; points `Q_b` on `C_0` with `u = 1`, `X^{3/2}H -> 1` | stdout |
| `dpoints.py` | `D`-points from convergents of the real cubic slope: `u = 0`, `x -> 1/2`, `nu_+ -> 1/4`, `Gamma_3` inequality satisfied and asymptotically tight | stdout |
| `usearch.py` | all points with `max(|a|,|b|) <= 400`, `0 < |8u| <= 24`: four, all with `max <= 11` | `usearch_400_24.out` |
| `toric_order.py` | max order at `e` of sections of `k(-K)`: 3, 6, 10; of `k(2(-K) - B_{1,4})`: 5, 10, 15 (each exact by a full-rank check one order higher) | stdout |

Working notes: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/dp2_notes.md`.

**Open points for a referee.** (1) The cone step of Theorem A uses the cone theorem for the rational elliptic surface and
Shioda's height formula; both are standard but should be checked against a reference. (2) Lemma 4.2's degree count is
verified for five `(k,a,b)` and argued in general; the general argument should be written out. (3) The `beta`-values in
§6.1 inside the boundary configuration assume no further negative curves; only the Hodge-index upper bound `2/3` is
proved. (4) `tau_e(2(-K) - B_{1,4}) < 6` is computed for `k <= 3` only. (5) The dimension counts in §2.4 are heuristic.
