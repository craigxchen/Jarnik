# Round 27, kill report on lane dp2

Status: complete. Killer: adversarial check of `round27/dp2.md`. All numbers below were recomputed with my own code in
`round27/kill_dp2_checks/`; I did not import or run the author's scripts. Unrefereed (killer-checked only).

## 0. Verdict

**The mathematics survives. The novelty claim does not.**

*What survives.* Every exact claim of the write-up reproduces with independent code:
* smoothness;
* `P` off the 56 lines;
* `Pic(S) = ZK`;
* `D`, `T`, `F_3`, `F_4'` and the resultant;
* the cubic branch field;
* `eps = 4/3`, `tau = 3/2`, `beta = 17/18`;
* the degree-6 census;
* uniqueness of the sections defining `Gamma_2`, `Gamma_3`, `Gamma_{-2}`;
* the D-point calibration;
* `S_0` and `C_0`;
* the bounded-`u` search;
* the toric orders 3, 6, 10 and 5, 10, 15.

Theorem A, Prop A1, Prop A3, Prop B1 and Theorem C (conditional) are correct as stated, up to the minor defects listed in
§§2-6. I found no false mathematical statement that the conclusions rest on. I also **close the author's open
point (4)**: `tau_e(2(-K) - B_{1,4}) = 5` exactly, for all multiples. Moreover the log canonical class of the
diagonal-integrality pair has degree `-1` on the covering family of lines through `e`, **for every k >= 3**. So the
negative transfer holds uniformly; see §6.

*What does not survive: the novelty.* Every item the write-up lists as "new mathematics" is one of three things:
* a known barrier class of the record;
* a known mechanism of the literature;
* a weaker conditional than one the record already has.

Item by item:
* **Theorem A** is Liouville with forced divisibility (T4 class, the `tau` cap). For `gcd(a,b) = 1` it is *literally*
  the statement `tau_P(-K) = 3/2`. Its only gain beyond `tau` comes from `p`-adic closeness to a second fixed rational
  point `P'`, which is forced divisibility again. The cone computation itself is routine Mordell-Weil lattice
  bookkeeping on a rational elliptic surface.
* **Prop A1** (pencil = sections) is the dimension-2 instance of the record's change-of-units/Lemma K1 phenomenon.
* **Prop A3** is McKinnon-Roth's cubic-surface mechanism transported to a dP2: MR-L Thm 4.5, where a cuspidal
  tangent-plane section gives `alpha = 3/2` via MR-I Thm 2.16. As a "fibre-local barrier" it is the
  dimension-2 analogue of the record's calibration curves (split Pell curves, round-24 real family), not a new
  barrier.
* **Theorem C** is a consequence of Vojta's conjecture (Lang-Vojta for integral points on `S \ D`). It is weaker
  than ERR Prop 0.1 / LMS Prop 1.4, which already give `alpha^ess >= 2` from Vojta on `S`. Its exceptional set
  `Z_K` depends on `K`, so it does not reach `alpha^ess` (the author says so). It does not transfer to `X`.
* **The transfer computation** is correct. But it follows from a one-line covering-curve argument (§6), and it restates
  the record's existing conclusion that "only Vojta's `D = E` term and the heuristic count see the reduction"
  (Lemma K1, expert note §5 item 7).

Nothing in the write-up gives `alpha^ess > tau` at any rational point of any variety, conditionally or not. That is
beyond what Vojta already predicts. The lane's honest headline ("no unconditional `alpha^ess > tau` for generic
sequences; an exact map of the window") is correct.

**Novelty verdict: repackaged known barrier.** The section/Liouville `tau` cap with forced divisibility (T4) and the
Vojta-conditional class. The new pieces are an explicit certified dP2 instance and the exact two-point cone. Both are
useful for the map but are not new mechanisms.

## 1. Exact recomputation of the certified data (§§1-2 of dp2.md)

| claim | my check | script / output | result |
|---|---|---|---|
| `S` smooth; good reduction at 3, 5 | grevlex Groebner basis of `(f_x,f_y,f_z)` over `Q`, `F_3`, `F_5` has pure powers `x^3,y^3,z^7` (`x^3,y^5,z^7` mod 5) among leading monomials | `k_geom.py` | confirmed |
| `f(0,0,1) = 1`, `P` off ramification | direct | `k_geom.py` | confirmed |
| 2-jet `T = z^2 + z(x - y/2) - x^2/2 - xy/2 + 3y^2/8`; `f - T^2 = zF_3 + F_4'`; `F_3 = 2w_3` | series of `sqrt(f(x,y,1))` | `k_geom.py` | confirmed |
| `Res(F_3, F_4') = -80353/1024`; direction `(0,1)`: `F_3 = 19/8`, `F_4' = 119/64` | dehomogenised resultant with full degrees 3, 4, plus the point at infinity | `k_geom.py` | confirmed; so `D` is irreducible and no bitangent passes through `p` |
| branch field `19t^3 - 10t^2 + 20t + 8`, disc `-2^10 * 1667`, one real root | sympy; hand computation of the discriminant `-1707008` | `k_geom.py` | confirmed (1667 prime) |
| tangent lines from `p`: disc `= t^2 * (irreducible degree 10)` | discriminant of `f(m, mt, 1)` in `m` | `k_geom.py` | confirmed (10 `I_1` + one type II: Euler number 12) |
| `#S(F_{3^n})`, traces `0,2,3,2,0,-1,0,2`, char poly `(X-1)(X+1)Phi_18` | my own `F_{3^n}` arithmetic (primitive polynomial, log tables), all of `P^2(F_q)` for `n <= 8`; Newton identities | `k_frob.py` / `k_frob.out` | confirmed. This is the class of a Coxeter element of `W(E_7)` (exponents 1,5,7,9,11,13,17), so eigenvalue 1 has multiplicity 1 and `Pic(S) = ZK` |
| `eps_P = 4/3`, `tau_P = 3/2`, `beta_P = 17/18`; Zariski decomposition on `[4/3,3/2]` | by hand: `L - (4/3)E = L/9 + (4/9)Dt`; positive part `(3-2t)(3L - 4E)` with square `2(3-2t)^2`; `152/81 + 1/81 = 17/9` | notes | confirmed |
| `alpha_{P,D} = 2` (MR-I Thm 2.16, `r = 2` real branch, `r = 0` complex) | formula check | -- | confirmed |
| degree-6 census: the net `|-3K - 4P|` has dimension 3 and contains `xD`, `yD`; its square-tangent-cone members are not over `Q` | my elimination: the members with tangent cone `q^2` correspond to the intersection of two conics in `P^2_q`; the eliminant is an irreducible quartic over `Q` (`182753a^4 + 558544a^3 - 647968a^2 - 87808a - 157440`); no solution at infinity | `k_deg6.py` / `.out` | confirmed (4 members, quartic field) |
| the census logic "`alpha < 2` in degree `<= 6` needs mult 4 with square tangent cone" | redone by hand from `Dt.Ct >= 0` (`m <= 4n/3`) and MR-I Thm 2.16; degree 4 admits nothing (`3m <= 4 < 4m` has no solution) | notes | confirmed |

Two remarks on the census. Odd degrees are excluded only because `Pic(S) = ZK`. The census bounds `alpha` (the
infimum), not `alpha^ess`, so it constrains only which curves an exceptional set must contain.

## 2. Theorem A and Prop A1 (elliptic pencil, two-point section cap)

**2.1 The geometry is right.**
* `Bl_{P,P'} S` is a rational elliptic surface with fibre `F = L - E_P - E_{P'}`.
* All fibres are irreducible: a fibre splits iff its line is a bitangent through `p`, and §1 excludes that.
* So `MW = E_8` over `Qbar`. Galois-invariant rank 3 gives MW rank 1 over `Q`.
* `sigma'` has the minimal norm 2, hence generates.
* Shioda's formula gives `[Gamma_k] = (k^2-k)L - (k^2-1)E_P - (k^2-2k)E_{P'}`.

I checked independently (`k_gamma3.py`) that the corresponding linear systems are single members:
* `(d,a,b) = (2,3,0)`: dimension 1, so `Gamma_2 = D`;
* `(2,0,3)`: dimension 1, so `Gamma_{-1} = D'`;
* `(6,8,3)`: dimension 1, so `Gamma_3`; `(6,9,3)` and `(6,8,4)` have dimension 0;
* `(6,3,8)`: dimension 1, so `Gamma_{-2}`.

Here `h^0(-6K) = 43`, and the 42 imposed conditions are independent. The fibre identity behind Prop A1 is also right. On
`E_l`, `-K|_{E_l} ~ P + P'`, so a section of `-kK|_{E_l}` with `ord_P = 2k-1` has residual zero `kP' - (k-1)P`, which is
`[k]P'` in the group law with origin `P`.

**2.2 The cone step is right, and the closedness gap is harmless.** The author cites "the cone theorem". The missing
point is that `NE(St)` intersected with `F^perp` equals `R_+ F`. Proof: write `z = aF + w`, with `w` in the negative
definite `E_8 (x) R`, and take the Zariski decomposition `z = P + N`. Since `F` is nef, `F.P = F.N = 0`. There are no
vertical negative curves, because all fibres are irreducible, so `N = 0`. Then `z` is nef with `z^2 = w^2 <= 0`,
which forces `w = 0`.

Galois averaging is the orthogonal projection onto `V`, because the average of isometries is self-adjoint and
idempotent. The projected generators satisfy `Gamma_k / k^2 -> F`, so the cone generated by `F` and the `Gamma_k` is
closed. The `tau = k + 1/2` case is right: `n + 1 >= k^2 + k + 1`.

So Theorem A is correct. I also checked `x_sec`:
* `x_sec(0, nu_-) = 2/3` iff `nu_- <= min_{k >= 2} (k-1)/(3k) = 2/9`, the minimum being at `k = 3`;
* the fibre term needs `nu_- <= 1/3`, and `k = -2` needs `nu_- <= 1/2`.

**2.3 What Theorem A is.** At `m = gcd(a,b) = 1` we have `nu_+ = nu_- = 0`, and Theorem A reads `x <= 2/3`. The
proof needs nothing but the inequality `ord_P G <= (3/2) deg G` for every section `G`, which is the definition of
`tau_P = 3/2`. So the "complete cap" adds nothing for the generic approximant; the author says this too (§4.3). Its only
non-trivial content is the `nu_- > 2/9` regime. That is forced divisibility at a second fixed rational point, which is
the record's T4 class (toric forced divisibility; Lemma K1: forced divisibility and gcd scale identically).

Prop A1 is the record's change of units: the fibre height cost `A'^{-mu}` is exactly the twist turning a fibre section
into a global one. The author classifies both correctly as T4. The novelty claim that "the pencil identifies the full
cone" is accurate but routine.

**2.4 The restricted theorem may be vacuous.** The restricted theorem is `alpha > 3/2` along sequences off `Gamma_3`
with `m_- >= H^{nu_-}`, `nu_- > 2/9`. The author does not show that its domain contains any window points. A
heuristic count (mine, **heuristic**) gives about `C^{1-2x-nu}` window points per dyadic height `C`, at exponent `x`
with `m ~ C^nu`. This uses `#(a',b') ~ A'^2` with `A' = C^{1-x-nu}`, about `m` values of `m`, `|U| <~ C^{2-3x}`, and a
square probability of about `(mA')^{-3}`. For `x > 1/2` this is `< C^{-nu}`, so the domain is heuristically *no richer*
than the generic window. The extra `P'`-side congruence (`m^3 | s + T`) only lowers the count.

Further evidence:
* My D-point recomputation (14 convergents, `log10 c` from 6.7 to 32.9) has `x + nu_+` rising from 0.664 to 0.733,
  below `3/4`. This is consistent with the `Gamma_3` inequality being asymptotically tight on `D`.
* `E_0` has `nu_- = 1/2` but `x <= 1/2`, outside the window.
* `iota(Q)` has `nu_- = 2/3` but `x/3`.

None of the known points with large `nu_-` lies in the window.

## 3. Prop A3 (calibration surface S_0)

**3.1 Reproduced** (`k_s0.py`):
* `S_0` is smooth (pure powers `x^3, y^5, z^7`).
* `F_3 = f_3`, `F_4' = x^4`, `Res = 1`, `F_3(0,1) = 2`, so `P_0` is on no line and `D_0` is irreducible.
* `C_0` lies on `S_0` identically.
* For `b = 11, 101, ..., 100001`: `u(Q_b) = 1`, `gcd = 1`, and `X^{3/2}H = 0.9737, ..., 0.9999975 -> 1`.

The Mordell-Weil independence argument is right:
* heights 6 and 2 would force `lambda^2 = 3`;
* `E_8` is torsion-free;
* Silverman's specialization theorem applies (injective for all but finitely many fibres).

**3.2 A detail the author omits.** On `S_0`, `F_3 = f_3 = (x + y)(3x^2 - xy + 2y^2)`, which `sympy` reports as not
irreducible. So `D_0` has a **rational** branch at `P_0`, with `alpha_{P_0,D_0} = 4` (the other two branches are
complex), whereas `D` on `S` has a real cubic branch with `alpha = 2`. `S_0` is therefore special in four ways:
* it has `Q`-lines (`z = 0` is a hyperflex line);
* `f_1 = f_2 = 0` in every fibre quartic;
* `P_0` lies on the parabolic curve of the cubic surface obtained by contracting `sigma(E_0)`, which is what makes
  `C_0` cuspidal;
* the tangent cone of `D_0` has a rational line.

The T1 sentence "`S_0` agrees with `S` in every fibre-local and averaged statistic one could feed a sieve" is therefore
an overstatement (minor). The barrier itself is unaffected: it uses only the shape `|f_j| <~ A'^j` and
`(A', c, m_+-)`, and `S_0` satisfies those.

**3.3 Novelty.** `C_0` is the tangent-plane section of the cubic surface obtained by contracting `sigma(E_0)`, and it
is cuspidal at `P_0`. The value `alpha = 3/(1*2) = 3/2` is exactly the cuspidal case of McKinnon-Roth MR-L Thm 4.5,
for split cubic surfaces: `alpha_x(-K) = 3/2` iff the tangent-plane section `C_x` is cuspidal, or nodal with
real-quadratic slopes. This is in the record's full-text reading `papers_check/read_mr.md` §3.4; the paper's
existence and venue were confirmed by WebSearch. The computation itself is MR-I Thm 2.16. So Prop A3 is
a known mechanism on a new instance. Pulled back to a dP2 it gives `alpha = tau = 3/2` because there `tau = 3/2`,
whereas on the cubic surface `tau = 2` (by hand: `L - tE = (2-t)(2L - 3E) + (2t - 3)Ct_x` on `[3/2, 2]`).

The *barrier* reading ("no fibre-local uniform bound can beat 3/2") is a fake-test-style meta-statement. It is correct
as stated: the points `Q_b` have `m = 1` and are fibre-locally generic. It is the dimension-2 twin of the record's
calibration curves, the split Pell/Thm 2 curves and the round-24 real family at exactly `15/4`.

**3.4 (UE) is a restatement, not a reduction.** For `m = 1`, `X = A'/c`, so "`log c <= (3 - eta) log A'` off a fixed
curve" *is* the statement "`x <= 2/3 - eta'` off a fixed curve for `m = 1` points". Calling it a reduction of
Question B(dP2) "to (UE)" (§4.5 verdict) reframes the question. It does not reduce it.

## 4. Prop B1, Section 5 (square condition)

* `(5.1)` is right: with `T` the 2-jet, `f - T^2 = c(f_3 - 2LQ) + (f_4 - Q^2) = cF_3 + F_4'`, and
  `u(2T + u) = f - T^2`.
* `(5.3)` is right after clearing denominators: `16U c^2 + (2UL_8 - 8G_3)c + (2UQ_8 + U^2 - G_4) = 0`, with
  `G_3 = 8F_3`, `G_4 = 64F_4'`. My first transcription had a wrong factor and found nothing. After the fix,
  `k_usearch.py 400 24` finds **exactly the author's four points**: `(-11,4,18)`, `(-11,4,246)`, `(-1,0,2)` (on the
  ramification curve, `s = 0`, `U = -12`), `(9,-3,43)`. The ratio `X^{3/2}H = 2.326` at `(-11,4,246)`, `x = 0.564`.
  This is evidence only, as the author says.
* Prop B1 (finitely many `D`-points in the window, by Roth at the real cubic slope) is right. It is MR-I Thm 2.16 in
  quantitative form.
* The sieve/descent verdict (averaged, or the factorisation of `D`, i.e. Liouville) agrees with the record's barrier
  items 2 and 3. Nothing new is claimed there, and nothing is wrong.

## 5. Theorem C (conditional, log general type pair)

**5.1 The proof is correct.** Each step checks:
* `K_Y + Dt + E = (-L + E) + (2L - 3E) + E = L - E`, with `(L - E)^2 = 1`, and `L - E` nef (the pencil of fibres
  through `P`).
* `Dt` is smooth rational (`p_a(Dt) = 1 + (-1 + (-1))/2 = 0`) and meets `E` transversally in three conjugate points.
* `Y \ (Dt + E) = S \ D`.
* A point with `8u = U` fixed reduces into `D` only at primes dividing `U` (and 2, and bad primes), so it is
  `Sigma_U`-integral.
* `X^{3/2}H <= K` gives `|u| <= C K^2`, hence finitely many `U`.

The K3 model also checks out. `Y_U` (with `v^2 = Delta_U`) is the twisted double cover `W_U : (U/8) w'^2 = s - T` of
`S` branched along `D`, read as a double cover of the `(a:b:w')`-plane with `zeta = w'`:
* `Delta_{lambda^2 U}(lambda a, lambda b) = lambda^6 Delta_U` holds (term of degree `6-2j` carries `U^j`), so `Y_U`
  depends on `U` modulo squares, through the involution `zeta -> -zeta`;
* at the three roots of `F_3` on `zeta = 0`, the quadratic part `64^2 (F_3'(theta)^2 eps^2 + U F_4'(theta) zeta^2)`
  is nondegenerate, so these are nodes;
* the five boundary curves `Lt_+, Lt_-, N_1, N_2, N_3` match the ramification curve together with the `D_4`
  resolution graph over the triple point of `D` on `W_U`.

**5.2 What Theorem C is.** The hypothesis is Lang-Vojta for the integral points of the log general type surface
`S \ D`, with `D` irreducible in `|-2K|`. That is a special case of Vojta's conjecture. The conclusion is weaker than
what Vojta already gives through ERR Prop 0.1 / LMS Prop 1.4 (`alpha^ess_P >= 2`) (record, papers_check §1.3).
`Z_K` depends on `K` and on `Sigma_U`, so the statement says nothing about sequences with `X^{3/2}H -> infinity`
however slowly. That is the handoff constraint "do not treat an `S`-dependent exceptional set as uniform", and the
author respects it (§6.1: "strictly weaker than `alpha^ess > 3/2`").

The reformulation "Liouville tight `<=>` the `tau`-computing section has bounded value `<=>` integral points on the
complement of its divisor" is the standard bridge between approximation and integral points. I found no source that
states it for dP2, but I would not call it new mathematics. It is the definition of the local height of `D` read
twice. Lang-Vojta for the complement of a single irreducible curve in `|-2K|` on a dP2 (equivalently, integral
points on a K3 surface minus a rational curve plus a `D_4` configuration) is open. The known unconditional methods
(Corvaja-Zannier, Levin, Autissier, Ru-Vojta) need several large boundary components. So Theorem C moves an open
problem to an open problem.

**5.3 A minor imprecision in the Subspace paragraph.** For near-tight points with `L = pi^*(-K)` one has
`h_E ~ (2/3) h_L`. So the Ru-Vojta criterion `beta(L,E) h_E > (1+eps) h_L` needs `beta(L,E) > 3/2`, not `> 1`.

For `L = pi^*(-K) - sE` the requirement is `beta(s) > 3/2 - s`. In fact `beta(s) = (3-2s)/6` on `[4/3, 3/2)`, and
`beta(1) = 2/9 < 1/2`. So the criterion fails on the whole family, and by a larger margin than the write-up states.
The conclusion is unchanged. On the K3 model the threshold is about 1 (the height sits on `L_+`), and
`beta(H, Lt_+) <= 2/3` by Hodge index is right.

**5.4 T2 check of the mechanism.** If the analogous integrality pair at `k = 4` (on `T_4`) were of log general type,
Lang-Vojta would contradict gauss Thm 4.1 (dense balanced near-tight 4-tuples). §6 proves that it is *not* of log
general type, for every `k`. So Theorem C's mechanism makes no false prediction at `k = 4`. On `P^2` the analogous pair
is also not of log general type: it is `F_1` with `K + sum l_i + E = (d-3)H - (d-2)E`, which has degree `-1` on the
ruling. Its conclusion holds there for elementary reasons, as the author says.

## 6. Transfer to X (Section 8)

**6.1 Reproduced.** My own toric code (`k_toric.py`, rank mod `2^61 - 1` certifies non-existence one order higher;
exact `fmpz` rank certifies existence) gives these maximal orders at `e` of sections of `tD` on `L_5`:

| divisor | t = 1 | t = 2 | t = 3 |
|---|---|---|---|
| `-K` | 3 | 6 | 10 |
| `L' = 2(-K) - B_{1,4}` | 5 | 10 | 15 |
| `K + sum Z_ab` (orientation representative) | 5 | 10 | 15 |

The `-K` row matches round 23's `floor(10t/3)`. The identity `div_{3(-K)}(V) = sum Z_ab + B_{1,4}`, with
`ord_{E_A} V = binom(|A|,2) - 2|A|`, is right. So is `K + sum Zt_ab + E = pi^*L' - 6E`.

**6.2 The open point (4) closes, for all multiples and all `k` (proved here).** Let `k >= 3`. Use Kapranov's model:
`Lbar_k` is the iterated blow-up of `P^{k-1}` (coordinates `zeta_1, ..., zeta_k`) along the coordinate linear
subspaces, the singleton rays being the hyperplanes `zeta_i = 0`.

*Lemma (killer).* Let `l` be a general line in `P^{k-1}` through `e = (1:...:1)` and `lt` its strict transform on
`Bl_e Lbar_k`. Then `(K + sum_{a<b} Zt_ab + E) . lt = -1`. Since such lines cover `Bl_e Lbar_k`,
`K + sum Zt_ab + E` is not pseudo-effective.

*Proof.* `l` avoids every coordinate subspace of codimension `>= 2`: a line through `e` meets a fixed codimension-2
subspace only if it lies in the hyperplane spanned by `e` and that subspace, a proper closed condition. So `l`
meets `E_A` only for `|A| = 1`, once each, and `-K . l = k`. `Z_ab` is the strict transform of the hyperplane
`{zeta_a = zeta_b}`, which contains `e`. So `Z_ab . l = 1` and `Zt_ab . lt = 0`. On `Bl_e` we have
`K = pi^*K + (k-2)E`, so `K . lt = -k + (k-2) = -2`, and `E . lt = 1`. Total: `-2 + 0 + 1 = -1`. `QED`

At `k = 5` this gives `tau_e(L') <= 5`. There is also an elementary version: for `F` in `H^0(tL')`, the product
`F(e + s v) * prod_j (1 + s v_j)^t` is a polynomial in `s` of degree `<= sum_j (m_j + t) = 5t`, because the singleton
coefficients of `L'` are 1. So `ord_e F <= 5t`, or `F` vanishes on every line through `e`. With the order-5 section
at `t = 1`, **`tau_e(L') = 5` exactly**, and `pi^*L' - 6E` is not pseudo-effective. This removes the author's caveat
"exact for `k <= 3`, not proved for all `k`".

At `k = 4` the bound is `binom(4,2) - 4 = 2`, and my computation gives orders 2, 4, 6. So it is attained.

**6.3 Consequence for the novelty claim.** The negative transfer is uniform in `k`, and it has a transparent reason:
a general line through `e`, minus one point, is an `A^1` in the complement of `sum Zt_ab + E`. So the diagonal
complement is log-uniruled. Integrality with respect to the diagonals can therefore never carry a Lang-Vojta
prediction at `e`, for any `k`. That includes `k = 4`, where the prediction would be false (gauss Thm 4.1). It makes
no reference to `k = 5` or to balance.

The author's conclusion "balance must enter on `X`" is correct. But the record already states it: Lemma K1 (expert note
§3) says only Vojta's `D = E` term and the heuristic count see the gcd reduction, i.e. balance. So §8.2 is a correct
negative computation that restates a mapped barrier.

**6.4 §8.1 is a sketch.** The claim that fibre-wise Liouville on any covering family of curves through `e` on `X` "is
inside the record's section cap" is argued by analogy with Prop A1. No proof is given for families parametrised by
`P^3` (where the direct image is a vector bundle, not a line bundle). It should be labelled heuristic. Nothing else
depends on it.

## 7. Pre-registered tests T1-T5 (run by the killer)

**T1 (fake test).** The write-up makes no growth claim and nothing about the circle problem, so T1 is vacuous for the
lane as a whole. No statement in it would bound `M`. Its dimension-2 analogue, `S_0`, passes. Theorem A is tight on
`C_0` (`x = 2/3`, `nu = 0`), and Theorem C's exceptional set may contain `C_0`. Caveat (§3.2): `S_0` is not
"statistically identical" to `S`. **Pass, with an overstatement corrected.**

**T2 (`k = 4` negative control).**
* Theorem A and Prop A1 are unconditional only up to `tau`, plus forced divisibility at a second point. Their `X`
  analogue is the record's section cap `10/3` (`15/4` balanced). That cap is a valid lower bound at `k = 4` too, so
  it excludes nothing at `k = 4` that exists.
* Theorem C's mechanism at `k = 4`: by the Lemma in §6.2, the diagonal-integrality pair has log canonical degree `-1`
  on lines through `e`. My toric orders give `tau_e(K + sum Z_ab) = 2 < 3 = binom(4,2) - (k-1)` at `k = 4`. So
  Lang-Vojta predicts nothing, consistent with the dense balanced near-tight 4-tuples of gauss Thm 4.1.
* `P^2` and Huang's `Y_4` (`alpha^ess = tau`): nothing in the lane would exclude Huang's Pell points. Theorem A-type
  inequalities are Liouville and hold on them. Theorem C-type conclusions concern bounded `K` only, and Huang's
  points with `d^2 H <= K` lie on finitely many `C_{a,b}`, as the author says.

**Pass.** The lane proves no false `k = 4` analogue, because it proves nothing beyond `tau`.

**T3 (split test).** The lane makes no statement on `X` that could imply split Question B or Erdős #887. The only `X`
statement is the negative transfer, which excludes nothing. So the 2026 five-divisor curve (`alpha = 15/4`) is
trivially not excluded. On the dP2 side, a split surface has `Q`-lines and `Q`-rational tangent-plane cubics with
`alpha = 3/2`, and Theorem A must be recomputed there, as the author notes. **Pass (vacuous).**

**T4 (criticality).** Each of the record's caps, against each mechanism:

| record cap | where it applies in the lane |
|---|---|
| section/`tau` cap | Theorem A **is** it. At `m = 1` it is exactly `tau_P = 3/2`. |
| Lemma K1 (forced divisibility scales like the gcd) | Prop A1: the fibre cost `A'^{-mu}` is the change of units to a global section, and forced divisibility at `P` cancels `m_+` exactly (`x <= 2/3 - nu_+`). |
| multi-point = one-point | Forced divisibility at the fixed second point `P'` is a one-point (two-fixed-point) Liouville bound on `Bl_{P,P'}`, not a Thue-Dyson argument. The author's tensor-product remark on Dyson is right but unproved in detail. |
| class L | not applicable (no torus). |
| fixed-`S` Subspace | Theorem C is conditional on Lang-Vojta with `Sigma_U` depending on `U`, and `Z_K` depends on `K`. Its unconditional Subspace route collapses to McKinnon-Roth `beta = 17/18`, with the threshold corrected to `3/2` in §5.3. |

**Verdict: every unconditional statement in the lane is inside the T4 class.** The author says so for Theorem A and
Prop A1.

**T5 (calibration).** Points that exist:
* `D`-points: recomputed, 14 convergents; they satisfy the `Gamma_3` inequality, with `x + nu_+` from 0.664 to
  0.733, below `3/4`.
* `E_0` points: hand check; with `f = e + 2j`, `nu_+ = delta/2` and `nu_- = 1/2`, and the fibre term is an equality.
* `C_0`: tight.
* The record's `X`-families (round-24 real family, split Pell curves, 8- and 10-point families) are not constrained by
  anything in the lane.

**Pass.**

## 8. Literature check

These are WebSearch snippet-level checks (October 2026); I did not read any full text.

| write-up claim | what WebSearch showed | status |
|---|---|---|
| Broustet 2006, arXiv math/0601225: `eps(-K) = 4/3` at a general point of a dP2 | paper found (arXiv, HAL); main theorem as quoted in arXiv 1806.05087: general-point value `1, 4/3, 3/2, 2, 3` for degrees `1, 2, 3, 4-8, 9` | confirmed (secondary quote) |
| Castañeda thesis, Waterloo 2019 (UWSpace 10012/14859), advisor McKinnon; per LMS, all split smooth rational surfaces of Picard rank `<= 5` | thesis record found; abstract: Hirzebruch-surface fibre method; the rank-5 statement is LMS's | confirmed |
| Huang 1809.02001: earlier computations "mostly consider rational surfaces, often (weak) del Pezzo surfaces of degree >= 3" | the same sentence appears in my search snippet | confirmed |
| no source computes `alpha` or `alpha^ess` on dP2/dP1 | my search returned only Huang, McKinnon 2007, MS21, LMS; nothing on degree 2 | consistent (absence at snippet level) |
| McKinnon-Roth, Eur. J. Math. 2 (2016), cubic surfaces | found (arXiv 1306.2977, Springer). The value `3/2` for cuspidal or real-nodal tangent-plane sections is MR-L Thm 4.5, read in full text in `papers_check/read_mr.md` §3.4 | confirmed; **not cited by the write-up as the source of Prop A3's mechanism** |
| Silverman specialization; Shioda height `2 + 2 P.O` | found (secondary statements in arXiv 0811.3109, 2408.02419, 2601.17182) | confirmed |
| Masser `|x| << H^3`; David / Hirata-Kohno elliptic-logarithm bounds | not checked by me; the author labels them memory/snippet | unverified, labelled |

There is no literature misstatement that matters. The omission that matters for novelty is MR-L Thm 4.5 (§3.3).

**Other minor defects**
* §6.3, "4-sections of genus `<= 9` ... (finitely many rational points each, Faltings)". An upper bound on the genus
  does not license Faltings, which needs genus `>= 2`. My rough Riemann-Hurwitz count gives genus 3 for the
  2-division curve of `sigma'`: 10 `I_1` fibres contribute 1 each, the type II fibre contributes 2, so
  `2g - 2 = -8 + 12`. This is unverified, and irreducibility is not checked. The section is labelled "proves
  nothing", so nothing depends on it.
* Write-up §6.1, Ru-Vojta threshold `beta > 1` (should be `> 3/2` for `L = pi^*(-K)`): see §5.3 of this report.
* Write-up §8.1: a sketch stated as fact (§6.4 of this report).
* Write-up §4.5: "(UE)" is a restatement (§3.4 of this report).

## 9. Novelty decision and corrected statement

**Decision: repackaged known barrier.** Mapped against the record's barrier list:

| write-up item | barrier it repackages |
|---|---|
| Theorem A, Prop A1 | the `tau`/section cap with forced divisibility and Lemma K1's change of units (expert note §5 items 3 and 7; T4). At `m = 1`, Theorem A is literally `tau_P = 3/2`. |
| Prop A3 | McKinnon-Roth's tangent-plane-cusp mechanism (MR-L Thm 4.5 via MR-I Thm 2.16), on a dP2 instead of a cubic surface. As a barrier it is the dimension-2 twin of the record's calibration curves at `15/4` (T5 items). |
| Theorem C | a consequence of Vojta's conjecture, weaker than ERR Prop 0.1 / LMS Prop 1.4. The record already classifies predictions beyond Liouville as "only conjectural: Vojta". |
| §8.2 transfer | correct, and true for all `k` by the line argument of §6.2. It restates Lemma K1's "only Vojta's `D = E` term and the heuristic count see the reduction". |
| §5 sieves and descent | the record's items 2-3, as the author says. |

What is genuinely new, but small and not a mechanism:
* a certified explicit dP2 `(S, P)` with `Pic(S) = ZK`, `tau = 3/2`, real-cubic branch field, and no `Q`-curve of degree
  `<= 6` with `alpha < 2`;
* the exact Galois-invariant effective cone of `Bl_{P,P'} S` (generated by `F` and the `Gamma_k = k sigma'`);
* `tau_e(2(-K) - B_{1,4}) = 5` (now proved for all multiples).

**Corrected statement.**
> Let `S : w^2 = f` (the write-up's `f`), `P = (0:0:1:1)`. Then `S` is smooth, `P` is on no line and off the
> ramification curve, and `Pic(S) = ZK` (Frobenius at 3 acts on `NS` as a Coxeter element of `W(E_7)`). Also
> `eps_P(-K) = 4/3`, `tau_P(-K) = 3/2`, `beta_P(-K) = 17/18`, `alpha_{P,D} = 2`, and no curve over `Q` of
> `-K`-degree `<= 6` through `P` has `alpha < 2`. Hence `3/2 <= alpha^ess_P(-K)`, Vojta predicts `>= 2`, and no
> improvement of `3/2` is known for any Zariski-dense sequence.
>
> (A) *[proved; T4 class]* The `Q`-effective classes of `Bl_{P,P'} S` (`P' = ` Geiser image) lie in the cone spanned
> by `F` and `Gamma_k = k sigma'`, `k in Z`. Hence Liouville with forced divisibility at `P` and `P'`, from any
> section of any multiple of `-K`, gives exactly `x <= x_sec(nu_+, nu_-)`. This equals `1/tau = 2/3` when
> `gcd(a,b) = 1`, and is `< 2/3` only if `nu_+ > 0` or `nu_- > 2/9`. Whether the latter set contains infinitely many
> window points is open; heuristically it is no richer than the generic window.
>
> (A3) *[proved; known mechanism, MR-L Thm 4.5]* On `S_0`, `alpha_{P_0}(-K) = tau = 3/2` is attained on a cuspidal
> tangent-plane cubic `C_0`, whose points are fibre-locally generic (independent of `P' - P`). So no bound uniform in
> fibre data alone exceeds `3/2`.
>
> (C) *[conditional on Lang-Vojta for `S \ D`]* For each `K`, `{X^{3/2} H <= K}` lies in a proper closed `Z_K`. This is
> weaker than `alpha^ess > 3/2`, and weaker than what Vojta's conjecture already gives.
>
> (Transfer) *[proved for every `k >= 3`]* On `Bl_e Lbar_k` the log canonical class of `(sum Zt_ab + E)` has degree
> `-1` on lines through `e`. So integrality with respect to the diagonals yields no Lang-Vojta prediction at `e`, and
> at `k = 5`, `tau_e(2(-K) - B_{1,4}) = 5` exactly.
>
> Nothing here proves, or reduces to anything weaker than Vojta, any case of `alpha^ess > tau` at a rational point
> with `tau < dim`.

## 10. Scripts

All in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/kill_dp2_checks/`.
They are my own code (python3, sympy 1.14, numpy, python-flint 0.9); none imports the author's scripts.

| script | output | what it checks |
|---|---|---|
| `k_geom.py` | `k_geom.out` | smoothness of `S` over `Q`, `F_3`, `F_5`; `T`, `F_3`, `F_4'`; resultant; cubic field and discriminant; tangent-line discriminant |
| `k_frob.py` | `k_frob.out` | `#S(F_{3^n})` for `n <= 8` by exhaustive count over `P^2(F_q)` (own `F_q` arithmetic); traces; char poly by Newton identities |
| `k_deg6.py` | `k_deg6.out` | the net `|-3K - 4P|` and its square-tangent-cone members (irreducible quartic eliminant) |
| `k_gamma3.py` | `k_gamma3.out` | dimensions of `|-dK - aP - bP'|` for 10 triples (`Gamma_2`, `Gamma_3`, `Gamma_{-2}` unique); D-points vs the `Gamma_3` bound |
| `k_s0.py` | `k_s0.out` | `S_0` smooth, `Res = 1`, `F_3` factors (rational branch of `D_0`), points `Q_b` with `u = 1`, `X^{3/2}H -> 1` |
| `k_usearch.py` | `k_usearch.out` | all points with `max(|a|,|b|) <= 400`, `0 < |8u| <= 24`: four |
| `k_toric.py` | `k_toric.out` (`t <= 2`), `k_toric3.out` (`t <= 3`) | maximal orders at `e` on `Lbar_5` and `Lbar_4` for `-K`, `L'`, `K + sum Z_ab` |
| `notes.md` | -- | working notes and hand derivations |

Housekeeping note: one script (`k_usearch.py`) was first written by mistake to the repository root
`/home/user/Jarnik/`. A backgrounded `cd` caused this. I moved it out at once, and `git status` there is clean.
