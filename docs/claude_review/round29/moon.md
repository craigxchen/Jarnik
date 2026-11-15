# Round 29, lane moon: moonshots against F1 (fake circles) and F2 (criticality)

*Status labels as in the record: proved, exact, evidence, heuristic, open. Everything here is unrefereed.*

## 0. Verdict

**No growth improvement and no uniform bound.** $M\le(\frac23+\varepsilon)\log R/\log\log R$ stands. Eighteen ideas from
far fields were screened (§§2–3). Fifteen fail the fake filter F1 or the criticality filter F2, each for a stated reason
(several by an exact one-line computation: Erdős–Turán is vacuous for every $N$; parallelogram completions are pair data;
the whole circle is a product set, so every symmetric function of it is fake data plus integrality). A machine look at
the five record clusters (exact) finds that their reduced chords carry exactly the pigeonhole-forced collisions at
$\ell=3$ and that every prime $\ell>4k$ divides at most one chord: no hidden identity. The best idea, the
**function-field analogue**, was developed (§4) and gives a calibration, not a route:
* over $\mathbb F_q(t)$ the analogue of the growth theorem holds with the same $\log R/\log\log R$ shape (Prop 4.3,
  proved under a coprimality hypothesis (H));
* over $\mathbb R(t)$ there are no finite residue fields, every local inequality admits Hadamard profiles with
  $\log R\asymp M$ (Obs 4.4), and the problem becomes a Littlewood–Offord question on the real moment curve (Lemma 4.5);
  its first case is exactly Sperner, $M=6$ at $D=4$ (Prop 4.6, proved), matching the record's polynomial templates;
* composition $t\mapsto h(t)$ is the function-field form of criticality (Lemma 4.1(iv)).
So the $\log M$ in Theorem 1.1 is entirely a finite-residue-field effect, and a function-field proof of the uniform
bound would need a derivation (Wronskian) argument, whose number-field shadow is abc-type. The uploaded Goldbach
paper's device (GRH plus a spectral operator) was screened as idea 14 and fails F1: the spectral side here is a
function of the angles, and GRH controls all primes, not the $\omega$ primes of a given $N$.

## 1. The two filters, made operational

Notation as in `expert_note/expert_note.md`: a normalised cluster is $z_x=\varepsilon_x\prod_p\pi_p^{a_{xp}}\bar\pi_p^{e_p-a_{xp}}$,
$x=1,\dots,M$, on $|z|^2=N$, in an arc of length $C\sqrt R$; $W=\log N$; $\phi_p=\arg\pi_p$; $c_{xy}$ the reduced chords.

**F1 (fakes).** A screen *fails F1* if the input it uses is (a) a function of $(p,e_p,a_{xp},\phi_p)$ in which $\phi_p$
enters only as a real number, or (b) one of the classes listed in expert note §5 items 1–2, 4–6 (pairs, cuts, gcds,
Liouville gaps, rational characters, residues at moduli $\le M^A$, averaged analytic inputs, single linear forms,
$|S|$-dependent $S$-unit bounds, abc on single relations). An idea *passes F1* only if it names a property of genuine
circles that some recorded fake at $W=O(M\log M)$ violates.

**F2 (criticality).** A screen *fails F2* if its only arithmetic input is one auxiliary polynomial (or one finite family)
with Liouville and toric forced divisibility, or a multi-point argument among subtuples of one cluster (one-point toric
Liouville on $T_M$; round 25), or any statement invariant under $z\mapsto gz$ in the sense of Lemma K1. A mechanism that
would also prove the false 4-point statement $(L_4^\circ)$ (gauss Thm 4.1) fails F2 automatically.

A useful consequence used repeatedly below:

> **Observation 1.1 (the whole circle is a product set; elementary).** Let $\mathcal E_R$ be all lattice points of
> $|z|^2=N$ (after removing the gcd-trivial factors $2$ and inert primes, which only multiply all points by a fixed
> integer). Then $\{z/|z|\}=\mu_4\cdot\prod_p\{e^{i(2a-e_p)\phi_p}:0\le a\le e_p\}$ (Minkowski product on the circle). Every
> symmetric function of $\mathcal E_R$ (power sums $\sum z^{m}$, the circle polynomial $F_N$, its discriminant, critical
> points, Coulomb energy, angular discrepancy) is therefore a function of $(R,(e_p,\phi_p)_p)$, i.e. of the data a fake
> keeps, except for integrality, which is the only extra information.

*Proof.* Unique factorisation in $\mathbb Z[i]$: every point is $\varepsilon\prod\pi_p^{a}\bar\pi_p^{e_p-a}$ with
$0\le a\le e_p$ independently. $\square$

So an idea that uses "the points outside the cluster" passes F1 only through integrality of some expression, and
then it is an auxiliary-polynomial argument and must face F2.

## 2. Eighteen ideas, one paragraph each

1. **Erdős–Turán / Lee–Yang rigidity of the circle polynomial.** $F_N(X)=\prod_{z\in\mathcal E_R}(X-z)\in\mathbb Z[i][X]$ is
self-inversive with all zeros on $|X|=R$. Erdős–Turán bounds the angular discrepancy of the zeros by
$c\sqrt{\log(L(F)/\sqrt{|a_0a_n|})/n}$; a cluster is a discrepancy spike of size $M/n$ in an interval of length $CR^{-1/2}$.
2. **Critical points and Coulomb energy.** By Walsh/Gauss–Lucas a cluster of $M$ zeros in a disc of radius $r$ forces
$M-1$ critical points of $F_N$ nearby; the Coulomb energy $-\sum\log|z-z'|=-\log|\mathrm{disc}F_N|\le0$ since the
discriminant is a nonzero integer. Do these see clusters?
3. **Power sums and Turán's power-sum method.** $s_m=\sum_{z\in\mathcal E_R}z^{4m}$ is an integer linear recurrence (a
Hecke-eigenform coefficient for the character $(z/|z|)^{4m}$); a cluster makes $M$ characteristic roots nearly coincide,
which Turán-type lower bounds for $\max_{m\in[K,K+n]}|s_m|$ might detect.
4. **The points outside the cluster.** Use all $r(N)$ points as auxiliary data: polynomial identities between cluster
points and the rest of the circle, e.g. integrality of $\prod_{z\notin I}(z_x-z)$ or resultants against $F_N$.
5. **Neighbouring circles.** Parallelogram completions $z_x+z_y-z_w$, Pick's theorem for the lattice polygon of the
cluster, lattice points of $N\pm k$ inside the circular segment, and abc on $N,k,N+k$.
6. **Non-archimedean twins.** The arc condition says $z_x/z_y\to1$ at $\infty$. Twins: (a) $z_x\equiv z_y\bmod\lambda$ for every
prime $\lambda\mid c_{xy}$ (the chord primes, unbounded moduli, outside the residue-design scope); (b) the $q$-adic arc
problem at $q\mid N$.
7. **Homogeneous dynamics.** Rotate so that $z_1$ is vertical and apply $\mathrm{diag}(R^{-1/4},R^{1/4})$: the cluster becomes
$M$ points of a unimodular lattice on a parabola of curvature $R^{-1/4}$, i.e. after rescaling, points of a lattice of
covolume $R^{-1/2}$ on a fixed parabola in a unit box; use equidistribution of geodesic translates.
8. **Entropy compression / Kolmogorov complexity.** A large cluster should give a short description of $N$, so few $N$ can
carry one; Moser–Tardos-style resampling for the profile.
9. **Codes: Plotkin stability.** Profiles are binary codes (cuts); near-extremal profiles are near-Hadamard by Plotkin
stability; then use one arithmetic fact about Hadamard designs of primes.
10. **Reciprocity / Rédei matrices.** The quartic residue symbols $[\pi_p/\pi_q]_4$ of the primes of $N$ are residues at
moduli that are primes of $N$ itself (far above $M^A$); quartic reciprocity ties the matrix together.
11. **Orbit parametrisations.** Parametrise $k$ points on the conic $x^2+y^2=N$ by a Bhargava-type space and use its
invariants.
12. **Tropical / matroidal structure** of the profile (Bergman fans, valuated matroids of the $2\times M$ point matrix).
13. **Finite free probability / sumsets.** The angle set is a Minkowski sum of arithmetic progressions
$\{(2a-e_p)\phi_p\}$, i.e. $F_N$ is a finite multiplicative convolution of local polynomials; use free-convolution root
rigidity or additive combinatorics of generalised arithmetic progressions in short intervals.
14. **L-functions of the individual $N$; GRH; a Goldbach-type spectral operator.** Following the uploaded Li–Liu paper
(GRH plus a spectral hypothesis for an operator built from characters mod $N$): count points in an arc by a
Beurling–Selberg majorant, $\#\le\sum_{|m|\le K}\hat B(m)e^{-4im\theta}\lambda_m(N)$, $K\asymp R^{1/2}/C$, and try to bound the
spectral side under GRH for Hecke $L$-functions.
15. **Sum–product.** The cluster is additively thin (all chords in a disc of radius $C\sqrt R$) and multiplicatively a box
(a subset of a rank-$\omega$ multiplicative progression); sum–product phenomena forbid both.
16. **Function fields.** Pose the problem over $\mathbb C(t)$ and over $\mathbb F_q(t)$ ($q\equiv3\bmod4$), where derivations
(Mason–Stothers, Wronskians), composition $t\mapsto h(t)$ and Frobenius are available, and see which structure
decides it there.
17. **Machine search for identities** among invariants (reduced chords, their prime factors, cut data) of the record
clusters.
18. **Approximating fakes by genuine circles.** Choose genuine primes with angles in prescribed sectors (Kubilius; under
GRH sectors of width $p^{-1/2+\varepsilon}$) to turn a fake into a genuine circle: a constructive (refutation) direction.

## 3. Screening against F1 and F2

Each screen states the outcome and the reason. "Fails F1(a)/(b)" refers to §1. Computations are in `moon_checks/`.

**1 (Erdős–Turán, Lee–Yang): vacuous; fails F1(a).** Normalise $G(X)=R^{-n}F_N(RX)$, $n=r(N)$: monic, zeros on $|X|=1$,
$|a_0|=1$, so $L(G)\le2^n$ and the Erdős–Turán discrepancy bound is $\le16\sqrt{\log 2}>1$ for **every** $N$: vacuous
(exact). Self-inversiveness (Lee–Yang) holds for any point set on a circle, fakes included.

**2 (critical points, Coulomb energy): fails F1(a) and F2.** Walsh's two-circle theorem is a statement about arbitrary
complex zeros. $|\mathrm{disc}\,F_N|\ge1$ is the product over all pairs of the pair inequalities $|c_{zz'}|\ge1$ inside the cut
identity for the full set, so it is implied by pair data; it is one integrality statement and gives no more than the
pair inequalities it multiplies.

**3 (power sums, Turán): fails F1(a).** By Observation 1.1, $s_m=R^{4m}\prod_p\sin(4m(e_p+1)\phi_p)/\sin(4m\phi_p)$ (with the
obvious reading at zeros of the sine), a function of fake data; Turán-type bounds hold for arbitrary complex roots.

**4 (points outside the cluster): fails F2 (screen, not a theorem).** Every lattice point of the circle is a monomial
$\varepsilon\prod\pi_p^{a}\bar\pi_p^{e_p-a}$, so a polynomial in all points is a Laurent polynomial in the $\pi_p,\bar\pi_p$ evaluated at
one point of the torus $R^1_{\mathbb Q(i)/\mathbb Q}\mathbb G_m^{\omega}$; the outside points add no new variables and no
archimedean smallness (they are at distance $\gg C\sqrt R$ from the arc, generically $\asymp R$). The argument is therefore an
auxiliary-polynomial/Liouville argument at one point with toric forced divisibility, the class capped by Lemma K1 and
the $\tau$-caps. I did not prove that every such argument reduces to the round-25 class on $T_M$; the screen only shows it
is of the F2 type.

**5 (neighbouring circles): fails F1(b); the geometric version is weak.** $|z_x+z_y-z_w|^2=N+n_{xw}+n_{yw}-n_{xy}$ with
$n_{xy}=|z_x-z_y|^2$ (law of cosines; exact), so the neighbouring circles reached by completions are pair data, at
$k\asymp C^2R$, not small $k$. Pick plus Andrews (a convex lattice $M$-gon has area $\gg M^3$) against the segment area
$\approx C^3\sqrt R/12$ gives only $M\ll CR^{1/6}$. abc for $N,k,N+k$ with $k\le2C^2\sqrt N$ is satisfied trivially ($\mathrm{rad}(k)\le N^{1/2}$).

**6 (non-archimedean twins): fails F1(b).** (a) $\lambda\mid c_{xy}\iff z_x\equiv z_y\bmod\lambda\gcd(z_x,z_y)$; the total mass
$\sum\log N(c_{xy})$ is already the last factor of the cut identity, and collisions are *forced* only at $\ell\lesssim M$
(pigeonhole), which is the proof of Theorem 1.1. Idea 17 below shows that actual record clusters have exactly the forced
collisions and nothing else. (b) The $q$-adic arc data at $q\mid N$ are the cut profile.

**7 (homogeneous dynamics): fails F1(b).** Equidistribution of translates gives counts averaged over the arc position
or over $N$, with error terms; exact incidence with a conic has measure zero and is invisible to it. Averaged inputs are
capped by grafts (expert note §5 item 2).

**8 (entropy, Kolmogorov): fails F1(b).** Description-length arguments bound how many $N\le X$ carry a cluster ("almost
all"); grafts cap them. Moser–Tardos proves existence, not non-existence.

**9 (Plotkin stability): fails F1.** The recorded fakes *are* Paley (skew-Hadamard) profiles (`round8/lw.md` Thm 6.1), so
stability pushes exactly into the fake region; the "one arithmetic fact" would have to be non-local, which is the
original problem.

**10 (reciprocity, Rédei matrices): no mechanism; heuristic F1 failure.** The symbols $[\pi_p/\pi_q]_4$ are residues at moduli
outside the residue-design scope, so they are not formally capped. But I found no implication from the arc condition to
the symbols: the arc constrains $z_x/z_y$ at $\infty$, the symbols constrain $\pi_p\bmod\pi_q$, and reduced chords are
not forced to be divisible by primes of $N$. Joint equidistribution of Gaussian primes in sectors and in ray classes
(Hecke characters times ray-class characters) predicts independence (**heuristic**). Outcome: no input to test.

**11 (orbit parametrisations), 12 (tropical), 13 (free probability, GAPs): fail F1(a)/(b).** A $k$-tuple on the conic is
parametrised by its Gaussian factorisation data, which fakes keep; the tropicalisation of the cluster point in $T_M$ is the
cut profile (round 18's chain/caterpillar conditions are its tropical content); the angle set and all its subsets are
functions of $\phi$, and free probability describes expected characteristic polynomials, not individual subsets.

**14 (L-functions of $N$, GRH, Goldbach-type operator): fails F1(a).** $\lambda_m(N)$ is a finite product of local factors in
$\phi_p$. GRH for Hecke $L$-functions controls sums over all primes, not the $\omega$ primes of a given $N$, which are
arbitrary. The Li–Liu device builds an operator whose quadratic form *is* the Goldbach error; here the analogous exact
identity (the Beurling–Selberg count) already has a spectral side that is a function of the angles.

**15 (sum–product): vacuous.** Sum–product theorems need $|A|\ge(\text{scale})^{\varepsilon}$; here $M\ll\log R$.

**17 (machine search on record clusters): negative evidence** (`cluster_invariants.py`, exact Gaussian arithmetic,
output `cluster_invariants_output.txt`). For the five record clusters ($R^2=537606725$, $32988780417125$,
$348984761607530165$ with 8 points, $235370525$ with 7, $683617125425$ with 5), after gcd normalisation the odd-prime
collisions of the reduced chords are exactly the pigeonhole minimum at $\ell=3$ (4, 4, 4, 3, 1 pairs), plus at most two
pairs at $\ell\in\{5,7,31\}$, and **every prime $\ell>4k$ divides at most one reduced chord**. No hidden identity: at the
level of chord arithmetic the record clusters look like the fake model (forced collisions plus one-off large primes).

**18 (genuine primes approximating fakes): heuristic barrier at 5 points.** Steering each $\phi_p$ into a sector of width
$p^{-1/2+\varepsilon}$ (GRH) is far too coarse: the cluster needs $M-1$ independent combinations accurate to
$N^{-1/4}=e^{-W/4}$. Joint choice of primes gives about $e^{W}$ choices against probability $e^{-W(M-1)/4}$, so it succeeds
only for $M\le5$: the record's heuristic count again (**heuristic**).

**16 (function fields): passes both filters in a degenerate sense, developed in §4.** Over $\mathbb C(t)$ there are no
fakes distinct from genuine circles (the analogue of $\phi_p$ is an arbitrary complex root, but all relations are exact),
and the analogue of Lemma K1 is the composition symmetry $t\mapsto h(t)$. The point of developing it is to learn which
structure decides the function-field problem.

| idea | F1 | F2 | outcome |
|---|---|---|---|
| 1 Erdős–Turán/Lee–Yang | fail (a) | — | vacuous (exact) |
| 2 critical points/energy | fail (a) | fail | implied by pair data |
| 3 power sums/Turán | fail (a) | — | function of angles |
| 4 outside points | (integrality only) | fail | Liouville type at one point (screen) |
| 5 neighbouring circles | fail (b) | — | pair data; Pick gives $R^{1/6}$ |
| 6 non-arch. twins | fail (b) | — | = cut identity / proof of Thm 1.1 |
| 7 dynamics | fail (b) | — | averaged |
| 8 entropy | fail (b) | — | averaged |
| 9 Plotkin stability | fail | — | lands on Paley fakes |
| 10 Rédei/reciprocity | not formally | — | no mechanism found |
| 11–13 | fail | — | profile/angle data |
| 14 L-functions/GRH | fail (a) | — | angles; GRH is about all primes |
| 15 sum–product | — | — | vacuous at $M\ll\log R$ |
| 16 function fields | see §4 | see §4 | developed |
| 17 data search | — | — | negative evidence |
| 18 approximating fakes | — | — | heuristic cap 5 |

## 4. Development: the circle problem over function fields

**4.1 Setting.** Let $k$ be a field of characteristic $\ne2$ in which $-1$ is not a square ($k=\mathbb R$, or
$k=\mathbb F_q$ with $q\equiv3\bmod4$), $A=k[t]$, $B=k(i)[t]$, and $\bar{\ }$ the conjugation of $k(i)/k$ on coefficients
(for $k=\mathbb F_q$ it is the $q$-Frobenius on coefficients). Let $N\in A$ be monic of degree $2D$. A *point* is a monic
$z\in B$ with $z\bar z=N$; $\deg z=D$. Use the absolute value $|f|=e^{\deg f}$ ($q^{\deg f}$ over $\mathbb F_q$), so
$R=|z|=e^{D}$, $\sqrt R=e^{D/2}$, and an arc of length $C\sqrt R$ at critical scale becomes
$$\deg(z_x-z_y)\le s:=\lfloor D/2+c\rfloor\qquad(\text{all }x,y),\qquad c\leftrightarrow\log C .$$
$M_k(c)$ denotes the supremum over $N$ of the size of such clusters. (Monic is the unit-class normalisation: points with
different leading coefficients are at distance $\deg=D$.)

> **Lemma 4.1 (dictionary; proved).** (i) *Arcs in any position are balls*: by the ultrametric inequality a set has all
> pairwise $\deg(z_x-z_y)\le s$ iff it lies in one ball $\{\deg(z-z_1)\le s\}$, iff its points share the coefficients of
> degrees $>s$. (ii) *Unique factorisation*: for a gcd-normalised cluster every prime of $A$ dividing $N$ splits,
> $P=\pi\bar\pi$, and $z_x=\prod_\pi\pi^{a_{x\pi}}\bar\pi^{e_\pi-a_{x\pi}}$. (iii) *FF Liouville*: with $g=\gcd(z_x,z_y)$ and
> $c_{xy}=(z_x-z_y)/g\in B\setminus\{0\}$, $\deg(z_x-z_y)=\deg g+\deg c_{xy}$, $\deg g=\sum_\pi\deg\pi\,(e_\pi-|a_{x\pi}-a_{y\pi}|)$.
> (iv) *Composition*: if $h\in k[t]$ has degree $\kappa\ge1$, $z\mapsto z\circ h$ maps a cluster for $(N,D,c)$ to one for
> $(N\circ h,\kappa D,\kappa c)$, of the same size.

*Proof.* (i) ultrametric. (ii) $B$ is a PID; $B/A$ is the constant-field extension, unramified as $\mathrm{char}\,k\ne2$; an
inert prime $Q$ divides every point to the same power and is removed by normalisation. (iii) $c_{xy}\ne0$ because the
points are distinct; degrees add. (iv) degrees multiply by $\kappa$ and conjugation commutes with composition by
$h\in k[t]$. $\square$

Item (iv) is the function-field face of F2: every degree inequality that holds must be homogeneous under
$(D,c)\mapsto(\kappa D,\kappa c)$. In particular $M_k(c)$ is the same for all $c<0$ (a cluster at $c=-1$ gives clusters at
every $c\le-1$), so the analogue of "$C\to0$" is "$c\to-\infty$", and $c=0$ is a fixed point.

> **Proposition 4.2 (FF cut bound; proved).** For a gcd-normalised cluster of $M$ points, with
> $S_{\pi\tau}=\#\{x:a_{x\pi}<\tau\}$,
> $$\sum_{x<y}\deg c_{xy}\ \le\ \binom M2c+\frac{MD}4-\sum_\pi\deg\pi\sum_{\tau=1}^{e_\pi}\frac{(M-2S_{\pi\tau})^2}4 .$$

*Proof.* By Lemma 4.1(iii) and the arc condition, $\sum_{x<y}\deg c_{xy}\le\binom M2 s-\sum_{x<y}\deg g_{xy}$, and
$\sum_{x<y}|a_{x\pi}-a_{y\pi}|=\sum_\tau S_{\pi\tau}(M-S_{\pi\tau})=\sum_\tau(M^2-(M-2S_{\pi\tau})^2)/4$. Insert $\sum_\pi e_\pi\deg\pi=D$
and $s\le D/2+c$. $\square$

*Exact check* (`ff_check.py`): over $\mathbb F_3$, all 23 238 ball clusters of size $\ge2$ at $c=0$ of all squarefree split $N$
with $\deg N\le12$ satisfy (iii) with exact division and Proposition 4.2 (0 violations).

> **Proposition 4.3 (FF analogue of the growth theorem over $\mathbb F_q$; proved under (H)).** Let $k=\mathbb F_q$,
> $q\equiv3\bmod4$, $M\ge q^2$, $d_0=\lfloor\log_qM\rfloor-1$, and assume
> (H) $N$ has no prime factor of degree $\le d_0$. Then
> $$\tfrac{M^2}2\big(\log_qM-K_q\big)\ \le\ \sum_{x<y}\deg c_{xy}\ \le\ \binom M2c+\frac{MD}4,\qquad K_q=2+\tfrac{3}{\sqrt q-1}+\tfrac1{q-1},$$
> hence $M\le\dfrac{D/2}{\log_qM-K_q-c}$, i.e. $M\le(\tfrac12+o(1))\,\log_qR/\log_q\log_qR$.

*Proof.* Upper bound: Proposition 4.2. Lower bound: let $\ell\in A$ be monic irreducible of degree $d\le d_0$; by (H)
$\ell\nmid N$, so $\ell$ is coprime to every $g_{xy}$. If $d$ is odd, $\ell$ is inert, $B/\ell=\mathbb F_{q^{2d}}$, conjugation
induces $u\mapsto u^{q^d}$, and $z^{1+q^d}\equiv N\not\equiv0$: the residues of the points lie in one coset of the norm-one
group, $h=q^d+1$ classes. If $d$ is even, $\ell=\lambda\bar\lambda$, $B/\lambda=\mathbb F_{q^d}$, $h=q^d-1$ classes, and
$z_x\equiv z_y\ (\lambda)$ implies $\bar z_x\equiv\bar z_y\ (\bar\lambda)$, hence $z_x=N/\bar z_x\equiv z_y\ (\bar\lambda)$. In both cases a
collision mod $\lambda$ gives $\ell\mid z_x-z_y$, hence $\ell\mid c_{xy}$, contributing $d$ to $\deg c_{xy}$. Pigeonhole gives at
least $\frac M2(\frac Mh-1)\ge\frac{M^2}{2(q^d+1)}-\frac M2$ colliding pairs. With $d\,I_q(d)\ge q^d-2q^{d/2}$ monic
irreducibles of degree $d$ (from $\sum_{e\mid d}eI_q(e)=q^d$) and $d\,I_q(d)\le q^d$,
$$\sum_{x<y}\deg c_{xy}\ge\sum_{d\le d_0}\Big[\frac{(q^d-2q^{d/2})M^2}{2(q^d+1)}-\frac{q^dM}2\Big]
\ge\frac{M^2}2\Big(d_0-\sum_{d\ge1}\frac{2q^{d/2}+1}{q^d+1}\Big)-\frac M2\cdot\frac{q^{d_0+1}}{q-1},$$
and $q^{d_0+1}\le M$, $\sum_d(2q^{d/2}+1)/(q^d+1)\le 2/(\sqrt q-1)+1/(q-1)$. So the sum is
$\ge\frac{M^2}2(\log_qM-2-\frac2{\sqrt q-1}-\frac2{q-1})$, which is the claim with the stated (slightly generous) $K_q$.
Rearranging the two bounds gives $M(\log_qM-K_q-c)\le D/2+O(1)$. $\square$

(H) is a simplification: dropping it needs the cut term of Proposition 4.2 to pay for primes of $N$ of small degree,
as in `two_thirds/proof.md`; I did not do this. Proposition 4.3 is only an illustration that the number-field mechanism
(small residue fields) transfers to $\mathbb F_q(t)$ with the same $\log R/\log\log R$ shape.

**4.2 Over $\mathbb R(t)$ the local mechanism disappears.** Over $k=\mathbb R$ every residue field of $A$ is $\mathbb R$ or
$\mathbb C$, so no pigeonhole collision is forced, and Proposition 4.2 is the only local constraint. It is satisfiable with
$M$ linear in $D$:

> **Observation 4.4 (FF "fakes" at $W\asymp M$; proved, combinatorial).** Let $D$ be a Hadamard order, $N$ squarefree with
> $D$ conjugate root pairs, and take the $M=2D$ sign vectors $\pm$(rows of a Hadamard matrix). All pairwise weighted
> distances are $D/2$ or $D$, so Lemma 4.1(iii) and Proposition 4.2 at $c=0$ are satisfied with $\deg c_{xy}=0$ for the
> distance-$D/2$ pairs and $\deg c_{xy}\le D/2$ for the $M/2$ complementary pairs (sum $\le M D/4$, equality budget).

So over $\mathbb R(t)$ the analogue of every *local* input (pairs, cuts, Liouville, collisions) admits growing profiles
with $\log R=D=M/2$, which is even better for the counterfeiter than the number-field fakes ($W\asymp M\log M$): the
$\log M$ of Theorem 1.1 is exactly the contribution of the finite residue fields. Whether such profiles are *realised* is
a purely algebraic question:

> **Lemma 4.5 (realisation over $\mathbb R(t)$; proved).** Let $N=\prod_{k=1}^D(t-a_k)(t-\bar a_k)$, $\mathrm{Im}\,a_k>0$, and
> $z_\varepsilon=\prod_k(t-a_k^{\varepsilon_k})$ ($a^+=a$, $a^-=\bar a$). Then $\deg(z_\varepsilon-z_{\varepsilon'})\le s$ iff
> $\sum_k(\varepsilon_k-\varepsilon'_k)\,\mathrm{Im}(a_k^r)=0$ for $r=1,\dots,m:=D-1-s$. A cluster is a fibre of
> $\varepsilon\mapsto\big(\sum_k\varepsilon_k\,\mathrm{Im}(a_k^r)\big)_{r\le m}$: a Littlewood–Offord problem for the vectors
> $w_k=(\mathrm{Im}\,a_k^r)_{r\le m}$ on the real moment curve.

*Proof.* $z_\varepsilon-z_{\varepsilon'}$ has degree $\le s$ iff the elementary symmetric functions $e_1,\dots,e_m$ of the two
root multisets agree, iff (Newton, characteristic 0) the power sums $p_1,\dots,p_m$ agree. $2\,\mathrm{Re}\,p_r(z_\varepsilon)=p_r(N)$
is the same for every point (the roots of $z$ and $\bar z$ make up those of $N$), so only the imaginary parts matter, and
$\mathrm{Im}\,p_r(z_\varepsilon)=\sum_k\varepsilon_k\mathrm{Im}(a_k^r)$. $\square$

This is `growth/families.md` Lemma P with complex instead of Gaussian-integer roots.

> **Proposition 4.6 (the first case is Sperner; proved).** Over $\mathbb R(t)$, at $D=4$, $c=0$: every gcd-normalised
> cluster has at most $6$ points, and $6$ is attained.

*Proof.* Here $s=2$, $m=1$. For non-squarefree $N=\prod(t-a_k)^{e_k}(t-\bar a_k)^{e_k}$ ($\sum e_k=4$; real roots are excluded
by normalisation since they divide every point), a point is $j\in\prod_k[0,e_k]$ and the condition is that
$\sum_k(2j_k-e_k)\,\mathrm{Im}\,a_k$ is constant. With $\mathrm{Im}\,a_k>0$ this functional is strictly monotone on the product
of chains, so a cluster is an antichain, of size at most the largest rank (de Bruijn–Tengbergen–Kruyswijk): $6,4,3,2,1$
for $(e_k)=(1,1,1,1),(2,1,1),(2,2),(3,1),(4)$. Equality: $\mathrm{Im}\,a_k=1$ for all $k$, distinct real parts, and the six
$\varepsilon$ with two plus signs. $\square$

This is exactly the record's $D=4$ polynomial-template record (6 points, "CG type, all $y=1$", `growth/families.md` §4),
now shown optimal over $\mathbb R$, not only over the searched integer box.

*Larger $D$ (evidence only, inconclusive).* The generic count is $2D-2$ real unknowns (real affine changes of $t$) against
$(M-1)m\approx(M-1)D/2$ equations, i.e. $M\le5$ for large $D$, the same 5 as the record's heuristic count and the dimension
of $T_5$. A numerical search (`ff_complex.py`: random designs with pairwise distance $\ge m+1$, Levenberg–Marquardt/trf
from random starts) found no non-degenerate solution even where the count allows one; random designs are almost never
fibres, so the method is unsuitable and **no conclusion is drawn** from it.

*Data over $\mathbb F_3$* (`ff_scan.py`, exhaustive over all $1.2\cdot10^6$ squarefree split $N$ with $\deg N\le16$; balls by
Lemma 4.1(i); **exact**): the maximal cluster size is $\le3$ at $c=-1$ and $\le4$ at $c=0$ for every degree; at $c=+1$ the
values $6$–$8$ occur only where the ball condition is nearly empty. Nothing grows with the degree in this range.

## 5. What survives, and exactly what structure it uses

**Nothing survives as a mechanism for the goal.** Of the 18 ideas, 15 fail F1 or F2 for a stated reason (§3), one
(Rédei matrices, idea 10) is not formally capped but has no input to test, idea 17 is negative evidence, and idea 16
(function fields) is a calibration, not a route. It does sharpen the diagnosis in three ways.

1. *Where the $\log M$ in Theorem 1.1 comes from.* Over $\mathbb F_q(t)$ the same small-residue-field collisions give the
   same $\log R/\log\log R$ shape (Proposition 4.3, under (H)). Over $\mathbb R(t)$ there are no finite residue fields and
   the local inequalities admit Hadamard profiles with $\log R\asymp M$ (Observation 4.4). So in the function-field world
   the *entire* growth bound is a finite-residue-field effect, and the uniform question is a purely algebraic one about
   Littlewood–Offord fibres of moment-curve vectors (Lemma 4.5). Its first case is exactly Sperner (Proposition 4.6:
   $M=6$ at $D=4$, matching the record's integer templates).
2. *F2 in function fields is composition.* $t\mapsto h(t)$ multiplies $(D,c)$ by $\deg h$ (Lemma 4.1(iv)). Any argument that
   is homogeneous in degrees is blind to the difference between $c=0$ and the criticality class it lives in, which is
   the function-field form of Lemma K1. A function-field proof of $M_{\mathbb R}(0)<\infty$ would have to use a
   non-homogeneous input; the natural one is the derivation (Wronskian, Mason–Stothers), whose number-field analogue
   is abc-type. abc on single relations is vacuous on squarefree profiles (expert note §5 item 6), so even a positive
   function-field theorem would point only to a multi-relation abc, i.e. back to Vojta-type statements.
3. *Which structure of genuine circles fakes lack.* Every idea that passed the formal part of F1 did so only through
   integrality of some expression in the $\pi_p$ (Observation 1.1), and then met F2. The only structures not of this kind
   that I could name are (a) the joint law of angles and residue symbols of the primes of $N$ (idea 10), for which
   no implication from the arc condition was found, and (b) derivations, which exist only in function fields. I have no
   candidate that uses a structure of genuine circles which fakes lack and that is not integrality.

**Open items this lane leaves (precise).**
* Is $M_{\mathbb R}(0)<\infty$ over $\mathbb R(t)$, uniformly in $D$? Equivalently: for vectors
  $w_k=(\mathrm{Im}\,a_k^r)_{r\le D-1-\lfloor D/2\rfloor}$ with distinct $a_k$ in the upper half-plane, is the largest fibre of
  $\varepsilon\mapsto\sum\varepsilon_kw_k$ bounded independently of $D$? The generic count says $5$; Proposition 4.6 gives $6$ at
  $D=4$. A construction with growing fibres would show that the conjecture is not "geometric"; a Wronskian proof would be
  the first complete function-field analogue of the uniform bound. **Open.**
* Remove (H) from Proposition 4.3 using the cut term (routine but not done).

## 6. Files

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round29/`:
* `moon.md` (this write-up), `moon_notes.md` (working notes).
* `moon_checks/cluster_invariants.py`, `cluster_invariants_output.txt`: idea 17 (exact Gaussian arithmetic, sympy factoring).
* `moon_checks/ff_scan.py`: exhaustive $\mathbb F_q(t)$ ball-cluster scan (run: `python3 ff_scan.py 3 16 4`, about 2 min;
  the $q=7$, $\deg N\le12$ run did not finish in 15 min and is not reported).
* `moon_checks/ff_check.py`: exact check of Lemma 4.1(iii), Proposition 4.2 and the forced collisions over $\mathbb F_3$.
* `moon_checks/ff_complex.py`: the inconclusive numerical search of §4.2 (no claim rests on it).

Literature: none newly cited. De Bruijn–Tengbergen–Kruyswijk (the divisor lattice is Sperner) and Erdős's
Littlewood–Offord lemma are quoted from memory (**memory-quote**). The Goldbach paper is the uploaded PDF (Li–Liu, "A
conditional proof of the Goldbach conjecture"), read at abstract level only.
