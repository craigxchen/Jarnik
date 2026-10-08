# A composite-block transfer from the Gaussian moat argument

The source is [*Bounded-Step Walks on Gaussian Primes*, Lemma 4.2, “A randomly signed product avoids a thin rectangle”](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Bounded-Step-Walks-on-Gaussian-Primes-September-26-2026/paper.pdf). Its [LaTeX proof](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Bounded-Step-Walks-on-Gaussian-Primes-September-26-2026/build/geometry.tex) treats one Gaussian prime factor over each of several comparable split rational primes. The lemma below is a **new transfer of that proof**, not a theorem stated in the source: a Gaussian factor may be replaced by a conjugate-primitive composite block. It supplies a sparse-sign conclusion, but does not prove a radius-uniform bound for lattice points in arcs of length `C sqrt(R)`.

## Composite-block signed-product lemma

Let `H_1,...,H_r` be nonunits in `Z[i]`. Assume

1. `gcd_G(H_j,bar(H_j))=1` for every `j`;
2. the rational-prime supports of `N(H_j)` and `N(H_k)` are disjoint when `j != k`;
3. for some fixed `kappa >= 1` and `L>0`, `L <= log N(H_j) <= kappa L` for every `j`.

The first condition rules out inert and ramified factors and forbids both orientations over one split rational prime within a block. Write `n_j=N(H_j)`, `P=prod_j n_j`, and `B_sigma=prod_j K_{j,sigma_j}`, where `K_{j,+}=H_j` and `K_{j,-}=bar(H_j)`. Let `Q` be any origin-centered rectangle, in any orientation, with half-lengths `R_1 >= W_1 >= 1`. Put `A=R_1 W_1`. There are positive constants `c_kappa,C_kappa`, depending only on `kappa`, such that if

```text
0 < d < 1/4,    A <= P^(1-d),    d^3 r >= C_kappa,
```

then, for independent fair signs,

```text
Pr_sigma[ Q intersect B_sigma Z[i] contains a nonzero point ]
    <= C_kappa exp(-c_kappa d^3 r).                 (1)
```

The following weighted statements hold before using comparability. Set `S=log P`, `V_2=sum_j (log n_j)^2`, and `L_max=max_j log n_j`. For a fixed primitive lattice direction `v`, the probability that `Q intersect B_sigma Z[i]` contains a nonzero point on `R v` is at most

```text
exp(-2 (d S - log sqrt(2))^2 / V_2)               (2)
```

when `d S > log sqrt(2)`. Further, witnesses for two sign vectors `sigma,tau` are necessarily collinear whenever

```text
sum_{j: sigma_j != tau_j} log n_j < d S - log 2. (3)
```

Thus (2) and (3) are directly usable even when the block norms differ; the balanced hypothesis gives the uniform cube-enlargement estimate (1).

### Proof

Every integer lattice point on the line `R v`, with `v` primitive as a vector in `Z^2`, is `a v` for some `a in Z`. Define

```text
G_sigma(v) = N(gcd_G(v,B_sigma)).
```

Each `B_sigma` contains at most one Gaussian orientation over every rational prime. If its valuation there is `e` and the valuation of `v` at that orientation is `t`, divisibility `B_sigma | a v` asks that the rational integer `a` contain `p^max(e-t,0)`. Consequently the least positive such `a` is exactly

```text
a_min = P / G_sigma(v).                           (4)
```

The Gaussian factors `gcd_G(v,H_j)` and `gcd_G(v,bar(H_j))`, over all `j`, are pairwise coprime and their product divides `v`. Taking norms and averaging the independent signs gives

```text
E_sigma log G_sigma(v)
  = (1/2) sum_j [log N(gcd_G(v,H_j))
                 + log N(gcd_G(v,bar(H_j)))]
  <= log |v|.                                      (5)
```

If `a v` is a witness in `Q`, then `|a v| <= sqrt(2) R_1`. Since `R_1 <= A <= P^(1-d)`, (4) forces

```text
log G_sigma(v) >= log |v| + d S - log sqrt(2).
```

The independent summands in `log G_sigma(v)` have ranges of length at most `log n_j`. Hoeffding's inequality applied to (5) proves (2).

If `sigma,tau` differ on the block set `J`, their selected products share a Gaussian factor `c` of norm `P/prod_{j in J} n_j`. Corresponding witnesses `u,u' in Q` satisfy `N(c) | det(u,u')`. The rectangle gives `|det(u,u')| <= 2A <= 2P^(1-d)`. Under (3), `N(c)>2A`, so the determinant is zero. This also makes the witness line unique for each bad sign vector when `dS>log2`.

For completeness, the rest of the source's Boolean-cube argument survives unchanged. Under comparability, `S>=rL`, `V_2<=r kappa^2 L^2`, and `L_max<=kappa L`; hence (2) bounds each fixed-line sign class by `exp(-c d^2 r/kappa^2)`, after enlarging `C_kappa`. By (3), different line classes have Hamming distance at least `(dS-log2)/L_max`. Their Hamming neighborhoods of radius `h=floor((dS-log2)/(4L_max))` are disjoint. The cube vertex-boundary inequality used in the source enlarges each class by at least `exp(c d^2 h/kappa^2)` before it reaches relative size `exp(-c d^2 r/kappa^2/2)`. Here `h>=c_kappa d r`, so summing the disjoint neighborhoods proves (1), with adjusted constants. Every block has norm at least five, so `L>=log(5)/kappa`. The threshold `d^3 r>=C_kappa` therefore absorbs integer rounding and the constant logarithms uniformly over the blocks.

## Application to an endpoint circle arc

Let `N=R^2` and `W=log N`. Center angular coordinates at the arc midpoint. Every point has angle `|theta|<=C/(2 sqrt(R))`, so its tangent coordinate has absolute value at most `C sqrt(R)/2` and its normal coordinate lies between `-C^2/8` and `0`. Its entire chord-difference set therefore lies, for large `R`, in an origin-centered rectangle with

```text
A = R_1 W_1 <= K_C sqrt(R) = exp(W/4+O_C(1)).   (6)
```

For a balanced `m`-row Boolean profile from the endpoint extraction, there are `2^m` blocks including the empty block, each with `log N(H_T)=W/2^m+o(W)` for fixed `m`. One pair of rows has equal orientation on `2^(m-1)` blocks, of total log norm `W/2+o(W)`. The corresponding signed block product divides their chord. When `m` is large enough for the lemma's block-count hypothesis, a selection of about `0.4*2^m` of those blocks meets its area hypothesis with a fixed `d>0`; the actual chord merely identifies one of the exponentially sparse exceptional sign vectors. Nothing in the profile makes that vector random.

To compare two chords sharing a row, the blocks with the **same signed orientation in both** are exactly those on which all three rows agree. They are `2^(m-2)` of the full `2^m` blocks, of total log norm `W/4+o(W)`. After omitting the empty block, as in the primitive numerator profile, only `2^(m-2)-1` such blocks remain, with log norm `(1/4-2^(-m))W+o(W)`. For each fixed `m`, this retained product fails to dominate the worst-case determinant ceiling `exp(W/4+O_C(1))` from (6). Including the empty block places the weight at the critical `W/4+o(W)` scale. The fixed-`d` rectangle argument needs a positive proportion of `W` as guaranteed slack and supplies none here; a smaller rectangle or lower-order factors could still decide the strict determinant comparison (3) for a particular configuration. The same issue arises for two disjoint pairs: their shared signed-orientation blocks are at most one quarter of all blocks. Thus the source lemma's fixed-`d` nearby-sign argument does not force collinearity of the actual endpoint chords. A new arithmetic excess above this quarter-weight barrier, or another source of independent sign samples, would be needed.

The comparison does not change if one uses the actual extracted points `z_i=B_i corez_i`: a core block dividing two `corez_i` also divides the corresponding original `z_i`, but the correcting factors supply no additional guaranteed common signed block weight. Pairwise determinant residuals and the finite correction factors remain part of the endpoint problem.

## Displacement entropy on one short circle arc

There is a separate obstruction to transferring the source's coverage induction. Let `S` be a set of `M` points on one circle in an arc shorter than a semicircle, and let `Y,Z_1,...,Z_n` be any jointly distributed variables taking values in `S`. Put `D=(Z_1-Y,...,Z_n-Y)` and `p_0=Pr(D=0)`. For any nonzero displacement `h=Z_j-Y`, the equations `|Y|=|Y+h|=R` intersect a circle with the line `2 Re(bar(h)Y)=-|h|^2`. Their two possible ordered-pair solutions are `(Y,Z_j)` and `(-Z_j,-Y)`. The second is excluded because `S` contains `Z_j` and cannot also contain its antipode `-Z_j`. Therefore `D != 0` determines `Y`, and

```text
H(Y|D) <= p_0 log M,
H(D) >= H(Y) - p_0 log M.                      (7)
```

The sharper residue form avoids comparing `log M` with a selected residue alphabet. Let `F_U` be any homomorphism-valued residue vector selected by auxiliary randomness `U` independent of `(Y,D)`, and let its alphabet have size `a_U`. Since `D != 0` determines `Y`, for every fixed `U` one has

```text
H(F_U(Y)|D) <= p_0 log a_U,
H(D) >= H(F_U(Y)) - p_0 log a_U.              (8)
```

In the source's notation, averaging over independent signed subsets of `s` prime coordinates gives `E_U log a_U=s L_*`. Its Lemma 5.1 and Theorem 5.2 ask for `E_U H(F_U(Y)) >= (1-epsilon)s L_*` and `H(D) <= kappa s L_*`, with `epsilon+kappa` small. Equation (8) instead forces

```text
epsilon + kappa >= 1 - p_0.                   (9)
```

For nontrivial shared displacements, `p_0` is small, so the source's cheap shared-data hypothesis cannot be supplied by points on one short circle arc. As an extreme exact case, choose `(Y,Z)` uniformly among the `M(M-1)` ordered distinct pairs of `S`. Then the displacement map is injective, `D=Z-Y` is uniform on its `M(M-1)` values, and `H(D)=log(M(M-1))`.

The source's `(tau,beta)` coverage condition also requires positive probability on more than `p-p^(1-beta)` residues of a coordinate group of order `p`. Any residue distribution induced from `S` has support at most `M`, so coverage of even one such coordinate requires `M>p-p^(1-beta)`. This support limit is independent of the sieve obstruction below.

## Why the moat sieve itself is constant on a norm fiber

The source defines `A(P_0)` by excluding zero modulo **both** Gaussian factors over every selected rational prime `p in P_0`. For the fixed norm fiber `S_N={z in Z[i]:N(z)=N}`, unique factorization gives the exact dichotomy

```text
A(P_0) intersect S_N = S_N    if gcd(N,prod_{p in P_0} p)=1,
A(P_0) intersect S_N = empty  otherwise.        (10)
```

Indeed `p | N(z)` iff one of the two factors over `p` divides `z`. Hence this finite periodic sieve cannot distinguish members of one equal-norm cluster. Its uniform component theorem concerns bounded *Euclidean* steps between Gaussian **primes**, whereas chords within an endpoint arc may have length `C sqrt(R)`, which grows with `R`, and the lattice points are generally composite. The source's information argument also samples a long self-avoiding bounded-step walk and uses repeated short increment words; a finite same-norm cluster supplies neither that walk nor the zero-avoidance condition. The transferable result at present is (1)--(3), with the precise quarter-weight obstruction above.

The arithmetic identities (4)--(5) have been checked on composite-block examples by [the exact checker](check_gaussian_moat_composite_block_transfer.py): 5,376 least-multiplier identities and 448 averaged-content inequalities pass. This finite check supports the bookkeeping; the proof above establishes the general statement.

## Compressed chord data: lengths and primitive residues

Squared chord lengths retain nearly the same entropy cost. Assume now that
`S` lies in an arc of angular width less than `pi` and length at most
`C sqrt(R)`, with `C<2 sqrt(2)`. The
[multiplicative-rectangle theorem](multiplicative_rectangle_separation.md)
makes its unordered pair products distinct. Equality of two nonzero chord
lengths gives equality of the absolute angular differences, because
`2R sin(theta/2)` is strictly increasing for `0<=theta<pi`.
Orient both pairs in increasing angular order. Equality of their angular
differences gives a multiplicative pair-product identity. The Sidon
conclusion identifies the oriented pairs, since neither is a repeated
point. Thus squared chord length is injective on unordered distinct pairs.

For jointly distributed `Y,Z_1,...,Z_n in S`, set

```text
L = (|Z_1-Y|^2,...,|Z_n-Y|^2),
p_0 = Pr(all entries of L are zero),
p_1 = Pr(exactly one distinct nonzero value occurs in L).
```

On the second event, the nonzero length identifies one unordered edge,
so there are at most two possible values of `Y`. If at least two distinct
nonzero lengths occur, they identify two different unordered edges sharing
`Y`; their intersection has exactly one vertex. Both events are determined
by `L`, so no additional entropy of an event indicator is needed. Hence

```text
H(Y|L) <= p_0 log M + p_1 log 2.
```

For the independent residue vectors in (8), the sharper version is
`H(F_U(Y)|L)<=p_0 log a_U+p_1 log 2`. Therefore high joint entropy
`E_U H(F_U(Y))>=(1-epsilon)s L_*` and the proposed shared-data budget
`H(L)<=kappa s L_*` necessarily imply

```text
epsilon+kappa >= 1-p_0-p_1 log(2)/(s L_*).       (11)
```

Replacing displacements by lengths thus saves at most one bit when a
nonzero continuation occurs, and saves none when two different nonzero
lengths occur. The strict constant restriction is substantive: the
multiplicative-rectangle theorem has sharp endpoint constant `2 sqrt(2)`.

### Primitive residuals lose the additive update

Work within one common Gaussian-unit class, with squared radius `N=R^2`.
The [primitive chord identity](balanced_bonus_phase_audit.md) writes each
distinct pair using a conjugate-primitive Gaussian factor `h=x+iy`, with

```text
t=|y|>=1,    q=N(h)=x^2+t^2 divides N,
|z-z'|^2 = 4t^2 N/q.                              (12)
```

The residual `t` alone need not identify a chord, even on arcs whose
endpoint constant tends to zero. Use the exact
[two-row family](small_imaginary_row_lattice_transference.md#4-a-balanced-full-incidence-two-row-family-at-the-boundary),
whose coprimality identities are covered by its
[checker](check_small_imaginary_row_lattice_transference.py). For positive
`a` divisible by `2210`, put

```text
H=a+i,       Q_1=2a+1-2i,       Q_2=4a+1-4i,
G=H Q_1 Q_2,
z_0=bar(G),  z_1=H Q_1 bar(Q_2),  z_2=H Q_2 bar(Q_1).
```

These three distinct Gaussian integers have radius `R=|G|~8a^3`.
The primitive anchor numerators are
`H Q_1=2(a^2+1)+a+i` and `H Q_2=4(a^2+1)+a+i`.
Both anchor residuals are therefore exactly one. Their positive angular
differences from `z_0` are `2 arctan(1/(2(a^2+1)+a))` and
`2 arctan(1/(4(a^2+1)+a))`; the containing arc has endpoint constant
`O(a^(-1/2))`. Nevertheless,

```text
z_1-z_0=2i bar(Q_2),
z_2-z_0=2i bar(Q_1),
(z_1-z_0)-(z_2-z_0)=4iH.                         (13)
```

The last quantity is nonzero modulo every Gaussian prime dividing `Q_1`,
by the fixture's disjoint norm supports and odd norms. Thus even the
starting residue together with `t=1` does not determine the terminal
residue in these actual block coordinates. There is no additive update
map depending only on those data, as required by the source's candidate
tests. Indeed, choose `Y` uniformly among these three points and move
`z_0` to `z_1`, and each of `z_1,z_2` to `z_0`. Every move has residual
one, so its residual data have entropy zero while `H(Y)=log 3`.

### An exact divisor-fiber bound

There is still a finite arithmetic bound on what the compression loses.
Under the same common-unit and short-arc hypotheses, define

```text
nu_N(t) = #{q | N : q-t^2 is a positive integer square,
                    q >= (4t^2/C^2) sqrt(N)}.
```

The square condition and lower bound follow from (12); positivity follows
because an imaginary primitive numerator would give antipodal points.
For fixed `t` and `q`, (12) determines the squared chord length, and hence
at most one unordered edge. Consequently

```text
#{unordered pairs with primitive residual t} <= nu_N(t) <= tau(N). (14)
```

For any random distinct ordered pair `(Y,Z)` with residual `T`, this gives
the candidate bound

```text
H(Y|T) <= E log(2 nu_N(T)).                         (15)
```

For a uniform unordered pair, with residual again denoted `T`, the exact
fiber decomposition gives the stronger counting statement

```text
H(T) >= log binom(M,2) - E log nu_N(T).            (16)
```

No independence is required in (15). These inequalities quantify the
missing information, but `tau(N)` is not bounded independently of `N`.
Restoring `q` to the observation restores the length and the two-candidate
entropy floor. Retaining `t` alone requires a new bound on its divisor
fibers, or a candidate-testing argument that tolerates the genuinely
different updates in (13). Neither such bound nor such a test is proved
here; compressed residuals presently give no radius-uniform count.

The new checker also verifies sixteen instances of (13), including the
primitive anchor residuals, exact common-factor norms and nonvanishing
of the update difference modulo every prime dividing `Q_1`. The family
and its all-parameter proof were independently checked; the finite
fixtures are supplementary.

### Ordinary divisor fibers can hide an unbounded orientation cost

The condition defining `nu_N(t)` is necessary, not sufficient for the
selected primitive numerators to fit the same circle. For an anchored
tuple `h_j/bar(h_j)`, its least squared radius is
`Norm(lcm_G(h_j))`. The ordinary condition `Norm(h_j)|N` need not retain
both Gaussian orientations over a rational prime. The following family
shows that the resulting error has no absolute bound, even for two
positive primitive numerators with imaginary coordinate one.

Take any even integer `r>=2` and integer `s>=1`, and put

```text
D=r^2+1,        a=2rDs,
x_+=a+r,       x_-=a-r,
h_+=x_++i,     h_-=x_-+i,
q_+=Norm(h_+), q_-=Norm(h_-),
F_+=4r^2 D s^2+4r^2 s+1,
F_-=4r^2 D s^2-4r^2 s+1.
```

Expansion gives `q_+=D F_+` and `q_-=D F_-`. The two `F` are odd,
are congruent to one modulo `r^2 s`, and differ by `8r^2 s`.
Consequently `gcd(F_+,F_-)=1` and `gcd(q_+,q_-)=D`. In contrast,
a Gaussian common divisor of `h_+,h_-` divides `2r`. Both norms are
odd and congruent to one modulo `r`, so that common divisor is a unit.
Each `h` is conjugate-primitive since its real coordinate is even and
its imaginary coordinate is one. Thus

```text
N_ord=lcm_Z(q_+,q_-)=q_+ q_-/D,
Norm(lcm_G(h_+,h_-))=q_+ q_-=D N_ord.             (17)
```

Both ordinary divisors meet the threshold for `nu_(N_ord)(1)` at `C=2`.
Indeed `q_-^2>=N_ord` is equivalent to `D F_- >= F_+`, and

```text
D F_- - F_+
  =r^2[4r^2 D s^2-4s(D+1)+1]>0.                 (18)
```

Here `r^2 D>=D+1` and `s>=1`. Hence `nu_(N_ord)(1)>=2`.
Nevertheless no Gaussian integer of norm `N_ord` is divisible by both
primitive numerators. This construction makes the discrepancy `D`
arbitrarily large while keeping both numerators odd-norm and primitive.

The actual primitive circle realization is

```text
z_0=bar(h_+ h_-),
z_+=h_+ bar(h_-),       z_-=h_- bar(h_+),
R=sqrt(q_+ q_-).
```

Its three points are distinct and have unit common Gaussian gcd. Their
angles relative to `z_0` range from zero to `2 arctan(1/x_-)`, so the
exact normalized arc length is

```text
C_*=2(q_+ q_-)^(1/4) arctan(1/x_-).             (19)
```

For every fixed even `r`, as `s` tends to infinity, `C_* -> 2`.
Using the ordinary radius instead would predict `D^(-1/4) C_*`, tending
to `2D^(-1/4)`. Choosing large `r` makes this false prediction arbitrarily
small. The family does not disprove an upper bound on `nu_N`, nor the
uniform endpoint conjecture. It shows why its qualifying divisors cannot
be promoted to simultaneous circle points at radius `sqrt(N)`.

Root derived the family and Luna independently checked its all-parameter
proof. The [signed-block checker](check_signed_block_flip_witness_height.py)
also checks 120 exact instances, including the ordinary and Gaussian
lcm norms and the primitive actual circle. The all-edge radius formula
in [the cotangent note](integer_cotangent_lcm_height_target.md) retains
the missing orientations and is consistent with (17).
