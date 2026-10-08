# Recent literature check and the short-difference corollary

Checked 12 September 2026. Target: the number of integer points on an arbitrary arc of length at most `C sqrt(R)` of the origin-centered circle of radius `R`, with a bound depending only on `C`.

## Primary-source findings

No endpoint resolution or uniform primitive-Pythagorean short-arc theorem was located in this bounded search. This is a search conclusion, not a claim that absence from the search proves the literature status.

1. **Explicit recent status evidence.** Germain–Moyano–Zhu, *On the vanishing of eigenfunctions of the Laplacian on tori*, version 3 dated 21 September 2025, Remark 6.6, still describes uniform boundedness on arcs of length `lambda^(1-epsilon)` on the radius-`lambda` circle as the Cilleruelo–Granville conjecture. Appendix C.2 uses the classical bound of `m` points for lengths less than `sqrt(2) lambda^(1/2-delta(m))`, where `delta(m) = 1/(4 floor(m/2)+2)`. Here `lambda = R`, not `R^2`. [Full text](https://arxiv.org/html/2406.19925v3); [version history](https://arxiv.org/abs/2406.19925).

2. **A relevant technical theorem.** Gabdullin, *Trigonometric Polynomials with Frequencies in the Set of Squares and Divisors in a Short Interval*, JFAA **30:2 (2024)**, Theorem 1.6: let `gamma_0 = (sqrt(5)-1)/2` and `k = N^(gamma_0-epsilon)`. For every positive integer `m <= 3Nk`,

   `#{d | m : 2N <= d <= 2N+2k} << epsilon^(-1)`.

   Its Theorem 1.4 gives `||f||_4 << epsilon^(-1/4) ||f||_2` for square frequencies `n^2` with `N <= n <= N+N^(gamma_0-epsilon)`. This is an additive-energy theorem, not a pointwise sum-of-two-squares theorem. [Full text, Theorems 1.4 and 1.6](https://arxiv.org/pdf/2205.13611); [journal record](https://link.springer.com/article/10.1007/s00041-023-10064-w).

3. **A 2026 near-match does not improve this endpoint.** *Trace and observability inequalities for Laplace eigenfunctions on the torus*, Forum Math. Sigma **14:e47 (2026)**, Section 3.3, uses classical cap estimates for dimension at least three; its two-dimensional clusters have diameter `O(lambda^(1/3))` and at most two points. No `R^(1/2)` circle bound is supplied. [Primary article, Lemmas 3.4–3.6](https://doi.org/10.1017/fms.2026.10199).

Oganesyan's *Lattice points on small arcs* is not evidence of a disproof: the current arXiv record explicitly states withdrawal on 22 August 2021 because of a mistake in Lemma 2. The obsolete disproof abstract still appears in search results. [Withdrawal record](https://arxiv.org/abs/2107.09991).

## Exact nonzero-difference corollary

Let `N >= 1`, `0 <= k <= N`, and

`r_I(m) = #{(x,x') in ([N,N+k] intersect Z)^2 : x^2-x'^2 = m}`.

For nonzero integer `m`, replacing `(x,x')` by `(x',x)` reduces to `m > 0`. Factor

`m = (x-x')(x+x')`.

The divisor `d = x+x'` lies in `[2N,2N+2k]`, and `m <= 2Nk+k^2 <= 3Nk`. A divisor `d` determines at most one ordered pair, since

`x = (d+m/d)/2`, `x' = (d-m/d)/2`.

Parity or interval restrictions only discard candidates. Therefore Gabdullin's theorem bounds `r_I(m)` uniformly when `k <= N^(gamma_0-epsilon)`, for fixed `0 < epsilon < gamma_0-1/2`, after enlarging the interval to the theorem's length if necessary. In particular this applies to `k <= A sqrt(N)` once `N` is sufficiently large depending on `A` and `epsilon`. The zero difference must be excluded: `r_I(0) = #([N,N+k] intersect Z)`.

At the exact endpoint there is also a simpler, explicit elementary bound:

`r_I(m) <= 1 + floor(k^2/N)` for every nonzero `m`.

Indeed, `q = x-x'` is a positive integer in

`[m/(2N+2k), m/(2N)]`.

This interval has length `mk/(2N(N+k)) <= k^2/N`, using `m <= 2Nk+k^2`; each `q` determines at most one pair. Consequently `k <= A sqrt(N)` gives `r_I(m) <= 1+floor(A^2)`. Thus endpoint difference multiplicity alone is already elementary. Gabdullin extends this kind of control to longer intervals.

## Application to arbitrary short arcs, including arcs near axes

One must not assume both original coordinates have magnitude comparable with `R`. This can be arranged at constant cost using an injective integer linear map.

Identify a point with the Gaussian integer `z = x+iy`. Choose either `z -> z` or `z -> (1+i)z` so that the midpoint direction of the resulting arc is at angular distance at least `pi/8` from every coordinate axis. Such a choice always exists: the two rotations differ by `pi/4`. For sufficiently large `R` depending on `C`, the entire arc then stays at distance at least `pi/16` from the axes. Coordinate sign changes make both coordinates positive. Both maps preserve integrality and are injective.

Writing `q` for the scale factor, `q` is either `1` or `sqrt(2)`, the new radius is `R' = qR`, and the arc length is

`L' <= q C sqrt(R) = sqrt(q) C sqrt(R') <= 2^(1/4) C sqrt(R')`.

Both coordinate sets lie in integer intervals `[N_x,N_x+k_x]` and `[N_y,N_y+k_y]` with `N_x,N_y` comparable with `R'` and `k_x,k_y = O_C(sqrt(R'))`. The preceding corollary therefore gives a bound depending only on `C` for every nonzero squared-coordinate difference, separately in each coordinate. The remaining bounded range of radii can be bounded directly. This normalization does not assert that `N_x` and `N_y` are equal or close to each other.

## Why this does not bound the number of circle points

Let the normalized points be `(x_i,y_i)`, with

`x_i^2+y_i^2 = n'`, where `n' = (R')^2`.

They lie in the positive quadrant, so distinct circle points have distinct `x_i` and distinct `y_i`. Put `a_i=x_i^2`, `b_i=y_i^2`. Circle compatibility gives

`a_i-a_j = -(b_i-b_j)`.

The short-difference corollary limits how often one fixed nonzero difference can occur. It does not limit the number of different differences generated by the points. With `M` points, there can be order `M^2` distinct differences, each with bounded multiplicity. Counting the displayed compatibilities therefore gives an upper bound of order `M^2`, which is consistent with the `M(M-1)` compatibilities already present.

Likewise, on the union of the two square-coordinate sets, the `M` representations `a_i+b_i=n'` contribute only order `M^2` to additive energy. A uniform Lambda-4 estimate allows precisely this order of energy for a set of order `M` elements.

There is an explicit abstract obstruction to deducing bounded sum multiplicity from these facts alone. For any `M`, choose `T > 5^(M-1)` and set

`A = {T+5^j : 0 <= j < M}`, `B = {T-5^j : 0 <= j < M}`, `S = A union B`.

Each of `A` and `B` has nonzero difference multiplicity at most one; `S` has nonzero difference multiplicity at most two. Uniqueness follows by comparing the balanced base-five coefficients of a difference; these coefficients are between `-2` and `2`. Yet `2T` has `2M` ordered sum representations in `S`. Moreover,

`E(S) = |S|^2 + sum_(d != 0) r_(S-S)(d)^2 <= 3|S|^2`.

Even the weighted Lambda-4 inequality holds with constant at most `3^(1/4)`, by Cauchy–Schwarz on the at most two terms representing each nonzero difference. Thus uniform difference multiplicity and a uniform weighted energy inequality are compatible with arbitrarily many matched sums.

This example is not a circle counterexample: its entries need not be squares. It identifies the exact gap. A successful continuation must exploit additional arithmetic of the simultaneous square conditions beyond fixed-difference multiplicity or Lambda-4 control. No valid further implication to a radius-independent circle-point bound was obtained in this check.
