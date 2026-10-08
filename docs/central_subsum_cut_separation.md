# Prime cuts exclude vanishing subsums of the central vector

The nonvanishing-subsum assumption in the
[six-point radius formula](six_point_circle_cover_radius_content.md)
has a direct arithmetic test. In particular it holds eventually on a full
fair core whose independent cut factors dominate the normalization and
private-factor errors. This removes a degeneracy from that reduction;
it does not bound the height of its moving isotropic frame or resultant.

## 1. A valuation criterion

Let `n=2h`, and let `u` be the ordinary primitive central vector of `n`
distinct rational circle nodes in an integer cotangent normalization `L`.
Its entries are nonzero, and `sum u_i=0`. At a split prime `p`, let `e_i`
be the Gaussian allocation exponents of the primitive circle points.
Write

```text
B_i(e) = (sum_j |e_i-e_j| - H_h(e))/2,
H_h(e) = sum of the largest h entries - sum of the smallest h entries.
```

The [primitive-column calculation](integer_cotangent_eigenvector_heights.md)
gives the all-prime error estimate

```text
|v_p(u_i)-B_i(e)| <= (n-1)v_p(2L).                 (1)
```

Indeed each of the `n-1` barycentric denominator differences has excess
valuation between zero and `v_pi(2L)`. Subtracting the smallest raw column
valuation changes a coordinate's model error by at most that total.
Average over the two conjugate Gaussian primes; the normalized vector is
ordinary integral. No assertion that `p` avoids the normalization is needed.

Suppose for a subset `T` of size `s<h` that, after a common translation and
possibly conjugating the Gaussian prime, the allocation vector is within
`epsilon_p` in maximum norm of the two-level vector

```text
e_i^* = a if i belongs to T, and 0 otherwise,       a>0.
```

For this vector `B_i(e^*)=(h-s)a` on `T` and zero off `T`.
Also

```text
|B_i(e)-B_i(e^*)| <= 2n epsilon_p.                 (2)
```

For proof, the sum of the `n` absolute differences changes by at most
`2n epsilon_p`. The function `H_h` is the minimum over real `t` of
`sum_j |e_j-t|`, so it changes by at most `n epsilon_p`. Dividing their
combined error by two even gives the slightly smaller constant `3n/2`.
The integer constant in (2) is convenient.

Set `E_p=(n-1)v_p(2L)+2n epsilon_p`. If

```text
(h-s)a > 2E_p,                                    (3)
```

then every coefficient indexed by `T` has strictly larger `p`-valuation
than every coefficient outside `T`. Consequently, for any single `k`
outside `T`,

```text
sum_(i in T union {k}) u_i != 0.                  (4)
```

In that sum `u_k` is the unique term of smallest valuation, so its valuation
cannot cancel. For an exact two-level cut at `p` not dividing `2L`, this
applies with `epsilon_p=E_p=0` and every positive height `a`.

## 2. A full collection of dominant cuts removes every proper subsum

It suffices, for each subset `T` with `1<=|T|<h`, to have one split prime
with the corresponding dominant cut (3). For any subset `S` of size
`2<=|S|<=h`, choose `k` in `S` and put `T=S minus {k}`. Equation (4)
excludes its vanishing sum. Singleton sums are already nonzero; sums over
larger proper subsets are negatives of sums over their complements.

The same argument applies when a cut's mass is spread over many primes.
For a fixed `T`, suppose those primes have heights `a_p`, approximation
errors `epsilon_p`, and

```text
W_T = sum_p a_p log p,
sum_p [(n-1)v_p(2L)+2n epsilon_p] log p = o(W_T).
```

Then some prime satisfies (3) for all sufficiently large members of the
sequence. Otherwise summing the reverse inequality would give
`(h-|T|)W_T <= o(W_T)`, a contradiction.

This is the precise error condition needed for an application to a fair
core. In the [independent-block reconstruction](independent_block_reconstruction.md),
distinct blocks and their conjugates are coprime. Thus the prime factors of
one cut have exactly its two-level allocation until the private factors
are included. If `r_p=max_i v_p(Norm(K_i))`, those extra allocations, after
common primitive normalization, differ from the model by at most `r_p`
up to a common translation. Moreover

```text
sum_p r_p log p <= sum_i log Norm(K_i) = 2 sum_i log |K_i|.
```

For fixed `n`, `log L=o(w)`, total private log height `o(w)`, and every cut mass
`W_T=w+o(w)` imply the stated error condition. Factors of two do not matter
here; using Gaussian log modulus instead of log norm rescales `w` by two.
This application retains prime powers, including cut factors supported on
primes dividing `L`; it does not discard their valuations by prime support.
If the number of labels also grows, the sufficient error assumption is
instead `n(log(2L)+sum_i log |K_i|)=o(w)`. The fixed six-label consequence
below needs no such additional uniformity in `n`.

## 3. Consequence and limit for six labels

For six labels it is enough to have dominant singleton and pair cuts.
The primitive central vector then has no vanishing proper subsum. Its
isotropic pencils have no colliding columns, and the degree-twelve form
`Delta` and degree-twenty form `K` in the radius calculation have no common
projective zero. The exact radius and its fixed-frame comparisons therefore
apply to every sufficiently developed full fair six-label core.

This does not make the comparisons uniform. Both `u` and its actual
isotropic frame still move with the core; the Bezout constants and the
minimum of `|K|` on the real circle locus can vary. Nor does this criterion
show that `Delta` is squarefree: nonvanishing proper subsums excludes its
common zeros with `K`, which is a different condition.

The checker `check_central_subsum_cut_separation.py` verifies (1), exact
two-level witnesses for (4), and the Lipschitz estimate (2) on finite exact
data. These checks supplement the valuation proof, not the asymptotic
existence of any fair core.
