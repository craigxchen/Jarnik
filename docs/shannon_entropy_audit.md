# Shannon entropy audit for the Gaussian divisor-angle measure

## Status

This note studies whether multiscale Shannon entropy can control endpoint-scale lattice-point concentration on circles.

The conclusion is negative for ordinary and branchwise Shannon entropy, for a precise reason: the divisor-angle measure has only `R_N = product_j (e_j+1) = N^{o(1)}` atoms, so its total information content is `log R_N = o(log N)`, whereas the endpoint angular scale has depth `(1/4) log N`.

This does not rule out all entropy methods. It isolates what a successful replacement would need: an exceptional-set-sensitive, non-averaged quantity such as high-order Renyi/min-entropy together with an inverse theorem uniform in the order, or an arithmetic conditional-entropy theorem that retains the product structure after conditioning.

---

## 1. The measure

Write

```text
N = product_j p_j^(e_j)
```

and let

```text
R_N = product_j (e_j+1)
```

be the number of divisor-angle choices, ignoring the fixed unit factor.

The normalized angular counting measure is

```text
mu_N = (1/R_N) sum_a delta_(alpha(a)),
```

where

```text
alpha(a) = sum_j (2 a_j-e_j) theta_j mod 2 pi.
```

It factors exactly as a convolution of the block measures

```text
mu_j = (1/(e_j+1)) sum_(a=0)^e_j delta_((2a-e_j)theta_j).
```

Thus convolution is natural. The question is whether Shannon entropy at angular scale

```text
delta = C N^(-1/4)
```

can detect an interval containing many atoms.

---

## 2. Total entropy ceiling

For any finite or countable measurable partition `P`,

```text
H_P(mu_N) <= log R_N.
```

This is simply because `mu_N` is supported on at most `R_N` atoms.

The standard divisor bound gives

```text
log R_N = o(log N).
```

At the endpoint scale, a dyadic partition has depth

```text
n = (1/(4 log 2)) log N + O_C(1).
```

Therefore

```text
H_n(mu_N)/n -> 0.
```

In the language of multiscale entropy, the measures have asymptotic entropy dimension zero for a trivial support-size reason, regardless of whether any endpoint arc has multiplicity one or multiplicity tending to infinity.

This is a fundamental mismatch with the positive-dimension self-similar setting in which Hochman-type entropy growth theorems are strongest.

---

## 3. A single exceptional cell is almost invisible to Shannon entropy

Suppose one endpoint interval contains exactly `M` of the `R_N` atoms, and assume for simplicity that all remaining atoms occupy separate endpoint cells.

The endpoint-scale Shannon entropy is then

```text
H_delta(mu_N)
  = log R_N - (M/R_N) log M.
```

Thus the entropy deficit from the cluster is only

```text
(M/R_N) log M.
```

This can tend to zero even when `M -> infinity`, for example whenever `M = R_N^{o(1)}`.

Hence ordinary Shannon entropy is an average statistic and cannot detect an arbitrarily sparse exceptional cluster.

---

## 4. Branchwise entropy

Fix a nested dyadic chain of cells ending at an endpoint cell `I` containing `M` atoms.

Let `I_0` be the whole circle and

```text
I_0 superset I_1 superset ... superset I_n = I.
```

The sum of the conditional information increments along this branch is exactly

```text
-log mu_N(I) = log(R_N/M).
```

Therefore even the entire branchwise information down to angular scale `N^(-1/4)` is bounded by

```text
log R_N = o(log N).
```

The branch has geometric depth comparable to `(1/4) log N`, but carries only `o(log N)` Shannon information.

So branchwise Shannon entropy does not repair the scale mismatch.

---

## 5. Conditioning on the exceptional arc

Let `A=(A_j)` be the independent exponent vector, uniformly distributed on the conductor box, and let

```text
E = {alpha(A) lies in I}.
```

If `I` contains `M` points, then

```text
P(E) = M/R_N.
```

The conditional distribution of `A` given `E` is uniform on those `M` exponent vectors. Its relative entropy with respect to the product distribution is

```text
D(P_(A|E) || P_A) = log(R_N/M).
```

By the chain rule, this information can be decomposed among coordinates or packets. However, its total size is still at most `log R_N = o(log N)`.

More importantly, conditioning destroys independence. Standard convolution entropy-growth theorems apply before conditioning, while the measure carrying the exceptional cluster after conditioning is no longer a convolution of independent block measures.

Thus conditional Shannon entropy identifies a global dependence pattern but does not by itself preserve the factorization needed for an inverse convolution theorem.

---

## 6. Smoothing and differential entropy

Convolve `mu_N` with a uniform kernel of width `delta`.

For well-separated atoms, the differential entropy is approximately

```text
log delta + log R_N.
```

Collisions among `M` atoms alter only the discrete part, by at most order `log M`, and in the global averaged entropy by approximately `(M/R_N) log M`.

Hence smoothing does not create new Shannon information; it merely adds the common term `log delta` to every atom. The useful variable part remains bounded by `log R_N`.

---

## 7. Why conductor-weighting is not Shannon entropy

One might try to assign the block `j` an entropy weight `log p_j`. But a random variable with only `e_j+1` outcomes has Shannon entropy at most

```text
log(e_j+1).
```

No choice of probability distribution on that block can give it entropy `e_j log p_j` when `p_j` is large.

Therefore a functional such as

```text
sum_j (log p_j) Delta_j
```

may be a useful energy, but it is not the Shannon entropy of the divisor-angle probability space and does not automatically satisfy the chain rule, data processing, or the hypotheses of known inverse entropy theorems.

Artificially splitting one prime block into `log p_j` independent random steps changes the measure and introduces randomness that the circle does not possess.

---

## 8. Relation to Hochman and Shmerkin machinery

Hochman-type inverse theorems classify failure of entropy growth for measures carrying entropy at a positive proportion of the available scales. Here

```text
H_n(mu_N)/n -> 0
```

before any exceptional clustering is considered.

Thus the generic structured conclusion is simply that the measure is sparse or atomic, which is already known and is not the arithmetic rigidity required for the endpoint theorem.

Shmerkin's inverse `L^q` theorems are more sensitive to concentration, but any fixed `q` still averages over the complete measure. A cluster of size `M` among `R_N` atoms changes the `q`-energy by a relative amount of order

```text
M^q / R_N.
```

For an unbounded family with `M = R_N^{o(1)}`, this is invisible for every fixed `q`. Detecting it requires `q` growing with `R_N`, while available inverse theorems are not uniform in such a regime.

---

## 9. Needle-in-a-haystack obstruction

The common issue is now precise.

The conjecture asks for control of the maximum endpoint-cell multiplicity:

```text
max_I # {a : alpha(a) in I}.
```

Shannon entropy, fixed-order Renyi entropy, Fourier moments, and global additive energy are averaged statistics. They can all ignore a cluster whose cardinality tends to infinity but occupies a vanishing fraction of the complete divisor set.

Thus a successful global tool must be simultaneously:

1. sensitive to the worst cell rather than a typical atom;
2. compatible with the convolution/product structure;
3. uniform when the exceptional mass fraction tends to zero;
4. arithmetic enough to exploit the actual Gaussian prime angles.

Ordinary Shannon entropy satisfies the second property but not the first or third.

---

## 10. What entropy direction remains plausible

Two entropy-adjacent directions remain logically possible.

### A. Uniform high-order Renyi inverse theorem

Develop an inverse theorem for `q` growing with the conductor, strong enough to detect a cell with `M -> infinity` even when `M = R_N^{o(1)}`, while retaining arithmetic information about the block angles.

No such theorem is currently available in the needed form, and its constants would have to be controlled very sharply.

### B. Product-event information theorem

Work directly with the rare event

```text
E = {sum_j X_j lies in I}
```

for independent block variables `X_j`, and classify the conditional dependence structure when `E` contains many atoms despite exponentially small geometric width.

This would be an inverse small-ball theorem phrased in information-theoretic language. It must use more than the scalar relative entropy `log(R_N/M)`; it would need to retain the arithmetic locations and prove that the event is generated by a bounded collection of coordinates or by an exact common resonance.

This is closer to inverse Littlewood--Offord theory than to standard Shannon entropy growth.

---

## Verdict

Shannon entropy was worth exploring, but it does not solve the scale mismatch:

```text
available Shannon information = log R_N = o(log N),
endpoint geometric depth      = (1/4) log N.
```

Moreover, one unbounded exceptional cluster can have vanishing global entropy effect.

The most useful lesson is not merely that entropy fails. It identifies the required new feature: a global invariant must retain convolution structure while being worst-case sensitive. Any averaged statistic of fixed order is vulnerable to the same sparse-exception obstruction.
