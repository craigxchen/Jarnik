# Fourier/Riesz-product route for the uniform endpoint problem

## Status

This is a new active route. It is harmonic rather than determinant-, cube-, or exceptional-set based. The uniform bound is not proved here.

The core advantage is that the complete angular measure of lattice points on one circle factors exactly over Gaussian prime-power blocks. Arc counting can therefore be attacked through a short-frequency majorant and a product of one-dimensional Dirichlet kernels.

## 1. Angular measure and exact factorization

Write the split part of the common norm as

```text
N = product_j p_j^(e_j),
```

choose one Gaussian prime `pi_j` above `p_j`, and write

```text
theta_j = arg(pi_j).
```

Ignoring a fixed global rotation, every lattice-point angle is

```text
alpha(a) = sum_j (2 a_j - e_j) theta_j mod 2 pi,
0 <= a_j <= e_j.
```

Let `mu_N` be the normalized counting measure on these angles. Then

```text
mu_N = convolution_j mu_j,
```

where

```text
mu_j = (1/(e_j+1)) sum_(a=0)^e_j delta_((2a-e_j)theta_j).
```

Hence the Fourier coefficients factor exactly:

```text
mu_hat_N(k)
  = product_j mu_hat_j(k),

mu_hat_j(k)
  = (1/(e_j+1))
      sum_(a=0)^e_j exp(i k (2a-e_j) theta_j).
```

Equivalently,

```text
|mu_hat_j(k)|
  = |sin((e_j+1) k theta_j)|
      / ((e_j+1)|sin(k theta_j)|),
```

with the usual limiting value `1` when the denominator vanishes.

Thus every nonzero frequency is damped multiplicatively unless many prime blocks are simultaneously resonant.

## 2. Arc majorization

Let `I` be an angular interval of length

```text
delta = C N^(-1/4).
```

Choose `T` comparable to `delta^(-1)`, hence `T ~ N^(1/4)`. A Fejer/Selberg majorant gives schematically

```text
mu_N(I)
  <= O(delta)
     + O(1/T) sum_(1 <= |k| < T)
         |mu_hat_N(k)|.
```

Multiplying by the total number of lattice points

```text
R_N = product_j (e_j+1)
```

reduces the uniform endpoint conjecture to a short-frequency product estimate:

```text
R_N/T
  + (R_N/T) sum_(1 <= |k| < T)
      product_j |mu_hat_j(k)|
  = O_C(1).
```

The first term is not automatically bounded, since `R_N` can grow. The Fourier damping must therefore cancel both the zero-frequency mass and the nonzero-frequency contribution at the endpoint scale.

A smoother majorant or a block-adapted Riesz product may be preferable to the raw Fejer kernel because it can build cancellation directly from selected prime blocks.

## 3. Resonant blocks

For a fixed frequency `k`, call block `j` resonant when

```text
|| k theta_j / pi ||
```

is small enough that `|mu_hat_j(k)|` is close to `1`.

A crude inequality is

```text
|mu_hat_j(k)|
  <= min(1, 1 / ((e_j+1)|sin(k theta_j)|)).
```

Thus a large product coefficient requires a large portion of conductor weight to lie on blocks satisfying near-root-of-unity conditions

```text
(pi_j / conjugate(pi_j))^k approximately 1.
```

The new target is a conductor-weighted resonance theorem:

> Uniformly for `1 <= k <= N^(1/4)`, either a bounded set of prime blocks carries a positive proportion of `log N`, or the total conductor weight of blocks with `|mu_hat_j(k)|` close to `1` is bounded away from `log N`.

If such a statement holds with a quantitative gap, the product of the remaining block factors decays exponentially in conductor height and the Fourier majorant gives an absolute arc bound.

## 4. Why this is genuinely different

Previous methods selected a finite subset of points and reduced local arithmetic to cut metrics. The present method uses the complete divisor-angle measure. It retains all prime-block independence through convolution and does not encode the configuration as pairwise distances.

The endpoint obstruction is now simultaneous frequency resonance, not Plotkin/max-cut equality.

The route can therefore succeed even when the local point set is additive-Sidon-like or cube-free.

## 5. Immediate obstacles

### 5.1 One-block resonances can be extremely accurate

The elementary algebraic lower bound for

```text
|(pi/conjugate(pi))^k - 1|
```

is exponentially weak in `k log p`. Since `k` may be as large as `N^(1/4)`, treating each block separately is hopeless.

The argument must control simultaneous resonance across many distinct Gaussian primes, or average over frequencies.

### 5.2 The zero-frequency term

A standard positive majorant contributes roughly `delta R_N`. This is not uniformly bounded. A successful proof needs either:

- enough cancellation in a signed extremal majorant;
- a multiscale decomposition in which prime-block smoothing reduces the effective mass before the final arc test; or
- an induction/descent step for circles with unusually large total divisor count.

### 5.3 Average versus worst frequency

It may be easier to prove

```text
sum_(k <= T) product_j |mu_hat_j(k)|^2
```

is small than to control every frequency. This suggests using an `L^2` or higher-moment majorant rather than `L^1` Fourier inversion.

## 6. Three concrete subroutes

### A. Conductor-weighted large sieve

Prove a short-frequency estimate of the form

```text
sum_(|k| <= T) |mu_hat_N(k)|^2
  <= C T / R_N^2
```

or a sufficiently strong variant. Combined with Cauchy-Schwarz and a band-limited interval majorant, this would bound local mass.

Expanding the square converts the left side into a kernel sum over pairs of divisor angles. Unlike earlier pairwise arguments, all pairs are averaged with oscillating signs through the frequency kernel.

### B. Block-adapted Riesz product

Select prime blocks greedily and construct a nonnegative trigonometric polynomial

```text
P(alpha) = product_j (1 + rho_j cos(k_j(alpha-alpha_0)))
```

that majorizes the target arc but has very small integral against `mu_N` because each selected factor is nearly orthogonal to one block measure. The frequencies must remain below `N^(1/4)` after multiplication.

This is an uncertainty-principle formulation: too many independent prime-angle convolutions should prevent the measure from concentrating at endpoint scale.

### C. Resonance graph across frequencies

For each block `j`, record the set

```text
R_j = {1 <= k <= T : |mu_hat_j(k)| >= 1-epsilon}.
```

A large arc count forces many frequencies lying in the intersection of many heavy-block resonance sets. Show that such a common intersection implies a low-denominator rational approximation shared by many Gaussian prime angles. Converting that approximation into Gaussian integer geometry may force all those primes into one narrow ray or yield an impossible multiplicative relation.

## 7. First theorem to prove

The most concrete first target is:

> **Short-frequency damping dichotomy.** There exist constants `eta,c>0` such that for every circle conductor and `T=N^(1/4)`, either a set of at most `O_C(1)` prime blocks carries at least `eta log N`, or
>
> ```text
> #{1 <= k <= T : |mu_hat_N(k)| >= exp(-c log N)}
>   <= T exp(-c log N).
> ```

The heavy-block alternative is compatible with existing descent. The damping alternative should imply a uniform arc bound through a suitable Fourier majorant.

## 8. Falsification tests

Before investing in a proof, test the route against:

1. the known four-point endpoint circle;
2. squarefree conductors built from Gaussian primes in a very narrow sector;
3. repeated prime powers with one dominant block;
4. artificial angle sets realizing simultaneous rational resonances;
5. random Gaussian-prime conductors, where strong Fourier damping should be visible.

The route should be abandoned if explicit actual Gaussian conductors produce `R_N`-sized low-frequency coefficients throughout a positive proportion of `1 <= k <= N^(1/4)` without a heavy-block descent structure.
