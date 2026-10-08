# Fourier route: moment stress test and phase barrier

## Status

This note stress-tests the Fourier/Riesz-product route for the uniform endpoint
problem. The uniform bound is not proved.

The exact product formula for the complete divisor-angle measure is valid and
useful. However, every estimate using only the magnitudes of the Fourier
coefficients fails by an exponential factor even in the ideal independent
squarefree model. A successful harmonic proof would have to exploit coherent
phase cancellation in the short Fourier inversion sum. At that point the
problem is no longer reduced by the product formula alone.

## 1. Setup

For

```text
N = product_j p_j^(e_j),
R_N = product_j (e_j+1),
H = log N,
T asymp exp(H/4),
```

let `mu_N` be the normalized angular measure on the circle points. Then

```text
mu_hat_N(k) = product_j mu_hat_j(k),

|mu_hat_j(k)|
  = |sin((e_j+1)k theta_j)|
      / ((e_j+1)|sin(k theta_j)|).
```

For an interval `I` of length `delta asymp 1/T`, a Selberg/Fejer majorant gives

```text
#(points in I)
  <= O(R_N/T)
     + (R_N/T) sum_(1 <= |k| < T)
         c_k exp(-ik x_I) mu_hat_N(k),
```

with bounded coefficients `c_k` of triangular type.

The zero-frequency term is not the real obstruction. Since `R_N <= N^epsilon`
for every fixed `epsilon>0` and sufficiently large `N`,

```text
R_N/T = R_N N^(-1/4) = o(1).
```

The problem is the nonzero-frequency sum.

## 2. L1 magnitude estimates fail in the independent model

Consider the squarefree model `e_j=1` with `r` blocks, so

```text
R_N = 2^r,
mu_hat_N(k) = product_(j=1)^r cos(k theta_j).
```

If the phases `k theta_j` behave independently and uniformly, then

```text
average_k |mu_hat_N(k)|
  approximately (2/pi)^r
  = R_N^(-gamma_1),

gamma_1 = log_2(pi/2) approximately 0.6515.
```

Therefore absolute-value Fourier inversion yields only

```text
#I <= R_N^(1-gamma_1+o(1)),
```

which is exponentially large in `r`, not uniformly bounded.

Thus even the ideal generic model does not support the required estimate

```text
average_k |mu_hat_N(k)| << 1/R_N.
```

## 3. L2 estimates fail

A Fejer-kernel pair count gives

```text
M^2
  <= C R_N^2
       (1/T) sum_(|k|<T) |mu_hat_N(k)|^2,
```

for `M` points in one interval of length `asymp 1/T`.

In the same independent squarefree model,

```text
average_k |mu_hat_N(k)|^2
  approximately (1/2)^r
  = 1/R_N.
```

Hence

```text
M <= C sqrt(R_N),
```

again far from an absolute bound.

This shows that a conductor-weighted large sieve with the natural random-model
strength cannot finish the endpoint theorem.

## 4. Arbitrary Fourier moments still fail

Using `q`-fold sums of the points in the target interval and a Fejer kernel at
scale `T/q` gives schematically

```text
M^(2q)
  <= C_q R_N^(2q)
       (q/T) sum_(|k|<T/q) |mu_hat_N(k)|^(2q).
```

For one squarefree block,

```text
average |cos x|^(2q)
  = binom(2q,q)/4^q
  asymp 1/sqrt(pi q).
```

Therefore, for `r` independent blocks,

```text
average |mu_hat_N(k)|^(2q)
  approximately (pi q)^(-r/2).
```

Taking `2q`-th roots gives

```text
M
  <= C_q 2^r q^(-r/(4q)).
```

The function `(log q)/q` is bounded, so no choice of fixed or growing `q`
turns the right side into an absolute constant. The exponential factor `2^r`
remains.

Thus all direct magnitude-moment arguments are structurally insufficient, not
merely quantitatively unoptimized.

## 5. Probabilistic concentration inequalities give the same loss

The measure `mu_N` is a convolution of independent block measures, so one may
apply Kolmogorov--Rogozin or Littlewood--Offord concentration inequalities.
For squarefree diffuse support these naturally give normalized concentration of
order at best

```text
Q(mu_N,delta) << 1/sqrt(r)
```

without additional arithmetic hypotheses. Multiplying by `R_N=2^r` again leaves
an exponentially large absolute count.

Stronger inverse concentration results return to additive structure of the
angle coefficients, which is the program already audited through affine cubes
and sparse inverse problems.

## 6. What remains of the harmonic route

A harmonic proof must exploit the signed, center-dependent sum

```text
sum_(1 <= |k| < T)
  c_k exp(-ik x_I) mu_hat_N(k),
```

rather than its absolute values or moments.

The exact factors `mu_hat_j(k)` are real for the symmetric block measures, but
their signs vary with `k`; the interval center supplies an additional phase.
A majorant can succeed only through substantial coherent cancellation across
frequencies.

However, proving such cancellation uniformly in the adversarial interval center
is essentially a short exponential-sum theorem for the complete divisor-angle
set. The product formula by itself does not control this phase cancellation.
A block-adapted Riesz product has the same issue: positivity forces a nonnegative
majorant, while its integral against `mu_N` still contains center-dependent
Fourier correlations that cannot be bounded through coefficient magnitudes.

## 7. Verdict

The exact convolution structure remains mathematically correct, but the proposed
large-sieve, resonance-counting, and Fourier-moment implementations do not yield
a uniform arc bound.

The failure is concrete:

```text
magnitude-only Fourier information is too weak even for an ideal independent
squarefree conductor.
```

The only surviving harmonic subroute is a phase-sensitive short exponential-sum
estimate. Unless a new arithmetic mechanism controls those phases, that subroute
is essentially a reformulation of the original local-counting problem rather
than a reduction.
