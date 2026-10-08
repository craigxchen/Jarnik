# Nonlinear consecutive Ptolemy test

This note tests whether geometric order can retain more information than
the symmetrized four-point Ptolemy sum. The candidate is the exact
log-sum-exp term in each consecutive window. It gives a nonlinear
strengthening of the usual `-2 log 2` estimate, but ordered gap geometry
does not make that strengthening grow with the conductor weight.

## 1. The candidate identity

Take points `z_1,...,z_M` in geometric order on an arc, and consider the
consecutive window `(i,i+1,i+2,i+3)`. Set

```text
A_i = |z_i-z_(i+1)| |z_(i+2)-z_(i+3)|,
B_i = |z_i-z_(i+3)| |z_(i+1)-z_(i+2)|,
C_i = |z_i-z_(i+2)| |z_(i+1)-z_(i+3)|.
```

Ptolemy gives `C_i=A_i+B_i`, with all three terms positive. Define

```text
rho_i = log(A_i/B_i),
Gamma_i = log(C_i^2/(A_i B_i)).
```

Then exactly

```text
Gamma_i = 2 log 2 + 2 log cosh(rho_i/2) >= 2 log 2.  (1)
```

In [the common-unit primitive chord notation](ordered_residue_growth.md), after canceling the common
positive factor of the three matching products, put

```text
Ahat_i = X_(1,i)t_i,i+1 t_i+2,i+3,
Bhat_i = X_(3,i)t_i,i+3 t_i+1,i+2,
Chat_i = X_(2,i)t_i,i+2 t_i+1,i+3.
```

The hats differ from the corresponding geometric products by the same
factor, so `rho_i` and `Gamma_i` are unchanged. Thus (1) is the exact
arithmetic identity

```text
2 log X_(2,i)-log X_(1,i)-log X_(3,i)
 +2u_i,i+2+2u_i+1,i+3
 -u_i,i+1-u_i+2,i+3-u_i,i+3-u_i+1,i+2
 = 2 log 2 + 2 log cosh(rho_i/2),                 (2)
```

where `u_ab=log t_ab`. The earlier ordered Ptolemy inequality retains
only `Gamma_i >= 2 log 2`; the new term is the nonnegative cosh term.
Summing (2) over consecutive windows preserves order, unlike the sum over
all seventy four-subsets.

## 2. What angular gaps say about the cosh term

Let `alpha_j=(theta_(j+1)-theta_j)/2` be the positive half-angle gaps.
For a window, direct chord formulas give

```text
B_i/A_i
 = sin(alpha_i+alpha_i+1+alpha_i+2) sin(alpha_i+1)
   / (sin(alpha_i) sin(alpha_i+2)).                  (3)
```

Consequently the ordered geometry determines `rho_i`, but (3) has no
radius factor. For equal small gaps, `B_i/A_i` tends to `3`, so

```text
Gamma_i tends to log(16/3)=1.673976...               (4)
```

This is a fixed positive constant. Unequal gaps can make the ratio large
or small. For an individual window, the lower bound in (1) is sharp as
the positive gaps vary. This pointwise observation does not establish an
optimal lower bound for a sum over many overlapping windows.

In particular, the extra term in (1) can only yield a conductor-sized
gain if one proves an additional arithmetic statement forcing many of the
ratios in (3), or the corresponding primitive residue ratios, to be
exponentially imbalanced. The chord formulas alone do not provide this.

## 3. An actual integer family with bounded nonlinear surplus

The absence of a radius factor is realized by actual Gaussian points. Let

```text
H_j=2(n+j)+i,                 0 <= j <= 3,
z_j=H_j product_(l != j) conjugate(H_l),
R_n=product_j |H_j|.
```

These are distinct equal-radius Gaussian points in geometric order (up to
reversal). Write `h_j=|H_j|`. Since

```text
|z_i-z_j| = 4 R_n |j-i|/(h_i h_j),                 (5)
```

the two non-diagonal products in the unique four-point window satisfy

```text
B_0/A_0 = 3,
C_0/A_0 = 4,
Gamma_0 = log(16/3),                               (6)
```

for every `n`. On the other hand,

```text
log(R_n^2) = 2 sum_j log h_j = 8 log n + O(1),      (7)
```

so the logarithmic radius tends to infinity while the nonlinear Ptolemy
surplus stays constant. The endpoint normalization of this particular
four-point family is `Delta_n sqrt(R_n) -> 12`; this is a test of the
scale-free ordered identity, not an endpoint counterexample at `C=1/2`.
It rules out a positive multiple of `W=log(R^2)` as a universal lower
bound for this surplus on arbitrary integer-circle configurations.

## 4. Scope of the test

The candidate identity (2) is strictly stronger than keeping only the
constant arithmetic-geometric-mean term, and it retains consecutive order.
The actual family above shows that its additional term is not forced to
grow with `W`, even on genuine integer circles. To turn it into the
desired balanced-layer inequality would require a new result linking
near-uniform conductor layers to exponential imbalance of the consecutive
Ptolemy terms or their primitive residues. No such link is proved here.
