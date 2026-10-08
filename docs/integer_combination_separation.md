# Integer-combination separation and its exact limitation

This note proves a finite height bound for every nontrivial integer
combination of equal-unit circle points. It then constructs a relaxation
with arbitrarily many points satisfying all the corresponding imaginary
part lower bounds simultaneously, including their actual phases.
The relaxation uses complex factors, not Gaussian integers, and is not a
counterexample to the endpoint conjecture.

## 1. A finite inequality for arbitrary integer combinations

Use the split-prime representation

\[
z_i=\epsilon\prod_p\pi_p^{a_i(p)}
                 \overline{\pi_p}^{e_p-a_i(p)},\qquad
|z_i|=R,\qquad W=\log(R^2).
\]

Let c be an integer vector with sum_i c_i = 0, and put

\[
\gamma_p=\sum_i c_i a_i(p),\qquad
H(c)=\sum_p|\gamma_p|\log p,
\]

\[
h_c=\prod_{\gamma_p>0}\pi_p^{\gamma_p}
       \prod_{\gamma_p<0}\overline{\pi_p}^{-\gamma_p}.
\tag{1}
\]

Assume at least one gamma_p is nonzero. Then h_c is nonreal: it and its
conjugate are coprime and it is not a unit. Thus
t_c = |Im h_c| is a positive integer. Exactly,

\[
\prod_i z_i^{c_i}=\frac{h_c}{\overline{h_c}},\qquad
|h_c|=e^{H(c)/2}.
\tag{2}
\]

Choose argument lifts phi_i in an arc of angular width Delta. The zero
coefficient sum gives

\[
\left|\sum_i c_i\phi_i\right|
\leq\frac{\|c\|_1\Delta}{2}.
\]

Taking the sine of half the phase in (2) is valid without choosing a
small-angle branch. Consequently

\[
t_c=e^{H(c)/2}
\left|\sin\left(\frac12\sum_i c_i\phi_i\right)\right|
\leq e^{H(c)/2}\frac{\|c\|_1\Delta}{4}.
\tag{3}
\]

At Delta <= C/sqrt(R), this proves

\[
\boxed{H(c)\geq\frac W2+
 2\log\frac{4t_c}{C\|c\|_1}
 \geq\frac W2+2\log\frac4{C\|c\|_1}.}
\tag{4}
\]

This works with arbitrary prime powers. The absolute value is taken
after summing a prime's exponent contributions; replacing it by a sum
over threshold-layer absolute values would in general change H(c).
When all gamma_p vanish, (4) is not asserted.

For c = e_i-e_j, (3) and (4) recover the primitive chord bound. Thus
the following obstruction addresses even the extension to every integer
combination, not just the pairwise inequalities.

## 2. Complete half-cuts have a uniform lower bound

Fix any even M >= 2. Index B = binomial(M,M/2) columns by the subsets
S of size M/2. Give each column weight w = W/B and valuation entries
a_i(S) = 1 if i belongs to S and zero otherwise. Define

\[
H(c)=\frac WB\sum_{|S|=M/2}\left|\sum_{i\in S}c_i\right|.
\tag{5}
\]

For a nonzero integer c with sum c_i = 0, let
D = max_i c_i - min_i c_i. Then D >= 2 and

\[
\boxed{H(c)\geq\frac{WM D}{4(M-1)}.}
\tag{6}
\]

Choose indices attaining the maximum and minimum. A uniform half-subset
contains exactly one of them with probability M/(2(M-1)). Conditional
on this event and the other chosen indices, the two possible subset sums
differ by D. Their mean absolute value is at least D/2. All other terms
are nonnegative, proving (6).

In particular H(c) >= WM/(2(M-1)); equality is attained by a difference
of two coordinate vectors. Thus even all the arc-width consequences (4)
with t_c replaced by one are simultaneously compatible with arbitrary M
once W >= 4(M-1) max(0,log(2/C)).

## 3. Actual phases can satisfy all integer-combination bounds

Here is the stronger statement. Fix 0 < C <= 1 and choose

\[
\boxed{W\geq4(M-1)\left(M\log5+\log\frac{56}{C}\right).}
\tag{7}
\]

Put R = exp(W/2) and Delta = C exp(-W/4). There exist M distinct real
angles in an interval of width Delta such that, for every nonzero
integer c of coefficient sum zero,

\[
\boxed{e^{H(c)/2}
 \left|\sin\left(\frac12\sum_i c_i\phi_i\right)\right|\geq1.}
\tag{8}
\]

### Proof by a countable union estimate

Choose independent uniform X_i in [0,Delta]. For a fixed c and
epsilon = exp(-H(c)/2), condition on every coordinate except one with
c_i != 0. The phase sum then varies uniformly on an interval of length
|c_i| Delta. In each period of length 2 pi, the set where
|sin(t/2)| < epsilon has length 4 arcsin(epsilon) <= 2 pi epsilon.
The conditional probability of failure of (8) is at most

\[
\epsilon\left(1+\frac{4\pi}{|c_i|\Delta}\right)
\leq\frac{14\epsilon}{\Delta}.
\tag{9}
\]

Here Delta <= 1. There are at most (2D+1)^M integer vectors of range D
and sum zero: every coordinate lies between -D and D. For D >= 2,
2D+1 <= 5^(D-1). Set alpha = M/(8(M-1)). By (6), the union of all
bad events has probability at most

\[
\frac{14}{C}\sum_{D=2}^{\infty}
 5^{M(D-1)}e^{W/4-\alpha WD}
=\frac{14}{C}
 \frac{5^M e^{-W/(4(M-1))}}
      {1-5^M e^{-\alpha W}}.
\tag{10}
\]

Under (7), the numerator including 14/C is at most 1/4, and
5^M exp(-alpha W) <= C/56 <= 1/56. Thus (10) is at most 14/55 < 1.
The angles are distinct with probability one. An admissible choice
therefore exists. Subtracting their common mean gives centered angles
phi_i with the same width and unchanged zero-sum combinations.

## 4. The phases are compatible with the formal factor products

This relaxation can respect factor-phase identities as well. Let
T_{iS} = 2 a_i(S)-1. Its columns are balanced, and

\[
TT^{\mathsf T}=\frac{BM}{M-1}
 \left(I-\frac1M\mathbf1\mathbf1^{\mathsf T}\right).
\]

For the centered angle vector phi, set

\[
\theta=\frac{M-1}{BM}T^{\mathsf T}\phi,\qquad
g_S=\exp(w/2+i\theta_S).
\]

Then T theta = phi. The points

\[
Z_i=\prod_S g_S^{a_i(S)}\overline{g_S}^{1-a_i(S)}
\]

have radius R and arguments phi_i. They occupy an arc of length at most
C sqrt(R). For gamma_S = sum_i c_i a_i(S), form the reduced monomial

\[
G_c=\prod_{\gamma_S>0}g_S^{\gamma_S}
     \prod_{\gamma_S<0}\overline{g_S}^{-\gamma_S}.
\]

Its modulus is exp(H(c)/2), and its argument is
(sum_i c_i phi_i)/2 modulo 2 pi. Therefore (8) says

\[
|\operatorname{Im}G_c|\geq1
\quad\hbox{for every nonzero integer c with sum zero.}
\tag{11}
\]

Both the actual-angle identities and all these lower bounds hold
simultaneously. The construction applies for arbitrarily large even M.

### 4.1. The canonical factors cannot all be Gaussian integers

This particular realization has a stronger failure of integrality than
merely an unverified hypothesis. Suppose `M>=4`, `0<C<=1`, and `W>0`.
Put `rho=exp(W/(2B))` and `Delta=C exp(-W/4)`.
Each column of `T` has equally many positive and negative entries.
Since the `phi_i` occupy an interval of width at most `Delta`,

```text
|(T^T phi)_S| <= M Delta/2,
|theta_S| <= (M-1)Delta/(2B).
```

Consequently

```text
|Im g_S| <= rho |theta_S|
 <= C(M-1)/(2B) exp(-(B-2)W/(4B)) < 1.             (12)
```

Here `B=binom(M,M/2)>=M`, and `B>2`. Also `|theta_S|<1/2`.
If all factors were Gaussian integers, their imaginary parts would
vanish. The small argument interval would force every `theta_S=0`,
and `T theta=phi` would contradict the distinct point phases. Thus
at least one of these factors is provably not a Gaussian integer.

Even a common rotation cannot repair them. The factors
`exp(i alpha) g_S` have one common modulus, and their pairwise
distances satisfy

```text
|exp(i alpha)g_S-exp(i alpha)g_U|
 <= C(M-1)/B exp(-(B-2)W/(4B)) < 1 < sqrt(2).      (13)
```

If all rotated factors were Gaussian integers, equal-norm separation
would make them all equal. The small range of the original arguments
would then make all `theta_S` equal. But `T 1=0`, so this again
forces `phi=0`. Common rotation leaves the point products unchanged,
since each row of `T` also sums to zero.

### 4.2. The remaining phase freedom is substantial

The preceding calculation excludes the canonical choice and its
common rotations, not every compatible choice of factors. In the
zero-winding branch `T psi=phi`, another real phase vector has the
form `psi=theta+eta` with `T eta=0`. If all `rho exp(i psi_S)`
were Gaussian integers producing distinct row points, at least two
factors would differ. Equal-norm separation implies

```text
osc(eta) >= sqrt(2) exp(-W/(2B)) - (M-1)Delta/B.    (14)
```

Indeed, their chord is at least `sqrt(2)` and at most
`rho (osc(theta)+osc(eta))`. For fixed `M>2`, the leading scale in
(14) is much larger than `Delta` as `W` grows. Thus a Gaussian
realization cannot come from a perturbation of size comparable to
the short point phases. A nonconstant kernel component, or different
integer winding lifts, remains possible and has not been excluded.
Unequal factor moduli and prime logarithmic weights are additional
freedoms in the genuine arithmetic problem.

### 4.3. Complementary-pair rotations cannot repair the construction

The matrix contains both `S` and its complement `S^c`, and their
columns are negatives. Consequently every vector satisfying
`eta_S=eta_(S^c)` lies in the kernel of `T`. These are independent
common rotations of complementary factor pairs, a space of dimension
`B/2`, not merely one common rotation of all factors.

The canonical phases satisfy `theta_(S^c)=-theta_S`. Under any such
pair rotation, the two factors in a complementary pair retain distance

```text
2 rho |sin(theta_S)| <= 2 rho |theta_S| < 1.
```

If both were Gaussian integers, equal-norm separation would force
them to coincide. If this happened in every complementary pair,
every row product would equal `rho^B`: each pair contributes one
factor and the conjugate of that same factor, regardless of the row.
This contradicts distinct point phases. Thus no vector in this entire
complement-symmetric kernel subspace repairs Gaussian integrality.

For a general kernel vector, some complementary pair must instead
have distinct factors. Its separation gives the more specific bound

```text
max_S |eta_S-eta_(S^c)|
 >= sqrt(2) exp(-W/(2B)) - (M-1)Delta/B.           (15)
```

Only the complement-antisymmetric part of the kernel can provide
this difference. Its kernel dimension is exactly `B/2-M+1`: the
full antisymmetric coordinate space has dimension `B/2`, and the
restriction of `T` to it still has rank `M-1`, because `T` annihilates
the symmetric part. In particular this residual dimension is zero
for `M=4`, but is already five for `M=6` and grows rapidly.

This identifies the remaining zero-winding freedom more precisely.
It does not bound that freedom or exclude the other winding branches.

## 5. What this rules out, and what it does not

The canonical g_S cannot all be Gaussian integers, as Section 4.1
proves; the Z_i are not asserted to be Gaussian integers either.
The imaginary parts and the quantities in (11) need not be integers.
The weights are formal real weights, not asserted to be logarithms of
distinct rational primes. No lattice-point counterexample follows.

What is proved is that the Gaussian separation bound applied to all
integer combinations of the valuation rows, even with their compatible
actual phases, does not force a uniform endpoint count in this
relaxation. Merely adding more such lower bounds cannot close that route.
One must use arithmetic not retained here, such as integrality of the
original factors and other factor products, or stronger relations among
the integer residues. This is consistent with the earlier balanced
rational-phase obstructions, which use additional integral projections.

The unrestricted uniform exponent-1/2 theorem remains unproved.
