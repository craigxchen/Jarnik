# Fixed Jensen derivatives cannot detect a reciprocal-square cloud

The finite approximate-coefficient statement in
[the uniform reciprocal-tail note](critical_moment_uniform_reciprocal_tail.md)
is qualitative. Here is an exact obstruction to obtaining its tail bound
from a fixed-degree real-rooted Jensen derivative alone. For every fixed
`J>=1`, there are arbitrarily high-degree real-rooted products with

```text
max_j |y_j|<=1,     |sum_j y_j|<=2,     sum_j y_j^2<=4,
e_(2r)(y)=0          for 1<=r<=J,
```

whose reciprocal lists have a positive squared tail after **every** fixed
index `K`. The first `2J+1` coefficients, and thus the corresponding
degree-`2J+1` derivative polynomial after any fixed normalization, are
identical across products with different escaped tail masses.

This does not conflict with the qualitative theorem: its number of
coefficient conditions `J` must depend on the requested tail tolerance.
The examples do not satisfy the *full* critical odd-moment identities.
The construction does not rule out a bound `K(J)` at a tail tolerance
`epsilon(J)` that decreases sufficiently quickly as `J` grows; it only
rules out forcing a vanishing tail at each fixed `J` from the displayed
identities alone.

## Construction

Fix `d=2J+1`. Choose any degree-`d` real polynomial `E_0(t)` with
constant term one, simple nonzero real roots, and zero coefficients at
every positive even degree. Such a polynomial is obtained by taking

```text
E_0(t)=1-2 T_d(lambda t),
```

where `T_d` is the odd Chebyshev polynomial, and `lambda>0` is small.
The equation `T_d(x)=1/2` has `d` simple real nonzero solutions. Reduce
`lambda` until all reciprocal roots `y_i^0` have absolute value below
`1/2`, `|sum_i y_i^0|<1`, and `sum_i (y_i^0)^2<1`.

For even `n`, and a parameter `a>0`, define

```text
C_(n,a)(t)=(1-a t^2/n)^(n/2),
B_(n,a)(t)=[E_0(t)/C_(n,a)(t)]_(<=d),
E_(n,a)(t)=C_(n,a)(t) B_(n,a)(t),              (1)
```

where brackets mean truncation of a formal power series after degree
`d`. Thus `E_(n,a)(t)-E_0(t)` is divisible by `t^(d+1)`. In particular
all its even coefficients through degree `2J` vanish **exactly**, and
its first and second power sums are exactly those of `E_0`.

As `n` tends to infinity, coefficientwise,

```text
B_(n,a)(t) -> B_(infinity,a)(t)
               =[E_0(t) exp(a t^2/2)]_(<=d).          (2)
```

At `a=0`, the right side is `E_0`, with `d` simple nonzero real roots
and nonzero leading coefficient. These properties are open in the
real coefficient topology. Hence there is an `a_0>0` such that, for
each fixed `0<a<a_0` and all sufficiently large even `n`, the
polynomial `B_(n,a)` also has `d` simple nonzero real roots. Shrinking
`a_0` keeps their reciprocal magnitudes below one. The factor
`C_(n,a)` supplies `n/2` copies each of the reciprocal roots
`+sqrt(a/n)` and `-sqrt(a/n)`. Thus `E_(n,a)` is a real-rooted product
of degree `n+d` satisfying all displayed bounds. The bounds on its
first two power sums follow exactly from the coefficient matching,
rather than from an asymptotic estimate.

For large `n`, the `d` reciprocals from `B_(n,a)` are bounded away
from zero, while all `n` cloud reciprocals have magnitude
`sqrt(a/n)`. After sorting by decreasing magnitude, for any fixed
`K>=d`,

```text
sum_(j>K) y_j^2 = [n-(K-d)] a/n -> a.              (3)
```

For `K<d`, the lower bound is at least the cloud mass `a`. Therefore
no index `K` controls the squared tail below, say, `a/2` under this
fixed list of `J` exact coefficient conditions.

## Why derivative compression loses the tail

Let `N=n+d` and let `P_(n,a)(x)=x^N E_(n,a)(1/x)` have roots equal to
the reciprocal entries `y_j`. Its normalized `(N-d)`th derivative is
a degree-`d` real-rooted polynomial. In reversed form its elementary
coefficients are

```text
e_k^(derivative)=e_k(y) (d)_k/(N)_k,       0<=k<=d.   (4)
```

Equation (1) fixes every `e_k(y)` for `k<=d` to the corresponding
coefficient of `E_0`, independently of `a`. Therefore, for the same
`n`, the entire derivative polynomial in (4) is **identical** for every
`0<a<a_0` for which the construction is real-rooted. In particular
rescaling its roots by `N/d` and using their preserved first and
second power sums cannot recover the original tail mass in (3) or
infer a vanishing tail at this fixed `J`. Interlacing does not change
this information loss: many tiny
original roots can be compressed into a fixed number of macroscopic
derivative roots.

The construction uses repeated cloud magnitudes, as allowed in the
approximate-coefficient theorem. They can be made pairwise distinct
without changing the conclusion. Index the `n` cloud entries with
`j=1,...,n`, take `n/2` positive and `n/2` negative base values
`+/-sqrt(a/n)`, and add `j/n^3` to entry `j`. For large `n`, all
positive magnitudes are distinct and exceed `sqrt(a/n)`, while all
negative magnitudes are distinct and lie below it. The perturbations
sum to `O(1/n)`, so the cloud still has vanishing first power sum,
second power sum tending to `a`, and maximum magnitude `O(n^(-1/2))`.
Replace `C_(n,a)` by this distinct-magnitude cloud product and redefine
`B_(n,a)` by the same degree-`d` inverse truncation. The convergence
in (2) and exact low-coefficient matching persist. Its macro roots
remain simple and have distinct absolute values: those of `E_0` do,
because an odd polynomial minus a nonzero constant has no opposite
root pair. The cloud magnitudes are eventually smaller than every
macro reciprocal magnitude. Thus the full product can also have
simple real roots with pairwise distinct absolute reciprocal values.

An [exact rational degree-three fixture](check_critical_moment_finite_coefficient_gaussian_cloud.py)
checks the coefficient cancellation, the first two power sums, and
real-root intervals for one finite cloud. It is only a finite algebra
check; the open-root and limiting-cloud arguments above establish the
arbitrary-degree statement.
