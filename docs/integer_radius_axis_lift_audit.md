# Integer-radius axis lifts at the `R^(3/4)` scale

This note audits the proposed lift

```text
Z=(N/P) A^2,
```

where `P=Norm(A)` divides `N`. Since `|A|^2=P`, one has `|Z|=N`; if
`A=a+ib`, then `arg Z=2 arg A`. Thus a near-real `A` from a pair in a
fixed-unit endpoint cluster produces an integer lattice point `Z=X+iY` on
the circle `X^2+Y^2=N^2` inside an axis-centered arc at the `N^(3/4)` scale.
This is a useful reformulation, but I found no published
theorem giving an `O(log N)` bound for this special near-axis arc uniformly
in `N`.

For actual pair factors this lift has a simpler exact description:
`A^2/P=z_a/z_b`, and hence `Z=z_a conjugate(z_b)`. It is the pair
autocorrelation point, so its integrality is already present in the
source cluster. Conjugating selects the positive-angle representative.
The `binom(M,2)` pair classes give that many distinct points in the
one-sided arc of length at most `C N^(3/4)` from the positive axis.

## Exact endpoint equation

On the right-hand axis, write `X=N-h`, with `h>=1`. Then

```text
Y^2=h(2N-h).                                          (1)
```

If the original pair angular bound is `|arg Z|<=delta` with
`delta=C N^(-1/4)`, then the lifted arc has total length at most
`2N delta=2C N^(3/4)`. (If `C N^(3/4)` denotes the total arc length,
replace `C` below by `C/2`.) The exact relation between the pair factor and
the lifted endpoint is

```text
h=2N b^2/P,
b^2/P = sin^2(arg Z/2) <= sin^2(C N^(-1/4)/2).       (2)
```

The linear estimate `sin(u)<=u` is only the weaker necessary relaxation
`b^2<=C^2 P/(4 sqrt(N))`. It gives

```text
h <= (C^2/2) sqrt(N).                                (3)
```

Direct geometry gives exactly the same bound:
`h=N(1-cos(arg Z))<=N delta^2/2`. There is no stronger information in
the linear pair estimate than in this geometric consequence of the arc.

Thus the question is an integer factorization problem with
`h=O_C(sqrt(N))`,
namely whether `h(2N-h)` is a square. Factoring
`(N-X)(N+X)=Y^2` does not immediately give `O(log N)`: the gcd
`gcd(h,2N-h)=gcd(h,2N)` varies with `h`, and after removing it one obtains
primitive square representations whose count is controlled in general by
divisor/representation functions rather than a logarithm.

For odd N, these points have X odd and Y even. Conversely, every point
on the radius-N circle in that parity branch and in the positive-axis
first-quadrant arc has a unique parametrization of this form with
`a>b>0`, `gcd(a,b)=1`, a and b of opposite parity, and `P=a^2+b^2|N`.
To see this, let `g=gcd(X,Y)`. Then `g|N`, and the primitive Pythagorean
triple `(X/g,Y/g,N/g)` has its odd leg first. Its standard parametrization
is `(a^2-b^2,2ab,a^2+b^2)`, giving `g=N/P`.
Thus `P|N` and the sine condition (2) are intrinsic to this parity branch
of the lifted arc; they are not additional restrictions beyond it.

The selected near-critical pairs do have an extra restriction:

```text
N^(1/2-1/M) <= gcd(X,Y)=N/P <= (C^2/4)sqrt(N).
```

Moreover `gcd(h,2N-h)=2g` and the two factors after division by `2g`
are exactly `b^2` and `a^2`. These identities recover the original
divisor coordinates; they do not by themselves give a new spacing bound.
The pair points also retain their mutual source-cluster compatibility,
which is absent from an arbitrary collection in this arc.

## What the literature supplies

The relevant primary results are at different scales:

* Bourgain--Rudnick, [*On the Geometry of the Nodal Lines of Eigenfunctions of the Two-Dimensional Torus*](https://arxiv.org/pdf/1012.3843), Section 2, Lemma 4, recalls the Cilleruelo--Córdoba bound of at most m points on arcs shorter than `sqrt(2) R^(1/2-1/(4floor(m/2)+2))`. Choosing m proportional to log R and subdividing a fixed-C arc gives `O_C(log R)` at scale `R^(1/2)`. This argument does not reach the lifted `R^(3/4)` scale. Root checked the displayed theorem in the primary PDF; an earlier draft's citation to “Lemma 2.1” was incorrect.
* Oganesyan, [*Lattice points on small arcs*](https://arxiv.org/abs/2107.09991), is withdrawn (the arXiv record cites a crucial error in Lemma 2). Its former abstract claimed unbounded counts for every `alpha` in `(1/2,1)`, but that claim is not an available theorem and is not used here.
* The global circle count is `r_2(R^2)=4 prod_(p=1 mod 4)(2 v_p(R)+1)`, so the standard `R^epsilon` estimate is weaker than `O(log R)` and cannot be used here. The fact that the global count can exceed every fixed power of `log R` does not imply the same for the near-axis subset.

These sources supply no applicable `O(log N)` bound at the lifted scale.
This audit does not claim an exhaustive classification of what is known
about this specialized integer-radius, near-axis counting problem.

A later [audit of Chan's almost-square theorems](chan_almost_square_axis_audit.md)
found a valid constant bound for integer-radius points only in the much
shorter strip `|Y|<N^(1/2)(2log N)^(1/7)`. The lifted target has
`|Y|<=C N^(3/4)`. Its small primitive imaginary coordinate does not remove
the large squarefree coefficient contributed by `N/P` in the Pell equations.

## Pell sanity check

The exact checker [check_gaussian_near_real_divisor_reduction.py](/Users/cxc/Github/Jarnik/docs/check_gaussian_near_real_divisor_reduction.py) constructs the first Pell four-point block and its six pair factors. For Pell indices `1,11,21`, all six ordinary norms are distinct in each block, not only the five/six near-critical factors selected by the Plotkin threshold. At index `1` the six norms are

```text
1540145, 4032065, 103685, 20737, 806413, 308029.
```

The same checker verifies the nested layer-cake identity, conjugate
coprimality, quotient identities, and the near-real norm uniqueness test
through `d<=2000`. This is finite evidence for injectivity of the ordinary-
norm lift, not an `O(log R)` theorem.

The safe consequence for the main argument is thus the exact lift and the
distinct ordinary norms. If an `O(log N)` bound for the lifted points were
proved, the `binom(M,2)` pair count would imply `M=O(sqrt(log N))`. Obtaining
that improvement still requires a new upper bound for (1), potentially
using the near-critical gcd restriction or compatibility between pair
points. The exact sine condition (2) alone is simply the arc condition.
