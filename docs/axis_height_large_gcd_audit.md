# Large-gcd audit for the axis-height lift

This note audits the suggested strengthening of the integer-radius lift. The
notation is as follows. `N=R^2` is odd, `A=a+ib` is primitive with
`d=a^2+b^2 | N`, and `g=N/d`. The lifted point is

```text
Z = g A^2 = X+iY,
X=N-h,
h=2g b^2,
Y=2g ab.
```

For the good pairs in the endpoint cluster, `d` is distinct and lies in the
prescribed interval, while

```text
b^2 <= C^2 d/(4 sqrt(N)),   h <= (C^2/2) sqrt(N),
g = N/d >= N^(1/2-1/M).
```

The last inequality is a genuinely large common divisor, but by itself it
does not give a new counting mechanism.

## Exact gcd identity

Because `N` is odd, `d` is odd and a primitive pair `(a,b)` has opposite
parity. Hence

```text
gcd(a^2-b^2, 2ab)=1.
```

Consequently

```text
gcd(X,Y)=g,
gcd(h,2N-h)=gcd(2g b^2,2g a^2)=2g.                 (1)
```

The second equality uses `gcd(a,b)=1`. Dividing the endpoint equation by
the gcd in (1) gives

```text
h/(2g)=b^2,             (2N-h)/(2g)=a^2.             (2)
```

Thus a point with a large gcd is exactly a factorization

```text
N=g(a^2+b^2),  gcd(a,b)=1,
```

and the normalized height is not merely controlled but is exactly
`h/g=2b^2`. In particular, the proposed “gap in `h/g`” is just the gap
between twice consecutive integer squares. It supplies at most the already
known bound on `b`.

## What can be counted immediately

Let `G` be a lower bound for the gcd and `H` a height bound. From (2), for
each fixed `g`, the possible heights are indexed by integers

```text
1 <= b <= sqrt(H/(2g)),
```

and each `(g,b)` determines at most one `a`, since
`a^2=N/g-b^2`. Therefore

```text
#points <= sum_{g | N, g >= G} sqrt(H/(2g)).          (3)
```

This is an exact elementary bound. Inserting
`G ~ sqrt(N)` and `H ~ sqrt(N)`, each summand is `O(1)`, so (3) is
`O(tau(N))=N^{o(1)}`. This is still not an `O(log N)` estimate: the
standard divisor bound is larger than every fixed power of `log N` along
some integers. If the condition `g|N` is discarded entirely, the resulting
count is only the much cruder `O(sqrt(N))` bound. Restricting to the
good-pair range improves the individual `b` bound to
`b <= (C/2)N^(1/(2M))`, but (3) still contains the number of divisors `g`
in a short interval. No logarithmic bound follows from the size of `g`
alone.

For a fixed `b`, the remaining condition is

```text
d=a^2+b^2 | N,       g=N/d in the short interval.    (4)
```

This is precisely the original near-real divisor problem. The difference of
two candidate norms with the same `b` is
`a_1^2-a_2^2=(a_1-a_2)(a_1+a_2)`, but the divisibility of each norm by `N`
does not force a useful lower bound on this difference. Since `N` is odd,
all admissible `g` are odd, but this parity spacing still does not create a
useful multiplicative gap: distinct admissible divisors can be close in the
relevant interval, and no stronger spacing is supplied by the hypotheses.

## A useful reformulation, and its limitation

Equation (1) says that every lifted point has a canonical primitive endpoint
factorization

```text
(h/(2g), (2N-h)/(2g)) = (b^2,a^2).
```

Hence counting lifted points with `g >= G` and `h <= H` is equivalent to
counting primitive representations `N=g(a^2+b^2)` subject to
`g >= G` and `gb^2 <= H/2`. The lift has not converted the problem into a
general square-product equation with an extra gcd saving; it has exposed the
same Gaussian divisors in endpoint coordinates.

The only unconditional improvements available from this route are therefore
the exact gcd identity, the normalized-square identity (2), and bounds such
as (3). An `O(log N)` conclusion would require a new theorem controlling the
number of divisors `d | N` for which `d=a^2+b^2` with small `b` in the
narrow near-`sqrt(N)` interval. Standard divisor bounds or the fact
`gcd(X,Y) >= N^(1/2-1/M)` do not provide that theorem.

## Status of the proposed strengthening

The large-gcd observation is valid and useful for bookkeeping, but the
strengthened claim “large gcd plus small height forces an `O(log N)` count”
is not established by the stated hypotheses. Any proof must exploit
additional correlation among the distinct divisors `d` (for example a new
spacing or factorization lemma); treating the gcds as independently spaced
quantities is invalid. Conversely, the calculation above does not construct
a counterexample to an eventual specialized `O(log N)` theorem. It
identifies the exact remaining problem and prevents claiming that the gcd
alone solves it.
