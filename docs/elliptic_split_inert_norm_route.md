# Elliptic denominator splitting: an exact CM norm model and the rational gap

The restriction to Gaussian split primes is substantive for a rational
elliptic boundary model. It cannot, however, be imposed through a
field-insensitive assertion that elliptic denominator ideals must have
a positive proportion of inert prime weight. The explicit construction
below has quadratically growing, mutually independent denominator
ideals over `Q(i)` with no odd inert prime factors. Its failure to
descend together with the prescribed rational boundary points is exact.

This note does not decide whether the ten rational sequences attached
to `y^2=x^3-2`, `P=(3,5)` can simultaneously have inert parts of
height `o(N^2)` along a subsequence. It identifies that question
precisely rather than assuming a prime-distribution statement.

## 1. The exact rational question

On the proposed fixed rational curve write

```text
E_0: y^2=x^3-2,       P_0=(3,5),
nP_0=(A_n/B_n^2,C_n/B_n^3),       B_n>0,
gcd(A_n,B_n)=1,
K={0,1,3,7,-1,-3,-7,-4,-8,-10}.
```

At every good odd prime `p`, let `r_p` be the order of the reduction
of `P_0`, and put `a_p=v_p(B_(r_p))`. The exact denominator formula is

```text
v_p(B_(N-k))=0                           if r_p does not divide N-k,
v_p(B_(N-k))=a_p+v_p((N-k)/r_p)          otherwise.     (1)
```

This follows from reduction to the identity and the formal group;
the odd-prime equality, including arbitrary powers, is
[Silverman–Stange, Lemmas 4–5](https://arxiv.org/pdf/1001.5303).
We exclude the finitely many values `N=k`.

Remove the fixed primes dividing the bad reduction or any
`B_(|k-l|)` for distinct `k,l` in `K`. At every remaining prime at
most one shifted denominator is divisible by `p`: simultaneous
contacts would imply `r_p | k-l`. Consequently the total inert mass
in the ten disjoint contacts is exactly

```text
I(N)=sum_(p=3 mod4, outside the fixed set)
       sum_(k in K, r_p | N-k)
          (a_p+v_p((N-k)/r_p)) log p.                  (2)
```

Each fixed prime contributes only `O(log N)`. This includes the bad
prime three, without assigning it a good-reduction order: the formal
subgroup is open of finite index in the compact group `E_0(Q_p)`, so
some fixed positive multiple `rP_0` first enters it. For indices
divisible by `r`, the formal parameter valuation is
`v_p(B_r)+v_p(n/r)+O(1)`; other indices have no denominator there.
The fixed exceptional correction at two is absorbed in `O(1)`.
Thus removing that fixed set changes the total mass by `O(log N)`. The desired
obstruction for this fixed model would be the nontrivial assertion

```text
liminf_(N -> infinity) I(N)/N^2 > 0.                   (3)
```

If the liminf is zero, a subsequence makes the inert mass of every
one of the ten contacts `o(N^2)`, since the terms are nonnegative.
Ordinary primitive-divisor existence does not prove (3): it supplies
neither the required splitting type nor a quadratic amount of
logarithmic mass. The CM primitive-divisor theorem in
[Streng's paper](https://pub.math.leidenuniv.nl/~strengtc/cmeds.pdf)
has explicit split-prime consequences and does not supply this
inert-height lower bound. No prime-distribution assertion resolving
(3) is invoked here.

## 2. A simultaneous construction over Q(i)

Use instead the fixed CM curve and rational point

```text
E: y^2=x^3-2x,       P=(2,2),
[i](x,y)=(-x,iy),
Q_n=[n+i]P=nP+[i]P.                                   (4)
```

The point has infinite order. The good reductions have
`#E(F_3)=4`, `#E(F_5)=10`; prime-to-residue-characteristic torsion
injection forces the rational torsion order to divide two, while
`P` is not two-torsion. The exact counts are checked below.

Let `b_n` be the Gaussian denominator ideal of `Q_n`: at a Gaussian
prime `pfrak`, its exponent is

```text
v_pfrak(b_n)=max(0,-v_pfrak(x(Q_n))/2).                 (5)
```

Pole orders in an integral Weierstrass equation make these integers.
The denominator and formal-group conventions agree with
[Streng, Section 2](https://pub.math.leidenuniv.nl/~strengtc/cmeds.pdf).

For every odd inert rational prime `p=3 mod4`,

```text
v_(p)(b_n)=0       for every integer n.                (6)
```

Indeed, if `Q_n` reduced to the identity over `F_(p^2)`, its
Frobenius conjugate `[n-i]P` would also reduce to the identity.
Subtracting yields `[2i]P=O` in the reduction, hence `[2]P=O`.
But `P=(2,2)` has nonzero second coordinate at every odd prime,
and the curve has good reduction there. This is impossible.
The same proof works simultaneously for all `Q_(N-k)`, `k in K`.

In fact `b_n` is coprime to its conjugate away from two. If both
conjugate prime ideals divided `b_n`, then at one of them both
`Q_n` and `bar(Q_n)` would reduce to the identity. Their difference
is the fixed point `[2i]P`, whose coordinates have denominators
supported only at two:

```text
2P=(9/4,-21/8),       [2i]P=(-9/4,-21i/8).             (7)
```

At the ramified prime `1+i`, conjugation preserves valuation. If
`Q_n` is in the formal group, so is its conjugate, and the formal
depth of their difference is at least their common depth. The
parameter `z=-x/y` of `[2i]P` is `6i/7`, of `(1+i)`-valuation two.
Therefore

```text
v_(1+i)(b_n)<=2.                                      (8)
```

Thus even the ramified denominator part is uniformly bounded.

## 3. Exact quadratic norm heights and bounded common factors

For any real point `R=(x,y)` on `E`, direct addition with `[i]P`
gives

```text
x(R+[i]P)=[-2(x-2)(x+1)-4iy]/(x+2)^2,
|x(R+[i]P)|^2=4(x-1)^2/(x+2)^2.                       (9)
```

Real points satisfy `x>=sqrt(2)` or `-sqrt(2)<=x<=0`.
In particular `x+2>0`, and (9) implies the uniform bound

```text
|x(Q_n)| <= 4+3sqrt(2).                               (10)
```

It also covers `nP=O` separately, where `x(Q_n)=-2`.
Both complex embeddings satisfy the same bound. Use absolute
logarithmic height `h_x(Q)=h(x(Q))` and define the canonical normalization
explicitly by `hhat(Q)=(1/2) lim_(r->infinity)4^(-r)h_x([2^r]Q)`.
Then

```text
log Norm(b_n)=h(x(Q_n))+O(1)
            =2(n^2+1) hhat(P)+O(1).                  (11)
```

For the last equality, canonical height is invariant under `[i]`
and under conjugation. Its bilinear pairing obeys
`<P,[i]P>=<P,-[i]P>=-<P,[i]P>`, so the cross term is zero.
Equivalently `[n+i]` has degree `n^2+1`. The height normalization
and isogeny degree law are stated in
[Streng, Section 3 following (3.10)](https://pub.math.leidenuniv.nl/~strengtc/cmeds.pdf).
The uniform archimedean bound (10) makes the error `O(1)` here;
no elliptic-logarithm estimate is needed.

The ten denominator ideals also have bounded mutual common factors,
including opposite Gaussian orientations. At any prime where two
points are in the formal group, the smaller of their depths is at
most the depth of their difference. Apply this to

```text
Q_(N-k)-Q_(N-l)=(l-k)P,
Q_(N-k)-bar(Q_(N-l))=[l-k+2i]P.                        (12)
```

The right sides are fixed nonzero points. Thus both
`gcd(b_(N-k),b_(N-l))` and
`gcd(b_(N-k),bar(b_(N-l)))` divide fixed denominator ideals.
This bounds their full prime-power norms independently of `N`,
not just their sets of prime divisors.

Remove the bounded ramified and common divisors. For example,
divide each `b_j` by

```text
gcd(b_j,bar(b_j) product_(l!=j) b_l bar(b_l)).
```

The removed ideal has bounded norm by (7)–(12). The resulting
ideals are pairwise coprime and coprime to every conjugate, including
their own. Since `Z[i]` is principal, choose Gaussian generators
`H_(N,k)`. These actual integers have only odd split-prime support
and satisfy, uniformly for the ten fixed labels,

```text
log N(H_(N,k))=2((N-k)^2+1) hhat(P)+O(1)
              =2N^2 hhat(P)+O(N).                    (13)
```

This realizes independent Gaussian factors of equal leading norm
height from genuine elliptic denominator arithmetic. It does not
realize the full lattice-arc equations or make the moving elliptic
point rational.

## 4. Rational trace and descent lose the protection

The trace is rational, but

```text
Q_n+bar(Q_n)=2nP                                      (14)
```

can have inert denominators even when both summands are integral
at those primes. Already

```text
Q_2=((-26+168i)/17^2, (4452-3646i)/17^3),
4P=(12769/84^2,900271/84^3).                           (15)
```

Both `Q_2` and its conjugate are integral at three and seven, while
their sum has both inert primes in its denominator.

One can make the `Q_n` rational for a translated descent datum:

```text
c'(Q)=bar(Q)+[2i]P,
c'(Q_n)=Q_n,       c'([i]P)=[i]P.                     (16)
```

But the rational boundary points for this datum are `[i]P+kP`,
not the original `kP`. Subtracting them from `Q_N` recovers
`(N-k)P` and its ordinary rational denominator sequence. The old
boundary points `kP` are not fixed by (16). Thus this descent does
not preserve both the moving rational point and the protected
boundary contacts.

The conclusion is limited and concrete: split support and equal
elliptic norm heights alone do not force residue growth. A successful
obstruction must retain the rationality of the moving moduli point
and its boundary divisors, or some other actual lattice metric
condition. For the fixed rational test, the weighted rank-of-apparition
question (3) remains unresolved by these arguments. Uniformity over
varying curves, generators, and boundary sets would be a further
requirement even if that fixed test were settled.

The exact arithmetic in (4), (7)–(10), (12), and (15) is tested by
[check_elliptic_split_inert_norm_route.py](check_elliptic_split_inert_norm_route.py).

Independent audit: the root checked the inert-prime argument, ramified
formal depth, full prime-power trimming, height normalization, and
rational trace/descent limitation, and reran the exact checker.
The `fresh_algebraic` agent independently repeated the complete proof
audit, including the height normalization and bad-prime clarification,
and reran the checker successfully.
