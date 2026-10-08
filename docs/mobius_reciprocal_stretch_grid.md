# Clipped conductor gives a finite divisor set for reciprocal stretches

For a nonconformal triangular integer map at actual source and target
anchors, the number of corresponding lattice points satisfies

```text
H = 4 delta^2 N N' / Nclip^2 is a positive integer,
m <= 2 tau(H),                                      (1)
```

Here `delta>0` is the determinant of that triangular representative;
`N,N'` are the exact least squared circle radii. In particular full
clipping and `N'<=N` give `N'=N` and `m<=2 tau(4delta^2)`;
at determinant one this is `m<=6`. Here `tau` counts positive divisors.
This is a conditional map
obstruction, not a map-existence theorem or a uniform endpoint bound.

The proof uses integer projections in the two actual configurations.
It arose from the isometric-anchor question, but needs no virtual point,
rational isometric direction, or increase of either least radius.

## 1. Frames and the physical projection identity

Let `z_i,w_i` be primitive Gaussian integral physical tuples of norms
`N,N'`, with corresponding actual anchors `z_0,w_0`. Use rational phases
`q_i=z_i/z_0`, represented by half-angle rows `H_i/conjugate(H_i)`, with
`H_0=1`. Assume the map has integral representative

```text
M=((a,b),(0,d)),       a>0, d>0, delta=ad,
```

so `M H_0` also represents phase one in the actual target frame. The
inverse representative `J=adj(M)=((d,-b),(0,a))` has determinant `delta`
and likewise preserves this anchor. Ordinary primitive reduction is
allowed; all occurrences of `delta` refer to the chosen representative.
If the original triangular determinant is negative, conjugating the
whole target tuple changes `d` to `-d`, preserves its least radius and
all arc lengths, and gives the stated positive-determinant setup with
`delta=|ad|`.

Write the Gram entries of `M` as `P,Q,S`, and define

```text
T=P+S,       W=(P-S)+2iQ,       K=Norm(W)>0,
C=W z_0,    g=gcd_Z(Re C,Im C).
```

For each actual point, its normalized squared half-angle stretch is

```text
x_i = Norm(MH_i)/(delta Norm(H_i))
    = T/(2delta) + Re(conjugate(C) z_i)/(2delta N).   (2)
```

Thus the `x_i` belong to an affine arithmetic lattice of spacing
`Delta=g/(2delta N)`. Each value occurs at most twice: a fixed
nonzero linear projection intersects a circle in at most two points.

Apply the same construction to `J` at the actual target anchor `w_0`,
obtaining `W',C'=W'w_0,g'` and `Delta'=g'/(2delta N')`. Its normalized
squared stretch at the corresponding point is exactly `1/x_i`.
Consequently reciprocals lie in a second affine lattice, also with
multiplicity at most two. Both statements concern the original sampled
points, not intermediate points of their arcs.

The actual-anchor hypothesis is material. For a general raw matrix whose
first column is not horizontal, the inverse formula uses an output phase
reference which need not be a Gaussian integer. One may first triangularize,
but must use the determinant of that resulting integral representative.

## 2. Reciprocal lattice counting

If positive numbers `x_i` occupy an affine lattice of spacing `Delta`,
with multiplicity at most two, and their reciprocals occupy a lattice
of spacing `Delta'` with the same multiplicity bound, then for every `X>0`,

```text
#{i:x_i<=X} <= 2+2X/Delta,
#{i:x_i> X} <= 2+2/(X Delta').
```

Choose `X=sqrt(Delta/Delta')`. This gives the exact useful estimate

```text
m <= 4+4/sqrt(Delta Delta')
  = 4+8delta sqrt(N N'/(g g')).                     (3)
```

No bound on the location of either arc relative to a rational normal is
needed. Large stretch forces small inverse stretch, so the two projection
counts control opposite ends of the same set.

## 3. The clipped conductor divides each projection content

At an odd split prime `p=pi bar(pi)`, put

```text
U=a+d-ib,       V=a-d+ib,       W=V conjugate(U),
u=v_pi(U), u'=v_bar(pi)(U), v=v_pi(V), v'=v_bar(pi)(V).
```

The source coefficient interval has endpoints `v-u` and `u'-v'`.
Both endpoints lie in the larger interval

```text
[-(v'+u), v+u']=[-v_bar(pi)(W),v_pi(W)].
```

Let `e=v_p(N)` and `e_0=v_pi(z_0)`. The source phase interval is
`[-e_0,e-e_0]`. If its clipped width is `c_p`, the preceding containment
gives

```text
c_p <= min(v_pi(W)+e_0, v_bar(pi)(W)+e-e_0)
    = v_p(g).
```

This uses both Gaussian orientations, including primes dividing the
determinant or common coefficient contents. There is no source conductor
at two. Hence

```text
Nclip | g.                                          (4)
```

The inverse has exactly the same clipped conductor. Here is the precise
endpoint statement establishing this symmetry. Its coefficient interval
has endpoints `v-u'` and `u-v'`. Let `I,I'` denote the forward and inverse
intervals, and let `t,t'` be corresponding phase valuations. Then

```text
clip_I'(t') = F(clip_I(t)),                          (5)
```

where `F` is the affine isometry carrying the forward coefficient
endpoints to the inverse endpoints in the order prescribed by the map.
If `v-u >= u'-v'`, then `F(t)=u-u'+t`; otherwise
`F(t)=v-v'-t`. Outside `I`, the target valuation equals the corresponding
endpoint of `I'`. At an endpoint, any numerator or denominator cancellation
moves the valuation outward; clipping to `I'` removes it. If `I` is a
point, both clipped widths are zero. Thus (5) includes all endpoint
cancellations and both orientations.

Taking widths in (5) proves `Nclip_inverse=Nclip`. Applying (4) to the
actual-target-anchored inverse gives `Nclip|g'`. Substitution in (3)
gives the elementary bound

```text
m <= 4+8delta sqrt(N N')/Nclip.                     (6)
```

There is a stronger integer consequence. Put `Q=Nclip`, `E=N/Q` and
`E'=N'/Q`, all integers by the clipped-conductor divisibility theorem.
Equation (2), `Q|g`, and its inverse give positive integers

```text
n_i=2delta E x_i,       n_i'=2delta E'/x_i,
n_i n_i'=4delta^2 E E'=H.                           (7)
```

Indeed `n_i=TE+Re(conjugate(C)z_i)/Q` is integral, and likewise for
the inverse. Thus each `n_i` is a positive divisor of the same integer
`H`; each occurs at most twice by the physical projection argument.
This proves (1). The factor two in the denominator is retained; no
fixed half-angle parity is assumed.

## 4. Scope and relation to isometric anchors

If both isometric physical points happen to be integral on the existing
source circle and their coordinates share a divisor `G`, their difference
divided by `G` gives a short integer chord. For parabolic shear ratio
`q=B/A`, a primitive chord normal has length at most
`4sqrt(N)/(G sqrt(q^2+4))`. The stretch projection step is
`q sqrt(q^2+4)/(2sqrt(N) times normal_length)`. The analogous inverse
estimate holds on the actual target circle. This recovers a conditional
many-point obstruction from high-content virtual pairs, without any
assumption about the sampled arc's angular offset.

Equations (1)--(7) are stronger as a route: they directly use the clipped
conductor and do not require either virtual point to be adjoined. They
also apply when isometric directions are quadratic irrational.

In the critical triangular reduction, the currently available
`delta=N^O(1/m)` and `N/Nclip=N^O(1/m)` give `H=N^O(1/m)` up to
known constants. Equation (1) is stronger than the elementary bound (6)
and allows divisor estimates to improve conditional growth bounds.
It does not make `H` uniformly bounded. Small endpoint constants enter
the existing estimates of these quantities; no automatic map existence
or full endpoint theorem follows.

Using the imaginary projection coordinate as well gives integer square
equations. The [squareclass continuation](reciprocal_squareclass_runge_growth.md)
combines those equations with the small kernel to prove the conditional
bound `m=O(sqrt(log R/loglog R))`.

The [exact checker](check_mobius_reciprocal_stretch_grid.py) verifies the
physical projection formulas, reciprocal stretches, both Gaussian prime
orientations, equality of forward and inverse clipped conductors, the two
content divisibilities, the positive integer product (7), and the
resulting count inequalities.
