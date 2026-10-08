# The exact conductor cost of a secant involution seeded by four points

Two intersecting chords through four actual lattice points give a rational
involution of the original circle. This guarantees four integral partners,
not integral partners for the remaining points. The formulas below identify
the exact source conductor retained, the extra Gaussian denominator needed
for all partners, and the actual-anchor determinant. They provide a concrete
test for this proposed bridge to reciprocal divisor counting; they do not
prove that a suitable center exists for a large endpoint tuple.

## 1. Rational center and exact partner integrality

Let `z_i` be a primitive Gaussian integral tuple with `Norm(z_i)=N`.
Write the secant center in lowest ordinary rational coordinates as

```text
P=A/s,       A in Z[i], s>0, gcd(Re A,Im A,s)=1,
K=Norm(A)-Ns^2 != 0.
```

The second circle intersection on the line through `P,z_i` is

```text
w_i=P+[K/Norm(sz_i-A)](z_i-P)
   =(A conjugate(z_i)-Ns)/(s conjugate(z_i)-conjugate(A)).   (1)
```

This is an involution and `Norm(w_i)=N` exactly. Its denominator never
vanishes at a source point because `K!=0`. It is a rational projective
map in half-angle coordinates, rather than a new higher-degree
correspondence escaping the Möbius framework.

Put `h_i=sz_i-A=g_i v_i`, where `g_i>0` is ordinary coordinate content
and `v_i` is integer primitive, and set `S_i=Norm(v_i)`. Then

```text
w_i=(A g_i S_i+K v_i)/(s g_i S_i).
```

Its necessary and sufficient integrality conditions are

```text
g_i S_i | K,
A+[K/(g_i S_i)]v_i == 0 mod s, coordinatewise.        (2)
```

Necessity of the first condition follows by reducing modulo `g_i S_i`
and using the coprimality of the coordinates of `v_i`; the second then
gives exact sufficiency. The least positive rational integer clearing
this one partner is also exact:

```text
d_i = s g_i S_i /
 gcd(s g_i S_i, Re(A g_i S_i+K v_i), Im(A g_i S_i+K v_i)).    (3)
```

Intersecting chords through pairs `(z_0,z_1)` and `(z_2,z_3)` makes
`w_0=z_1,w_1=z_0,w_2=z_3,w_3=z_2`, so (2) holds for these four rows.
It imposes no immediate version of (2) on the other rows.

## 2. Source clipping is an explicit interval at every prime

Assume `A!=0`; the center zero gives the conformal antipodal map.
At an odd split source prime `p=pi bar(pi)`, write

```text
e=v_p(N), r=v_pi(A), r'=v_bar(pi)(A), j=v_p(s).
```

Use absolute physical allocation coordinates `t_i=v_pi(z_i)`. Because
the tuple is Gaussian primitive, their interval is exactly `[0,e]`.
Rewriting (1) on the circle gives

```text
w=N(A-sz)/(Ns-conjugate(A)z).
```

The two coefficient-interval endpoints in these coordinates are
`r-j` and `e+j-r'`. Actual source and target anchoring translate the
appropriate intervals and do not change clipped widths. Consequently

```text
v_p(Q)=length([0,e] intersect
             [min(r-j,e+j-r'),max(r-j,e+j-r')]),       (4)
Q=Nclip.
```

Both Gaussian orientations give the same width. In particular every
source prime not dividing `Norm(A)` is fully retained. More precisely,
let `A_red=A/gcd_G(A,s)`. Then

```text
N/Q | gcd(N,Norm(A_red)).                            (5)
```

Indeed the erased width is at most
`max(r-j,0)+max(r'-j,0)=v_p(Norm(A_red))`. This remains true when the
coefficient endpoints reverse order: then this bound is at least `e`
and the assertion is immediate. No squarefree hypothesis is used.

Thus the four-chord construction has a measurable retention benefit:
source primes outside the center numerator automatically survive. It
does not imply that the common divisor on the right side of (5) is small.

## 3. Exact least norms use a Gaussian lcm, without a squared denominator loss

Reduce the Gaussian fraction (1) separately for each row, writing
`w_i=n_i/b_i` with coprime Gaussian numerator and denominator. Let

```text
L=lcm_G(b_0,...,b_(m-1)),
G=gcd_G(Lw_0,...,Lw_(m-1)).
```

Units in these choices do not affect any norm. The exact least squared
circle norms of the union and of the target alone are

```text
N_union = N Norm(L),
N_target = N Norm(L)/Norm(G).                        (6)
```

For the union, source Gaussian primitiveness forces every common
similarity multiplier clearing the source to be a Gaussian integer.
Such a multiplier also clears the partners exactly when divisible by
`L`. Moreover the resulting full union has Gaussian gcd one: at every
Gaussian prime dividing `L`, a row with maximal reduced denominator
has `Lw_i` a unit there. This proves the first formula.

The same maximal-denominator argument proves `gcd_G(G,L)=1`. If even
one partner `w_0` is integral, as it is for a four-chord seed, then
`G|Lw_0` implies `G|w_0`. Therefore

```text
Norm(G)|N,
N_target=[N/Norm(G)] Norm(L),       Norm(L)|N_target. (7)
```

The target can discard only a factor of the old source norm after the
Gaussian lcm is paid. Formula (6) retains the exact orientation costs;
clearing the rational denominators (3) and charging their square would
generally overestimate them.

At a split prime `p` not dividing `NsK`, the cost is especially explicit.
Here each `h_i=sz_i-A` has at most one nonzero Gaussian orientation:
a rational factor of `h_i` would divide `K`. No numerator cancellation
in (1) is possible at an orientation dividing its denominator, since
the corresponding resultant is `K` up to a unit. Hence

```text
v_p(N_target)=max_i v_pi(h_i)+max_i v_bar(pi)(h_i).    (8)
```

The ordinary lcm of the `S_i` has at this prime the maximum of these
two maxima, rather than their sum. It divides the new target conductor
but can miss the cost of the other conjugate orientation. Inert primes
cannot divide a primitive `S_i`; the good-prime case at two has no such
denominator. Equation (8) is an exact split-prime assertion, not an
assumption that new primes cancel.

## 4. Actual-anchor matrix and the price of a chord intersection

Take the guaranteed paired point `z_0` as source anchor and `w_0=z_1`
as target anchor. Put

```text
alpha=Re(A conjugate(z_0)), beta=Im(A conjugate(z_0)),
D_0=Norm(sz_0-A), c=gcd(D_0,2s beta,K).
```

An exact primitive integral triangular representative is

```text
M = (1/c) ((D_0,-2s beta),(0,-K)),
|det M|=|K|D_0/c^2.                                 (9)
```

One verification starts with the half-angle involution
`((-beta,alpha+Ns),(alpha-Ns,beta))`. Multiplication on the output by
the conjugate of its first column gives
`N ((D_0,-2s beta),(0,-K))`; removal of `N` and the ordinary content
gives (9). The target reference in (9) is the actual integral partner,
so its determinant is valid for the reciprocal-stretch divisor theorem.

For intersecting chords, put

```text
u=z_1-z_0, v=z_3-z_2,
D=det(u,v) != 0, t=det(z_2-z_0,v),
A_raw=D z_0+t u, r0=gcd(Re A_raw,Im A_raw,D).
```

After a common sign adjustment, `s=|D|/r0` and `A=A_raw/r0`. Exact
power-of-a-point identities give

```text
D_0=t^2 Norm(u)/r0^2,
K=t(t-D) Norm(u)/r0^2.                              (10)
```

For an endpoint arc of length `L_arc<=C N^(1/4)` in the small-angle
range, all its chord directions differ by at most `L_arc/sqrt(N)`.
Thus

```text
|D|, |t|, |t-D| <= C^3 N^(1/4),
s <= C^3 N^(1/4)/r0,
D_0, |K| <= C^8 N/r0^2,
|det M| <= C^16 N^2/(r0^4 c^2).                     (11)
```

The small area bounds the rational center denominator, but (11) alone
does not make the actual determinant subpower or bound the Gaussian lcm
in (6) by a fixed small factor. That lcm still ranges over all primitive
distances `sz_i-A`, with its uncancelled new-prime part given by (8).
No claim that the lcm is unavoidably large on every endpoint tuple is
made; controlling it is the remaining construction problem.

## 5. Comparison with the existing geometric transforms

The [primitive-chord strip transform](endpoint_primitive_chord_transform.md)
keeps every original point integral in a dot-determinant sublattice and
obtains another positive-definite conic. The
[near-tangent Bezout reduction](near_tangent_bezout_reduction.md) identifies
that sublattice precisely and retains a large boundary term. Neither
construction pairs each original point with the other intersection of
a fixed secant, nor computes the target denominator ideal (6).

Here the source circle and all partner norms remain `N` as rational
points; integrality, rather than Euclidean circle membership, is the
missing condition. Equations (4)--(9) connect that gap to the source
conductor and to the determinant/least-radius parameters required by the
reciprocal divisor bound. They do not supply the missing global choice
of center or target endpoint estimate.

## 6. Every three-against-one seed cut is completely erased

Suppose the full source allocation at a prime `p` is binary, with values
`0,e`, and the four paired seeds have a `3|1` split between these values.
Then every one of the three finite-center pairings has

```text
v_p(Q)=0.                                           (12)
```

No binary or integral hypothesis is required for the remaining target
partners. This is stronger than applying a binary-source/binary-target
fiber lemma to a full image already assumed to be integral.

Here is the exact involution argument. In the absolute physical
allocation coordinates used in (4), the secant map and its inverse are
the same rational map. Their coefficient interval is consequently the
same interval `I`. The inverse-clipping identity says

```text
clip_I(t_i')=F(clip_I(t_i)),
```

where `F` is an affine isometry on `I`; endpoint cancellations disappear
under clipping. In any pairing of a `3|1` seed pattern, two members of
the majority source fiber have targets at opposite values `0,e`.
Their source clipped values agree, so their target clipped values agree.
Thus `clip_I(0)=clip_I(e)`, which is precisely (12).

If one pairing has parallel chords, its secant center is at infinity
and the pairing is a common circle reflection. Such a reflection acts
affinely with slope `-1` on physical allocations; it cannot permute a
`3|1` binary seed pattern. The center zero gives the antipodal map,
which preserves every allocation and likewise cannot pair a `3|1`
pattern. Hence whenever the following divisor is nontrivial, all three
pairings have finite nonzero centers and (12) applies.

Define the source-supported seed divisor

```text
B_seed=product p^e over binary source primes with a 3|1 seed split.
```

For each pairing, the exact consequence is

```text
B_seed | E=N/Q.                                    (13)
```

If that pairing were also an admissible nonconformal endpoint map with
`C,C'<=2` and actual `N'<=N`, the erased-conductor estimate would give

```text
B_seed <= E <= 2^(m q_m) N^q_m.
```

Thus `B_seed>=N^eta` requires

```text
eta <= q_m + m q_m log(2)/log(N).                   (14)
```

This is a simultaneous obstruction to all three maps from the same
four seeds. For fixed `eta>0`, their admissibility fails as both `m`
and `N` grow, since `q_m<=8/m`.

For illustration only, consider the **formal weighted allocation model**
giving equal total weight to every nonconstant binary sign pattern on
`m>4` labels. For any four selected labels, the fraction of its total
weight in `3|1` cuts is exactly

```text
2^(m-1)/(2^m-2) -> 1/2.                             (15)
```

Indeed exactly half of all binary patterns have an odd number of ones
on those four labels, and the two excluded constant patterns contribute
none. Thus such a model violates (14) for every choice of four seeds
once `m,N` are large. This formal allocation calculation does not
construct Gaussian primes with those weights or an actual endpoint
tuple. For actual tuples the usable hypothesis is the measurable
lower bound on `B_seed` itself. No other projective-map construction is
excluded by (13)--(15).

The [checker](check_fixed_secant_conductor_bridge.py) verifies paired
anchors, all remaining partner denominators, exact least union/target
norms, the triangular phase gauge, both prime orientations in (4),
and the retained and new-conductor formulas.
