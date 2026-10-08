# One-coordinate cotangent mutations: exact costs and a scoped obstruction

This note tests arbitrary replacements and the natural circle monomial
mutation. The monomial has an exact valuation obstruction whenever the full
threshold-cut profile is available. The obstruction applies directly to
adding a point. Replacing a point can delete a prime supported only on that
point, so it needs a separate argument.

## 1. Exact replacement set and lcm update

For (L>0), write (Q_L(x,y)=(xy+L^2)/(x-y)). For a retained set (Y),
a possible replacement (x') is characterized exactly by

```text
x'-y | y^2+L^2       for every y in Y.                 (1)
```

This follows by reducing the numerator modulo (x'-y). Thus the candidate
set is finite: after choosing (y_0\in Y), every candidate has the form

```text
x'=y_0+d,   d | y_0^2+L^2,
```

followed by the checks (1) for the other retained coordinates.

If (N_Y) is the lcm of the anchor edges and retained-retained edges, then
the full conductor update is exactly

```text
N(Y union {x})  = lcm(N_Y, n(x), {n(Q_L(x,y)): y in Y}),
N(Y union {x'}) = lcm(N_Y, n(x'), {n(Q_L(x',y)): y in Y}).   (2)
```

There is no general monotonicity in (2).

## 2. A valuation lemma for affine monomials

Let (e_0,\ldots,e_m) be Gaussian-prime valuations after a common
projective normalization, and let

```text
e_* = sum_i c_i e_i,       c_i in Z,       sum_i c_i=1.     (3)
```

Suppose (c) is not a unit vector. Its positive support (S=\{i:c_i>0\})
has positive mass (P=\sum_{i\in S}c_i\ge2), while the negative mass is
(-(P-1)). Assume that one actual split-prime threshold layer separates
(S):

```text
g=min_(i in S)e_i - max_(j not in S)e_j > 0.              (4)
```

Let (E=\max_i e_i), attained at (r\in S). Since (c_r\ge1),

```text
e_* >= c_r E + (P-c_r) min_S(e) - (P-1) max_notS(e),
e_*-E >= (c_r-1)(E-min_S(e)) + (P-1)g >= (P-1)g > 0.    (5)
```

Therefore an added monomial point has valuation strictly beyond the old
range. If the old range is (s=E-\min_i e_i), the new all-edge lcm pays a
factor at least `p^((P-1)g)` at this split prime. The reverse threshold gives
the analogous valuation below the old minimum. This uses total valuation at
one prime and allows nested layers of that prime.

More generally, if (m=\min_i e_i), (M=\max_i e_i), and a new point has
valuation (e_*), the exact exponent change in the all-edge lcm at this
prime is

```text
dist(e_*, [m,M]) = max(0, m-e_*, e_*-M).                 (6)
```

The valuation width is invariant under common projective normalization, so
this conductor cost cannot be removed by a common Gaussian factor. The
upper-cut estimate gives `dist(e_*,[m,M]) >= (P-1)g`.

Summing over the unordered threshold cut `S|S^c`, define
`W_S=sum_p gamma_p(S) log p`, where

```text
gamma_p(S)=max(0,min_S e-max_(S^c) e)
          +max(0,min_(S^c) e-max_S e).
```

At most one summand is positive. This includes both orientations and is
independent of the choice of Gaussian prime above `p`. The prime-power argument
then gives the global exact-cost inequality

```text
log(N_added/N_old) >= (P-1) W_S.                         (7)
```

Nested layers of one prime are retained in `gamma_p(S)`.

Consequently, in a profile containing every threshold cut, every genuinely
new affine integer monomial has an upper or lower cut on its positive support
and cannot be added while preserving or reducing (N). The conclusion is
scoped to profiles where those threshold cuts are actually present; an
arbitrary finite clique need not contain them.

## 3. Replacement can delete a private prime

For `L=6`, exact positive cliques give

```text
C0=(8,9,10,12),       N(C0)=27625,
C1=(9,10,12,18),      N(C1)=5525,
C2=(10,12,18,27),     N(C2)=425.                         (9)
```

The first move replaces `8` by `18`, raises the minimum from `8` to
`9`, and lowers `N` by `5`. The second replaces `9` by `27`, raises
the minimum from `9` to `10`, and lowers `N` by `13`. It is the exact
circle monomial

```text
q(27)=q(10)/q(18),   q(x)=(x+6i)/(x-6i).                 (10)
```

The old `13`-factor is carried by the anchor edge of `9` and disappears
when that coordinate is removed. Thus a replacement theorem must retain the
full update (2); pair integrality does not imply `N'|N`.

For a replacement that deletes source row `l`, the full cut-allocation
bookkeeping can first adjoin the new rational point and then delete `l`. The
deletion can recover at most the source singleton/private weight `W_l`, since
adding a point cannot create a new private layer for a row that is then
deleted. Thus

```text
log(N_replacement/N_source) >= (P-1) W_S - W_l.          (8)
```

For fixed `P`, a full fair profile with every cut weight `w+o(w)` gives
`(P-2)w+o(w)`. More uniformly, if every cut weight is between
`(1-epsilon)w` and `(1+epsilon)w`, the lower bound is
`(P-2-P epsilon)w`. For all `P>=3` this is at least
`(P-2)(1-3epsilon)w`, so it is positive once `epsilon<1/3`, even if
the integer coefficients vary with the configuration. Thus positive mass
`P>=3` cannot shrink the radius in this regime. Mass two is the critical
replacement case: its coefficient vectors have the form `2e_i-e_j` or
`e_i+e_j-e_k`, with the indicated positive and negative indices distinct.

The chain then terminates sharply. `C2` has no fifth positive coordinate
compatible with all four entries, and no replacement that both raises its
minimum and has (N'\le425). Exact candidate sets from (1) are

```text
remove 10: {10,22,24,30,42}; raising candidates have N'=5525,7225,5525,2125;
remove 12: {12,78};
remove 18: {18};
remove 27: {6,8,9,14,27}.                                (11)
```

This is a finite obstruction to treating a successful replacement as an
iterable uniform descent.

## 4. Angle test and the minimal-radius inequalities

For the minimality argument assume explicitly `0<C<=1` and
`N^(1/4) Delta<=C`. In particular `Delta<=1`, and the source arc is
a B2 set by the product-rigidity theorem.

Write the source angles in increasing order as

```text
0=alpha_0<alpha_1<...<alpha_(k-1)=Delta,
```

where `Delta<=1`, and delete `alpha_l`. For retained indices `i,j`, the
reflection has the unwrapped angle

```text
beta=alpha_i+alpha_j-alpha_l.                            (12)
```

There is no modulo ambiguity at this scale. The reflected point lies inside
the old span exactly when

```text
alpha_l <= alpha_i+alpha_j <= alpha_l+Delta.             (13)
```

For an endpoint deletion this specializes to

```text
l=0:       alpha_i+alpha_j <= Delta,
l=k-1:     alpha_i+alpha_j >= Delta.                    (14)
```

For an interior deletion both inequalities in (13) are required. The
endpoint cases are strict-span cases whenever (13) holds: deleting the left
endpoint leaves a positive minimum, and deleting the right endpoint leaves
a maximum below `Delta`. Deleting an interior point leaves both endpoints,
so an inside reflection has exactly the old span.

For every triple, including endpoint deletions, the new span is at most
`2 Delta`. For an interior `l`, use `-alpha_l<=beta<=2 Delta-alpha_l` and
split according as `beta` is nonnegative. If `beta>=0`, the span is at most
`max(Delta,2 Delta-alpha_l)<=2 Delta`; if `beta<0`, it is at most
`Delta+alpha_l<=2 Delta`. If `l=0`, then `beta>=0` and `beta<=2 Delta`;
if `l=k-1`, then `-Delta<=beta<Delta`. These give the same bound.

Now choose, among primitive `k`-point tuples with
`N^(1/4) Delta<=C`, one with least `N`, and use the secondary least-span
choice at equal `N`. The B2 product rigidity assumption makes every
reflection with retained `i,j` a genuinely new point. Since

```text
N_new/N = Norm(d)/Norm(G),
```

the universal span bound implies

```text
N_new/N > 1/16.                                         (15)
```

Otherwise the reflected tuple still satisfies the same `C` constraint and
has smaller `N`. If (13) holds, its span is no larger than `Delta`, so

```text
N_new/N >= 1,
```

with strict inequality when the span strictly decreases. Thus (13)--(15)
give exact angle-aware inequalities, but only for the selected minimal
tuple. This extremal normalization is not asserted to survive a later
fair-profile extraction; the inequalities do not yet bound the number of
points.

Separately, the finite chain above illustrates inside-span replacement at
larger normalized constants. The cotangent clique includes the infinity
anchor at angle zero. The step `C1=(9,10,12,18)` to
`C2=(10,12,18,27)` deletes 9 and uses `q(27)=q(10)/q(18)`.
This is a singular reconstruction from retained coordinates; it is
**not** reflection with the deleted coordinate as denominator, since
`q(27)!=q(10)q(18)/q(9)`. Its norm ratio is `1/13`, and its anchored
span strictly decreases.

There is also a genuine inside-span reflection with smaller norm in C1:
delete 18 and use the infinity anchor and 10. It again adds 27, now
giving rational directions represented by `(9,10,12,27)` at L=6,
primitive norm 1105, and norm ratio `1/5`. These directions no longer
form an integral cotangent clique at L=6: the 9--27 edge has value
`-31/2`. Using L=12 and finite coordinates `(18,20,24,54)` restores
an integral clique and gives the same norm 1105. The anchored span
stays `[0,alpha(9)]`. Neither example is a tuple selected by the
minimality construction above.

The four-point Pell family is likewise compatible with these inequalities.
At its first audited index `(U,V)=(9,4)`, ordering by angle gives the anchor
followed by `X=7184,3723,3456`, with common cotangent denominator `L=24`.
Exact enumeration of all 24 choices of deleted row and unordered retained
anchor pair, including repeated anchors, gives eight inside-span
reflections. Their smallest norm ratio is **233**. The smallest ratio
over all 24 reflections is 89, attained only outside the old span.
The checker reconstructs the Pell points and verifies these counts and
ratios using rational Gaussian arithmetic, including the infinity anchor.
This is a finite illustration, not a claim that this first Pell index
satisfies `C<=1` or was chosen by the minimality construction. The
stronger whole-orbit radius floor for the family is proved in
[monomial_basis_valuation_content_floor.md](monomial_basis_valuation_content_floor.md).

## 5. Exact Gaussian reflection replacement identity

There is a useful exact formula for the special replacement

```text
w=z_i z_j/z_l,                                            (16)
```

inside a primitive tuple `z_0,...,z_m` of common Gaussian norm `N`, where
`i,j!=l`; the case `i=j` is allowed. Let
`G=gcd_G(z_k:k!=l)`, let `h=gcd_G(z_l,z_i z_j)`, and put

```text
d=z_l/h,       A=z_i z_j/h.
```

Primitivity gives `gcd_G(G,z_l)=1`. Hence `gcd_G(d,A)=1`, `G^2|A`, and
`gcd_G(dG,A)=G`. Clearing the denominator in (16) gives the raw tuple

```text
tilde z_k=d z_k (k!=l),       tilde z_l=A,
```

whose exact common Gaussian gcd is `G`. Therefore the new primitive tuple is

```text
z'_k=d z_k/G (k!=l),         z'_l=A/G,                    (17)
```

and its common norm is exactly

```text
N_new=N * Norm(d)/Norm(G).                                (18)
```

The reverse operation exposes the two costs. From (17),

```text
z'_i z'_j/z'_l=d^2 h/G.
```

Clearing this reverse denominator by `G` gives the tuple
`G z'_k=d z_k` and `G(z'_i z'_j/z'_l)=d z_l`. Its exact common gcd is `d`,
and dividing by `d` recovers the original primitive tuple. Thus the forward
reflection has private factor `G` and denominator `d`, while the reverse has
private factor `d` and denominator `G`, up to Gaussian units.

This identity retains every Gaussian prime power and is independent of any
archimedean estimate. It is the precise conductor bookkeeping for the
mass-two replacement `q_i q_j/q_l`.

## 6. The global reanchor reflection is different

For the least endpoint `A`, the exact reanchor map is

```text
F_A(x)=A+(A^2+L^2)/(x-A).                                (19)
```

It must be applied to the whole labelled tuple, including the infinity
anchor. It sends (q(x)) to (q(A)/q(x)), preserves the edge norms up to
orientation, and preserves (N) and the endpoint span. Applying it to only
one coordinate requires (Q_L(F_A(x),y)\in\mathbb Z) for every retained
(y), which is not implied by the old pair conditions: the old clique
controls angle differences, while the new conditions control
(\theta_A-\theta_x-\theta_y).

The values above were checked by exact divisor enumeration and the all-edge
reduced Gaussian norm formula. The valuation inequalities (3)--(5) are exact
prime-power statements and use no numerical asymptotics.
