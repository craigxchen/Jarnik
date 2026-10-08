# What projective rigidity still needs from a large cluster

The fixed-template theorem in `fresh_algebraic_parallel.md` is a new
restriction, but it does not yet improve the count's dependence on the
radius. This note tracks the remaining height issue. All arguments here
are prose proofs.

## 1. The upper bound supplied by the lattice itself

For four distinct Gaussian points on a common centered circle, put

```
P=(z_i-z_k)(z_j-z_l),
Q=(z_i-z_l)(z_j-z_k),
rho=P/Q.
```

The circle identity `conjugate(z)=R^2/z` shows that
`conjugate(rho)=rho`. Since `rho` belongs to `Q(i)`, it is rational.
Write `rho=a/b` in lowest terms, with `b>0`, and let
`H(rho)=max(|a|,b)`. The equation `bP=aQ` and integer Bezout identity
give a nonzero Gaussian integer `g` such that

```
P=a g,                  Q=b g.
```

Consequently, if the four points lie on an arc of length `L`,

```
H(rho) <= max(|P|,|Q|) <= L^2.                         (1)
```

For normalized arc width `C_arc=L/sqrt(R)`, this is

```
H(rho) <= C_arc^2 R.                                  (2)
```

The new projective theorem gives the complementary lower bound

```
H(rho) >= (1/(16 C_arc))^(2/9).                       (3)
```

Its five-point consequence similarly forces the maximum height of two
normalized cross-ratios to be at least
`(R/(2^56 C_arc^10))^(1/96)`.
These lower bounds leave a large interval below (2); combining them
with (2) does not give a uniform count.

## 2. Bounded-height extraction is a sufficient theorem, not compactness

Here is a precise sufficient statement:

> There exist fixed integers `m_0>=4` and `B>=1` such that every
> `m_0`-point lattice cluster on an arc of length at most
> `(1/2)sqrt(R)` contains four points with cross-ratio height at most `B`.

If this statement were proved, put

```
c_B=1/(16 B^(9/2)),
c=min(1/2,c_B/2).
```

An arc of length `c sqrt(R)` could not contain `m_0` points: the
extraction statement and the four-point projective lower bound would
give contradictory bounds on its normalized width. Subdivision would
then give, for every fixed `C>0`,

```
number of points <= (m_0-1) ceil(C/c).                 (4)
```

Thus bounded-height extraction would itself complete the requested
uniformity argument. It has not been proved.

Conversely, if uniformity fails, then for **every** fixed `m_0>=4` and
`B>=1` there is an actual endpoint cluster with at least `m_0` points
such that **every** four-point cross-ratio has height greater than `B`.
Indeed, start with arbitrarily large clusters at one fixed arc constant,
and select a consecutive `m_0`-point window. Its normalized span is at
most `(m_0-1) C/(M-m_0+1)`, which tends to zero. Once this is smaller
than both `1/2` and `c_B`, (3) applies to every quadruple in the window.

Accordingly, escape to large rational heights is not an optional case
that a hypothetical counterexample might avoid. It is forced by the
new projective theorem. A proof of bounded-height extraction must
exclude that escape using additional arithmetic.

## 3. Exact obstruction to a real-compactness extraction

There is an elementary family showing why bounded real cross-ratios
cannot supply the missing height bound. For any fixed number `m>=4`,
integer `t>=2`, and `1<=j<=m`, set

```
q_j(t)=j+1/(t+j).
```

These rational directions are distinct and ordered. Every four-point
cross-ratio converges to the corresponding finite, nondegenerate
cross-ratio of the integer grid `1,...,m`. Nevertheless every one of
their rational heights tends to infinity, uniformly over quadruples.

For a direct calculation, write `u_j=t+j` and

```
rho_0=((i-k)(j-l))/((i-l)(j-k))=a/b  in lowest terms,
A=(u_i u_k-1)(u_j u_l-1),
E=(u_i u_l-1)(u_j u_k-1).
```

The difference formula

```
q_i-q_j=(i-j)(1-1/(u_i u_j))
```

gives `rho(t)=(a/b) A/E`, while

```
A-E=(i-j)(l-k)=K != 0.                               (5)
```

If `g=gcd(aA,bE)`, then `g` divides the fixed nonzero integer `abK`.
The reduced denominator of `rho(t)` is therefore at least

```
|E|/(|aK|) >= t^4/(4 m^4).                           (6)
```

The last inequality is a deliberately loose bound: each factor in `E`
is at least `t^2-1`, while `|a|,|K|<m^2`. Thus no quadruple escapes
the height growth, although all the real cross-ratios converge.

These direction sets can be realized by actual equal-radius Gaussian
points. Put

```
A_j=j(t+j)+1,       B_j=t+j,       U_j=B_j+i A_j,
z_j=U_j product_(l!=j) conjugate(U_l).
```

The points have common radius `product_j |U_j|`, and their complex
cross-ratios are those of the rational parameters `q_j=A_j/B_j`.
A rational projective compression of the parameters can also make
their angular spans arbitrarily small, while preserving all
cross-ratios. This does **not** make an endpoint counterexample:
denominator clearing increases the radius, and this construction does
not keep `L/sqrt(R)` bounded. It shows precisely which hypotheses
cannot justify extraction: rationality, concyclicity, cyclic order,
small angular span, and real convergence of projective invariants
do not by themselves bound rational height.

## 4. Remaining gap

The fixed-projective-template theorem rules out five-point endpoint
families whose projective invariants remain fixed, and gives explicit
height growth when the invariants vary. It does not show that a large
cluster has any invariant of small enough height. Neither the general
lattice upper bound (2) nor compactness of real cross-ratio values fills
that gap.

The next useful result would need an arithmetic upper bound or an
extraction principle for rational projective heights that is uniform in
`R`, or a different argument that handles their necessary divergence.
No such result is claimed here.
