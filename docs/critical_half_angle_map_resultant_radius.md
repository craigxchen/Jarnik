# Fixed critical half-angle maps: exact content and new norm factors

A fixed rational map can contract an angle to order `m>=2` near an actual
anchor. A coprime homogeneous polynomial representation of degree `d`
has `m<=d`. Its within-row primitive content is bounded by a fixed
resultant, so an individual output half-angle norm grows to degree `d`.
The old Gaussian conductor primes, however, generally disappear and new
ones take their place. The resultant does **not** bound the lcm across
different output rows. These are exact limits on a fixed-map radius
argument; no endpoint construction or uniform point bound follows.
For centrally extracted many-row tuples, the stronger all-row bound in
[fixed-degree rational-map inflation](fixed_degree_rational_map_inflation.md)
already excludes every non-pure **fixed-degree**, negligible-height
map once the row count exceeds twice its degree. The calculations here
isolate within-row content and source-prime overlap; they do not
supersede that theorem. A separate
[growing-degree coefficient bound](growing_degree_half_angle_height_barrier.md)
addresses maps whose degree varies with the tuple.

## Primitive output and its radius exponent

Let `F,G in Z[A,B]` be coprime homogeneous forms of the same degree
`d>=2`, with `F(1,0)=f_0 !=0`. Suppose `G(1,u)` has exact order
`m>=2` at `u=0`. For a primitive integer half-angle pair `(a,b)`, put

```text
u=b/a,  g=gcd(F(a,b),G(a,b)),
A'=F(a,b)/g,  B'=G(a,b)/g,
K=(A'+i B')/(1+i)^epsilon,
epsilon=1 if A',B' are both odd, and 0 otherwise.
```

The homogeneous resultant `rho=Res(F,G)` is a nonzero integer, and

```text
g | |rho|.                                             (1)
```

Indeed, the homogeneous resultant Bezout identities show that any
common divisor of the two values divides both `rho a^j` and
`rho b^j` for a fixed `j`; `gcd(a,b)=1` removes the powers. Equivalently,
at a prime not dividing `rho`, the two reductions have no common
projective zero. The parity division in `K` removes the only remaining
possible common Gaussian factor of `A'+iB'` and `A'-iB'`.

The primitive anchor norm is exactly

```text
P'=Norm(K)=(F(a,b)^2+G(a,b)^2)/(2^epsilon g^2).       (2)
```

The actual circle phase defined by the raw output pair is
`(F+iG)/(F-iG)=i^epsilon K/bar K`. The unit `i^epsilon` can vary by row;
it is retained when measuring angles, while it does not change the
Gaussian denominator `bar K` or its least-radius lcm.

For sufficiently small fixed `|u|`, constants depending only on
`F,G` give

```text
c_F (a^2+b^2)^d/|rho|^2 <= P' <= C_F (a^2+b^2)^d,
G(a,b)/F(a,b) = (gamma/f_0)u^m+O_F(u^(m+1)),
gamma !=0.                                            (3)
```

Thus the output circle angle, after subtracting its image at the
anchor, has exact local order `m`, while the individual anchor norm
has degree `d`. If a source tuple has squared radius `N`, an output
tuple retaining this anchor must have squared radius `N'` divisible
by `P'`. In the special case that the farthest source row has
`a^2+b^2 >= cN` and source angular span `Delta` comparable to
`N^(-1/4)`, (3) forces its output normalized span to be at least a
constant times `N^((d-m)/4)`. Since `m<=d`, there is no improving
power of `N` in this **single-row-dominates-conductor** case. The
condition is essential: a multirow source conductor can be much
larger than every individual half-angle norm.

For several output rows, after the displayed Gaussian primitivization,
the exact least squared radius relative to an anchor at `1` is

```text
N'_min=Norm(lcm_(Z[i]) {bar K_i}).                    (4)
```

Equivalently, the least physical radius is
`R'_min=|lcm_(Z[i]) {bar K_i}|`; the displayed norm is its square.

An actual target circle may have a larger radius. Formula (1) controls
each row separately; it says nothing about common factors among the
different `bar K_i` in (4).

## Exact source-prime overlap criterion

Let `D=a-ib` be the source denominator before its possible parity
division, and let `T=F(a,b)-iG(a,b)` be the raw output denominator.
Since `a=ib` modulo `D` and `gcd_(Z[i])(D,b)=1`, homogeneity gives

```text
gcd_(Z[i])(D,T) ~ gcd_(Z[i])(D,C_-),
gcd_(Z[i])(D,F(a,b)+iG(a,b)) ~ gcd_(Z[i])(D,C_+),
C_-=F(i,1)-iG(i,1),   C_+=F(i,1)+iG(i,1).           (5)
```

Here `~` means equality up to a Gaussian unit. If the corresponding
`C_±` is nonzero, only primes
dividing the corresponding **fixed** Gaussian integer can pass from a
source row's denominator to that same output row's denominator or
numerator. For the source numerator `a+ib`, the two constants are
`bar(C_+)` and `bar(C_-)` respectively. If `C_-=0`, the
linear form `A-iB` divides `F(A,B)-iG(A,B)` over `Z[i]`, giving
structural same-orientation persistence; `C_+=0` gives the crossed
orientation. Neither case controls newly created common factors among
different output rows. A statement that **all** old source norm primes
are erased from both target orientations requires both `C_-` and
`C_+` to be units (or at least nonzero for all but their fixed prime
divisors).

For the concrete critical degree-two map

```text
F(A,B)=A^2+B^2,    G(A,B)=B^2,
u -> u^2/(1+u^2),   rho=1,   C_-=-i, C_+=i,
```

the output integer pair is always primitive and (5) says that **every**
source Gaussian denominator factor is erased in its own output row.
Its new norm form is

```text
F^2+G^2=(A^2+B^2)^2+B^4.                            (6)
```

This map is locally injective on a one-sided positive interval; it
shows exact source-prime erasure, not a favorable endpoint-radius lcm.
The erasure occurs on literal circles, not just formal half-angle rows:
for any sufficiently large odd `a`, the two source points
`a-2i,a+2i` have radius `sqrt(a^2+4)` and angular separation
`2 arctan(2/a)`. Their primitive relative half-angle denominator is
`D=a-2i`. The mapped pair can be represented by the actual Gaussian
points `a^2+4-4i,a^2+4+4i`, at radius
`sqrt((a^2+4)^2+16)` and separation
`2 arctan(4/(a^2+4))`. Equation (5) gives `gcd(D,T)=1` for their
denominators. Both pairs lie on short arcs, although this two-point
family makes no claim about a large endpoint cluster.

## Why a critical map must create a norm factor

Suppose, conversely, that the output norm form uses only the source
norm divisor:

```text
F(A,B)^2+G(A,B)^2=c(A^2+B^2)^d,    c>0.             (7)
```

Over `C[A,B]`, the coprime conjugate forms `F+iG` and `F-iG` must
distribute the entire powers of the two linear factors `A+iB` and
`A-iB`. Conjugation then forces, for some nonzero constant `alpha`,

```text
F+iG=alpha(A+iB)^d   or   alpha(A-iB)^d.
```

On circle phases the induced map is therefore a fixed rotation times
`z^d` or `z^(-d)`. Its angular derivative has absolute value `d`,
so it has no critical order at the anchor. Every fixed critical map
must introduce a genuinely new norm factor, as in (6). This does not
prove that the new factor always raises the **multirow** least radius:
specialized output values may share new primes, and (4) is the
unresolved all-row cost.

[The exact checker](check_critical_half_angle_map_resultant_radius.py)
tests primitive content, Gaussian overlap, the erasure example, and
the norm-form classification for small power-map fixtures.
