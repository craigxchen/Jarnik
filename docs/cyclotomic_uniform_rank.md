# A uniform rank bound for all odd cyclotomic indices

The truncated Möbius matrix has a uniform nullity bound under the
endpoint support budget. No squarefree or bounded-prime-support
assumption is needed. Consequently the entire one-class cyclotomic
orientation model has a uniform finite count. This is a theorem about
that model; no reduction of arbitrary circle configurations to it is
asserted.

For an odd positive integer `h` and a finite set `E` of odd positive
integers, define

```text
M_E[d,e]=mu(e/d) if d|e, and 0 otherwise,   d odd, d<h,
D(E)=2*[1 in E]+sum_(e in E, e>1) phi(e).
```

Then

```text
D(E)<=4h  implies  dim_Q ker M_E<=1125.                (1)
```

The constant is not optimized. The
[nullity-four examples](cyclotomic_nullity_four_counterexample.md)
show that the earlier proposed bound of three was false. The new proof
compresses arbitrary prime-power supports to divisor downsets, then
uses an elementary intersection graph.

## 1. A prime-power compression preserving nullity

Fix an odd prime `p`. Slice the column support into finite sets

```text
E_j={s:p^j s in E, p does not divide s},   j>=0.
```

Use zero slices beyond the largest exponent. Write a coefficient
vector in corresponding slices `f_j(s)`, zero off `E_j`. For a row
index `p^j d`, with `p` not dividing `d`, its matrix equation is

```text
sum_(s:d|s) mu(s/d) [f_j(s)-f_(j+1)(s)]=0,
                                  d<h/p^j.             (2)
```

This is exact because the only nonzero values of `mu(p^t)` for
`t>=0` are `1,-1` at `t=0,1`.

Process the exponent levels `a` from largest down to `1`. Suppose
all higher levels are already prefix closed; in particular,

```text
E_(a+1) subset B,        B=intersection_(j=0)^a E_j.
```

Replace the support by

```text
E'_j=E_j union E_a    for j<a,
E'_a=B,
E'_j=E_j             for j>a.                         (3)
```

Here all right-hand sides use the support before this one step.

Let `V` and `V'` be the old and new rational kernels. For `f in V`,
put `v=f_a-f_(a+1)`, and define the linear transformation

```text
(Tf)_j=f_j-v     for j<=a,
(Tf)_j=f_j       for j>a.                              (4)
```

Because `E_(a+1) subset E_a`, the vector `v` is supported on `E_a`.
The transformed lower slices therefore fit (3), while
`(Tf)_a=f_(a+1)` is supported on `B`.
In every row block `j<a`, the difference in (2) is unchanged.
At `j=a` it becomes zero, and at `j>a` it is unchanged. Thus
`T(V) subset V'`.

The kernel of `T` on `V` is explicit. Its elements have one common
slice `v` at all levels `0,...,a`, zero slices above `a`, and

```text
support(v) subset B,
sum_(s:d|s) mu(s/d) v(s)=0 for d<h/p^a.                (5)
```

Call this space `W`. Every such vector also belongs to `V'`: all
differences below level `a` are zero, the level-`a` equation is (5),
and higher slices vanish. Moreover

```text
T(V) intersection W={0}.
```

Indeed every vector in `T(V)` has equal slices at levels `a,a+1`,
whereas a vector in `W` has those slices `v,0`. Therefore

```text
dim V' >= dim T(V)+dim W = dim V.                     (6)
```

This is why simply moving `p^a` to `p^(a-1)` is unnecessary and can
be incorrect: compression (3) can add several lower columns and
retains the kernel lost by the map (4).

## 2. The exact charged cost and termination

For each `s in E_a\B`, step (3) removes the column `p^a s` and adds
only previously absent columns among `s,ps,...,p^(a-1)s`. The removed
ordinary totient weight is

```text
phi(s)(p-1)p^(a-1),
```

whereas the added ordinary weight is at most

```text
phi(s) sum_(j=0)^(a-1) phi(p^j)=phi(s)p^(a-1).
```

It saves at least `(p-2)p^(a-1)phi(s)>=1` for each such `s`.
Columns corresponding to `s in B` need no additions, since they were
already present at every lower level. The only possible extra charge
in `D(E)` occurs if the base index `1` is newly added. This requires
`s=1` among the removed columns, and its ordinary-weight saving is at
least one, covering precisely the extra base charge. Consequently

```text
D(E')<=D(E).                                          (7)
```

After this step, the level-`a` columns occur at every lower level,
and higher prefix closure remains intact. Descending through all
levels therefore makes the support closed under division by `p`.
No new prime is introduced. If the support was already closed under
division by another prime `q`, every `E_j` is `q`-closed; unions and
intersections in (3) preserve that property. Processing the finitely
many primes of the original support yields a divisor-closed set `F`
with

```text
D(F)<=D(E),
dim ker M_F>=dim ker M_E.                             (8)
```

All operations are finite. The empty-support case is immediate.

## 3. Divisor downsets and their large orders

For divisor-closed `F`, a row indexed by `d notin F` is zero, since
`d|e in F` would imply `d in F`. The square submatrix on low indices
`d,e in F`, `d,e<h`, is upper triangular in increasing order and has
diagonal entries one. Hence

```text
dim ker M_F=#{n in F:n>=h}.                           (9)
```

Set `U=sum_(e in F)phi(e)<=D(F)<=4h`. Let `S` be the set of complex
roots of unity whose exact orders lie in `F`; it has cardinality `U`.
For every `n in F`, its full group `H_n` of `n`th roots is contained
in `S`, since `F` contains all divisors of `n`. Therefore

```text
n<=U<=4h,
|H_n intersection H_m|=gcd(n,m).                      (10)
```

Make a graph on the `r` indices `n in F` with `n>=h`, joining a
distinct pair when `gcd(n,m)>=2h/15`. For a fixed `n`, put
`g=gcd(n,m)`, `n=ag`, `m=bg`. The integers `a,b` are positive, odd,
coprime and at most `30`, by (10). There are only fifteen odd
positive integers at most `30`, hence at most `225` pairs `(a,b)`,
including `(1,1)` for `m=n`. Each pair determines `m=bn/a` uniquely.
Thus this graph has maximum degree at most `224`.

If `r>=1126`, greedy selection finds six independent vertices: each
selection removes at most `225` vertices, and five selections remove
at most `1125`. For their six orders, the first two terms of
inclusion-exclusion give

```text
|union_(i=1)^6 H_(n_i)|
 >=sum_i n_i-sum_(i<j) gcd(n_i,n_j)
 >6h-15*(2h/15)=4h.
```

This contradicts containment in `S`. Therefore `r<=1125`.
Equations (8)--(9) prove (1).

## 4. Uniform counts for the complete one-class model

Take a family of distinct one-class cyclotomic orientation polynomials
with common actual widths `W_e`, even base-layer counts, positive
degree `L=sum_e phi(e)W_e`, and minimum pair contact `h`, where
`L<=4h`. Every active nonbase width is at least one, and an active
base width is at least two. Thus its active support obeys
`D(E)<=L<=4h`, and (1) gives `rho=dim ker M_E<=1125`.

The [affine/obtuse whole-fiber theorem](cyclotomic_affine_rank_reduction.md#5-an-affine-dimension-bound-for-every-whole-polynomial-fibre)
then gives

```text
number of polynomials <=2rho+1<=2251.                 (11)
```

This permits arbitrary finite odd indices and arbitrary multiplicities.
If `L<4h`, its affine embedding has strictly negative pairwise inner
products. The vectors are then affinely independent, so the stronger
bound is

```text
number of polynomials <=rho+1<=1126.                 (12)
```

For actual growing fixed templates in one resonant Fibonacci class,
the [exact affine allocation bridge](cyclotomic_affine_rank_reduction.md#3-low-jets-and-the-actual-small-arc-count)
also gives at most `1126` distinct eventual points on an arc of length
at most `sqrt(2)*sqrt(R)`. For general fixed `C`, subdivision gives
the alternative bound

```text
number of points <=1126 max(1,ceil(C/sqrt(2))).         (13)
```

The whole-fiber bound (11) remains separately available; take the
smaller bound where useful. These are model counts independent of its
degree, rates, prime-power exponents, or number of primes in an index.
The fixed template's asymptotic entry threshold may still depend on
the template. Nonproportional affine classes and arbitrary lattice
circles have not been reduced to this one-class system. No improvement
to the established general growth bound in `R` is claimed.

The [prime-power compression checker](check_cyclotomic_prime_power_compression.py)
uses exact rational ranks on finite support families, including the
case in which a single-column downward move would fail. It supplements
the dimension and cost proof above; it is not its justification.

## 5. Compression for an arbitrary selected set of rows

The same compression has a more general exact linear-algebraic form.
Let `R` be any set of odd positive row indices, with no downward-closure
or interval assumption, and let `E` be a finite odd column support.
Define

```text
M_(R,E)[d,e]=mu(e/d) if d|e, and 0 otherwise,  d in R, e in E.
rho_R(E)=dim_Q ker M_(R,E).
```

Only finitely many rows can be nonzero. There is a finite divisor
downset `F` such that

```text
sum_(e in F)phi(e) <= sum_(e in E)phi(e),
D(F)<=D(E),
rho_R(E)<=rho_R(F)=|F\R|.                            (14)
```

To verify the extension, fix a prime `p` and apply exactly the support
operation (3) and map (4). Write the selected rows of block `j` as
`p^j d in R`, with `p` not dividing `d`. Their equations are

```text
sum_(s:d|s)mu(s/d)[f_j(s)-f_(j+1)(s)]=0,
                                       p^j d in R.  (15)
```

The transformation preserves every row expression in blocks `j!=a`
and sets the entire block `a` to zero. Thus it maps the old selected-row
kernel into the new one regardless of how rows were selected. Its
kernel consists of the common-slice vectors supported on `B` that
satisfy the selected block-`a` equations. Those vectors also lie in
the new selected-row kernel: all their other blocks are zero.
The image still has equal slices at `a,a+1`, whereas these embedded
kernel vectors have slices `v,0`. The intersection is zero, proving
the same nullity inequality (6). Neither this argument nor the cost
and termination proof uses the former cutoff `d<h`.

After full compression, selected rows outside `F` vanish. The full
square matrix on rows and columns `F` is triangular with diagonal
one. Any subset of its rows is independent, so the selected-row rank
is `|R intersection F|`, which proves (14).

This statement may be applied independently to several colored
supports with different selected row sets. Multiplying the cost of
each color by any fixed positive weight preserves the summed cost
inequality. It does not give a uniform bound for arbitrary selected
row sets: the complement `F\R` may contain many small indices.

For example, an overlap problem may select all low rows except those
divisible by a specified odd integer. Formula (14) applies to that
selection, but it does **not** preserve the values on the omitted
rows. The map explicitly zeros a full block of row expressions.
Consequently (14) alone cannot transport equations coupling omitted
rows of different colors, or replace them by separate zero equations.
