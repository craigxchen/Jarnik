# A quantitative gap between endpoint parameters with fixed central weights

Fix a nonresonant primitive central vector `u`, write `U=max |u_i|`, and
choose one of its two controlled rational hyperbolic pencils. The
[moving-frame estimates](six_point_moving_frame_radius_height.md) give
integer binary forms `Delta_u` of degree twelve, squarefree over `Q`, with
coefficient height `U^O(1)`. For a primitive integer parameter `(s,t)`, put
`H=max(|s|,|t|)`. An actual six-node circle configuration of bounded endpoint
scale satisfies

```text
d^2=Delta_u(s,t)>0,       0<|d| <= C_end U^a H,       (1)
```

where `a` is absolute. Throughout this note the endpoint constant
`C_end>0` is fixed. Constants denoted `C(C_end)` are allowed to depend on it.
This note proves a gap for solutions of (1); it does not bound their first
height uniformly in the moving weights.

## 1. Very close approximation to one of twelve algebraic parameters

There are absolute exponents `b,c` and a constant `C(C_end)>=1` such that,
whenever `H>=C(C_end) U^b`, every solution of (1) belongs to a branch chart
with a real root `alpha` of the dehomogenized `Delta_u` and

```text
0<|r-alpha| <= C(C_end) U^c H^(-10).                (2)
```

Here the chart is either `r=s/t` or `r=t/s`, chosen with the denominator
of absolute value `H`; thus `|r|<=1` and the fraction is reduced.
Representatives differing by a common sign give the same parameter.

To prove (2), work in either bounded chart. Its squarefree integer
polynomial `f(r)` has degree at most twelve and coefficient height
`U^O(1)`. Its distinct finite roots have separation at least `U^(-O(1))`:
use the nonzero integer discriminant, Cauchy's upper bound on root size,
and the discriminant product formula, as in the moving-frame note.
Let `rho=U^(-O(1))<=1` be a sufficiently small common separation bound.
If a real `r` is at least `rho/4` from every root, factorization and the
nonzero integer leading coefficient give `|f(r)|>=U^(-O(1))`.
For the solutions in (1), however,

```text
|f(r)|=d^2/H^12 <= C_end^2 U^(2a) H^(-10).          (3)
```

Above the stated polynomial threshold, `r` must therefore lie within
`rho/4` of a unique root `alpha`. This root is real: otherwise its
conjugate is a distinct root also within `rho/4` of the real number `r`,
contradicting separation. Every other root is at least `3rho/4` from `r`.
Factoring `f` once more gives

```text
|r-alpha| <= U^O(1) |f(r)|,
```

which proves (2). The strict lower inequality uses `d!=0`. Each of the
two charts has at most twelve roots; using twenty-four chart-root pairs
is a harmless uniform overcount, including parameter infinity through
the other chart.

## 2. An elementary ninth-power gap

Write `A=C(C_end) U^c`, increasing it to be at least one. Suppose two
distinct projective parameters obey (2) for the same chart and the same
root, and their heights are `H_1<=H_2`. Reduced rational separation gives

```text
1/(H_1 H_2) <= |r_1-r_2|
             <= A(H_1^(-10)+H_2^(-10))
             <= 2A H_1^(-10).
```

Consequently

```text
H_2 >= H_1^9/(2A).                                 (4)
```

The two choices `d` and `-d` at one projective parameter are not distinct
parameters in (4); they give the conjugate normalized circle tuple.
The gap does not rely on Roth's theorem or a rational-point finiteness
theorem.

For example enlarge the polynomial lower threshold to `H_*>=max(2,2A)`.
Then successive distinct solutions attached to one chart-root pair satisfy
`H_{j+1}>=H_j^8`. Thus the number of projective parameters with

```text
H_* <= H <= B
```

is at most

```text
24 [1 + max(0, log(log B / log H_*)/log 8)]          (5)
```

when `B>=H_*`; it is zero when `B<H_*`. Taking integer parts would improve
the displayed harmless real upper bound. More generally one can replace
eight by any fixed number less than nine by increasing `H_*` by a fixed
power of `A`.

The finitely many parameters below `H_*` number at most `O(H_*^2)`, so
this yields an effective bound `U^O(1)+O(log log B)` for the parameters
in one controlled pencil up to height `B`, for fixed `C_end`. There are
two pencils. Neither count is a count of lattice points on one circle.

## 3. A scoped consequence when the relevant branch has low degree

Suppose the nearby real root `alpha` lies on a rational irreducible
homogeneous factor `g` of `Delta_u` of degree `m`. Normalize `g` to be
primitive integral. Its coefficients have size `U^O(1)`, since all degrees
are bounded and an integer factor bound follows from the root bound and
Gauss's lemma. At an integer parameter satisfying (1), `g(s,t)` is a
nonzero integer: otherwise `Delta_u(s,t)=0`. On the chosen normalized
chart,

```text
H^(-m) <= |g(r,1)| <= U^O(1) |r-alpha|.             (6)
```

The upper bound follows by factoring out `r-alpha`, with every other
factor bounded above on `|r|<=1`. The same argument applies with variables
interchanged. Combining (2) and (6) gives

```text
H^(10-m) <= C(C_end) U^O(1).                        (7)
```

If `m<=9`, the parameter height is therefore polynomially bounded in
`U`. This conclusion is conditional on the degree of the branch factor.
The explicit nonresonant witnesses in the project have irreducible
degree-twelve branch forms, so this case does not cover arbitrary central
vectors.

The moving weights remain essential in both (4) and (7). Their polynomial
threshold and constants cannot be dropped, and arbitrarily high isolated
parameters are not excluded by the gap alone. In particular these results
do not change the general arc-count growth rate or prove the desired bound
independent of the radius.
