# Global matching content, triangle content, and the retained circle center

Let `m>=5` distinct Gaussian integers `z_i` have common norm `N=R^2`
and Gaussian gcd one. Define

\[
g=\gcd_{i<j<k}|\det(z_j-z_i,z_k-z_i)|,\qquad
D=\gcd_{i<j}(z_j-z_i),
\]

and let `G` be the Gaussian gcd of all products of two disjoint
differences. Thus `G` is also the gcd of the quartet factors `gamma_Q`.
Gaussian gcds are understood up to units. Finally define

\[
H_i=\gcd_{j\ne i}z_j,\qquad
E=\operatorname{Norm}\left(\prod_i H_i\right).
\tag{1}
\]

The following are unconditional integer divisibilities:

\[
\boxed{E\mid N,\qquad \operatorname{Norm}D\mid g,\qquad
\operatorname{Norm}G\,\operatorname{Norm}D\mid4Eg^2,
\qquad \operatorname{Norm}G\mid2Eg^2.}\tag{2}
\]

The factor `E` is necessary in a local comparison: it records source
prime powers carried by all but one row. The exact local formulas below
separate those factors from accidental common differences.

There is also a strengthened triangle-product divisibility. Set

\[
K=\binom m3,\quad
b=\binom{\lfloor m/2\rfloor}3+\binom{\lceil m/2\rceil}3,
\quad d=\binom{m-1}3-b.
\]

Then

\[
\boxed{N^bE^d\mid\prod_{|I|=3}|D_I/g|.}\tag{3}
\]

On an arc of length `C sqrt(R)`, this implies, for every `m>=9`,

\[
\boxed{\operatorname{Norm}G
\le\frac{C^6}{32}R^{2\gamma_m},\qquad
\gamma_m=\begin{cases}
3/[2(m-1)],&m\text{ even},\\
3/(2m),&m\text{ odd}.
\end{cases}}\tag{4}
\]

Thus growing cardinality forces the full global matching content to be
subpower in the radius. This does not select one quartet with bounded
content and does not establish a uniform endpoint count.
The equality `Norm(G)=gcd_Q Gamma_Q` for circle tuples, proved in the
independent quartet audit, makes (4) also a bound on that integer gcd.

## 1. Exact source-prime formulas

Primitivity implies that `N` is odd with only split prime factors. Fix
`p^e || N`, choose `pi` above `p`, and sort the row levels:

\[
0=t_1\le a=t_2\le\cdots\le b'=t_{m-1}\le t_m=e.
\]

The endpoint levels follow from primitivity at both conjugate primes.
The `H_i` are pairwise Gaussian-coprime: a prime dividing `H_i` and
`H_j` would divide every row. Their total valuations at the two primes
are respectively `a` and `e-b'`. Therefore

\[
\boxed{v_p(E)=a+e-b'\le e.}\tag{5}
\]

For triangle content, there are two cases already proved in
`affine_shape_radius_divisibility.md`:

* With at least three distinct levels, `v_p(g)=0`.
* With exactly the two levels `0,e`, let `r` be the minimum extra
  valuation within either same-level class, beyond its common level.
  Then `v_p(g)=r`.

For completeness, the first case uses a triple at levels `0,t,e` with
`0<t<e`; its determinant has valuation zero. In the second, a mixed
triple has determinant valuation equal to its repeated-level pair's
extra valuation, and a minimizing pair attains `r`.

The matching gcd has these exact exponents:

\[
\begin{array}{c|c|c}
\text{row levels}&v_\pi(G)&v_{\bar\pi}(G)\\ \hline
a<e,\ b'>0&a&e-b'\\
0,e,\ldots,e&e+r&r\\
0,\ldots,0,e&r&e+r
\end{array}\tag{6}
\]

For the first row of the table, the quartet at levels `0,a,b',e` has
the matching `(0,b'),(a,e)`. Both edges join unequal levels, so the
matching attains the lower bounds `a,e-b'` simultaneously. The
remaining cases use a minimum-residue pair in the majority class,
a third majority point, and the outlier. In particular,

\[
v_p(\operatorname{Norm}G)=
\begin{cases}
v_p(E)+2v_p(g),&\text{an extreme class has one row},\\
v_p(E),&\text{otherwise}.
\end{cases}\tag{7}
\]

The first condition in (7) means the two-level singleton cases in (6).
In a configuration with three or more levels, (7)'s second line applies
even if the minimum or maximum is unique. Also `v_pi(D)=v_bar(pi)(D)=0`,
because a difference between a level-zero and a level-`e` row is a unit
at each relevant orientation. Hence `gcd(N,Norm(D))=1`.

## 2. Exact accidental-prime formulas, including two

Fix an odd prime `p` not dividing `N`. If it splits, choose a Gaussian
prime `pi` above it; if it is inert, use `pi=p`. Put

\[
\delta_{ij}=v_\pi(z_j-z_i),\qquad r=\min_{i<j}\delta_{ij}.
\]

Partition indices by `(z_i-z_1)/pi^r` modulo `pi`. There are at least
two classes. Define `tau` as the minimum sum of the three pair
valuations in a triangle, and `h` as the minimum sum over a pair of
disjoint edges. Their exact values are:

\[
\begin{array}{c|c|c}
\text{residue partition}&\tau&h\\ \hline
\text{at least three classes}&3r&2r\\
\text{two classes, both sizes at least two}&3r+u&2r\\
\text{one singleton and a majority}&3r+u&2r+u
\end{array}\tag{8}
\]

In the second row, `u>0` is the minimum excess `delta_ij-r` within
either class; in the third it is the minimum within the majority.
Three different classes attain the first triangle minimum. For two
classes, a minimum within-class pair and a point in the other class
attain the triangle minimum. The matching minima are the global-gcd
theorem proved in `quartet_matching_gcd_cut_budget.md`.

Since all rows are units at these primes, the circle determinant
identity gives

\[
v_p(g)=\tau,\qquad
v_p(\operatorname{Norm}D)=2r,\qquad
v_p(\operatorname{Norm}G)=2h.\tag{9}
\]

At split primes, equality of norms makes every difference's two
conjugate valuations equal. At inert primes the factor two in the
norm comes from `Norm(p)=p^2`. In all three rows of (8), `h+r<=tau`.
Thus these primes satisfy the stronger comparison

\[
v_p(\operatorname{Norm}G\operatorname{Norm}D)\le2v_p(g).
\tag{10}
\]

At the ramified prime `rho=1+i`, all rows are units because `N` is
odd. Their differences satisfy `r>=1`. Use the same definitions of
`r,tau,h`; the residue field has two elements, so only the last two
rows of (8) occur. The factor `2i` in the determinant identity now
has valuation two, while a rational integer's `rho` valuation is twice
its ordinary two-adic valuation. Therefore exactly

\[
\boxed{v_2(g)=(\tau-2)/2,\quad
v_2(\operatorname{Norm}D)=r,\quad
v_2(\operatorname{Norm}G)=h.}\tag{11}
\]

The same inequality `h+r<=tau` gives
`v_2(Norm(G)Norm(D))<=2v_2(g)+2`. Also `Norm(D)` is even because
`r>=1`. Combining (7), (10), and (11) proves the last two
divisibilities in (2). The divisibility `Norm(D)|g` follows directly
by dividing all differences by `D`: every triangle determinant is
multiplied by `Norm(D)` under Gaussian multiplication by `D`.

## 3. Retaining singleton layers in the triangle product

For the source-prime levels in Section 1, put `S_h={i:t_i>=h}`,
`1<=h<=e`. Equation (5) says exactly that

\[
v_p(E)=\#\{h:|S_h|=1\text{ or }m-1\}.\tag{12}
\]

The triangle valuation formula counts at least
`binom(|S_h|,3)+binom(m-|S_h|,3)` monochromatic triples at each layer.
Its minimum is `b`; at a singleton layer the count is exactly `b+d`.
With at least three distinct levels, `v_p(g)=0`, so immediately

\[
\sum_Iv_p(D_I/g)\ge eb+d\,v_p(E).
\tag{13}
\]

For exactly two levels, subtracting the common triangle valuation `r`
does not lose the layer contribution. In the notation
`s=|S_h|`, `E_s=binom(s,2)+binom(m-s,2)`, the exact earlier argument
gives

\[
\sum_Iv_p(D_I/g)\ge
e\left[\binom s3+\binom{m-s}3\right]
 +\bigl((m-2)E_s-K\bigr)r.
\]

The coefficient of `r` is nonnegative for `m>=5`. This proves (13)
also in that case. Multiplying the primewise divisibilities proves (3).

## 4. Endpoint consequence and its limit

Write `A=C^3/8`. Each triangle on the endpoint arc has absolute
determinant at most `A sqrt(R)`. Equation (3) gives the joint estimate

\[
\boxed{g^K E^d\le A^K R^{K/2-2b}=A^K R^{K\gamma_m}.}\tag{14}
\]

For `m>=9`, one has `2d>=K` (equality at `m=9`). For even `m=2n`,
`d=(n-1)^3`; for odd `m=2n+1`, `d=n(n-1)(2n-1)/2`. These formulas
verify the threshold directly. Since `E>=1`, raising (14) to `2/K`
gives

\[
Eg^2\le E^{2d/K}g^2\le A^2R^{2\gamma_m}.
\]

Combine this with (2) to obtain (4). More precisely,

\[
\operatorname{Norm}G\operatorname{Norm}D
\le4A^2R^{2\gamma_m}.
\]

For `5<=m<=8`, where `K>2d`, equation (14) and `g>=1` instead give
`Eg^2<=(A R^gamma_m)^(K/d)`. No asymptotic cardinality claim depends
on those four small cases.

The global gcd can be small even when every individual quartet content
is large, as the fair-cut example already demonstrates. Thus (4) does
not provide the quartet required for a uniform isolation argument.

## 5. A legitimate scaled circle, with its center denominator retained

The points

\[
u_i=(z_i-z_0)/D\in\mathbb Z[i]
\]

lie on the circle

\[
|u_i+c|^2=N/\operatorname{Norm}D,\qquad c=z_0/D.
\tag{15}
\]

This is a valid translated and scaled circle with radius `R/|D|`.
However `gcd(z_0,D)=1`: any common Gaussian prime would divide all
original rows. Thus the displayed Gaussian denominator is reduced.
An integer Gaussian multiplier `lambda` makes `lambda c` integral
if and only if `D|lambda`. Consequently every such multiplier has
`|lambda|>=|D|`; clearing the center denominator restores a radius at
least `R`.

There is no smaller nonintegral similarity multiplier hidden here:
the normalized differences `u_i` have Gaussian gcd one, so any complex
`lambda` with all `lambda u_i` integral must itself be a Gaussian
integer, by a Gaussian Bezout combination. A further translation
preserving the lattice points is integral because `u_0=0`, and cannot
remove the center's fractional part.

For an ordinary positive integer multiplier, write `D=a+ib` and
`c_D=gcd(|a|,|b|)`. The least possible multiplier is exactly

\[
\boxed{q=\operatorname{Norm}D/c_D\ge|D|.}\tag{16}
\]

Indeed `D|q` is equivalent to `Norm(D)|qa` and `Norm(D)|qb`.
The norm of `D` divides `g`, so this denominator is controlled by
triangle content, but it cannot be dropped. Also `G` divides matching
products, not necessarily individual differences; it supplies no
additional common similarity factor for all the points. These facts
prevent this normalization alone from becoming an origin-centered
radius descent.

The exact integer checker
`check_global_matching_triangle_content_normalization.py` verifies the
local formulas, both content divisibilities, the strengthened triangle
product, and the reduced center denominator on 1,200 primitive circle
tuples with squared radius at most 1000. These finite checks supplement
the proofs above.
