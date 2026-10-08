# A tensor-product multiplicity floor at each core prime

The coefficient floor from
[Weyl-related residuals](higher_rank_weyl_union_floor.md) can count
one core block more than once. The valid multiplicity is at least
the minimum of its full and empty color-block counts. The proof
uses tensor products of evaluation kernels over a Gaussian DVR;
it does not multiply repeated divisibility assertions about one
scalar. For rectangular width `d=2` the local multiplicity is sharp,
and the summed floor equals the determinant-height numerator in
every rank.

## 1. Setup and statement

Let `Q` be an integral multilinear `SL_n` invariant on `m=nd`
rows, of relative `GL_n` weight `det^d`, satisfying the actual
conformal-block annihilator condition from
[the higher-rank coherent map](higher_rank_terminal_fusion_grid.md).
Thus every binary residual for every coordinate pair vanishes at
the actual binary rows. A single scalar Veronese evaluation zero
would not suffice.

Fix a nonzero source coefficient `c`, whose coordinate assignment
partitions the labels into blocks `B_1,...,B_n`, each of size `d`.
Use the primitive odd split Gaussian core factorization

```text
P_i=X_i+iY_i=K_i product_(T containing i) H_T,
n_T=Norm(H_T).
```

Core blocks and their conjugates have pairwise disjoint prime
support. Each `P_i` is coprime to its conjugate at the odd core
primes. For a subset `T`, let

```text
f_T = number of a with B_a subset T,
e_T = number of a with B_a intersection T empty,
r_T = min(f_T,e_T).
```

Choose `r_T` disjoint ordered pairs `(a_j,b_j)` of full and empty
colors. At an oriented core prime `pi` in `H_T`, put
`a_pi=v_pi(H_T)` and

```text
L_(pi,j)=sum_(i in B_(b_j)) v_pi(K_i).
```

Then the following exact local estimate holds:

```text
v_pi(c) >= sum_(j=1)^r_T max(0,a_pi-L_(pi,j))
        >= r_T*a_pi - sum_(j=1)^r_T L_(pi,j).            (1)
```

In particular, without a correction overlap at this prime,
`pi^(r_T*a_pi)` divides `c`.

## 2. Restricting to independent coordinate pairs

For each selected pair, set `O_j=B_(a_j) union B_(b_j)`.
These row sets are disjoint. On a row in `O_j`, allow only the
two coordinate axes of that pair, denoting their entries by `Y_i`
in color `a_j` and `P_i` in color `b_j`. Fix every remaining row
to the coordinate basis vector prescribed by its color block.
This produces a polynomial tensor

```text
R in W_1 tensor ... tensor W_r,                       (2)
```

where `W_j` is the free module of binary row-multilinear
polynomials on `O_j` with total degrees `(d,d)`. Its monomial
basis is `Y_J P_(O_j\J)` for `|J|=d`.
The degrees follow from the source's coordinate weight `det^d`.
The coefficient of
`product_j Y_(B_(a_j)) P_(B_(b_j))` in `R` is exactly `c`.
Restriction neither combines nor rescales the surviving source
monomials.

Let `ell_j:W_j -> Z[i]_(pi)` evaluate at the actual rows. Every
partial evaluation of (2) is zero:

```text
(id tensor ... tensor ell_j tensor ... tensor id)(R)=0. (3)
```

To check (3), expand in the other factors' monomials. Each
coefficient fixes a coordinate assignment outside `O_j`, and
the remaining coefficient is precisely a binary residual of `Q`.
It vanishes by the CB condition. This argument does not require
generic residues or invertibility of any outside bracket.

We deliberately use the full balanced coefficient modules `W_j`.
No integral invariant-basis conversion or saturation theorem for
an invariant lattice is needed.

## 3. The tensor-kernel lemma adds valuations

Work over the DVR `A=Z[i]_(pi)`. Set `K_j=ker ell_j`. The image
of `ell_j` is a principal ideal, hence a free `A`-module. Therefore
the exact sequence splits and `W_j=K_j direct-sum L_j`, with
`L_j` free of rank one. Tensoring these decompositions proves

```text
intersection_j ker(id tensor ... tensor ell_j tensor ... tensor id)
               = K_1 tensor ... tensor K_r.             (4)
```

Indeed a summand containing `L_j` maps injectively under that
partial evaluation, and the other factors are free. Thus (3)
places `R` in the right side of (4).

Let `lambda_j` extract the coefficient of
`Y_(B_(a_j)) P_(B_(b_j))`. For any element of `K_j`, all other
balanced monomials have a `P_i` factor from the full block
`B_(a_j)`, and are divisible by `pi^(a_pi)` on evaluation.
The selected monomial's multiplier has valuation `L_(pi,j)`:
the full-block `Y_i` are units, and the empty-block `P_i` have
only their correction factors at this prime. Consequently

```text
lambda_j(K_j) subset pi^max(0,a_pi-L_(pi,j)) A.          (5)
```

Applying `lambda_1 tensor ... tensor lambda_r` to (4) adds the
orders in (5). Its value on `R` is `c`, proving (1).
The intermediate tensor decomposition is allowed over a Gaussian
DVR; integrality of its factors over `Z` is unnecessary.

## 4. Accumulation over core blocks and the exact count

The final coefficient `c` is an ordinary integer. Therefore its
valuation at the conjugate split prime equals its valuation at
`pi`, and `pi^a | c` gives the ordinary divisor `Norm(pi)^a`.
Different core blocks have disjoint rational prime support, so
(1) accumulates to

```text
log C_Q >= log|c|
 >= sum_T r_T log n_T - 2 sum_i log|K_i|.              (6)
```

Here the correction bound is independent of `n`: at each core
prime the selected empty blocks are disjoint, so each row's
correction valuation is charged at most once. Summing over
distinct core supports costs at most `sum_i log Norm(K_i)`.
Negative lower bounds in the second line of (1) cause no issue;
the first line is the exact nonnegative valuation bound.

Under `log n_T >= (1-eta)w`, define

```text
M(n,d)=sum_(T subset [nd]) min(f_T,e_T).
```

Then (6) yields the strengthened coefficient floor

```text
log C_Q >= M(n,d)*(1-eta)*w - 2 sum_i log|K_i|.          (7)
```

Writing `q=2^d` gives an explicit finite formula:

```text
M(n,d)=sum_(f,e>=0,f+e<=n)
 n!/[f!e!(n-f-e)!] * (q-2)^(n-f-e) * min(f,e).         (8)
```

The previous union count `U(n,d)` replaces `min(f,e)` by
the indicator of `f,e>=1`, so `M>=U`. For `n<=3` there is
no change. At `n=4`, `M(4,d)=U(4,d)+6` for every `d`.

## 5. Sharp local model at rectangular width two

For `d=2`, define an invariant directly from actual binary rows:

```text
Q_P(v)=det [ (P_i v_i^a)_(a=1)^n | (Y_i v_i^a)_(a=1)^n ]
                   with rows i=1,...,2n.               (9)
```

It is row-multilinear and has relative weight `det^2`.
Subtracting `i` times each second-block column from its corresponding
first-block column replaces `P_i` by `X_i`. Thus (9) has ordinary
integer coefficients despite its Gaussian display.
After setting any two coordinate families to `(P_i,Y_i)`,
two columns coincide. Thus all its coherent binary residuals
vanish. A coefficient indexed by color pairs is, up to sign,
the product of the binary brackets within those pairs. Every
perfect matching occurs this way, so the content of (9) is
exactly the gcd of these matching products.

Suppose at one core prime the binary bracket valuations are

```text
v_p([i,j]) = a  if both i,j belong to T,
            0  otherwise.                             (10)
```

Actual primitive Gaussian rows with (10) exist: multiply the
rows in `T` by a fixed oriented Gaussian prime to the power `a`,
choose the remaining primitive row factors generically modulo
that prime, and keep the row directions distinct. The attached
checker gives explicit integer-coordinate examples.

A perfect matching has at least `max(0,|T|-n)` internal pairs
in `T`, and a matching attaining that number exists. Thus the
content of (9) has exactly this multiple of `a` as its `p`-order.
For the chosen color pairs, the unnormalized coefficient has
order `a*f_T`. Since `|T|-n=f_T-e_T` at `d=2`, the corresponding
coefficient of the primitive invariant has exact order

```text
a*(f_T-max(0,f_T-e_T))=a*min(f_T,e_T).                 (11)
```

This attains (1). Consequently no universally larger local
multiplicity, including one obtained solely from partial-block
counts, is possible at `d=2`. These are actual one-prime
arithmetic fixtures, not short-arc or near-uniform-profile models.

The [general local sharpness construction](local_tensor_coefficient_floor_sharpness.md)
extends this determinant by disjoint rainbow factors. It attains
`min(f_T,e_T)` for every occupancy pattern and every `n,d>=2`
over the formal collision DVR. It does not assert that the local
minima can all be attained by one global integer section.

## 6. The entire width-two family reaches equality

For `d=2`, the signed full-minus-empty count has generating
function

```text
(x+x^(-1)+2)^n=x^(-n)*(1+x)^(2n).
```

Using `min(f,e)=(f+e-|f-e|)/2` and the same central-binomial
identity as in the transpose note gives

```text
M(n,2)=n*4^(n-1)-n*binom(2n,n)/2=B(n,2).              (12)
```

Since `h(n,2)=1`, the stronger arithmetic floor meets the
determinant first-minimum exponent for every `n`. It produces
no strict exponent gain in that family. For `n=2`, `M(2,d)=2`
for every `d`, so the strict boundary comparison at `d>=3` in
[the transpose note](terminal_fusion_transpose_and_boundary_comparison.md)
is unchanged. An all-parameter comparison of `M` with `B/h`
for `n,d>=3` is not proved here.

The [checker](check_higher_rank_tensor_multiplicity_floor.py)
verifies occupancy counts, (12), and exact one-core valuation
sharpness for all subsets on four, six, and eight actual Gaussian
rows at two prime-power depths. Its finite checks supplement the
module proof; they do not assume endpoint geometry.
