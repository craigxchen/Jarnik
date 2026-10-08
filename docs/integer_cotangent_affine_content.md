# Affine triangle content in an integer cotangent clique

This note gives an exact bound for the common content of the primitive
triangle determinants attached to an integer cotangent clique.  It also
records two direct bounds for the Gaussian content of the universal rows.
These bounds do not prove the endpoint height target.

Let

```text
H_i=X_i+iL,                     n_i=X_i^2+L^2,
W_0=product_i conjugate(H_i),
W_i=H_i product_(j!=i) conjugate(H_j),
```

where `L>0`, there are at least two distinct integers `X_i`, and every

```text
Q_ij=(X_i X_j+L^2)/(X_i-X_j)
```

is an integer.  Let `G` be a Gaussian gcd of the rows and put
`z_i=W_i/G`.  The `z_i` form the primitive integral circle tuple.  Write

```text
N=|z_i|^2,       D_abc=det(z_b-z_a,z_c-z_a),
g=gcd_(a<b<c) |D_abc|.
```

Thus `N` is the squared primitive radius.  The main conclusion is

```text
g divides 4L^3,                 gcd(g,N) divides L.       (1)
```

The second assertion is substantially sharper at primes which divide
the radius.  The power `L^3` in the first assertion allows primes absent
from the radius, where all three edges of a triangle can contribute.

The later [intrinsic-scale identity](intrinsic_cotangent_scale_triangle_divisibility.md)
separates these contributions exactly: every triangle half-determinant
is its Gaussian content norm times the product of its three reduced
cotangent denominators, with an explicit factor one or two from parity.
It identifies the minimal clearing scale without proving a new endpoint
height exponent.

## 1. Equal-allocation cancellation is bounded by the edge residue

Fix a Gaussian prime `pi`.  Put `f_a=v_pi(z_a)`.  For any edge `a,b`,
after choosing the edge orientation, the integer cotangent condition gives

```text
z_b/z_a=(Q_ab+iL)/(Q_ab-iL),                         (2)
```

where `Q_0i=X_i`; reversing orientation inverts the displayed quotient.
There is no additional independent unit factor. Reduce the fraction in
`Z[i]` by the same Gaussian gcd in numerator and denominator. If `f_a!=f_b`, the
two summands in `z_b-z_a` have unequal `pi`-orders, and hence

```text
v_pi(z_b-z_a)=min(f_a,f_b).                          (3)
```

If `f_a=f_b`, the reduced numerator and denominator in (2) are both
`pi`-units.  Their difference is obtained from `2iL` by division by
their Gaussian gcd.  Consequently

```text
0 <= v_pi(z_b-z_a)-f_a <= v_pi(2iL).                 (4)
```

The upper bound in (4) is where integrality of `Q_ab` is used.  It is
false for an unrestricted collection of rational circle parameters
before their pair-cotangent denominators are cleared.

## 2. Proof of the triangle-content bound

The primitive norm `N` is odd and has no inert prime divisor.  Indeed,
a ramified or inert prime dividing the common norm would divide every
`z_a`, contrary to Gaussian primitivity.

First let `p` be an odd split prime and write `p=pi conjugate(pi)` and
`e=v_p(N)`.  Then

```text
0<=f_a<=e,       min_a f_a=0,       max_a f_a=e.
```

For every triangle the standard circle identity is

```text
2i D_abc
 = N (z_b-z_a)(z_c-z_a)(z_b-z_c)/(z_a z_b z_c),      (5)
```

up to the orientation sign.  If `e>0`, choose a triangle containing a
global minimum-allocation row and a global maximum-allocation row.  In
(5) the baseline order is

```text
e-(maximum allocation-minimum allocation)=0.
```

If the third allocation is strictly intermediate, all three difference
orders are the minima from (3).  If it equals an endpoint, exactly one
pair has equal allocation.  Equation (4) therefore gives the sharper
bound

```text
v_p(g)<=v_p(L)                  when p divides N.      (6)
```

If `e=0`, all allocations are zero.  Applying (4) to the three edges in
(5) gives

```text
v_p(g)<=3v_p(L).                                       (7)
```

The same proof gives (7) for an odd inert prime: it is absent from `N`,
and every `z_a` is a unit at that prime.

Finally take `pi=1+i`.  Again `N` is odd, so all rows are `pi`-units.
Here

```text
v_pi(2iL)=2+2v_2(L),       v_pi(2iD)=2+2v_2(D).
```

The three applications of (4) in (5) yield

```text
v_2(g)<=3v_2(L)+2.                                      (8)
```

Equations (6)--(8) prove (1).  They retain all split, inert, and
ramified primes and all possible equal-allocation cancellations.

## 3. Exact relation with the universal-row content

For anchored triangles, direct expansion gives

```text
|det(W_i-W_0,W_j-W_0)|
 =4L^3 |X_i-X_j| product_(ell!=i,j) n_ell.             (9)
```

Every raw determinant is `N(G)` times its primitive counterpart.
Moreover the anchored triangle determinants generate the same integer
ideal as all triangle determinants, since every unanchored determinant
is an integral signed sum of three anchored ones.  Therefore (9) gives
the exact identity

```text
gcd_(i<j) [4L^3 |X_i-X_j| product_(ell!=i,j)n_ell]
 = N(G) g.                                             (10)
```

In particular `N(G)` divides the explicit gcd on the left.  By (1),
the quotient of that scalar gcd by the actual Gaussian-content norm is
a divisor of `4L^3`; its common part with the squared primitive radius
is already a divisor of `L`.

There is also a useful bound from one selected pair.  At `pi`, put
`alpha_i=v_pi(conjugate(H_i))` and `t=v_pi(2iL)`.  The row ideal gives
the exact formula

```text
v_pi(G)=sum_i alpha_i-max(0,max_i alpha_i-t).           (11)
```

For every pair `i,j`, (11) implies

```text
v_pi(G)<=t+min(alpha_i,alpha_j)
             +sum_(ell!=i,j)alpha_ell.
```

Hence

```text
G divides (2iL) gcd_G(conjugate(H_i),conjugate(H_j))
                  product_(ell!=i,j)conjugate(H_ell).  (12)
```

The exact pair-gcd formula gives

```text
N(gcd_G(H_i,H_j))=|X_i-X_j| d_ij,
d_ij=gcd((X_i^2+L^2)/(X_i-X_j),X_i,L,X_i-X_j),
d_ij divides L.
```

Taking norms in (12),

```text
N(G)<=4L^2 |X_i-X_j| d_ij
                 product_(ell!=i,j)(X_ell^2+L^2).      (13)
```

If `L<=X_1<X_2<...<X_k`, the pair `(1,2)` gives

```text
N(G)<=2^k L^3 X_2 product_(ell>=3)X_ell^2.             (14)
```

Thus the desired absolute-content target follows under the additional
condition `|X_2-X_1|d_12=O(L^2)`.  Without that condition, (14) gives
only

```text
N >= 2^(-k) X_1^2 X_2/L^3 >= 2^(-k)(X_1/L)^3,
```

the classical cubic height bound.

## 4. Relation to the existing fixed-point determinant bound

Combining the already proved affine-shape bound in
`affine_shape_radius_divisibility.md` with the angle
`2 arctan(L/X_1)<=2L/X_1` gives

```text
N >= (X_1/L)^(10/3)       for k=5,
N >= (X_1/L)^(7/2)        for k=7.                    (15)
```

These are translations of that existing determinant estimate, not new
growth exponents.  Both remain below the endpoint exponent four.  The
new information here is the exact content control (1), (10), and (13).

Run `python3 docs/check_integer_cotangent_affine_content.py` for exact
tests of the row gcd, triangle identity, pair divisibility, and the
prime-exponent consequences on small cliques, cleared random rational
parameter tuples, and the Pell counterfamily.  The checker illustrates
the identities; the valuation argument above is the proof.
