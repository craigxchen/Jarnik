# Conjugation and ordered slopes in the coherent kernel

This note isolates what the actual Gaussian rows add to the metric of the
integral coherent kernel. Conjugation gives an exact symmetry of the same
generator matrix, not a second independent matrix. Ordered near-real slopes,
even together with conjugate-primitivity, do not bound the first minimum at
any fixed invariant rank. The balanced full-Boolean norm profile is a separate
arithmetic constraint, and the argument below does not remove it.

## 1. Conjugation is a common shear

Write the actual primitive anchor numerators as

```text
P_i=X_i+iY_i,       gcd(P_i,bar(P_i))=1.
```

For a `2n`-set `I`, the unnormalized local column is the determinant

```text
L_I(P,Y)=det[(P_i V_i | Y_i V_i)_(i in I)].
```

Since `bar(P_i)=P_i-2iY_i`, subtracting `2i` times each column in the
second `n`-column block from the corresponding column in the first gives

```text
L_I(bar(P),Y)=L_I(P,Y)=L_I(X,Y).                    (1)
```

Complementary determinant factors do not depend on the binary rows, so (1)
holds for every local generator before and after coefficient primitivization.
It also follows directly from `SL_n` invariance. If `E_(X,Y)(Q)` is the
coherent residual polynomial of an invariant `Q`, a determinant-one shear of
the first two source-coordinate families gives

```text
E_(X+tY,Y)(Q)=E_(X,Y)(Q)       for every complex t.  (2)
```

In particular `E_(P,Y)=E_(bar(P),Y)=E_(X,Y)`. The matrix has ordinary
integer coefficients in the `X,Y` chart. Stacking it with its conjugate
does not increase its rank; its Hermitian Gram is just twice the usual real
Gram. This symmetry supplies no new positive definite form on the kernel
beyond restricting the ambient Euclidean form.

At an odd split core prime, conjugation exchanges the two Gaussian local
charts, but every kernel coefficient is an ordinary integer. Its valuations
at the two primes above `p` are equal. Passing from its absolute value to
its Gaussian norm squares both any divisibility floor and the upper height;
there is no extra exponent from counting the conjugate chart twice.

## 2. Ordered primitive Gaussian rows can have arbitrarily large first minimum

Fix `n,d>=2` and use any fixed integral basis of the multilinear `SL_n`
invariant lattice `M_Z`. For every `B>0` and `epsilon>0`, there are distinct
positive even integers `X_1<...<X_m`, `m=nd`, with `Y_i=1`, such that

```text
max_i Y_i/X_i < epsilon,
gcd(X_i+i, X_i-i)=1,
lambda_1(ker_Z E_(X,Y)) > B.                          (3)
```

For `(n,d)=(3,4)`, this is a rank-341 saturated kernel in the same
rank-462 ambient invariant lattice used in the terminal calculation.

To prove (3), take a nonzero `Q in M_Z`. The coherent residual
`E_(X,Y)(Q)` is obtained by leaving the other `n-2` source coordinates
as independent variables and substituting `(X_i,Y_i)` in the first two
coordinates of every row. It is a nonzero polynomial: row multilinearity
means distinct source-coordinate assignments remain distinct monomials.
The same is true after setting every `Y_i=1`; the rows in which the second
coordinate was selected are exactly the rows where neither `X_i` nor a free
coordinate appears. Thus some coefficient of `E_(X,1)(Q)` is a nonzero
polynomial in `X_1,...,X_m`.

There are only finitely many nonzero `Q` with coefficient sup norm at most
`B`. Choose one such nonzero coefficient polynomial for each and take their
product. The product is not zero. The lattice of strictly ordered positive
even integer tuples is Zariski dense in `R^m` (choose the variables
successively from infinite admissible integer tails), so one tuple avoids
the zero set of this product. No kernel vector has coefficient sup norm at
most `B`, and therefore its Euclidean first minimum exceeds `B`. The
distinct-configuration theorem gives the stated kernel rank for the chosen
tuple.

Because every `X_i` is even and `Y_i=1`, the Gaussian integer `X_i+i`
has odd norm and is coprime to its conjugate: a common Gaussian prime would
divide `2i`, hence lie above 2, contradicting the odd norm. Applying a
common even shear `X_i -> X_i+2t` preserves all brackets and the entire
kernel by (2), while making every ratio `1/X_i` as small as desired.
The ordering of these ratios is the reverse of the ordering of the `X_i`.

As a small exact illustration, at `(n,d)=(2,2)` take

```text
(X_1,X_2,X_3,X_4)=2(S,S+1,S+q,S+q+1),
Y_i=1,       S>=1, q>=2.
```

In the invariant basis `f_1=[12][34]`, `f_2=[14][23]`, their coherent
evaluations are `4` and `4(q^2-1)`. The primitive kernel vector is
`(q^2-1,-1)`, so its first minimum grows like `q^2` despite strict
ordering, positive bracket products, conjugate-primitivity, and arbitrarily
small slopes as `S` grows.

This theorem rules out a gain from real order, near-reality, the elementary
conjugation involution, and primitive parity alone. The rows in (3) are not
asserted to have the balanced full-Boolean Gaussian factorization with
subpower correction factors. Clearing their norms to place them on one
circle can enlarge the radius enough to lose the endpoint arc scale.

## 3. The extra Hermitian arithmetic is a simultaneous near-square system

The full-Boolean profile does give a genuine identity beyond bracket
valuations. Write

```text
P_i=K_i product_(T containing i) H_T,
n_T=Norm(H_T),
G_ij=product_(T containing i,j) n_T,
Z_ij=bar(P_i) P_j/G_ij=s_ij+i t_ij.
```

Each `Z_ij` is a Gaussian integer: every shared `H_T bar(H_T)=n_T`
divides `bar(P_i)P_j`. Its imaginary part is the small bracket residue
`t_ij=Delta_ij/G_ij`, while its real part is the integer
`s_ij=(X_iX_j+Y_iY_j)/G_ij`. The prescribed norm is

```text
s_ij^2+t_ij^2
 = Norm(K_i) Norm(K_j)
   product_(T containing exactly one of i,j) n_T.     (4)
```

For three distinct labels, cancellation in the rank-one Hermitian matrix
`(bar(P_i)P_j)` gives the exact cocycle

```text
Z_ij Z_jk=C_ijk Z_ik,
C_ijk=Norm(K_j)
  product_(T: T intersect {i,j,k} is {j} or {i,k}) n_T. (5)
```

Hence `s_ij t_jk+t_ij s_jk=C_ijk t_ik`. On a short ordered arc with
positive real coordinates, the `s_ij` and consistently oriented `t_ij`
are positive. Under equal full-Boolean block weights, both `s_ij` and
`C_ijk` have logarithmic size `2^(m-2)w+o(w)`, so (5) balances at the
expected scale. Equations (4)-(5) couple all small imaginary residues to
the oriented Gaussian norm factors. They do not, by themselves, yield a
smaller covolume or a short vector in the coherent kernel; a strict gain
would require a bound on their simultaneous global integer solutions.
