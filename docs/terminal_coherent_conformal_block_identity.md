# The terminal coherent kernel as a highest-root annihilator

Let `m=3d` with `d=2k>=2`. The terminal constituent `K_k` is the
space of multilinear `SL_3` invariant polynomials in `m` ternary rows
`v_i=(a_i,b_i,c_i)`. It is spanned by products of `d` ternary
determinants on disjoint triples, so every member has total degree
`d` in each of the three coordinate families and satisfies

```text
Q((a_i,b_i,c_i) g)_i = det(g)^d Q((a_i,b_i,c_i))_i
```

for `g in GL_3`. Fix scalar directions `z_i` and define

```text
D_z = sum_i z_i a_i partial_(c_i),
E_z Q(x) = Q((1,z_i,x_i)_i).
```

The following is an identity of polynomials, with no genericity or
distinctness assumption on `z`:

```text
D_z^d Q(a,b,c)
  = (-1)^d d! (product_i a_i)
      E_z Q(b_1/a_1,...,b_m/a_m).                       (1)
```

The right side is interpreted by squarefree homogenization, so it is
polynomial even when some `a_i=0`. Consequently

```text
ker(E_z:K_k -> squarefree degree-d polynomials)
  = ker(D_z^d:K_k -> multilinear polynomials).          (2)
```

For the terminal cases `d=2k`, the sign in (1) is positive.

## Direct proof

As a polynomial in the `c_i`, `Q` is homogeneous of total degree `d`.
The top directional derivative therefore gives

```text
D_z^d Q(a,b,c) = d! Q((a_i,b_i,z_i a_i)_i).            (3)
```

Swap the second and third coordinates globally. The swap has
determinant `-1`, so the relative-invariant law gives

```text
Q((a_i,b_i,z_i a_i)_i)
  = (-1)^d Q((a_i,z_i a_i,b_i)_i).
```

Finally, factor `a_i` out of each row by multilinearity. This proves
(1) on the dense locus where every `a_i` is nonzero, hence as a
polynomial identity everywhere. Since squarefree homogenization is an
injective linear map, (2) follows. The exact finite
[checker](check_terminal_coherent_conformal_block_identity.py) audits
the sign and coefficient indexing for explicit products of disjoint
triangles at `d=1,2,4`; the proof applies to every `d` and every
linear combination.

## Conformal-block convention

Let `V=(F^3)^(tensor m)` and let `A=V/(sl_3 V)` be its invariant
coinvariant quotient. A multilinear invariant polynomial `Q` is an
element of `A^*`. Let `E_31` send the first standard basis vector to
the third and set

```text
T_z = sum_i z_i E_31^(i) on V.
```

Precomposition of a functional with `T_z` is exactly the differential
operator `D_z` above. Therefore (2) can be written

```text
ker E_z = {Q in A^*: Q(T_z^d V)=0}
        = (A / image(T_z^d -> A))^*.                    (4)
```

The [genus-zero conformal-block criterion of
Belkale–Brosnan–Mukhopadhyay, Section 9.1](https://math.umd.edu/~pbrosnan/Papers/BBMArxiv.pdf)
describes the level-`ell` coinvariant quotient by the image of
`(sum_i z_i e_theta^(i))^(ell+1)`, where `e_theta` is a highest-root
operator. Our `E_31` is the opposite-root convention. A common Weyl
conjugation exchanges the two operators and acts trivially on invariant
functionals, so their annihilators agree. Taking `ell=d-1` and all
marked weights equal to the defining fundamental weight identifies
`ker E_z` with the **dual of that coinvariant quotient**, conventionally
the conformal-block functional space. This distinction matters because
some authors call the quotient itself the dual conformal-block space.

The identification is pointwise for every distinct configuration of
marked points. It supplies the representation-theoretic factorization
and connection language for studying collision divisors, but it does
not by itself quantify the integral content of a kernel Pluecker
vector at the Gaussian full-profile rows. In particular, a conformal
block rank or a fusion-path count alone is not a coefficient-height
floor.
