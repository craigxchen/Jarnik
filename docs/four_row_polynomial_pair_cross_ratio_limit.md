# A six-pair-block limit on four-row polynomial exclusion

This note isolates a necessary cross-ratio equation for a direct
four-row polynomial full-Boolean profile. It also gives six distinct
Gaussian linear factors whose norms satisfy that equation with
coefficients arising from four ordered real leading terms. Thus the
pair-block cross ratio and positivity alone cannot exclude a
four-row profile. The example does **not** construct the fifteen
factors or constant-imaginary rows, and proves nothing about integer
circle endpoints.

## 1. Exact pair differences and the cross ratio

Suppose the fifteen nonempty subsets `T` of `[4]` index Gaussian
polynomials `H_T in Q(i)[t]` of one positive degree `e`. Assume
they and their conjugates are pairwise coprime, and put
`n_T=H_T bar(H_T)`. Rescale leading coefficients into the constant
`K_i`, so every `H_T` is monic. Suppose

```text
P_i=K_i product_(T contains i)H_T=X_i+i c_i,
X_i in Q[t],      c_i in Q\{0}.                         (1)
```

Each row has degree `d=8e`. Write `f_i=X_i/c_i`. For any pair,
the four blocks with `T` containing both labels divide
`Im(bar(P_i)P_j)=c_i c_j(f_i-f_j)`. Their norm product has degree
`4*(2e)=d`. The difference is nonzero: otherwise two rows in
(1) would have the same Gaussian zeros, contrary to the prescribed
distinct incident blocks. Hence

```text
f_i-f_j=lambda_ij product_(T contains i,j)n_T,
lambda_ij in Q\{0}.                                   (2)
```

For real `t`, every `n_T(t)>0`; a real zero of `H_T` would also
be a zero of its conjugate. Consequently every `f_i-f_j` has a
constant sign on the real line. The four `f_i` are globally ordered,
and their degree-`d` leading coefficients `a_i` are distinct, with
`lambda_ij=a_i-a_j`.

The elementary difference identity

```text
(f_1-f_2)(f_3-f_4)
 -(f_1-f_3)(f_2-f_4)
 +(f_1-f_4)(f_2-f_3)=0
```

and (2) have a large common factor: every triple block occurs once
in each term, and the full block occurs twice. Cancelling them gives
the exact six-pair-block condition

```text
lambda_12 lambda_34 n_12 n_34
 -lambda_13 lambda_24 n_13 n_24
 +lambda_14 lambda_23 n_14 n_23=0.                 (3)
```

The leading coefficients of (3) satisfy the same identity because
the `lambda_ij` are differences of four numbers. This is a genuine
necessary condition of the full polynomial profile; it contains no
singleton, triple, or full-cut factor after cancellation.

## 2. Exact independent Gaussian-norm solution of (3)

Take six pair norms

```text
n_12=t^2+1,     n_34=t^2+9,
n_13=t^2+36,    n_24=t^2+49,
n_14=t^2+25,    n_23=t^2+225.                        (4)
```

They are norms of the six monic Gaussian linear factors `t+i b`
with `b=1,3,6,7,5,15`, respectively. These factors and all their
conjugates have distinct roots in `Q(i)[t]`. Direct expansion gives

```text
11 n_12 n_34-16 n_13 n_24+5 n_14 n_23=0.        (5)
```

The coefficient ratio is compatible with four globally ordered
leading coefficients. For example set

```text
(a_1,a_2,a_3,a_4)=(1,38/27,3/2,2).
```

With positive differences `mu_ij=a_j-a_i` for `i<j`, one gets

```text
(mu_12 mu_34,mu_13 mu_24,mu_14 mu_23)
                         =(11,16,5)/54.
```

Using `lambda_ij=-mu_ij` in (3) recovers (5). All six norms are
positive on the real line, so this also meets the positivity part
of the pair-block condition. The triple additive equations, row
factorizations (1), and the remaining nine Gaussian blocks have
not been supplied. In particular (5) is a countermodel only to an
exclusion based solely on (3), even when its six factors have the
required degree and polynomial coprimality.

The [exact checker](check_four_row_polynomial_pair_cross_ratio_limit.py)
verifies the coefficient identity, polynomial norm factorization,
root disjointness, and leading-coefficient ratios.
