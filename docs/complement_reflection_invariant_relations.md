# Complement reflection doubles the invariant-relation modulus

The smaller-cut congruences can be applied a second time after an exact
common reflection of the actual Gaussian directions. This preserves
every invariant zero with its original coefficients and preserves
primitive pair residues up to sign. After a common trimming of the
core blocks, the second modulus is supported on the complementary
block. Its coprimality with the original block doubles the leading
coefficient threshold from `w` to `2w`.

The finite constants below are deliberately not optimized. The exact
Gaussian construction, including a zero imaginary coordinate after
reflection, is checked by
[check_complement_reflection_invariant_relations.py](check_complement_reflection_invariant_relations.py).
This is a restriction on exact invariant relations, not a proof of
the uniform lattice-arc bound.

## 1. An integral reflected configuration

Use `m` actual conjugate-primitive Gaussian rows

```text
P_i=K_i A_i,       A_i=product_(T containing i)H_T,
gcd_G(P_i,bar(P_i))=1,       log|K_i|<=sigma.
```

All `H_T` and their conjugates have disjoint support at odd split
primes. Define

```text
K_all=product_i K_i,       H_all=product_T H_T,
Z=K_all H_all,
R_i=Z/P_i=(product_(j!=i)K_j) product_(T not containing i)H_T.
```

These are actual Gaussian integers. Every original `P_i` has odd
norm: otherwise `1+i` would divide both `P_i` and its conjugate.
Hence all `K_i,H_T,Z,R_i` have odd norm. Put

```text
d_i=gcd(|Re R_i|,|Im R_i|)>0,       P_i'=R_i/d_i.          (1)
```

The division is by a positive ordinary integer. Since the norm is
odd, the resulting row is conjugate-primitive. In particular there
is no division by `1+i` rotating individual directions.

Let

```text
E_T=gcd_G(H_T,bar(K_all)),       H_T^*=H_T/E_T.           (2)
```

Choose the associates in this expression consistently; their units
do not affect any norm or divisibility assertion. For every row,

```text
P_i'=K_i^* product_(T not containing i)H_T^*,
K_i^* in Z[i],       log|K_i^*| <= (2m-1)sigma.          (3)
```

To prove the divisibility in (3), take `pi^e || H_T`. The opposite
orientation does not occur in the core product of any `R_i`. Thus

```text
v_p(d_i)=min(v_pi(R_i),v_barpi(R_i)) <= v_barpi(K_all).
```

If `i` is outside `T`, the row `R_i` contains `pi^e`; after division
by `d_i`, it still contains
`pi^max(0,e-v_barpi(K_all))`, exactly the power in `H_T^*`.
This works simultaneously for all active rows and all blocks.

Because the original blocks are pairwise coprime,
`product_T E_T` divides `bar(K_all)`. It follows that

```text
sum_T log N(E_T) <= 2 log|K_all| <= 2m sigma,             (4)
log N(H_T^*) >= log N(H_T)-2m sigma.                    (5)
```

For (3), compare norms in (1): before division the correction is
`K_all/K_i`, of modulus at most `exp((m-1)sigma)`; replacing the
row's old core by the trimmed core increases that correction by
a factor of modulus at most `product_T|E_T|<=exp(m sigma)`.
The positive divisor `d_i` only decreases it.

The new blocks, indexed by their sets of dividing rows, are

```text
H_S'=H_(S^c)^*.                                        (6)
```

They retain the independent conjugate-coprime support property.

## 2. Exact zero relations and primitive residues are preserved

The transformation of the actual binary rows is

```text
P_i' = (1/(N(P_i)d_i)) Z bar(P_i).                       (7)
```

The map `z -> Z bar(z)` is one common real linear map of determinant
`-N(Z)`. The remaining factors in (7) are positive rational row
scalings. Therefore every simultaneous `SL_2` invariant of degree
two in each row obeys its usual `GL_2` covariance and row
homogeneity. In particular

```text
Q(P_1,...,P_m)=0  implies  Q(P_1',...,P_m')=0,            (8)
```

with exactly the same polynomial and coefficients. No permutation
of the rows or change of invariant basis is needed.

For a pair, `P_i' bar(P_j')` is a positive rational multiple of
`bar(P_i)P_j`. Divide each pair product by the norm of the Gaussian
gcd of its two row numerators. The resulting conjugate-primitive
Gaussian integers are still positive rational multiples of one
another, and so are equal: ordinary primitive integer coordinate
pairs on the same ray coincide. Consequently

```text
t_ij'=-t_ij,       |t_ij'|<=T.                          (9)
```

The reflected configuration is not asserted to preserve a chosen
small imaginary-coordinate bound. That bound is not used in the
smaller-cut congruence proof.

A reflected imaginary coordinate can be zero. This causes no gap:
a real conjugate-primitive Gaussian integer is `+1` or `-1`.
Any dividing core block is then a unit. More directly, for an
inside row the argument using `2iY_i'=P_i'-bar(P_i')` still gives
`gcd_G(H,Y_i')=1`, even if `Y_i'=0`. Outside imaginary coordinates
may also vanish; the distinguished matching value used in the
minor proof is nonzero because the directions remain distinct.
Thus the local proof applies without requiring every `Y_i'` to be
nonzero.

## 3. The doubled maximal-minor theorem

Let `m=2q>=6`, `n=q+1`, and `t=n(n-3)/2`. Fix an inside set `S`
of size `q-1`. In the integral restriction basis of
[smaller_cut_relation_rank.md](smaller_cut_relation_rank.md), every
maximal `t` by `t` coefficient minor `delta_S` of actual numerical
zero relations in the balanced kernel obeys

```text
M_S | delta_S,
log M_S >= log N(H_S)-(4n-4)sigma-2log T.                (10)
```

Apply the same theorem to (3), (6), (8), and (9). Its coefficient
minor is literally the same integer `delta_S`, since the invariant
polynomials and the labelled restriction are unchanged. It yields

```text
M_S' | delta_S,
log M_S' >= log N(H_(S^c))-2m sigma
            -(4n-4)(2m-1)sigma-2log T.                 (11)
```

The first modulus divides `N(H_S)` and the second divides
`N(H_(S^c)^*)`. These supports are disjoint, so the product divides
the same integer. Thus

```text
M_S M_S' | delta_S,
log(M_S M_S') >= log N(H_S)+log N(H_(S^c))
                  -2m(4n-3)sigma-4log T.               (12)
```

In a full profile this is at least
`2(1-eta)w-2m(4n-3)sigma-4log T`.
No prime-power exponent has been discarded except for the explicitly
charged common trimming (2).

Since the determinant has absolute value at most the product of
the coefficient sums of its relations, all relations of coefficient
height at most `h` have a proper restriction image at every such
cut whenever

```text
t h+2m(4n-3)sigma+4log T < 2(1-eta)w.                  (13)
```

## 4. The six-row conclusion

For six rows, `n=4`, `t=2`, and `2m(4n-3)=156`. The exact
exterior-product injectivity from
[two_small_invariant_relations.md](two_small_invariant_relations.md)
therefore gives the stronger bound

```text
log C_Q+log C_R >= 2(1-eta)w-156sigma-4log T              (14)
```

for independent numerical zero relations in the balanced kernel.
Every pair whose coefficient-height sum is strictly below this
right side is automatically in that kernel: each individual height
is below the earlier balanced threshold `2(1-eta)w-24sigma`.
Thus all exact degree-two-each numerical zeros of coefficient
height below

```text
(1-eta)w-78sigma-2log T                                  (15)
```

span a space of dimension at most one. This doubles the leading
threshold of the one-pass proof. It still does not supply the two
independent relations needed for a contradiction on six rows.
