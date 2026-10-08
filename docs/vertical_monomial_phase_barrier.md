# First-order monomial phase separation does not remove the four-factor deficit

This is an exact obstruction in the formal ten-block model from
[the four-cofactor note](four_vertical_lcm_five_thirds.md). It is not a
construction of Gaussian integers with a common real part, and it does
not disprove a quadratic four-cofactor bound.

The model satisfies not only the pairwise gcd bounds but **every fixed
integer monomial first-order phase-separation bound** described below.
Thus taking further products and quotients of the four cofactors and
applying off-axis Gaussian integrality alone does not repair that model.
The statement does not cover extra cancellation of the leading reciprocal
moment, sums of monomials, or identities with coefficients growing with
the heights.

## 1. The actual first-order phase inequality

Let `w_j=-d+i t_j`, where `d` is a positive integer and `t_j>=T>0`.
Fix a nonzero vector `r=(r_1,...,r_4)` of ordinary integers and put

```text
s = sum_j r_j,
U = product_(r_j>0) w_j^r_j,
V = product_(r_j<0) w_j^(-r_j),
g = gcd_G(U,V),       N=U/g,       D=V/g,
Z = N conjugate(D) in Z[i].
```

The unit chosen for `g` does not affect `Z`. Since
`arg(w_j)=pi/2+arctan(d/t_j)`, the argument of `i^(-s) Z`, modulo
`2pi`, is

```text
eta = sum_j r_j arctan(d/t_j),
|eta| <= d ||r||_1/T.
```

If `Im(i^(-s) Z)` is nonzero, it is a nonzero ordinary integer.
Consequently

```text
1 <= |Im(i^(-s) Z)|
  = |Z| |sin eta|
 <= |Z| d ||r||_1/T,

|Z| >= T/(d ||r||_1).                                  (1)
```

If that imaginary part vanishes, this argument yields no lower bound.
It is also permissible to remove the ordinary integer content of `Z`
before applying (1), giving the same inequality for the smaller modulus.
For fixed `r` and `d`, the exponent furnished by (1) is one.

The exact primewise modulus cost before ordinary-content removal is

```text
log |Z| = sum_pi |sum_j r_j v_pi(w_j)| log |pi|.         (2)
```

Thus this phase test charges the absolute signed incidence sum at each
Gaussian prime, rather than just the pair incidences.

## 2. All integer signed incidences in the ten-block model

Take one formal prime block `b_S` for each two-element or three-element
subset `S` of `{1,2,3,4}`. Give every block modulus `T^(1/6)` and put

```text
w_j = product_(S containing j) b_S.
```

Use distinct nonconjugate prime symbols, so cancellation leaves no
ordinary rational content. Each cofactor has formal modulus `T`, each
pair gcd has modulus `T^(1/2)`, and the lcm has modulus `T^(5/3)`.
For the arbitrary signed vector `r` in Section 1, (2) gives

```text
log_T |Z| = h(r)
 = [sum_(i<j) |r_i+r_j|
    + sum_(i<j<k) |r_i+r_j+r_k|]/6.                    (3)
```

For every nonzero integer vector `r`,

```text
h(r) >= 1.                                            (4)
```

Here is a proof covering all such vectors, rather than a finite search.
If `s=sum r_j` is nonzero, group the six pair sums into the three
complementary pairs. Each complementary pair has total `s`, so the
sum of the six absolute values is at least `3|s|`. The four triple
sums are `s-r_j` and have total `3s`, so their absolute values also
sum to at least `3|s|`. Since `s` is a nonzero integer, (4) follows.

If `s=0`, the triple contribution is `sum |r_j|>=2`. The pair
contribution is

```text
2 (|r_1+r_2|+|r_1+r_3|+|r_1+r_4|).
```

The sum in parentheses is a nonnegative even integer: replacing its
absolute values by the signed arguments does not change its parity,
and those arguments sum to `2r_1`. It cannot vanish for a nonzero
`r`, because then `r_2=r_3=r_4=-r_1` and `s=0` forces every entry
to vanish. The pair contribution is therefore at least four, proving
(4) in the remaining case.

Equality occurs for `r=(1,0,0,0)` and `r=(1,-1,0,0)` and their
permutations and negatives. In particular, higher signed monomials
do not expose a hidden exponent below one in this model.

## 3. Precise remaining gap

For `T>=1`, (4) gives `T^h(r)>=T>=T/(d ||r||_1)`. Hence the formal
`T^(5/3)` lcm model satisfies every necessary inequality (1), even
with arbitrary fixed integer exponents in the monomial. This includes
the possible removal of ordinary content because the formal model
has none.

This is only a barrier for the specified inequality family. Actual
common-real-part factors must satisfy more than incidence and (1).
For example, a separately established cancellation of
`sum_j r_j/t_j` could make `eta` substantially smaller; the sine bound
used here deliberately retains only its universal first-order size.
Positive partial fractions or additive identities could also introduce
constraints absent from (2). No such constraint sufficient for a
quadratic lcm estimate is established in this note.

Run `python3 docs/check_vertical_monomial_phase_barrier.py` for all 6,560
nonzero signed vectors in `[-4,4]^4`, and exact Gaussian quotient checks
of the squared inequality in Section 1, both before and after ordinary
content removal. The checker includes zero-projection cases and makes
no inference from them. The proofs above establish the infinite claims.
