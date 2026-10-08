# Simultaneous five-row Vieta cut bookkeeping

This note records the exact combinatorics obtained by applying the
five-row Vieta selected-core rule in every row. It is a bookkeeping
statement about inherited cut factors. It does not assert that the
actual swapped roots retain every displayed factor: the correction
integers `E_i` can absorb such factors.

Let `U` range over the nonempty cuts of `{1,2,3,4,5}`, and write
`s=|U|`. For a fixed swapped row `i`, the selected-core rule is

```text
i not in U, s >= 2,       or       i in U, s >= 4.       (1)
```

Thus the selected set `S_i` has `6+4+1+4+1=16` cuts: the six pairs
and four triples avoiding `i`, the six cuts of size at least four, and
no singleton. For every retained row `j != i`, exactly `3+3+4+1=11`
selected cuts contain `j`, split by sizes `2,3,4,5`.

Applying (1) simultaneously to all five rows gives the exact incidence
table.

| old cut size `s` | new rows in which `U` is selected | multiplicity over the five new rows |
|---:|---|---:|
| 1 | none | 0 |
| 2 | the three rows outside `U` | 3 |
| 3 | the two rows outside `U` | 2 |
| 4 | all five rows | 5 |
| 5 | all five rows | 5 |

Consequently the common transferred factor (at the level of the
divisibilities `H_U | E_i P_i'`) is

```text
G_common = product_{|U|=4} H_U * H_{12345}.              (2)
```

It has six blocks. Dividing this formal transferred pattern leaves ten
blocks in each row: six pair cuts avoiding that row and four triple cuts
avoiding it. Equation (2) must not be read as a divisor of every actual
`P_i'`, because the exact statement is `H_U | E_i P_i'`.

The content lower bound in the one-row swap has exponent

```text
h_i(U) = max(0, 2(s-1)-5),  if i in U,
          max(0, 2s-4),      if i not in U.              (3)
```

For fixed `i`, this gives

```text
F_i = product_{|U|=3, i not in U} n_U^2
      * product_{|U|=4, i in U} n_U
      * product_{|U|=4, i not in U} n_U^4
      * n_{12345}^3,                                    (4)
```

and the exponent sum is `4*2+4*1+1*4+3=19`. Multiplying (4) over all
five swaps gives the exact inherited-content product

```text
product_i F_i
 = product_{|U|=3} n_U^4
   * product_{|U|=4} n_U^8
   * n_{12345}^15.                                      (5)
```

Since `E_i=g_i Norm(K_i)/F_i`, the coupled correction product is exactly

```text
product_i E_i
 = (product_i g_i)(product_i Norm(K_i))
   / ((product_{|U|=3} n_U^4)
      (product_{|U|=4} n_U^8) n_{12345}^15).             (6)
```

There is no upper bound on the numerator in (6) from the selected-core
count, so simultaneous swapping does not justify treating the five
`E_i` as small or independent.

For anchor bookkeeping, let `A` contain one conjugate copy of every old
nonempty cut block, as forced by the divisibilities `conj(P_i) | z_0`.
The product of the five transferred cores has cut exponents

```text
product_i Gamma_i
 = product_{|U|=2} H_U^3
   * product_{|U|=3} H_U^2
   * product_{|U|=4} H_U^5
   * H_{12345}^5.                                       (7)
```

The least common multiple of the five core denominators uses one copy
of each size `2,3,4,5` block, so one copy of `A` cancels all inherited
selected cores; singleton blocks were never selected. A product of five
denominator corrections cannot cancel `A` five times: relative to `A^5`,
the product (7) is short by two copies of every pair block and three
copies of every triple block. After the one-anchor cancellation, the
uncontrolled contribution is the residual correction support coming
through the `E_i` (and any non-core denominator factors).

The checker
[check_five_row_simultaneous_vieta_cut_patterns.py](check_five_row_simultaneous_vieta_cut_patterns.py)
verifies (1)--(7), including the five-row incidence table and the exact
integer exponent vectors. No endpoint, radius, or small-`E_i` claim is
made here.
