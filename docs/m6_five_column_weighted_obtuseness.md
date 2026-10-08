# The seven balanced cuts for the `(1,1,1,1,2,2)` dependence

This is a bounded exact check for one unequal positive dependence in a
possible maximal binary profile. There are six rows,
four light and two heavy, in the notation below; the dependence is

```text
c = (1,1,1,1,2,2),       sum(c_i)=8.
```

Consequently a balanced cut has `c`-mass `4`.  Quotienting a cut by its
complement gives exactly seven classes.  Numbering the light rows `0,1,2,3`
and the heavy rows `4,5`, representatives are

```text
{0,1,2,3},
{0,1,4}, {0,2,4}, {0,3,4},
{0,1,5}, {0,2,5}, {0,3,5}.
```

The first is the four-light-versus-two-heavy star class.  The other six are
one heavy plus two lights, with the complementary representative chosen to
contain row `0` among the lights.  (Complementing a column does not change
any Hamming distance.)

For a selected set of columns `K` and nonnegative column weights `w`,
normalized by `sum_T w_T=1`, put

```text
m_ij(w) = sum_(T in K) sigma_ij(T) w_T,
sigma_ij(T) = +1 if T separates i,j, and -1 otherwise.
```

The strict weighted obtuseness condition is `m_ij(w)>0` for all fifteen row
pairs.  An exact Farkas certificate for its failure is a nonnegative integer
combination `y` of pair inequalities for which

```text
sum_(i<j) y_ij m_ij(w) < 0
```

for every positive `w`.  The checker finds such a certificate for every one
of the `binom(7,5)=21` five-column subsets.  Twelve subsets have a
certificate whose five coefficients are all `-1`, so they admit no
normalized nonnegative solution even at the weak boundary. The other nine have four
zero coefficients and one coefficient `-2`; their only possible weak
solutions must give weight zero to one selected column, and therefore cannot
be strict: the same nonzero combination of strictly positive margins
would be positive, whereas its coefficient form is nonpositive.
Each of these nine systems does have a weak solution: four selected
regular columns receive equal weight and the fifth receives zero.
The checker verifies those boundary witnesses separately from the
infeasibility certificates.

The full seven-column system is feasible: take every column weight equal to
`1`.  The margins are `1` for light-light pairs, `1` for light-heavy pairs,
and `5` for the heavy-heavy pair.  Thus deletion of any two classes destroys
strict obtuseness, even though the complete seven-class profile itself is
feasible.  This is a concrete finite obstruction for this prescribed
dependence; it gives no parity theorem for arbitrary odd `M` or for
multilevel allocations.

Run the standard-library checker with:

```text
python3 -B docs/check_m6_five_column_weighted_obtuseness.py
```
