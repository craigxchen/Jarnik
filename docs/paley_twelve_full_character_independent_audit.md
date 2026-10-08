# Independent audit of the order-twelve Paley character obstruction

The [order-twelve proof](paley_twelve_full_character_height_obstruction.md)
and its [checker](check_paley_twelve_full_character_height_obstruction.py)
agree. The checker passes and compares the actual set of 187 norm-one
ternary certificates with the generated row-pair, half-column-pair, and
single-column set. The support table gives all nonextremal ternary cases.
For amplitude at least two, Hadamard inversion and the triangle inequality
give `B(c) >= 18 max_x |c_x| >= 36`. For amplitude one, the table and the
single-flip perturbation bound give `B(c) >= 36` away from norm one.
The classified norm-one cases give `B(c) >= 28`, with equality attained
by the canonical assignment. Thus the stated bound applies to **all**
nonzero zero-sum integer `c`, rather than only the enumerated vectors.

The saturation argument uses the at least three unflipped physical copies
of each label: comparing any flipped copy with an unflipped copy forces
`2 lambda_x` integral. Since the row sum vanishes, the converse integral
coefficient claim follows by parity. Distinct split rational primes keep
the positive and negative Gaussian prime orientations disjoint, so the
composite height is exactly `sum_j |v_j| log p_j`. The fixed-modulus prime
count and pigeonhole argument give actual prime logarithms within `1/55`
of a common `ell`; the weighted off-diagonal Gram is at most
`-ell + 1 < 0`. With `B(c) >= 28`, its primitive half-height excess is
greater than `(ell-1)/2`.

For any normalized Hadamard matrix of order `M`, with five copies per
nonconstant column and at most two assigned flips per label, a ternary
norm-one certificate has the following exact behavior. Write
`h=|supp(c)|`, and let `k` count support rows whose selected flip label
has a nonzero transform coefficient. Equality in Hadamard inversion
forces every active transform sign to align at every support row; hence
every such selected flip lowers a coefficient magnitude by one. Therefore

```text
B(c) = 5M/2 + h - 2k.
B(c) > 5(M-1)/2  if and only if  2k-h <= 2.
```

If at least two transform labels are active, their signed columns agree
on every support row. Orthogonality bounds `h <= M/2`. This is a useful
general equality-case criterion, but it does not classify the certificates
or bound `k` at arbitrary orders. In particular, the Walsh four-row
example can have `2k-h=4`.

As a **bounded sparse diagnostic only**, the
[order-twenty checker](check_paley_twenty_sparse_diagnostic.py) verifies
the normalized Paley matrix for `q=19` and exhausts ternary zero-sum
certificates of supports four and six, up to global sign. It checks
14,535 and 387,600 certificates respectively and finds minimum raw
transform sums 24 and 28, both above the norm-one value `M=20`.
It says nothing about supports eight through twenty or arbitrary integer
amplitudes, and it does not establish an unbounded Paley obstruction.
