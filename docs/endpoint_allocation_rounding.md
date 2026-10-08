# Rounding prime-power allocations with negligible-height row twists

## Status

This is an independently audited prose reduction for a hypothetical
sequence of endpoint circle clusters whose cardinalities tend to infinity.
It does not establish the uniform endpoint bound. It removes the nested
prime-power structure from the conductor at the cost of moving, small-height
linear forms. The original points still supply the short-arc equations;
the rounded points themselves need not form a short arc.

The stronger [central-truncation theorem](endpoint_central_truncation.md)
now gives exact Gaussian integer row factors of one common modulus,
and permits the number of selected rows to grow like `c log_2 M`
for any fixed `c<1/2`. Its coefficients and primitive residues are
small even relative to each individual retained block.

Throughout, `h` means absolute logarithmic Weil height. All logarithms
are natural. In an asymptotic statement the number of selected rows is
fixed while the original cluster cardinality tends to infinity.

## 1. Pair distances bound the total interior allocation

Remove the actual common Gaussian divisor of an `M`-point circle cluster.
This can only decrease its endpoint arc constant. The remaining squared
radius `N=R^2>1` is odd and supported on split primes. Write

```text
z_i = epsilon_i product_p pi_p^(a_ip) bar(pi_p)^(e_p-a_ip),
0 <= a_ip <= e_p,
W = log R^2 = sum_p e_p log p,
d_ij = sum_p |a_ip-a_jp| log p.
```

Let `g_ij` be a Gaussian gcd and `c_ij=(z_i-z_j)/g_ij`.
After dividing by the gcd the two points still have equal norm, so
their nonzero difference has even squared norm. Thus `|c_ij|^2>=2`.
Also `N(g_ij)=exp(W-d_ij)`. An arc of length at most `C sqrt(R)` gives

```text
d_ij >= W/2 + log(2/C^2).                         (1)
```

This holds across all four Gaussian unit classes. Assume `C<=sqrt(2)`;
for the uniformity problem it suffices to fix any such positive `C`
and subdivide longer arcs into finitely many pieces.

Put `w_p=e_p log p`, `x_ip=a_ip/e_p`, and

```text
r_p(t) = #{i : x_ip >= t},
V_p = integral_0^1 (r_p(t)-M/2)^2 dt,
lambda = log(sqrt(2)/C) >= 0.
```

The coarea identity is

```text
sum_(i<j) |x_ip-x_jp| = M^2/4 - V_p.
```

Summing (1) over pairs therefore gives

```text
sum_p w_p V_p <= M W/4 - M(M-1)lambda.            (2)
```

The interior allocation of row `i` is

```text
E_i = sum_p min(a_ip,e_p-a_ip) log p.
```

For each prime, the following identity and Cauchy--Schwarz inequality
control it:

```text
sum_i min(x_ip,1-x_ip)
 = integral_0^(1/2) [r_p(t)-r_p(1-t)] dt
 <= sqrt(V_p).
```

Indeed, subtract `M/2` from both terms in the integrand and view the
result as the inner product with the function that is `+1` on the
first half of the interval and `-1` on the second. Its squared norm
is one. Weighted Cauchy--Schwarz and (2) now prove

```text
sum_i E_i <= sqrt[W (M W/4-M(M-1)lambda)]
          <= W sqrt(M)/2.                         (3)
```

In particular, all but at most `M^(3/4)/2` rows satisfy

```text
E_i <= W M^(-1/4).                                (4)
```

The `sqrt(M)` scale in (3) cannot be improved using the pair-distance
relaxation alone. For square `M`, symmetrize a vector with `sqrt(M)`
entries equal to `1/2` and `(M-sqrt(M))/2` entries at each endpoint.
Its average distance between any two rows is `1/2`, and its total
interior allocation is `sqrt(M)/2` per unit weight. This is a real
allocation example, not a claimed circle configuration.

## 2. Exact rounding and the height of the error

Round `a_ip` to its nearer endpoint `b_ip in {0,e_p}`, resolving ties
arbitrarily. Define actual Gaussian integers

```text
ztilde_i = epsilon_i product_p pi_p^(b_ip) bar(pi_p)^(e_p-b_ip),
U_i = product_p pi_p^((a_ip-b_ip)_+)
                  bar(pi_p)^((b_ip-a_ip)_+).
```

Then

```text
|ztilde_i|=R,
z_i = (U_i/bar(U_i)) ztilde_i,
log N(U_i)=E_i,
h(U_i/bar(U_i))=E_i/2.                            (5)
```

The last equality uses coprimality of `U_i` and its conjugate.
All their prime factors are split, and each is present in just one
orientation. The logarithmic height of the ordinary rational
rotation denominator, in real and imaginary coordinates, is `E_i`;
it should not be confused with the Weil height in (5).

For the rows in (4), the error height is `o(W)`. This does **not**
mean that the error angle tends to zero: a fixed rational rotation
already has negligible height on this scale. No assertion that the
rounded points cluster is used below.

Two good rounded allocation vectors cannot coincide once
`2 M^(-1/4)<1/2`. Otherwise
`d_ij<=E_i+E_j<W/2`, contradicting (1). Thus the good rounded rows
remain distinct, even before using their units.

## 3. Near-uniform profile extraction survives rounding

Keep the good rows and, if desired, take their largest common unit
class. Their number still tends to infinity. The sign vectors of
their original threshold allocations have pairwise nonpositive
inner products by (1). The extraction theorem in
[uniform_profile_extraction.md](uniform_profile_extraction.md)
therefore selects any fixed `k` rows with inherited oriented pattern
weights

```text
(1+error) W/2^k,   |error|<=eta,   eta -> 0.
```

For a selected row, the total weight on which its threshold bit
differs from its rounded endpoint bit is exactly `E_i`. Consequently
rounding changes the joint pattern distribution in total variation
by at most `sum_i E_i/W`. Each individual pattern weight changes by
at most `sum_i E_i`. If `E_i<=tau W` for the selected rows, the
rounded relative error is at most

```text
eta_end = eta + 2^k k tau.                        (6)
```

For fixed `k`, taking `tau=M^(-1/4)` makes `eta_end` tend to zero.
Discarding every prime that has an interior allocation would not
give this conclusion: all primes can have an interior row in the
sharpness example above. The estimate is specifically rowwise.

## 4. Independent Gaussian blocks and moving linear equations

Use common-unit selected rows numbered `0,...,m`, with `m=k-1`.
At each rounded prime, orient its entire power opposite to the
allocation of row zero, and let `S subset {1,...,m}` consist of the
rows whose rounded allocation differs from row zero. The product
of these oriented powers is `H_S`. Include `S=empty`.

Every prime belongs to one block only. The blocks are pairwise
coprime, each is coprime to its conjugate, and distinct blocks are
also coprime to each other's conjugates. The two complementary
oriented patterns give

```text
log N(H_S) = (1+error_S) W/2^m,
|error_S| <= eta_end.
```

Set

```text
A_i = product_(S contains i) H_S,
V_i = U_i bar(U_0).
```

The exact ratio identities and size estimates are

```text
ztilde_i/ztilde_0 = A_i/bar(A_i),
z_i/z_0 = (V_i A_i)/bar(V_i A_i),
(1-eta_end)W/4 <= log|A_i| <= (1+eta_end)W/4,
log|V_i| <= tau W.
```

The angular diameter of the **original** cluster is at most
`C exp(-W/4)`. Its distinctness and the last ratio identity imply

```text
0 < |Im(V_i A_i)|
  <= (C/2) exp[(tau+eta_end/4)W].                 (7)
```

Write `V_i=p_i+i q_i`, `A_i=X_i+i Y_i`, and
`t_i=Im(V_i A_i)`. Then (7) is the exact integer equation

```text
q_i X_i + p_i Y_i = t_i != 0,                    (8)
log max(1,|p_i|,|q_i|,|t_i|) = o(W).
```

Thus any unbounded sequence of endpoint clusters produces, for
every fixed number of rows, independent Gaussian cut blocks of
asymptotically equal weight whose row products obey simultaneous
moving integer linear equations with negligible-height coefficients
and nonzero values. This reduction permits arbitrary initial prime
exponents. It neither assumes a fixed prime support nor supplies one.

Selecting a common unit class is only a convenience. Any ratio of
Gaussian units is `v/bar(v)` for one of
`v in {1,i,1+i,1-i}`. Absorbing this factor into `V_i` costs at most
`sqrt(2)` and gives the analogous equations without that selection.

## 5. Normalization and scope

The empty block is common to all selected rounded points, up to
conjugation. Removing it changes the intrinsic weight to
`W'=(1-2^(1-k)+o(1))W`; all nonempty blocks remain uniformly
weighted. The ratio equations are unchanged.

One must retain the original angular precision `C exp(-W/4)`.
Replacing it by a fresh endpoint estimate at the smaller intrinsic
radius loses information. Equations (7)--(8) remain valid with
negligible height measured relative to `W'`, since `W'/W` tends to a
fixed positive constant.

Small coefficient height does not allow deletion of every conductor
prime dividing a coefficient. It controls the coefficient's own
valuation mass, not the entire conductor power above that prime.
This distinction matters when prime exponents grow.

The remaining arithmetic task is to rule out the simultaneous
system (8) for some fixed number of rows and the stated uniform
block weights, or to prove another consequence incompatible with
an original endpoint cluster. No such impossibility theorem is
claimed here. The existing fixed-coefficient polynomial divisibility
barrier and fixed-prime-support limitations require separate checks
when applied to these moving equations.

## Verification

The pair-distance constant, integral identity, weighted
Cauchy--Schwarz estimate, exact Gaussian multiplier, and preservation
of pattern weights were independently audited. Exact arithmetic
checks covered 15,624 allocation vectors, including exhaustive
small cases and random larger integer exponents. The symmetrized
sharpness identities were checked through `M=1024`. A separate
derivation checked the block orientation, unit factors, and
nonzero linear-form values in (7)--(8). These checks supplement the
proofs above; no Lean formalization of this note is claimed.
