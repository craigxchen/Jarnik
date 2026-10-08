# Bounded heights and one permutation recover every pair phase

The necessary critical/successor moment conditions admit a finite state
description independent of the degree. After the final row-height groups
and the ordering of one spectral coefficient are fixed, a state consists
of eight bounded integer heights, one of them zero, and one permutation
of the rows. This recovers every pair's pattern position modulo four,
including the exceptional initial part of successor B.

The graph retains necessary conditions, not sufficient conditions for
real moment coefficients. Exact pair distances and the number of columns
are **edge weights**, not information recoverable from the bounded state.
Unweighted reachability of a terminal state does not settle the problem.

## Setup and normalization

Use the hypotheses and notation of
[the spectral orientation note](critical_moment_spectral_orientation.md).
There are `n=4s+2` retained columns, four rows in each of `A,B`, and
pair distances

```
D_ij = 2s+1 across A,B;       D_ij = 2s+2 within a group.
```

After interchanging the groups if necessary, all final signed row sums
`H_a=sum_j S_aj sign(a_j)` on `A` are equal. Write the critical
orientation as `epsilon_ab=epsilon_b`. Choose row zero in `A` and set

```
delta_a=0 (a in A),       delta_b=-epsilon_b (b in B).
```

The total number of positive signed roots in row `i` is then
`p_i=p_0+delta_i`. Fix the strict coefficient ordering of the eight
`kappa_i` provided by the spectral lemma. Its allowed three-block form
and the parity of `s` determine every critical/successor first sign.
Only the parity of `s` matters for the graph; all distances below can
be reduced modulo four when computing states.

Process the columns by increasing `|a_j|`, absorbing `sign(a_j)` into
each column. At a prefix let `u_i` count positive entries already
processed in row `i`. Define

```
r_i=u_i-u_0,       r_0=0,
tau=(-1)^(p_0-u_0).
```

For positive `t` in the gap after that prefix, order the normalized
real values `tau P_i(t)` in **increasing** order. Let `pi` be this
permutation, and let

```
sigma_ij=sign(tau P_i(t)-tau P_j(t)),
g_ij=sign(kappa_i-kappa_j).
```

Thus `sigma_ij=+1` means row `i` occurs after row `j` in `pi`.
No separate bit for `tau` is stored. The signs of the normalized values
are determined by the height vector alone:

```
w_i=sign(tau P_i(t))=(-1)^(delta_i-r_i).                (1)
```

In particular row zero's normalized value is always positive. Every
negative row must precede every positive row in `pi`; the zero gap
therefore lies before row zero.

## Exact reconstruction of all pair phases

Let `k_ij` be the number of processed differing entries for a pair.
The spectral difference formula gives

```
sigma_ij = g_ij (-1)^[(delta_i+delta_j-D_ij-r_i-r_j+k_ij)/2].
```

The exponent is an integer. If `e_ij` is zero when `sigma_ij=g_ij`
and one otherwise, this recovers

```
k_ij == D_ij+r_i+r_j-delta_i-delta_j+2e_ij (mod 4).     (2)
```

Let `eta_ij` be the fixed first sign of the pair pattern and put
`h_ij=eta_ij(r_i-r_j)`. For `k=k_ij mod 4`, the complete local
height conditions are

```
critical:    h_ij=(0,1,2,1)[k],
successor A: h_ij=(0,1,0,-1)[k],
successor B: (h_ij=0 and k=0), or h_ij=(2,1,2,3)[k].    (3)
```

For successor B the zero-height case is precisely an untouched pair.
After its first differing column, its oriented partial sum is always
one, two, or three. It can never return to zero: a transition into
phase zero comes from phase three at height three and ends at height
two. Thus no hidden initial-occurrence flag is needed.

The next required oriented sign is therefore reconstructed exactly:

```
critical:    eta_ij*(+,+,-,-)[k],
successor A: eta_ij*(+,-,-,+)[k],
successor B: eta_ij if untouched;
             eta_ij*(-,+,+,-)[k] otherwise.            (4)
```

These formulas recover the full finite pattern state. They do not
recover the integers `floor(k_ij/4)`.

## Signed columns act by interval reversals

For a proposed signed column let `I` be its plus set and `b_i=1_(i in I)`.
The height and comparison updates are literal:

```
r_i' = r_i+b_i-b_0,
sigma_ij' = sigma_ij (-1)^(b_0+b_i b_j).               (5)
```

Indeed a common positive root flips the sign of `P_i-P_j` exactly
when both rows lie in `I`; normalization also flips all comparisons
when `u_0` increases. Substitution in (2) yields

```
k_ij' == k_ij + 1_(b_i != b_j) (mod 4),               (6)
```

because `b_i+b_j+2b_i b_j` is congruent to the separation indicator
modulo four.

At the new positive root, the rows in `I` vanish together. They must
form a contiguous block in the current value ordering, with that block
containing the zero gap. Its endpoints may coincide with the gap: all
selected values may initially have the same sign. Each shared-root
difference has a simple zero, so the order inside this block reverses.
All comparisons across the block remain unchanged. The permutation
update is consequently:

1. Reverse the interval occupied by `I`.
2. Reverse the entire permutation if row zero belongs to `I`.

The second reversal is the change of normalization in (5). It is what
keeps row zero's normalized value positive. Algebraically, flipping all
internal comparisons of `I` preserves a total order precisely when
`I` is contiguous; compatibility of old and new signs in (1) then
requires that interval to contain the zero gap. Thus the interval rule
is also sufficient for these comparison and sign updates.

If the zero gap follows `z` rows, there are
`(z+1)(9-z)-1<=24` nonempty intervals containing it, including the
full interval. In addition include the empty plus set. The empty and
full plus sets are distinct constant-column labels; both are self-loops
on the compressed state and have zero pair-distance increments.

Retain a transition only when the new state satisfies (1)--(3).
The height relations and (6) then force exactly the next signs in (4)
on all separated pairs. This is an exact reconstruction of the pair
rules together with the necessary spectral-order tests, with no
unbounded counters hidden in the vertex.

## A uniform finite vertex bound

The successor-A rules inside `A` give

```
r_a in {-1,0,1},       max_(a in A) r_a-min_(a in A) r_a <=1.
```

With `r_0=0`, the other three `A` coordinates have only
`2^3+2^3-1=15` possibilities. The critical rule against row zero gives

```
r_b in {0,delta_b,2delta_b}
```

for each of the four `B` rows. There are therefore at most
`15*3^4=1215` height vectors, before enforcing the other pair rules.
Each has at most `8!` value orderings. For each fixed parity and labeled
cross-orientation class the graph has at most

```
1215*8! = 48,988,800 vertices.                         (7)
```

This is a coarse bound. Sign compatibility, pair phases, and actual
reachability can reduce it substantially. Fixing the number `q` of
positive `epsilon_b`, with `0<=q<=4`, and relabeling within the groups
gives ten normalized cases across the two parities of `s`. The checker
uses all ten; no reduction under overall sign reversal is needed.

The subsequent [exact state census](critical_moment_state_census.md)
reduces the number of valid states to 322,560, 325,584 or 326,496,
depending on `q`, by counting compatible value orders as linear
extensions. This counts all valid states, not only reachable ones.

## Exact weighted-path target

Initially `r=0` and every `k_ij=0`. The initial comparison is

```
sigma_ij = g_ij (-1)^[(delta_i+delta_j-D_ij)/2],
```

which determines the initial permutation. At the end, all positive
roots have been processed and `u_i=p_i`, so the terminal state is

```
r=delta,       pi=the increasing kappa order.           (8)
```

Every column edge has a length weight one and a vector of 28 distance
weights `1_(b_i!=b_j)`. Any actual moment solution produces a path from
the initial state to (8) with exactly

```
length = 4s+2,
total distance on every cross pair = 2s+1,
total distance on every within pair = 2s+2.             (9)
```

These exact weights remain essential. Reaching (8) with different
counts or a different length neither proves existence nor refutes an
exclusion at the required degree. Conversely, even a path satisfying
(9) certifies only these necessary discrete conditions; it does not
construct distinct real magnitudes solving the moment equations.

The subsequent [closed-graph computation and integer-flow construction](critical_moment_compressed_graph_even_q2.md)
produce exact paths satisfying (9) in the even-parity `q=2` case for every
`s=2+20805120K`, `K>=1`. Thus these discrete conditions do admit arbitrarily
large degrees; an all-degree obstruction must use further information.

## Implementation and exact audit

[The reusable checker/module](check_critical_moment_compressed_automaton.py)
provides `Automaton(s_parity,q)`, its initial and terminal states, phase
reconstruction, state validation, and transitions labeled by the plus
mask and 28 distance increments. Exact cumulative distances and length
must be tracked by its caller.

The completed audit covers **45,890 states through depth five** over all
ten normalized cases and **413,010 retained transitions**. For every
audited state it exhausts all 256 signed columns, comparing the direct
pair-comparison update with the zero-gap interval reversal. It also
checks (6), the explicit next-pattern signs, the positivity of normalized
row zero, and both constant-column self-loops. These are bounded exact
checks of the implementation. The unrestricted compression and vertex
bound are proved above; no all-degree reachability exclusion is asserted.
