# Exact weighted paths in infinitely many compressed degrees

The even-parity `q=2` automaton admits paths with all required length and
pair-distance labels for every `s=2+20805120K`, `K>=1`. This rules out an
all-degree path exclusion using these necessary ordering conditions.
It does not construct real magnitudes satisfying the moment equations.

For the even-`s` parity and `q=2` case of the
[compressed critical-moment automaton](check_critical_moment_compressed_automaton.py),
the full graph reachable from its initial state is finite and closes.
The [C++ checker](check_critical_moment_compressed_graph_even_q2.cpp)
reimplements the Python state's heights, normalized value order, validity
rules, zero-gap interval transitions, and both constant-column self-loops.
It retains each transition's exact 28-bit pair-distance increment label.

Under explicit caps of one million states, thirty million edges, and sixty
seconds, the search closes at **322,564 states and 2,903,076 labeled
edges**. Every state has nine outgoing transitions: both constants and
one column of each nonconstant plus-set size `1,...,7`. Iterative exact
Tarjan search finds five strongly connected components. Four are transient
singletons; the only closed component has 322,560 states. The initial
state is transient, while the terminal state belongs to the closed
component. The terminal first appears at path length eleven.

One shortest path, encoded as the eight-bit masks of rows carrying a
positive signed column in increasing coefficient magnitude order, is

```text
128, 64, 206, 205, 35, 19, 240, 8, 4, 60, 195.
```

The Python reference checks that path's terminal state and labels.
All sixteen cross-pair distances are five; eight within-group distances
are six, and the four remaining within-group distances are two. Thus
state reachability alone does not give a punctured orthogonal matrix;
individual distance labels are essential.

The closed component has a stronger exact symmetry. For each plus-set
size `k=1,...,7`, its size-`k` transition is a permutation of the
322,560 states: every state has exactly one incoming and one outgoing
edge of that size. Every particular mask of size `k` occurs on exactly
`322560/binom(8,k)` edges. These statements were checked for all nine
sizes `0,...,8` and all 256 masks, including constants.

Consequently an exact integer circulation assigns weight `binom(8,k)`
to every size-`k` edge, with weight one on each constant edge. Flow is
conserved at every vertex, and every one of the 256 signed-column masks
has total weight 322,560. The total weighted length is 82,575,360;
each row pair is separated by exactly 41,287,680 weighted columns.
This is an **isotropic recession circulation** for the linear
weighted-flow relaxation. It defeats a homogeneous separation bound
based only on the compressed graph and pair labels. The circulation
alone does not supply the affine start-to-terminal flow with the required
within-pair offset; the integer correction below supplies that step.

The same conclusion holds if both constant-column masks are forbidden.
Assign weights by cut size
`[0,10,28,56,70,56,28,10,0]`. Every active edge has positive weight.
The uniform label counts and direct 28-pair sum give total length
83,220,480 and every pair distance 41,610,240. Removing constant
self-loops leaves the closed component strongly connected.

There is also no rational linear conservation obstruction in the
28 pair labels plus length. Choose the terminal vertex as root. For
each active edge `u→v`, concatenate fixed active paths from root to
`u` and from `v` back to root. The checker records 29 such closed
walks, including their exact mask paths and label vectors. Those 29
integer vectors have rank 29 modulo 1,000,000,007, hence full rank
over the rationals. Every arbitrary rational length-and-distance
correction to the eleven-column start-to-terminal path can therefore
be expressed as a signed rational combination of active closed walks.
Adding a sufficiently large multiple of the positive active
circulation makes all giant-component edge weights nonnegative. It
increases the target degree `s` by 20,805,120 for each integer multiple.
Rational multiples prove rational weighted-flow feasibility for every
sufficiently large even `s` in this parity case.

An exact integer correction is stronger. Number the 29 recorded cycles
in JSON order and take coefficients

```text
1, 0, -2, 1, 0, 2, -3, 4, -3, 0, 0, 0, -1, 1, 0, ..., 0.
```

The correction changes the eleven-column base path's length by `-1`,
adds four separations to each of the four pairs `(0,1)`, `(2,3)`,
`(4,5)`, `(6,7)`, and changes no other pair distance. The resulting
signed flow has exactly the target labels for `s=2`: length ten,
cross-pair distance five, and within-pair distance six. The cycle
correction alone uses 110 distinct active directed edges; 35 have
negative multiplicity, with minimum `-6`. Including the base path gives
117 touched directed edges, 33 negative multiplicities, and minimum `-5`.

Add **one** copy of the positive active circulation. Every giant-component
active edge then has strictly positive integral multiplicity: the
minimum on giant-component edges touched by the correction is seven,
and every untouched active edge has weight at least ten. The base
path's transient edges remain positive. Flow balance is from the
initial vertex to the terminal vertex, and the support is connected
because the giant component is strongly connected. The directed Euler
trail theorem therefore gives an actual finite **signed-column mask
word** with the prescribed pair distances. More generally every
integer `K≥1` yields such a mask word for
`s=2+20,805,120 K`. At `K=1` the word has length 83,220,490,
cross-pair distance 41,610,245, and within-pair distance 41,610,246.
This is a combinatorial path result only: it gives no pairwise-distinct
coefficient magnitudes solving the odd-moment equations.

[The machine-readable summary](critical_moment_compressed_graph_even_q2_summary.json)
records the graph, path, SCC, circulation counts, the 29 replayable
closed-walk witnesses, and the integral correction.

To reproduce the graph closure and exact assertions:

```text
c++ -O3 -std=c++17 docs/check_critical_moment_compressed_graph_even_q2.cpp -o /tmp/check_compressed_q2
/tmp/check_compressed_q2
python3 docs/check_critical_moment_compressed_graph_even_q2_reference.py /tmp/check_compressed_q2
```

The C++ checker also has `--sample`, which emits sixteen states spread
across the BFS order and all 144 of their labeled transitions. Those
transitions were compared exactly against the Python reference, including
both constant loops and distance masks. The reference checker also
replays the terminal path and all 29 closed walks, checks the
modular rank, and verifies the integral correction and multiplier
`K=1` independently. The separate Python automaton checker
proves its formulas and checks a larger depth-five prefix set against
all 256 signed columns per state.
