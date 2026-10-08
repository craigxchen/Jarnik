# Joint-cluster auxiliary tests with the exact twist content

Two concrete joint tests retain the opposite-twist data. The product of
two chord-line factors pulls back to an ordinary source polynomial and
has only one universally forced gap factor in its values. A complementary
four-row determinant has a further exact divisibility from mixed selected
and retained prime layers. The latter is sharp at a literal equal split
of the logarithmic norm, but its additional factor vanishes for primes
having only one allocation layer. Neither test proves a growth improvement
on arbitrary endpoint arcs.

The earlier [common-twist determinant notes](common_twist_triple_determinant_factorization.md)
concern products of scalar row characters. The determinant here instead
uses the actual two Gaussian factors of each original point, while the
line product uses the two actual descended clusters.

## Two source chord lines after one descent

Let a primitive source tuple of norm `N` have a cut `S|T`, with oriented
factor `alpha`, `K=Norm(alpha)`, and descended points

```text
z_i=alpha x_i (i in S),       z_j=conjugate(alpha) y_j (j in T),
Norm(x_i)=Norm(y_j)=M=N/K.
```

Choose two anchors in each group. Let `ell_X(w)` and `ell_Y(w)` be their
integer signed chord-line determinants on the descended circle, and
`L_S(z),L_T(z)` the corresponding source chord-line determinants. Similarity
of oriented area gives the exact polynomial identities

```text
L_S(z)=K ell_X(z/alpha),
L_T(z)=K ell_Y(z/conjugate(alpha)),
ell_X(z/alpha) ell_Y(z/conjugate(alpha))=L_S(z)L_T(z)/K^2. (1)
```

No approximation of the relative cluster phase is made: the two distinct
arguments on the left retain its full rational rotation
`alpha/conjugate(alpha)`. Using `ell_X(w)ell_Y(w)` at one common descended
argument is a different polynomial, generally containing a large
evaluation against the other cluster's chord line. Replacing that
evaluation by a small within-cluster determinant is not justified.

At a source row in `S`, `K|L_S(z_i)`; at a source row in `T`,
`K|L_T(z_i)`. Thus every value of `L_S L_T` on the tuple is divisible by
`K`. There is no second forced copy from the opposite group. If `c_S,c_T`
are the integer coefficient contents of the two lines, their primitive
product is `F=L_S L_T/(c_S c_T)`. The guaranteed divisor of its values is
only

```text
K/gcd(K,c_S c_T).                                      (2)
```

Gauss's lemma supplies the exact coefficient content `c_S c_T`; extra
evaluated content can occur, but is not implicit in (1).

An actual maximal-gap `2+3` fixture is

```text
N=1105,
z=((4,-33),(9,-32),(12,-31),(24,-23),(31,-12)),
S={0,2},             T={1,3,4},
alpha=4+i,           K=17=max_(cuts) K_cut,
w=((-1,-8),(4,-7),(1,-8),(7,-4),(8,-1)),              (3)
```

where the descended tuple is primitive, distinct, and has norm `65`.
Its rational relative rotation is exactly `(15+8i)/17`. With the two
`S` anchors and the first two `T` anchors, the primitive source lines are

```text
L_S/2=-X+4Y+136,          L_T/3=-3X+5Y+187.
```

At the final row their product is
`57*34=1938=17*114`, with exactly one factor of seventeen. This rules out
a universal second gap charge for this candidate even after exact
coefficient normalization and even for a maximal cut. All points are in
their displayed minor-arc order, but their endpoint chord already excludes
the `C=1/2` arc class. This fixture does not exclude an extra divisibility
theorem restricted to that class.

## A stacked determinant for selected and retained layers

Now partition any set of the original labelled cuts into selected and
retained cuts, as in the
[sequential descent](general_bipartition_opposite_twist_descent.md). Write

```text
z_i=beta_i w_i,       Norm(beta_i)=Q,       Norm(w_i)=M=N/Q.
```

The full labelled tuples `beta` and `w` are Gaussian-primitive: their
allocation layers are respectively the selected and retained layers of
the original primitive tuple. Original row units may be placed in `w`.
For any four distinct source labels define the integer

```text
D=det [Re(beta_i), Im(beta_i), Re(w_i), Im(w_i)]_(four rows). (4)
```

At a split prime let `e=v_p(N)`, `a=v_p(Q)`, `b=v_p(M)`, so `a+b=e`.
Sort the four original allocations as `t_1<=t_2<=t_3<=t_4`. The selected
and retained counts on the same four rows satisfy

```text
u_i=v_pi(beta_i),       v_i=v_pi(w_i),       t_i=u_i+v_i,
u_1<=...<=u_4,          v_1<=...<=v_4.
```

Monotonicity follows because each count is the number of its chosen
thresholds below the original allocation. Put

```text
f_p=min(u_2-u_1,v_2-v_1)+min(u_4-u_3,v_4-v_3).
```

Then the evaluated determinant obeys

```text
v_p(D) >= e-(t_4-t_1)+f_p.                            (5)
```

Divisibility includes the possibility `D=0`. If `G_4` is the Gaussian gcd
of the four original points and `F=product_p p^f_p`, the global statement
is

```text
Norm(G_4) F | D.                                       (6)
```

Thus the extra charge detects selected and retained layers simultaneously
present in either of the two outer allocation gaps of the sampled four
rows. Those gaps may contain other original rows; they need not be single
cut blocks of the full tuple. If every prime has one varying layer, each
minimum is zero, so this extra factor is one.

For proof, convert the real columns to
`(beta,conjugate(beta),w,conjugate(w))`; their determinant is `-4D`.
At the odd split prime their four row valuations are
`(u_i,a-u_i,v_i,b-v_i)`. Expansion gives a lower valuation equal to the
minimum assignment sum, which is

```text
e-max_(four distinct labels r,s,j,k) (u_r-u_s+v_j-v_k).
```

Monotonicity permits the two positive positions to be placed on the two
largest rows and the two negative positions on the two smallest. Exchanging
a lower positive position with a higher negative one cannot decrease the
sum. The remaining independent choices give the maximum

```text
(t_3-t_2)+max(u_2-u_1,v_2-v_1)+max(u_4-u_3,v_4-v_3).
```

Its difference from `e` is exactly (5). The term `e-(t_4-t_1)` is the
prime exponent of `Norm(G_4)`. Full source primitivity removes inert and
ramified norm primes, and `-4` is a unit at every remaining prime. This
proves (6) without assuming generic absence of cancellation.

## Exact sharpness at Q=sqrt(N)

Take `pi=4+i`, `p=17`, and an ambient five-tuple with allocations
`(0,1,2,3,6)` and all row units one:

```text
z(t)=pi^t conjugate(pi)^(6-t).
```

Select the first three threshold cuts. Each is a complete labelled cut
block of this ambient tuple. Then `N=p^6`, `Q=M=p^3`. On the four labels
with allocations `(0,1,2,6)`,

```text
u=(0,1,2,3),             v=(0,0,0,3),
beta=(conjugate(pi)^3,p conjugate(pi),p pi,pi^3),
w=(conjugate(pi)^3,conjugate(pi)^3,conjugate(pi)^3,pi^3).
```

Here `Norm(G_4)=1`, `f_p=1`, and

```text
D=16*4*p*Re(pi^3)*Im(pi^3)=2659072,
v_17(D)=1.                                            (7)
```

This is actual equal-norm Gaussian arithmetic at exactly complementary
half norms, so (5) can be sharp. The five source points lie in one minor
arc: their lifted phases are `(-6,-4,-2,0,6)arg(pi)` and
`12 arg(pi)<pi`, since `Re(pi^3)=52>0, Im(pi^3)=47<52` gives
`3arg(pi)<pi/4`. Their arc constant is large; this is not a target-scale
counterexample. A source chord explicitly excludes `C=1/2` in the checker.

## The size estimate still needs a further arithmetic input

Choose one of the four source rows as `z_0` and replace the last two real
columns in (4) by the real coordinates of

```text
epsilon_i=w_i-(z_0/Q)conjugate(beta_i)
         =(z_i-z_0)/beta_i.
```

This is a real column operation preserving `D`. Its anchor error is zero.
If the source arc has physical length `L`, all other errors have modulus
at most `L/sqrt(Q)`. Laplace expansion on the error columns has three
possibly nonzero terms. Every beta minor is at most `Q`, and every error
minor is at most `L^2/Q`. Hence

```text
|D|<=3L^2<=3C^2 sqrt(N).                              (8)
```

For a nonzero determinant, (6) and (8) give
`Norm(G_4)F<=3C^2 sqrt(N)`. The exact nested-layer factor `F` is an
additional arithmetic restriction in this test; it is not present on
one-layer-per-prime profiles. Nonvanishing is also an additional
requirement, not an automatic consequence of four distinct source points.
Neither the sharp fixture nor the assignment calculation rules out useful
special-value cancellation or additional divisibility on actual small arcs.

The [checker](check_joint_cluster_auxiliary_content_audit.py) verifies
the assignment minimum exhaustively, checks actual multi-prime evaluated
divisibility, and verifies both fixtures with exact contents and norms.
