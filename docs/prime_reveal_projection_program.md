# Order-free prime reveal via projections

## 1. Successful fiber

Let

\[
A=\left\{\varepsilon\in\{\pm1\}^r:
\sum_{j=1}^r \varepsilon_j\theta_j\in I\right\},
\qquad |A|=M,
\]

where the angles arise from Gaussian prime factors of the circle norm and
`I` has endpoint angular width.

For a subset `S` of coordinates, let `pi_S A` denote the projection of `A`
to the coordinates in `S`. This is canonical and independent of any reveal
order.

## 2. Projection as a union of predecessor windows

If `T` is the forgotten coordinate set, a projected sign vector on `S` lies
in `pi_S A` exactly when at least one orientation assignment on `T` repairs
its partial angle into `I`. Hence `pi_S A` is contained in a union of at most
`2^|T|` translated windows for the smaller Gaussian product supported on
`S`.

For one forgotten prime `p_i`,

\[
\pi_{-i}A\subseteq
\left\{\eta:\Theta_{-i}(\eta)\in I-\theta_i\right\}
\cup
\left\{\eta:\Theta_{-i}(\eta)\in I+\theta_i\right\}.
\]

The smaller circle has radius `R/sqrt(p_i)`, while each predecessor window
has physical length `C sqrt(R)/sqrt(p_i)`. Relative to the endpoint scale of
the smaller circle, this is smaller by the factor `p_i^{-1/4}`.

Thus forgetting a prime converts the final endpoint cluster into at most two
strictly subcritical clusters on the smaller circle.

## 3. Shearer inequality

For every finite `A subset prod_i X_i`, the coordinate projections satisfy

\[
|A|^{r-1}\le \prod_{i=1}^r |\pi_{-i}A|.
\]

More generally, any fractional cover of the coordinate set yields a weighted
projection inequality. This is the natural order-free replacement for
following one reveal path.

Combining the one-coordinate predecessor description with a subcritical arc
bound gives schematically

\[
|\pi_{-i}A|\ll_C 1+\frac{\log R}{\log p_i}.
\]

Therefore

\[
(r-1)\log M
\lesssim_C
\sum_i \log\left(1+\frac{\log R}{\log p_i}\right).
\]

This estimate alone is not expected to prove a constant bound, but it is a
real global inequality involving all split primes simultaneously.

## 4. Multi-coordinate projections

For a forgotten set `T`, the projection is contained in at most `2^|T|`
translated windows on the circle obtained by dividing out

\[
P_T=\prod_{i\in T}p_i.
\]

Each window is subcritical by the normalized factor `P_T^{-1/4}`. A fixed-k
subcritical theorem should then bound each window by roughly

\[
O\left(1+\frac{\log R}{\log P_T}\right),
\]

leading to

\[
|\pi_{T^c}A|
\lesssim_C
2^{|T|}\left(1+\frac{\log R}{\log P_T}\right).
\]

The factor `2^|T|` is the repair-choice cost. The conductor gain is through
`log P_T`. The central optimization problem is to choose a fractional family
of forgotten sets `T` so that Shearer's inequality amplifies the conductor
gain faster than the repair-choice cost.

## 5. Why this matches the growing-circle intuition

A final point has many paths through the lattice of partial prime products.
Projection forgets the path and records all possible ancestors at a chosen
smaller circle. The projection inequalities then force the same final fiber
to be simultaneously realizable through many different smaller-circle
ancestor systems.

This is stronger than analyzing one reveal order, but unlike raw order
commutativity it produces a quantitative inequality.

## 6. Concrete theorem target

Prove a weighted Gaussian projection theorem of the form

\[
\log |\pi_{T^c}A|
\le
|T|\log 2
+\Phi\left(\frac{\log R}{\log P_T}\right)
-C_0\,\mathsf{Overlap}(T),
\]

where the overlap term measures the impossibility of all `2^|T|` predecessor
windows being simultaneously populated by Gaussian divisor points.

Without the negative overlap term, the method likely yields only a
polylogarithmic bound. A uniform bound would follow if one can choose a
fractional cover for which the accumulated overlap penalty cancels the
binary repair cost.

## 7. Formalization plan

1. Define coordinate projection of a successful fiber.
2. Prove the one-prime projection is contained in the union of two shifted
   predecessor fibers.
3. Formalize finite Shearer/Loomis-Whitney for Boolean fibers.
4. Connect predecessor fibers to the existing fixed-k subcritical theorem.
5. Optimize weighted families of forgotten prime sets.

The additive face curvature remains useful as the two-coordinate infinitesimal
version of this projection theory, but the projection inequalities are the
more global invariant.
