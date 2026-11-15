## 0. Summary of results

**Status: obstruction plus calibration. No record is improved and no rank bound is proved. The growth rate
$(\frac23+\varepsilon)\log R/\log\log R$ is unchanged.**

1. **Template theorem (proved; Theorem 1.2).** A primitive Q(i)(t)-cluster with $D$ even at $c=0$ specialises at
   rational $t_0\to\infty$ (bounded denominators) to primitive integer clusters with
   $C(t_0)\to\Gamma/\sqrt{|g(t_0)|}$. Here $\Gamma=\max|\text{coeff}_{t^{D/2}}(z_i-z_j)|$ and $g$ is the fractional
   Gaussian gcd of the values. The best constant $C_\infty=\Gamma/\sqrt{\prod_\ell g_\ell^*}$ is a product of exactly
   computable local maxima. It is an invariant under affine changes and rescaling, it is never improved by
   composition, and it is attained along a progression. Clusters with $K<D/2$ (e.g. odd $D$) give $C\to0$; $K>D/2$
   gives nothing.
2. **Cube family as templates (exact).** The smallest constant found is $C_{\rm cube}=145877.658\ldots$, at
   $y=(-\frac94,-\frac23,\frac14,\frac83)$ with $t_0\in\frac12\mathbb Z$. This is $2.2\cdot10^4$ times better than
   the round-29 figure $3.18\cdot10^9$, and still $2.5\cdot10^4$ times worse than the golden $5.8977$. The other 37
   parameter classes of height $\le25$ are all worse, by $\ge10^2$. Sub-clusters: 5, 6, 7 points at $30.9$, $433$,
   $7946$. An explicit 74-digit integer example is exactly checked. **No record improved.**
3. **Records versus function fields (exact codes).**
   * Conjugation symmetry forces the cube code (Lemma 3.2). The record 8-point cluster with $C=5.754$ has
     **exactly the cube code**, with unequal weights.
   * The golden family is a specialisation, but along units, of a level-$D/4$ function-field cluster. Its asymptotic
     code is a sub-code of the cube.
   * The 5.674 record (non-squarefree, rank 7) and the 10-point cluster (rank 7) are not specialisations of any
     known $c=0$ function-field cluster. The 10-point cluster would need a $\ge10$-point Q(i)(t)-cluster at $c=0$,
     which would itself be new.
4. **Question B (exact plus numerical).** All 56 five-subsets of a cube cluster have gcd 1, contact 4,
   $-K\cdot C=16$, so **$\alpha=4$ and $\delta=0$**. They are not A-points and impose nothing on $Z$.
   * They lie on none of the 15 quadric cones and none of the 15 cubic threefolds (exact, via $\sigma$-invariance),
     and off $W$ and $\iota W$ (numerical).
   * 24 of them are Zariski dense in the subtori $\{u_au_b=u_cu_d\}$; the other 32 sweep hypersurfaces.
   * $(L_5^\circ)$, $A_5$ and the reduction $A_5\Rightarrow$ uniform bound are untouched.
   * General dictionary: $\alpha<4$ iff the gcd-reduced 5-point cluster has negative level.
5. **Calibration (proved conditionally on the named conjecture).** The function-field analogue FF-$A_5$ follows from
   the truncated function-field Vojta conjecture for $(\mathrm{Bl}_eX,\pi^*\partial+E)$.
   * **(S4) is corrected:** with $E$ in the divisor the inequality is $r\mu\le|S|-1+\varepsilon h$. It forbids rank
     $\ge5$ at $c=0$ and allows rank 4; Theorem D is its $\varepsilon=0$ free-curve case. So (R$_4$) has exactly the
     Vojta shape.
   * FF-$A_5$ sees only 5-subsets with large gcd, which a cluster has once $M\ge9$ (then $\alpha\to15/4$). It does
     not see the rank.
   * A function-field proof would be weak evidence for Question B. A function-field counterexample over
     $\mathbb Q(i)$ would refute $A_5$ and $V_E(5)$ arithmetically.
6. **Lower bounds.**
   * The cube gives $M_\infty(C)\ge8$ for $C>145877.66$, against the known $C>5.89771$.
   * Unbounded Q(i)(t)-clusters would give $\sup_CM_\infty(C)=\infty$, which refutes $A_5(\delta)$ for all
     $\delta<1/30$ and $V_E(5)$, but not $(L_5^\circ)$. It refutes Cilleruelo–Granville only if the constants stay
     bounded, and nothing in the model controls them.

**Stop rule.** This lane found no rank-5 cluster and proved no rank bound on a positive-dimensional locus. By the
pre-registered rule this lane's contribution does not keep the function-field line open.
