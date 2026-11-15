## 6. Lower bounds for $M_\infty(C)$, and what unbounded function-field clusters would imply

**Proposition 6.1 (proved; Theorem 1.2).**
(a) A primitive Q(i)(t)-cluster of size $M$ with $D$ even, $c=0$ and constant $C_\infty$ gives $M_\infty(C)\ge M$ for
every $C>C_\infty$. Sequences realising it: $t_0=r+mn$, $n\to\infty$.
(b) A Q(i)(t)-cluster of size $M$ with $K<D/2$ (for instance odd $D$ at the model's $c=0$) gives $M_\infty(C)\ge M$
for every $C>0$.
(c) Composition $z\mapsto z\circ h$ ($h\in\mathbb Q[t]$) never lowers $C_\infty$: its values are a subset of the
original's values. Rescalings $z\mapsto\lambda z$ and rational affine changes of $t$ leave $C_\infty$ unchanged.

**The cube family alone.** $M_\infty(C)\ge8$ for $C>C_{\rm cube}$, where the best value found is
$C_{\rm cube}\le145877.66$ (at $y=(-\frac94,-\frac23,\frac14,\frac83)$, $t_0\in\frac12\mathbb Z$ in the optimal
class). This is an explicit upper bound, proved by the example in §2. That it is the infimum over the family is
**evidence only**: it is the minimum over all rational $y$ of height $\le25$, and $C_\infty$ grows with height. The
record's $M_\infty(C)\ge8$ for $C>5.89771$ (golden family) is better by a factor of about $2.5\cdot10^4$. **No record
is improved.**

**Unbounded function-field clusters.** Suppose $\sup_DM(D,0)=\infty$.
1. *Over $\mathbb C(t)$ only, i.e. with no $\mathbb Q(i)$-rational members:* no integer consequence (Theorem 1.2
   needs $\mathbb Q(i)$-coefficients). Function-field consequence: clusters with $M\ge9$ and $D$ large have
   5-subsets with $\alpha\le4(1-\phi(M))+o(1)\to15/4$ (§5). If these do not lie in one fixed $Z\subset T_5$, then
   FF-$A_5$, and with it the truncated function-field Vojta conjecture for $(\mathrm{Bl}_eX,\pi^*\partial+E)$, fails.
2. *Over $\mathbb Q(i)(t)$, with constants $C_M$:* $M_\infty(C_M)\ge M$ for all $M$, so $\sup_CM_\infty(C)=\infty$. By
   the absolute form of the A-branch (round 23 J4, refereed), this **refutes $A_5(\delta)$ for every
   $\delta\in(0,1/30)$, and refutes $V_E(5)$**. It does **not** refute $(L_5^\circ)$ (J4's scope remark), and it does
   **not** refute Cilleruelo–Granville, which fixes $C$, unless $\sup_MC_M<\infty$. The function-field picture gives
   no control of $C_M$: $C_M=\Gamma_M/\sqrt{g^*_M}$ is governed by coefficient heights and denominators, which nothing
   in the model bounds. In the only family known, $C$ grows with the height of the parameters (§2).
3. *Bounded function-field clusters ((R$_4$) true)* imply nothing for integers. Integer clusters need not come from
   Q(i)(t)-clusters at all; the records come from unit templates and sporadic points (§3).

So the function-field line can produce at most refutations of the Vojta-type routes ($A_5$, $V_E(5)$), never a
proof, and never by itself a refutation of the conjecture at fixed $C$.
