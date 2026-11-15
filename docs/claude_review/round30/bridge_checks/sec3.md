## 3. The record's 8- and 10-point clusters: specialisations or sporadic?

"Specialisation of a c=0 function-field cluster" means $W=\{z_j(t_0)/g(t_0)\}$ for a primitive Q(i)(t)-cluster with
$D\ge2$, $c=0$, and some rational $t_0$.

**Lemma 3.1 (proved).** If $W$ is such a specialisation, then $|W|\le M$, the size of the function-field cluster.
For every Gaussian prime $\pi$ outside the finite set of Theorem 1.2(2), the $\pi$-column of $W$'s exponent code
equals a column of the function-field code, up to complement.
*Proof.* A good $\pi$ divides $t_0-\rho$ for at most one root $\rho$, and then the points containing $\rho$ (resp.
$\bar\rho$) have $v_\pi=v_\pi(t_0-\rho)$ (resp. $0$). Distinct points of $W$ come from distinct $z_j$. ∎

**Lemma 3.2 (proved; conjugation symmetry forces the cube code).** Let $W$ be an 8-point cluster on a squarefree
normalised $n$, closed under $w\mapsto\bar w$ with no point on the axis. Then every layer cut separates each of the 4
conjugate pairs. So the cuts are among the $2^4/2=8$ "transversal" cuts of the 4 pairs, and the code has rank $\le4$.
The cube code (ff Thm C) is the code with all 8 transversal cuts. The cube cluster is itself closed under
$z\mapsto\bar z$ (pairs $0/1$, $A/A^c$, $B/B^c$, $P/P^c$), so all of its specialisations are axis-type clusters in the
sense of `families.md` §3.
*Proof.* $v_\pi(\bar w)=v_{\bar\pi}(w)=1-v_\pi(w)$ when $p=\pi\bar\pi\,\|\,n$. In the cube, $\bar z_F=z_{F^c}$ because
conjugating every chosen root gives the complementary transversal. ∎

So the cube templates are axis clusters: $C^2\to8d_{\max}+4\lim y_0^2/a$ (`families.md` §3.1), with integer depth set
$\{d_j\}$. Their huge constants mean depth sets of size about $C^2/8\approx2.7\cdot10^9$. Corollary S
($C\ge4\sqrt2-o(1)$ for every symmetric 8-point family) applies to them.

**The record clusters** (`record_codes.py`, exact factorisation in $\mathbb Z[i]$; output in `record_codes.txt`):

| cluster | $n$ | $C$ | rank | cuts | structure |
|---|---|---|---|---|---|
| 8 pts, $a=590749322$ | $5\cdot13\cdot17\cdot29\cdot53\cdot181\cdot457\cdot733\cdot3389$ | 5.7542 | 4 | 8 | **exactly the cube code**: the 8 transversal cuts of the conjugate pairs, unequal weights $\log65,\log17,\dots,\log3389$ |
| 8 pts, $a=5743586$ (record) | $5^3 13^2 17^2\cdot89\cdot109\cdot557$ | 5.6739 | 7 | 10 | non-squarefree; layers of $5,13,17$ break the pair separation |
| golden, $a=23186$ / $7465178$ | $5^2\cdot17\cdot61\cdot89\cdot233$ / $5^2\cdot13\cdots1597$ | 5.8978 / 5.8977 | 5 | 6 | 4 heavy cuts (the blocks $J_{n+s}$) are transversal cuts; the two layers of 5 ($F$ merged with a block) are not |
| 10 pts | $5^3 13^2\cdot17\cdot29\cdot113$ | 8.2393 | 7 | 8 | no perfect matching separated by all cuts |

**Answers.**
1. *10 points: not a $c=0$ specialisation of any known function-field cluster.* By Lemma 3.1 it would need a
   Q(i)(t)-cluster of size $\ge10$ at $c=0$. None is known: $\sup_D$ of known values is 8, and $M(D,0)\le6$ for
   $D\le6$ and $\le8$ for $D=7$ (exact or Plotkin). Under (R$_4$) the LP would give $\le8$ for large $D$. So the
   10-point cluster is either sporadic or a specialisation of a function-field cluster that would itself be a new
   record ($M(D,0)\ge10$). The referee's diagonal-symmetric search had already found it unique up to $n\approx5\cdot10^{19}$.
2. *The golden family is a specialisation, but of a different function-field object.* Put $X=\varphi^n$. Then
   $J_{n+s}=a\varphi^sX+\tau(a)\varphi'^s(\pm1/X)$, so $Z_w(X)=X^{-4}z_w(X)$ with $z_w\in\mathbb Q(\sqrt5,i)[X]$ of
   degree 8. The roots come in pairs $\pm\rho_s$, the same $t\mapsto-t$ symmetry as the cube. The family is evaluated
   at the units $X=\varphi^n$, not at integers. Because of the factor $X^{-4}$, $|Z|\asymp X^{4}=X^{D/2}$, so bounded
   $C$ needs $\deg(z_w-z_{w'})\le3D/4$: **level $c=D/4$** ($L=1$ at $D=8$), not $c=0$. (Proved by the computation in
   `families.md` §2: angular spread $\asymp\varphi^{-2n}=X^{-2}$, so $\deg(z_w-z_{w'})=6$.) Its asymptotic code is the
   cube code restricted to the 4 transversal cuts "one $e$-side". In general, Pell/unit templates, which specialise
   along a pair of conjugate real-quadratic places, are function-field clusters at level $K=3D/4$ plus an
   involution. They are cheaper than polynomial templates by a factor 2 in the contact. This is why they hold the
   records.
3. *The 5.754 cluster has exactly the cube code (exact), with unequal weights.* Combinatorially it is a cube, as is
   every squarefree symmetric 8-point cluster of rank 4 that uses all 8 transversal cuts (Lemma 3.2). It cannot lie
   in the asymptotic range of a cube template from the explored parameters, since those have $C_\infty\ge1.46\cdot10^5$
   (§2). Whether it equals a cube value $z_F(t_0)/g$ at a small $t_0$ for some rational $y$ of large height was
   **not decided**. A finite $t_0$ has no bearing on $M_\infty$ in any case.
4. *The 5.674 record* has non-squarefree $n$ and rank 7, so it is not even combinatorially a cube. Small primes with
   multiplicity carry its structure. The table indicates it is **sporadic** with respect to the function-field
   picture; this is evidence, not a theorem.

**Bottom line for (1).** The polynomial ($c=0$, rational-branch) templates are the *weakest* mechanism. The record's
8-point families arise from unit templates (level $D/4$) or from axis surfaces (del Pezzo $X_{0134}$). Their
combinatorics is the cube code (Lemma 3.2), which explains why the cube is the function-field optimum's shape, but
the arithmetic realisation is different.
