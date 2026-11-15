## 2. The cube family as templates: the smallest $C$, and explicit integer examples

**Method (exact up to the stated search ranges).** By Theorem 1.2,
$C_\infty(y)=\Gamma(y)/\sqrt{\prod_{\ell\in P}g^*_\ell(y)}$. Each local maximum $g^*_\ell$ is computed exactly
(`localgcd.py`):
* for odd $\ell$, by the nearest-candidate argument in the proof of Theorem 1.2(2): the maximum is attained at a
  local image of a root (split $\ell$) or at $\mathrm{Re}\,\rho$ (inert $\ell$). Only finitely many candidates;
* for $\ell=2$ (ramified), by enumerating $t\in2^{-j}\mathbb Z/2^{k}$ with $k\le14$, together with the candidates
  $\mathrm{Re}\rho\pm\mathrm{Im}\rho$. A truncated enumeration can only *underestimate* $g_2^*$, so every stated
  $C_\infty$ is a valid upper bound, realised along an explicit progression.

A brute-force scan over 12 denominators $q\mid72$ (`template.py`) reproduces the optimum at $y_{\rm J}$.

**Parameters.** The rational points of the cube surface $\sum y=0$, $\prod y=1$ of height $\le25$ were enumerated
(`ratpts.py`, 76 points; after $y\mapsto-y$ and ordering, 38 classes). The visible one-parameter family
$\{\pm a,\pm1/a\}$ is degenerate for every $\varepsilon$-assignment (repeated roots; proved by inspection:
$1/a+ia$ and $-1/a-ia$, or $1/a\mp ia$, coincide up to sign or conjugation). The three $\varepsilon$-assignments
of a point gave the same $C_\infty$ in every case tested.

**Results** (`run_family.py`, output `fam25.txt`):

| $y$ | height | $\Gamma$ | local factors $g^*_\ell\ne1$ | $C_\infty$ |
|---|---|---|---|---|
| $(-\frac94,-\frac23,\frac14,\frac83)$ | 9 | 118.31 | $2{:}\,1.73\cdot10^{-4}$, $3{:}\,1.52\cdot10^{-4}$, $5{:}\,25$ | **145877.66** |
| $(-\frac{25}6,\frac23,\frac{18}5,-\frac1{10})$ | 25 | 527.9 | $2{:}\,2$, $3{:}\,1.5\cdot10^{-4}$, $5{:}\,3\cdot10^{-6}$ | $1.89\cdot10^7$ |
| other 36 classes of height $\le49$ | | | | $\ge1.9\cdot10^7$ (all) |

The denominators of the roots, $r_k=\varepsilon_kd_k/n_k+in_k/d_k$ for $y_k=n_k/d_k$, cost a factor
$\prod_k(n_kd_k)^2=(\prod n_k)^4$ in $|g|$. The gcd gains at the primes dividing root differences recover only a
bounded part of it. That is why $C_\infty$ grows quickly with height. **Smallest value found:
$C_{\rm cube}=145877.658383\ldots$**, a factor $2.2\cdot10^4$ below the unreduced $3.18\cdot10^9$ of ff Cor C′. The
gain comes entirely from the gcd reduction and from sampling $t_0\in\frac12\mathbb Z$.

**Sub-clusters** (same $y$, all subsets, `example.txt`). The best $k$-point sub-cluster constants are
$C_\infty^{(5)}=30.88$ (subset $\{0,A^c,B^c,P,P^c\}$), $C_\infty^{(6)}=$ 432.87 (subset $\{0,1,B,B^c,P,P^c\}$), $C_\infty^{(7)}=$ 7946.5 (drop $A^c$). Subsets of
size $\le4$ that share a root have $C\to0$: their reduced level is negative (4 points with $C\to0$, already known).

**Explicit integer example (exact, `example.py`).** Take $y=(-\frac94,-\frac23,\frac14,\frac83)$ and $t_0=13999/2$,
in the optimal class. Dividing the values by their fractional Gaussian gcd gives 8 distinct Gaussian integers with
gcd 1 on $|w|^2=R^2$, where
$R^2=76715325307615615132388329160620427624067222429804201398681082902106198205$ (74 digits, squarefree, 19 primes,
all $\equiv1\bmod4$). The constant is $C=145877.668212\ldots$. The points are
$(-6193356331893854325752232243998600469,\ 6193356331893549048567755974257051838)$ and its image under
$(x,y)\mapsto(-y,-x)$, and likewise for
$(\dots785479872085784293288446,\ \dots617894447902436590554517)$,
$(\dots639625293685896839599381,\ \dots763749026302324556008638)$ and
$(\dots701743192057945485163746,\ \dots701631127930276532348117)$; the full digits are in `example.txt`. The
symmetry is the conjugation symmetry of Lemma 3.2, after a unit rotation. All norms are checked equal, distinctness
and primitivity are checked, and $C$ is computed to 60 digits. At $t_0=(10^{12}+3999)/2$ and
$(10^{30}+423)/2$ (200 and 488 digits) the constant is $145877.658383\ldots$, equal to $C_\infty$ to 9 digits.
