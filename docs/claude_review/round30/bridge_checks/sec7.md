## 7. Status of the coordinator's observations used, and of computations

* **(S1)** Correct; used in §4.1. For $w=D/2$ pairs of the cube, $z_i-z_j=A\cdot2i\,\mathrm{Im}B$ with
  $\mathrm{Im}B$ constant.
* **(S2)** Not used.
* **(S3)** Not used.
* **(S4)** **Corrected** (§5). The right divisor is $\pi^*\partial X+E$, which gives $r\mu\le|S|-1+\varepsilon h$.
  At $c=0$ this allows $r=4$ and forbids $r\ge5$ off $Z_\varepsilon$ ($\varepsilon<1/4$). So (R$_4$) has the Vojta
  shape, and the cube (rank 4, $16=|S|$) violates the $\varepsilon=0$ form by exactly 1, i.e. it is
  Vojta-exceptional. The statement about the gcd theorems is confirmed in order of magnitude (Corvaja–Zannier shape
  $(\deg u\deg v\chi)^{1/3}\asymp D$; seen through search results only).
* **(S5)** Not used.
* **(S6)** **Proved in corrected and sharpened form** (Theorem 1.2):
  - $C_\infty=\Gamma/\sqrt{\prod_\ell g_\ell^*}$ is an invariant of the cluster;
  - only $c=0$ with $D$ even gives a positive constant;
  - odd $D$, or $K<D/2$, gives $C\to0$;
  - "polynomial templates" are only the rational-branch case, and the record's families are not of this type (§3).

**Computations** (all in `bridge_checks/`):

| script | content | status |
|---|---|---|
| `glib.py` | Gaussian-rational polynomials, $\mathbb Z[i]$ gcd | library |
| `cube_basic.py` | cube structure; pair degrees; 5-subsets: $g$, $\mu_0$, $\alpha$, rank, profiles; 15 cones | exact |
| `template.py` | brute-force $C(t_0)$ scan over 12 denominators | exact values, mpmath output |
| `localgcd.py` | exact local maxima $g_\ell^*$, $C_\infty$ | exact (2-adic part: lower bound for $g_2^*$) |
| `ratpts.py`, `ratpts25.txt` | rational points of the cube surface | exact |
| `run_family.py`, `fam25*.txt` | $C_\infty$ over the parameter classes | as `localgcd.py` |
| `example.py`, `example.txt` | explicit integer clusters, sub-cluster constants | exact checks |
| `record_codes.py`, `.txt` | exponent codes of the record clusters | exact factorisations |
| `sweep_rank.py` | swept dimension 3 of each 5-subset family | numerical (40 digits, SVD) |
| `wtest.py`, `wtest.txt` | membership in $W$ | numerical evidence |

Notes: `bridge_notes.md`.
