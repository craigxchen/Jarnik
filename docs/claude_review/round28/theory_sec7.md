## 7. Reproducibility: scripts, outputs, labels

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks/`.
Heavy runs used `prlimit --as=3500000000`, at most two of my processes at a time. Boxes are built by
`tlib.build_box(n, slabs)` exactly as the round-27 killer's `klib.build_box` (patterns `syl(3)`, `sdc(7)`; random boxes
by the killer's `k8/k11` seeds and generator, re-implemented in `t6/t7`).

| script -> output | content | label |
| --- | --- | --- |
| `tlib.py` | library: standard-monomial least spaces mod p, shift sums, exact `J`-membership (fmpz), power obstructions | — |
| `t1_profiles.py` -> `t1_{syl8,sdc7}_L{1,2,3}_out.txt` | full defect profiles w.r.t. `L_pair`; graded `H_A`, `H_B`; linkage identity | mod p (`P1`) |
| `t2_syl8_invariant.py` -> `t2_syl8_L2_e6_out.txt` | `GL(3,2)` (168), 8 invariant sextic orbits, invariant parts of `J(A)_6`, `J(A_l)_6`; obstruction | **exact** |
| `t3_firstorder.py` -> `t3_syl8_L2_out.txt` | first-order closure `J1(A,l)` vs `J(A_l)`, `e <= 8` | mod p |
| `t4_firstorder_inv.py` -> `t4_out.txt` | obstruction lies in `J1(A,l)_6` (invariant parts) | mod p |
| `t5_identify.py` -> `t5_out.txt` | obstruction `== c sigma^6 mod J(A)^G` | **exact** (sympy) |
| `t6_powers.py` -> `t6_syl8_2_e6_out.txt`, `t6_misc_out.txt` | power obstructions `phi^e`, `phi` ternary | mod p |
| `t7_bco_scan.py` -> `t7_bco_scan_out.txt` | `h(sigma)`, rigidity of the top level, 6 models + 27 random boxes | **exact** (enumeration) |
| `t8_subset.py` -> `t8_syl8_L2_sig1_out.txt` | `A cap {|sigma| <= 2}` on M8: reach 7 | mod p |
| `t9_subset_greedy.py` -> `t9_syl8_L2_out.txt` | greedy removal of extreme levels; `e = 7` defect spanned by powers | mod p |
| `t10_bco_ip.py`, `t10b_h_only.py` -> `t10_k3_out.txt`, `t10_out.txt`, `t10c_out.txt`, `t10b_6_12_out.txt` | Sylvester `n <= 128`: `h`, rigidity, box-zonotope (IP, HiGHS) | exact when status 0 (IP optimality; floating-point solver) |
| `t11_zeta_pair.py` -> `t11_out.txt` | `zeta_pair` by exhaustive DFS over `L_pair` generators | **exact** search |
| `t12_small_search.py` -> `t12_out.txt` | ~400 small boxes (`n <= 6`): no `reach > zeta_pair` | mod p / exact search |
| `t13_staircase.py` -> `t13_*_out.txt` | staircase (standard-monomial) criterion | mod p |
| `t14_linkage.py` -> `t14_syl8_L2_out.txt` | linkage form of `(SHIFT_e)`, `e = 5,6,7`; linkage identity in `C_l` | mod p |
| `t16_random_rigidity.py` -> `t16_out.txt` | (R) on random (unbalanced) sign-column profiles, `n <= 128` (IP) | exact when status 0/2 |
| `t17_product_check.py` -> `t17_out.txt` | Theorem 2.1 sanity: reach and `H` of small products add/convolve | mod p |
| `t15_sigma_exact.py` -> `t15_out.txt` | `sigma^6 notin J(A)_6`, `sigma^7 in J(A)_7`, `sigma^6 in J(A_l)_6`; top level | **exact** |

Notes file: `theory_notes.md`. Everything here is **unrefereed**.
