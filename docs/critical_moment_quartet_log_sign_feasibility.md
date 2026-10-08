# Quartet log signs on exact compressed Euler words

The particular magnitudes constructed here fail the exact two-quartet
tradeoff in [the common-jet note](critical_moment_common_jet_products.md),
equation (19). The stronger [multiple-block construction](critical_moment_exact_crossproduct_relaxation.md)
realizes all exact scalar cross ratios and the first moment with different
positive real magnitudes. It still does not supply higher moments or
identify the auxiliary coefficients with the actual middle coefficients.

The 36 strict quartet log-product inequalities in [the common-jet note](critical_moment_common_jet_products.md), equation (16), are **simultaneously feasible** with the first-moment equations on some exact-distance mask words from [the even-`q=2` Euler construction](critical_moment_compressed_graph_even_q2.md), when their signs are read from the automaton's prescribed spectral order. The constructed magnitudes are not claimed to produce that order for their actual middle coefficients. This is a compatibility result for two necessary conditions, not a solution of the full odd-moment equations or the polynomial product identities.

For `i<j` in rows `0,...,3` and `k<l` in rows `4,...,7`, let `A_q,B_q` be the two quartet split classes in the common-jet note. For a prefix of `t` columns, put `h_q(t)=#(A_q in prefix)−#(B_q in prefix)`. Direct inspection of one signed column gives the exact identity

```
h_q(t)=[D_il(t)+D_jk(t)−D_ik(t)−D_jl(t)]/2,          (1)
```

where `D_uv(t)` is the prefix Hamming distance of rows `u,v`. Because the complete exact-distance word has `|A_q|=|B_q|`, summation by parts gives

```
log(prod_(B_q) a_h / prod_(A_q) a_h)
  =sum_(t=1)^(n−1) h_q(t)(log a_(t+1)−log a_t).          (2)
```

The prescribed coefficient order is `6<7<0<1<2<3<4<5`. For all 36 ordered quartets above, the coefficient sign `epsilon_q=sign(w/v)` is `+1`: the two factors in `w=(c_i−c_j)(c_k−c_l)` and those in `v=(c_i−c_l)(c_j−c_k)` have the same product sign, whether `k,l` lie both below, both above, or straddle the first group in this order.

There is an exact terminal-rooted closed-flow label with a strict positive quartet prefix. Take the 29 active closed walks stored in [the graph JSON](critical_moment_compressed_graph_even_q2_summary.json), and combine them with the integer coefficients in the [checker](check_critical_moment_quartet_log_sign_feasibility.py). Its length and all twelve within-group pair distances are zero; its cross-pair distance correction is

```
L_ik=−4i(k−4),       0<=i<=3, 4<=k<=7.              (3)
```

Adding four copies of the positive active circulation makes every active edge of this first closed flow strictly positive. The complementary closed flow is the earlier integer `s=2` correction minus (3), plus `K−4` active circulations; it too is strictly positive on every active giant-component edge for every integer `K>=7`. The checker replays the 29 walks, checks (3), and verifies positivity for `K=7` on every touched edge (untouched edges have positive circulation weight). Each flow has a terminal-rooted Euler circuit because the giant component is strongly connected and its entire active edge set has positive multiplicity.

Concatenate the eleven-column initial-to-terminal path, the first circuit, and the complementary circuit. The resulting word has the required pair distances at `s=2+20,805,120 K` for every integer `K>=7`, with length `n=4s+2`. At the proper prefix `t_0` just after the first circuit, the initial path contributes equal cross distances of five and the active circulation contributes the same distance to every pair. Thus (1) and (3) give

```
h_q(t_0)=2(j−i)(l−k)>0                              (4)
```

for **all 36 quartets at once**. The compressed state at that prefix is the terminal height `delta=(0,0,0,0,−1,−1,1,1)`.

Both circuits use every active giant-component edge. Hence all fourteen axial heights `delta±e_r` occur at proper prefixes, and [the barycentric construction](critical_moment_first_moment_barycentric.md) assigns strictly increasing positive integers `a_1<...<a_n<=T=16n` satisfying the eight first-moment equations. Now replace every magnitude after `t_0` by `a_h+H`. The first moments still agree: the change in row `r` minus row `0` is `2H[delta_r−r_r(t_0)]=0`. The magnitudes remain positive, integral, and strictly increasing.

As `H` grows, the product ratio in (2) grows as a positive constant times `H^h`, where `h=h_q(t_0)>=1`. An explicit simultaneous choice is any integer

```
H >= (2T)^n+T.                                        (5)
```

Indeed `H>=T` makes shifted suffix nodes lie in `[H,2H]`, while unshifted nodes lie in `[2,T]`. Since `B_tail−A_tail=h_q(t_0)`, the ratio of products in (2) is at least `H^h/(2^n T^n)>1`. Thus all 36 required strict log signs hold alongside the first moments on this exact-distance word in every stated degree. The numerical values `c_i` and the full quartet polynomial identity need not be realized by these magnitudes; higher odd moments remain unproved.
