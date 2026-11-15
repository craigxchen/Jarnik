# kill_prove notes
## initial review observations (reasoning)
- KEY: Prop 2.5 claim "every integer vector of span is (a,-a) so balance holds exactly" is wrong: balance (Lemma 1.2(3)) is |sum a|<=||a||_1-L,
  which needs ||a||_1>=L too. Theorem-B kernel contains 1_H (H codim-2 subspace with H^perp in A u 0), weight D/4 < L+1. So Theorem-B fakes
  and their doubles VIOLATE lattice Lemma 1.2(1). The lattice form DOES separate those fakes. Test computationally.
- Prop 2.5 "any code": fails if a unit vector lies in span(C-C) (columns k, D0+k coincide). Needs hypothesis.
- Plan: lattice filter on equidistant 6-codes at D=8 (could PROVE no 6-pt equidistant at D=8).
- (S4) "forbids r>=6 for eps<1/2" is normalisation dependent (needs h<=D).
- PROVED (kill_prove_checks/equi_lattice.py, equi_witness.py): all 3 classes of 6-pt equidistant codes at D=8 (20160 cliques, s1=1 frame)
  have a weight-2 integer vector in span of differences -> violates Lemma 1.2(1) -> NO 6-pt equidistant cluster at D=8,c=0 (proof, upgrades lane evidence).
  Witness class 0: e2-e3 = (d2+d3-d4-d5)/2. 5-pt: 4 classes, 2 killed (min l1=2), 2 survive with 6 sigmas: matches lane's "2 classes realised".
- Cube Thm 3.1 CONFIRMED independently (cube_special_mine.py resultants+flint isolation: one nondeg |y|-set {0.12379,0.91118,2.60974,3.39713},
  orbit 48; cube_fibres_mine.py: generic (8,4,4,4) rank 4, special (8,4,4,...) gaps>=4 vs err 1e-42; cube_comb.py combinatorial 8 unique, rank 4).
  Degenerate curves: d=(a,a,b,b) -> y_i+y_j=0; d=(-1,0,0,1)-> y0=y3. all degenerate.
- Prop 2.5 WRONG for lattice form (fake_lattice.py): Theorem-B fake D=8 has 1_H (weight 2<4) in span; double at D=16 has (1_H,-1_H) l1=4<8, a+=2<3.5.
  Also 'any code' needs no unit vector in span(C-C) (C2 example).
- 5-pt equidistant D=8 class 2 CERTIFIED (Krawczyk, krawczyk5.py): box radius 1e-30, sep 0.32, |Im|>=0.28. Fibre profile (6,5,5,4,...):
  6-pt rank-3 fibre exact (diffs in span of cluster diffs). Lane's '(5,5,4,4,...)' profile is wrong for this member.
