# Referee notes for L7.md (round 7, second attempt)

## Plan
1. Independent exact recomputation of Prop 5.4 (comp counts on 1260 relabelled V_A members + mu), own code.
2. Independent check of Theorem B class table.
3. Logic checks of Sec 4 (Lemma 4.2, Cor 4.3, Prop 4.4, Cor 4.5), Prop 5.2, Lemma 5.3.
4. Scope: Sing(F_A) avoidance hypothesis dropped in summary; "moving section" claim.

## Observations so far
- r7a Thm 3 includes a density hypothesis; its proof uses the density hypothesis, not openness of G1.
  L7 says "its proof uses that it is [open]" and "Thm 3 is not correct" -> mischaracterisation (minor);
  the false part is r7a Remark 1 (G1-G3 fail only on codim-2 locus).
- Prop 5.2 proof uses Lemma 4.2 for images of general t_j-lines; needs them to avoid Sing(F_{V_k}),
  and lambda(gamma_k)=2 needs leaves to avoid Sing. 2D analog: P^2, lines through A, product chart
  (slopes from A and B): lambda(gamma_1)=1 not 2. Gap.
- Cor 4.5 "as every density argument needs [moving]" unjustified: density can come from varying C.
  BUT for same-component sections sigma' ~ Sigma + mF: sigma'^2 = lambda(C)+2m >= m >0 automatically,
  so Prop 4.4-type additions are covered without 'moving'.
- Summary drops "C' avoids Sing(F_A)" hypothesis. Curves through Sing get Sigma^2 <= comp(C') (saturation degree),
  which can exceed lambda(C').

## Results (referee's own code reflib.py, no author imports)
- t1_arrangement.py: (A1) ok; coaxial triples exactly 3; 3 points on 5 planes; 8 simple triple points (cube); web dim 4. CONFIRMED.
- t2_member.py: members through x=(6,-6,-3,-5,2,0,1), (3,-2,4,0,2,-5,-6), (1,7,-4,2,-9,5,3): good, class table of Thm B exactly; comp=2 on leaf. CONFIRMED.
- t4_mu.py seed 11: comp_A(mu) = -10 on 4 random filtered mu curves. CONFIRMED.
- t3_bound.py x=(6,-6,-3,-5,2,0,1): 1260 distinct classes, sum comp = -9576, avg -38/5, min -12, max 2;
  forgetful validation (Lemma 5.3 count vs Lemma 4.1 formula) 0 mismatches over 8820 (label,class) pairs.
  => Prop 5.4 bound lambda_A(beta) <= (35/12)(-38/5)+(1/12)(-10) = -23 CONFIRMED independently.
- t5_degen.py: non-unique-web-quadric locus in G = 6 sigma_2(P_e) + sigma_2(v13) + sigma_11(H_a1) + sigma_11(H_c1) + sigma_11(Pi),
  class 7 s2 + 3 s11 (matches Porteous c1^2-c2). All in boundary.
- t6/t7: m = P1∩P2 for partition {a1b1c1,a1b2c2,a2b2c1}|{a1b1c2,a1b2c1,a2b1c1}: m meets l_2 only => x_m in D_{O b1 b2} only
  (NOT in ∂^(2)); web quadric through m unique = P1∪P2 (rank 2). Pencils of lines through m (good, same class as generic pencils:
  3 coaxial triples + 12 non-coaxial pairs, 1 each) give comp_A = -1 (10/10), generic pencils comp_A = -2 (25/25).
  By Lemma 5.3 => x_m ∈ Sing(F_A) \ ∂^(2). So 'avoid Sing(F_A)' is NOT automatic for good curves.
- Prop 5.2: 2D analog P^2, slopes from A and B: lambda_1(gamma_1)=1 (leaves through radial point A), not 2. Proof gap.
  Also round4 Thm 2(b) allows dominant (P^1)^R, any R; Prop 5.2 / table 9.1 only cover birational R=4.
- t3 with members x=(3,-2,4,0,2,-5,-6) and (1,7,-4,2,-9,5,3): all 1260 per-class values identical to first member (evidence for equality).
- t8_muT.py: delta_T list reproduced exactly (+1: 3 coaxial, 6 events, Oa1c2, Oa2c1; -1: a2b1c2,a2b2c2,a2b1b2,b1b2c2), 3comp(mu)+sum = -23.
- t9_misc.py: det(web quadric) = 4900*(linear)*(cubic) squarefree -> double cover irreducible (closes Thm B class-constancy gap);
  lambda_{f_i}(tau gamma_A) in {-4,-3,-2}; lambda_{f_i}(beta)=-9.
- t10_variants.py: variants exact over Z: -23, -31, -32, -32 (comp(mu) -10,-14,-16,-16). Reproduced (still evidence: negative-coefficient decomposition).
- Web search on Thm A novelty: nothing specific found (inconclusive).
