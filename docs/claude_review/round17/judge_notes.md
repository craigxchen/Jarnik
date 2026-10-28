# Judge notes, round 17 (working file)

## Reading status
- read: round16 README, master README, lit18 README, growth README, round15 README, threshold/L8 §5, round10 analytic §6
- read: cached results for all 7 agents + calib/referee:calib (from task text)
- read: hlight §2, §4; survive §2-§4; cap §0, §4-5; calib §0, §7

## Synthesis points
- All four attacks: CIRCLE(q) undecided. Under R16 §7.4 rule -> route 5 STOPS.
- H_light: PROVED (c_0 = 1/6 all q >= 7; c_0 >= 0.2349 at q >= 3251 via Weil; exact 0.2416 at q=3251). Two referee passes.
- survive Thm 3.1/3.4 (Galois realised by multipliers, gap principle): the quantitative test = exact membership
  modulo a height window [~q tau/2, m' q tau/20]. Most valuable new structural result.
- cap: Gauss-sum identity, fat-point form, single families fail; major correction: cap may come from TC_H' \ TC_V.
- calib: Prop 2.1 (= cap Prop 4.1/4.2 independently); kappa_V evidence for 'no'; extrapolation model-dependent.
- Duplicate: calib Prop 2.1 == cap Prop 4.1+4.2 (both proved independently). Credit once.
- Numbers in cap/calib/survive used c_0 = 0.184 (evidence) or the c_0-free floor. With H_light proved, recompute (J1).

## Next step thinking
- low-height window P_+(N) < M^{2-delta} uncapped (analytic §6) but restricted: not a growth route by itself.
- (L_7): expected false, refutation closes route, no progress. (L_8) Vojta-type.
- H^gen+L1: individual-N worst-case statements; Prop C kills averages.
