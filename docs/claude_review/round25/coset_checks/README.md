# coset_checks (round 25, coset.md)
All scripts: python3, numpy; run under `prlimit --as=3500000000`.
* common.py: monomial enumeration (up to sign), exact window half-widths.
* density_check.py -> density_check_out.json, density_check_log.txt: §8.1 (Lemma 2.1 evidence).
* line_cosets.py -> line_cosets_out.json, line_cosets_log.txt: §8.2(a).
* line_verify.py -> line_verify_out.json: cover margins and survivor locations, §8.2(a).
* line_batch.py -> line_batch_out.json; batch_margins.py -> batch_margins_out.json: §8.2(b).
* seeded_base.py -> seeded_base_out.json: §8.2(c).
* quant_bound.py -> quant_bound_out.json: Corollary 2.3 numbers, §8.3.
* barrier_M0.py -> barrier_M0_out.json: Corollary 5.2 thresholds.
* graph_formula_check.py -> graph_formula_check_out.txt: Corollary 4.2 exact check.
