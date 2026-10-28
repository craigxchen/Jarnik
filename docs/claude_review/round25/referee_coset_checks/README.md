# referee_coset_checks (referee of round25/coset.md)
python3 + numpy/mpmath, run under prlimit --as=3500000000.
* constants.py -> constants_out.json: Cor 2.3 n0/log-mu (incl. corrected tau=sqrt(41/37) for 13..73), Cor 5.2 M_0.
* graph_check.py -> graph_check_out.json: Cor 4.2 formula vs brute force; Thm 5.1 total step mass.
* cover_check.py -> cover_check_out.json: independent certification of the 61 line covers (§8.2).
* char_check.py -> char_check_out.json: character margins of the eta-line up to weight 40.
* density_ref.py -> density_ref_out.json: §8.1 density near the transplant vector.
* seed_check.py -> seed_check_out.json: Lemma 7.2 exact Gaussian test (r=4, m=2, all (pi/4)-shifts).
* notes.md: referee working notes.
