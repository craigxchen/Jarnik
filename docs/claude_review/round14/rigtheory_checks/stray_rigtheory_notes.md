- K=3 eps, exact levels 3:4,7:2 (joint2.py): q<=31 no indiv-feasible non-eps-class label; 37 INFEASIBLE; 41+ feasible.
  exact all-level witnesses (jointG, G={3,7}) outside every eps-class at q0=41,43,47,53,71,103. => Core_J^eps(3) = [3,37] (pending 59,61,67,83 witnesses).
- K=4 fixed-eps with G={3,7} genuine all levels: INFEASIBLE at q0=11..61 and 73 (joint rigid given G); feasible 67,71,79 but lift at 11 failed.

# appended copy at 14:28
- K=4 witnesses (G={3,7,11}) at 67, 71 verified. q_core^{J,eps}(4)=61 given 3,7.
- COHERENCE: K=3 witness for q0=41 uses eps (1,-1,1,-1)@19, (1,1,-1,1)@31, (1,1,1,-1)@37 -> incoherent on [19,37]; coh.py K=3 S=[11..37] feasible (incoherent).
  K=4 witness q0=67 is coherent (eps=1) on 19..61. coh.py K=4 S=[11..61], G={3,7} running.
