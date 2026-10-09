import cProfile, pstats
from combi import *
import iso, modla
orbs=shell_orbits(8,10,6)
cProfile.run("r=iso.lam_profile((4,3,1),orbs,modla.P1,dguess=28)","prof4.out")
st=pstats.Stats("prof4.out"); st.sort_stats('tottime').print_stats(12)
