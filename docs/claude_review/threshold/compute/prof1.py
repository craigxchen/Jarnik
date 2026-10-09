import cProfile, pstats, sys
from combi import *
from iso import *
orbs=shell_orbits(8,10,6)
cProfile.run("r=lam_profile((4,3,1),orbs,P1,verbose=True)","prof.out")
print(r['N'],r['reg'],r['time'])
st=pstats.Stats("prof.out"); st.sort_stats('cumulative').print_stats(15)
