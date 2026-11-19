from fp import *
from construct_q3 import *
Delta = trim([1,0,2]); C = trim([0,2,0,1])
for s in list(range(1,50)):
    o = build(Delta, C, s, [(1,0),(0,1)], [1,1])
    print(s, o['L'], 'g',o['g'],'D=',o['D2'][0]/2, 'maxdiff=',o['maxdiff2']/2, 'crit',o['critical'], 'eqN',o['equal_norms'],'dist',o['distinct'],
          'RinF3[t]',o['R_in_base'],'Rsf',o['R_squarefree'],'N*G=R^2',o['N0_times_normG_eq_R2'],'maxexp',o['N0_max_exponent'],o['N0_mult_degrees'])
