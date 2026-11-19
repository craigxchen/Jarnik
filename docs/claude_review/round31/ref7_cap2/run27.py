import sys
from cap2 import run
X,Y,XpY,XmY=(1,0),(0,1),(1,1),(1,2)
cfg=sys.argv[1]
if cfg=='m2':
    for L in (174,175): run(27,2,[X,Y],[1,1],L)
if cfg=='m3':
    for L in (377,378): run(27,2,[X,Y,XpY],[1,1,1],L)
if cfg=='m4':
    for L in (907,908): run(27,2,[X,Y,XpY,XmY],[1,1,1,1],L)
if cfg=='astra':
    run(27,2,[X,Y,XpY,XmY],[1,1,1,1],973,Tkind='astra')
