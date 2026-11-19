from sqf import check
import sys
for d in (3,4,5,6):
    m=1<<(d-1); L0=4*m*(d-2)
    for L in (L0, L0+1, L0+2, L0+7):
        for seed in (1,2):
            r=check(d,L=L,seed=seed); print(r,flush=True)
# control: no trace condition, random S
for d in (3,4,5):
    m=1<<(d-1); L0=4*m*(d-2)
    cnt=sum(check(d,L=L0,seed=s,trace_cond=False)['sqfree'] for s in range(200))
    print('control d',d,'random S squarefree fraction',cnt,'/200',flush=True)
