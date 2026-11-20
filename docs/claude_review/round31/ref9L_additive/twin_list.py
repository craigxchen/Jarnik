import sys
from twin import analyse
Cn,Cd=int(sys.argv[1]),int(sys.argv[2]); maxch=int(sys.argv[3])
for q in map(int,sys.argv[4:]):
    st=analyse(q,(Cn,Cd),maxch)
    st.pop('ex_aff') if not st['affine3'] else None
    print(st, flush=True)
