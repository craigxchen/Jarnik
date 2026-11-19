"""Projective elimination certificate for a fixed set of leading cofactor differences."""
from p3_h_plane import calc,elim
import sys
p=int(sys.argv[1]) if len(sys.argv)>1 else 101
for chart in range(3):
    ctx,x,ds,eqs,D,Delta=calc(p,chart=chart)
    R=elim(list(ds.values()),'v'); G=elim(eqs,'v')
    assert R and G and R.is_constant() and G.is_constant()
    print('p',p,'chart',chart,'rank-drop bound',R.total_degree(),'necessary bound',G.total_degree(),'number of necessary forms',len(eqs),flush=True)
