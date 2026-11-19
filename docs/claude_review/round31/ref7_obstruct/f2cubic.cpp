// F2: squarefree odd sign functions sigma = chi(Q(x)) with Q an anisotropic (irreducible) binary cubic form over F_q.
// marks U_Q = prod_{v in V\0, chi(Q(v))=1} (X - v); cluster search with one common divisor.
#include "gf.h"
typedef GFlog F; typedef F::E E; typedef Poly<F> P;
int main(int argc,char**argv){ int q=atoi(argv[1]), r=atoi(argv[2]); long cap=atol(argv[3]); vector<int> m; for(int a=4;a<argc;a++) m.push_back(atoi(argv[a]));
  F f; f.init(q,m); ll qr=1; for(int j=0;j<r;j++) qr*=q; E g=f.gen(); E I=f.pw_(g,(qr+1)/2), th=f.pw_(g,qr+1);
  ll n=q*q; vector<E> V(n); vector<array<int,2>> X(n); for(ll idx=0;idx<n;idx++){ X[idx]={(int)(idx%q),(int)(idx/q)}; V[idx]=f.add(f.mul(f.fromint(X[idx][0]),I),f.mul(f.fromint(X[idx][1]),f.mul(I,th))); }
  auto chi=[&](int y){ y%=q; if(y<0)y+=q; if(!y) return 0; ll e=1; for(int j=0;j<(q-1)/2;j++) e=e*y%q; return e==1?1:-1; };
  set<vector<int>> seen; vector<vector<int>> sets; 
  for(int a=1;a<q;a++) for(int b=0;b<q;b++) for(int c=0;c<q;c++) for(int d=0;d<q;d++){ bool aniso=true; vector<int> sel;
    for(ll idx=1;idx<n;idx++){ ll x=X[idx][0], y=X[idx][1]; ll v=(a*x*x*x + b*x*x*y + c*x*y*y + d*y*y*y)%q; if(v==0){aniso=false;break;} if(chi(v)==1) sel.push_back(idx);} 
    if(aniso && !seen.count(sel)){ seen.insert(sel); sets.push_back(sel);} }
  int nm=sets.size(); printf("q=%d distinct anisotropic-cubic marks=%d D=%lld\n",q,nm,(n-1)/2);
  vector<P> pol(nm); for(int t=0;t<nm;t++){ pol[t]={f.one()}; for(int idx:sets[t]) mullin(f,pol[t],V[idx]); }
  int D=(n-1)/2; vector<vector<short>> dd(nm,vector<short>(nm)); int minpair=1e9;
  for(int s=0;s<nm;s++) for(int t=s+1;t<nm;t++){ dd[s][t]=dd[t][s]=deg(f,psub(f,pol[s],pol[t])); }
  int best=1; long nodes=0; vector<int> bs;
  function<void(vector<int>&,vector<int>&,int,int)> dfs=[&](vector<int>&C,vector<int>&cp,int mx,int st){ nodes++; if((int)C.size()>best){best=C.size();bs=C;} if(nodes>cap) return;
    for(int v=st;v<nm;v++){ int m2=mx; for(int x:C) m2=max(m2,(int)dd[x][v]); if(2*m2>D+(int)cp.size()) continue; vector<int> nc; set_intersection(cp.begin(),cp.end(),sets[v].begin(),sets[v].end(),back_inserter(nc)); if(2*m2>D+(int)nc.size()) continue; C.push_back(v); dfs(C,nc,m2,v+1); C.pop_back(); if(nodes>cap) return;} };
  for(int s0=0;s0<nm&&nodes<=cap;s0++){ vector<int> C={s0}; vector<int> cp=sets[s0]; dfs(C,cp,0,s0+1);} 
  // best pairwise excess
  int bex=1e9; for(int s=0;s<nm;s++) for(int t=s+1;t<nm;t++){ vector<int> nc; set_intersection(sets[s].begin(),sets[s].end(),sets[t].begin(),sets[t].end(),back_inserter(nc)); bex=min(bex,2*dd[s][t]-(D+(int)nc.size())); }
  printf("max cluster (one common divisor)=%d nodes=%ld; best pair excess=%d (<=0 critical)\n",best,nodes,bex);
}
