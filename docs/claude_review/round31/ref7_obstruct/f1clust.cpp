// F1 cluster search: marks = (form l, S subset F_q) selecting fibres {v+tau: l(v) in S} u {v-tau: l(v) in -(F_q\S)}.
// Cluster validity uses ONE common divisor: max deg(U_i-U_j) <= (D + |cap S_j|)/2.
#include "gf.h"
typedef GFlog F; typedef F::E E; typedef Poly<F> P;
int main(int argc,char**argv){ int q=atoi(argv[1]), r=atoi(argv[2]); long cap=atol(argv[3]); vector<int> m; for(int a=4;a<argc;a++) m.push_back(atoi(argv[a]));
  F f; f.init(q,m); ll qr=1; for(int j=0;j<r;j++) qr*=q; E g=f.gen(); E I=f.pw_(g,(qr+1)/2), th=f.pw_(g,qr+1);
  int k=2; ll n=q*q; vector<E> V(n); vector<array<int,2>> X(n); E b0=I, b1=f.mul(I,th), tau=f.mul(b1,th);
  for(ll idx=0;idx<n;idx++){ X[idx]={(int)(idx%q),(int)(idx/q)}; V[idx]=f.add(f.mul(f.fromint(X[idx][0]),b0),f.mul(f.fromint(X[idx][1]),b1)); }
  vector<array<int,2>> forms; forms.push_back({0,1}); for(int c=0;c<q;c++) forms.push_back({1,c});
  struct Mk{ vector<int> ids; P poly; int fm; int S; }; vector<Mk> mk;
  for(int fi=0;fi<(int)forms.size();fi++) for(int S=1;S<(1<<q)-1;S++){ Mk M; M.fm=fi; M.S=S; vector<E> rts;
     for(ll idx=0;idx<n;idx++){ int lam=(forms[fi][0]*X[idx][0]+forms[fi][1]*X[idx][1])%q; if((S>>lam)&1){ rts.push_back(f.add(V[idx],tau)); M.ids.push_back(idx);} int ml=(q-lam)%q; if(!((S>>ml)&1)){ rts.push_back(f.sub(V[idx],tau)); M.ids.push_back(n+idx);} }
     sort(M.ids.begin(),M.ids.end()); M.poly={f.one()}; for(E x:rts) mullin(f,M.poly,x); mk.push_back(M);} 
  // dedupe identical root sets
  sort(mk.begin(),mk.end(),[](const Mk&a,const Mk&b){return a.ids<b.ids;}); vector<Mk> u; for(auto&M:mk) if(u.empty()||u.back().ids!=M.ids) u.push_back(M); mk=u;
  int nm=mk.size(); int D=n; printf("q=%d marks(distinct)=%d D=%d\n",q,nm,D);
  vector<vector<short>> dd(nm,vector<short>(nm)); for(int s=0;s<nm;s++) for(int t=s+1;t<nm;t++){ dd[s][t]=dd[t][s]=deg(f,psub(f,mk[s].poly,mk[t].poly)); }
  int best=1; vector<int> bs; long nodes=0;
  function<void(vector<int>&,vector<int>&,int,int)> dfs=[&](vector<int>&C,vector<int>&cp,int mx,int st){ nodes++; if((int)C.size()>best){best=C.size(); bs=C;} if(nodes>cap) return;
    for(int v=st;v<nm;v++){ int m2=mx; for(int x:C) m2=max(m2,(int)dd[x][v]); if(2*m2>D+(int)cp.size()) continue; vector<int> nc; set_intersection(cp.begin(),cp.end(),mk[v].ids.begin(),mk[v].ids.end(),back_inserter(nc)); if(2*m2>D+(int)nc.size()) continue; C.push_back(v); dfs(C,nc,m2,v+1); C.pop_back(); if(nodes>cap) return; } };
  for(int s0=0;s0<nm && nodes<=cap;s0++){ vector<int> C={s0}; vector<int> cp=mk[s0].ids; dfs(C,cp,0,s0+1); }
  printf("max cluster with one common divisor=%d (nodes %ld%s)\n",best,nodes,nodes>cap?", capped":"");
  for(int v:bs) printf("  form(%d,%d) S=%x\n",forms[mk[v].fm][0],forms[mk[v].fm][1],mk[v].S);
  // cap of the best cluster and reduced degree
  if(bs.size()>1){ vector<int> cp=mk[bs[0]].ids; for(int v:bs){ vector<int> nc; set_intersection(cp.begin(),cp.end(),mk[v].ids.begin(),mk[v].ids.end(),back_inserter(nc)); cp=nc;} int mx=0; for(int a:bs) for(int b:bs) if(a!=b) mx=max(mx,(int)dd[a][b]); printf("  common divisor degree=%zu reduced degree=%d max reduced diff degree=%zu\n",cp.size(),D-(int)cp.size(),mx-cp.size()); }
}
