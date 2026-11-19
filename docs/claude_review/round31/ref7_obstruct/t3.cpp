// Mechanism T3 (char 3): V an F_3-subspace of i*F_{3^N} (anti-invariant), tau outside V.
// Roots (V+tau) u (V-tau), squarefree. For each F_3-hyperplane H of V and a in F_3 the mark selects the coset
// H+tau+a*mu from V+tau and the two cosets H-tau-j*mu (j!=a) from V-tau: one root of each conjugate pair.
// Claim: if rho_{H,a}=P_H(tau+a*mu) agree for two marks then deg(U-U') <= 4n/9 < n/2 = D/2 (critical).
#include "gf.h"
typedef GFlog F; typedef F::E E; typedef Poly<F> P;
F f; int N; ll Nq;
E conj(E a){ ll e=1; for(int j=0;j<N;j++) e*=3; return f.pw_(a,e); }
int main(int argc,char**argv){
  N=atoi(argv[1]); int d=atoi(argv[2]); int seed=atoi(argv[3]); int randomV=atoi(argv[4]); vector<int> m; for(int a=5;a<argc;a++) m.push_back(atoi(argv[a]));
  f.init(3,m); ll qN=1; for(int j=0;j<N;j++) qN*=3; if(f.Q!=qN*qN){fprintf(stderr,"field\n");return 1;}
  E g=f.gen(); E I=f.pw_(g,(qN+1)/2); E th=f.pw_(g,qN+1); if(conj(I)!=f.neg(I)){fprintf(stderr,"I\n");return 1;}
  mt19937 rng(seed);
  // basis of anti-invariant space: I*th^j, j<N. Choose V basis b_1..b_d and tau
  vector<E> ab; { E t=I; for(int j=0;j<N;j++){ ab.push_back(t); t=f.mul(t,th);} }
  auto randanti=[&](){ E v=f.zero(); for(int j=0;j<N;j++) v=f.add(v,f.mul(f.fromint(rng()%3),ab[j])); return v; };
  vector<E> B; E tau;
  if(!randomV){ for(int j=0;j<d;j++) B.push_back(ab[j]); tau=ab[d]; }
  else { // random independent
    while(true){ B.clear(); for(int j=0;j<d;j++) B.push_back(randanti()); tau=randanti();
      // check independence of B+tau via span size
      set<E> span; span.insert(f.zero()); vector<E> all=B; all.push_back(tau); for(E b:all){ set<E> ns; for(E s:span) for(int c=0;c<3;c++) ns.insert(f.add(s,f.mul(f.fromint(c),b))); span=ns;} if((ll)span.size()==(ll)pow(3,d+1)) break; } }
  ll n=1; for(int j=0;j<d;j++) n*=3;
  vector<E> Vel(n); vector<vector<int>> co(n,vector<int>(d)); for(ll idx=0;idx<n;idx++){ ll s=idx; E v=f.zero(); for(int j=0;j<d;j++){co[idx][j]=s%3; s/=3; v=f.add(v,f.mul(f.fromint(co[idx][j]),B[j]));} Vel[idx]=v; }
  // hyperplanes: functionals c (first nonzero =1)
  struct Mark{ vector<E> roots; P poly; E rho; int hid,a; };
  vector<Mark> marks; int hid=0;
  for(ll cidx=1;cidx<n;cidx++){ vector<int> c=co[cidx]; int fn=-1; for(int j=0;j<d;j++) if(c[j]){fn=j;break;} if(c[fn]!=1) continue;
    auto fval=[&](ll idx){ int s=0; for(int j=0;j<d;j++) s+=c[j]*co[idx][j]; return s%3; };
    vector<E> H; ll muidx=-1; for(ll idx=0;idx<n;idx++){ int v=fval(idx); if(v==0) H.push_back(Vel[idx]); else if(v==1 && muidx<0) muidx=idx; }
    E mu=Vel[muidx];
    for(int a=0;a<3;a++){ Mark M; M.hid=hid; M.a=a; E base=f.add(tau,f.mul(f.fromint(a),mu));
      for(E h:H) M.roots.push_back(f.add(base,h));
      for(int j=0;j<3;j++) if(j!=a){ E b2=f.neg(f.add(tau,f.mul(f.fromint(j),mu))); for(E h:H) M.roots.push_back(f.add(b2,h)); }
      E rho=f.one(); for(E h:H) rho=f.mul(rho,f.sub(base,h)); M.rho=rho;
      P pl={f.one()}; for(E x:M.roots) mullin(f,pl,x); M.poly=pl; marks.push_back(M);} hid++; }
  int nm=marks.size(); int D=n;
  // checks: norms equal & squarefree & in F_3^N[X]
  P N0=pmul(f,marks[0].poly,[&]{P c=marks[0].poly; for(auto&x:c) x=conj(x); return c;}());
  int normbad=0; for(int t=1;t<nm;t++){ P c=marks[t].poly; for(auto&x:c) x=conj(x); P Nt=pmul(f,marks[t].poly,c); if(!psub(f,Nt,N0).empty()) normbad++; }
  P h=pgcd(f,N0,pderiv(f,N0));
  printf("N=%d d=%d n=%lld marks=%d hyperplanes=%d normsdiffer=%d sqfree=%d degN=%d\n",N,d,n,nm,hid,normbad,(int)(deg(f,h)==0),deg(f,N0));
  // pairwise
  map<E,vector<int>> byrho; for(int t=0;t<nm;t++) byrho[marks[t].rho].push_back(t);
  int maxmult=0; for(auto&kv:byrho) maxmult=max(maxmult,(int)kv.second.size());
  printf("distinct rho values=%zu  max multiplicity=%d\n",byrho.size(),maxmult);
  ll crit=0, critEq=0, eqpairs=0, viol=0; int maxdeg_eq=0;
  vector<vector<char>> adj(nm,vector<char>(nm,0));
  for(int s=0;s<nm;s++) for(int t=s+1;t<nm;t++){ P df=psub(f,marks[s].poly,marks[t].poly); if(df.empty()){ printf("identical marks %d %d\n",s,t); continue;}
    int dd=deg(f,df); vector<E> A=marks[s].roots,Bv=marks[t].roots; sort(A.begin(),A.end()); sort(Bv.begin(),Bv.end()); vector<E> in; set_intersection(A.begin(),A.end(),Bv.begin(),Bv.end(),back_inserter(in)); int g=in.size();
    bool c= 2*dd <= D+g; if(c){crit++; adj[s][t]=adj[t][s]=1;}
    if(marks[s].rho==marks[t].rho){ eqpairs++; maxdeg_eq=max(maxdeg_eq,dd); if(c) critEq++; if(9*dd>4*n) viol++; } }
  printf("critical pairs=%lld ; equal-rho pairs=%lld (critical %lld, deg>4n/9: %lld, max deg %d, 4n/9=%.1f)\n",crit,eqpairs,critEq,viol,maxdeg_eq,4.0*n/9);
  // clusters with ONE common divisor: set C valid iff max_{i,j in C} deg(U_i-U_j) <= (D + |cap_{j in C} S_j|)/2 (downward closed)
  vector<vector<int>> dd(nm,vector<int>(nm,-1)); for(int s=0;s<nm;s++) for(int t=s+1;t<nm;t++){ P df=psub(f,marks[s].poly,marks[t].poly); dd[s][t]=dd[t][s]=deg(f,df);} 
  vector<vector<E>> S(nm); for(int t=0;t<nm;t++){ S[t]=marks[t].roots; sort(S[t].begin(),S[t].end()); }
  int best=1; vector<int> bestset={0}; long nodes=0;
  function<void(vector<int>&,vector<E>&,int,int)> dfs=[&](vector<int>&C,vector<E>&cap,int mx,int start){ nodes++; if((int)C.size()>best){best=C.size();bestset=C;} if(nodes>20000000) return;
    for(int v=start;v<nm;v++){ int m2=mx; bool ok=true; for(int u:C){ m2=max(m2,dd[u][v]); } vector<E> nc; set_intersection(cap.begin(),cap.end(),S[v].begin(),S[v].end(),back_inserter(nc)); if(2*m2> D+(int)nc.size()) continue; C.push_back(v); dfs(C,nc,m2,v+1); C.pop_back(); } };
  for(int s0=0;s0<nm;s0++){ vector<int> C={s0}; vector<E> cap=S[s0]; dfs(C,cap,0,s0+1); }
  printf("max cluster (single common divisor)=%d nodes=%ld :",best,nodes); for(int v:bestset) printf(" (H%d,a%d)",marks[v].hid,marks[v].a); printf("\n");
}
