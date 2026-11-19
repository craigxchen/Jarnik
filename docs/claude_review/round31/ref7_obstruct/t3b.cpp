// T3: scan all translates tau (mod V) for coincidences rho_{H,a} = P_H(tau)+a*nu_H; verify contact 4n/9 for equal-rho pairs.
#include "gf.h"
typedef GFlog F; typedef F::E E; typedef Poly<F> P;
F f; int N;
E conjE(E a){ ll e=1; for(int j=0;j<N;j++) e*=3; return f.pw_(a,e); }
E evalP(const P&a,E x){ E r=f.zero(); for(int i=a.size()-1;i>=0;i--) r=f.add(f.mul(r,x),a[i]); return r; }
int main(int argc,char**argv){
  N=atoi(argv[1]); int d=atoi(argv[2]); vector<int> m; for(int a=3;a<argc;a++) m.push_back(atoi(argv[a]));
  f.init(3,m); ll qN=1; for(int j=0;j<N;j++) qN*=3;
  E g=f.gen(); E I=f.pw_(g,(qN+1)/2); E th=f.pw_(g,qN+1);
  vector<E> ab; { E t=I; for(int j=0;j<N;j++){ ab.push_back(t); t=f.mul(t,th);} }
  ll n=1; for(int j=0;j<d;j++) n*=3;
  vector<E> Vel(n); vector<vector<int>> co(n,vector<int>(d)); for(ll idx=0;idx<n;idx++){ ll s=idx; E v=f.zero(); for(int j=0;j<d;j++){co[idx][j]=s%3; s/=3; v=f.add(v,f.mul(f.fromint(co[idx][j]),ab[j]));} Vel[idx]=v; }
  struct HP{ vector<E> H; E mu; P PH; E nu; };
  vector<HP> hps;
  for(ll cidx=1;cidx<n;cidx++){ vector<int> c=co[cidx]; int fn=-1; for(int j=0;j<d;j++) if(c[j]){fn=j;break;} if(c[fn]!=1) continue;
    HP hp; ll mi=-1; for(ll idx=0;idx<n;idx++){ int s=0; for(int j=0;j<d;j++) s+=c[j]*co[idx][j]; s%=3; if(!s) hp.H.push_back(Vel[idx]); else if(s==1&&mi<0) mi=idx; }
    hp.mu=Vel[mi]; hp.PH={f.one()}; for(E h:hp.H) mullin(f,hp.PH,h); hp.nu=evalP(hp.PH,hp.mu); hps.push_back(hp); }
  int nh=hps.size();
  { // identity check: U_{H,a} = P_V + rho*Y^2 - rho^2*Y - rho^3 + rho*nu^2, rho=P_H(tau+a mu), Y=P_H(X)
    P PV={f.one()}; for(E v:Vel) mullin(f,PV,v); E tau=ab[d]; int bad=0,tot=0;
    for(auto&hp:hps) for(int a=0;a<3;a++){ P pl={f.one()}; E b1=f.add(tau,f.mul(f.fromint(a),hp.mu)); for(E h:hp.H) mullin(f,pl,f.add(b1,h)); for(int j=0;j<3;j++) if(j!=a){ E b2=f.neg(f.add(tau,f.mul(f.fromint(j),hp.mu))); for(E h:hp.H) mullin(f,pl,f.add(b2,h)); }
      E rho=evalP(hp.PH,b1); P Y=hp.PH; P Y2=pmul(f,Y,Y); P rhs=PV; auto addsc=[&](const P&A,E c){ if(rhs.size()<A.size()) rhs.resize(A.size(),f.zero()); for(size_t i=0;i<A.size();i++) rhs[i]=f.add(rhs[i],f.mul(c,A[i])); };
      addsc(Y2,rho); addsc(Y,f.neg(f.mul(rho,rho))); E c0=f.add(f.neg(f.pw_(rho,3)),f.mul(rho,f.mul(hp.nu,hp.nu))); addsc(P{c0},f.one()); trim(f,rhs); tot++; if(!psub(f,pl,rhs).empty()) bad++; }
    printf("identity check: %d/%d marks fail\n",bad,tot); } ll ntau=1; for(int j=d;j<N;j++) ntau*=3;
  int globalmax=0; ll eqpairs=0, eqcrit=0; int worst=-1;
  for(ll ti=1;ti<ntau;ti++){ ll s=ti; E tau=f.zero(); for(int j=d;j<N;j++){ tau=f.add(tau,f.mul(f.fromint(s%3),ab[j])); s/=3; }
    map<E,vector<pair<int,int>>> by; for(int hI=0;hI<nh;hI++){ E base=evalP(hps[hI].PH,tau); for(int a=0;a<3;a++) by[f.add(base,f.mul(f.fromint(a),hps[hI].nu))].push_back({hI,a}); }
    int mx=0; for(auto&kv:by) mx=max(mx,(int)kv.second.size()); globalmax=max(globalmax,mx);
    for(auto&kv:by) if(kv.second.size()>=2){ // verify contact for these marks
      vector<P> polys; for(auto pr:kv.second){ auto&hp=hps[pr.first]; int a=pr.second; P pl={f.one()}; E b1=f.add(tau,f.mul(f.fromint(a),hp.mu)); for(E h:hp.H) mullin(f,pl,f.add(b1,h)); for(int j=0;j<3;j++) if(j!=a){ E b2=f.neg(f.add(tau,f.mul(f.fromint(j),hp.mu))); for(E h:hp.H) mullin(f,pl,f.add(b2,h)); } polys.push_back(pl);} 
      for(size_t x=0;x<polys.size();x++) for(size_t y=x+1;y<polys.size();y++){ P df=psub(f,polys[x],polys[y]); int dg=deg(f,df); eqpairs++; worst=max(worst,dg); if(2*dg<= (int)n) eqcrit++; } }
  }
  printf("N=%d d=%d n=%lld hyperplanes=%d taus=%lld : max rho-multiplicity over all tau=%d ; equal-rho pairs=%lld critical-without-division=%lld max deg diff=%d (4n/9=%.1f, n/2=%.1f)\n",N,d,n,nh,ntau-1,globalmax,eqpairs,eqcrit,worst,4.0*n/9,n/2.0);
}
