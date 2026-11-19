// Own finite-field arithmetic. Two implementations with the same interface.
// Elements are ints. GFlog: log representation (ZERO=-1), Zech addition, Q <= ~5e6.
// GFdir: digit representation (index base p), direct multiplication mod a primitive poly.
#pragma once
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
struct GFdir {
  int p,N; ll Q; vector<int> mod; // monic, low->high, size N+1
  vector<ll> pw;
  void init(int p_, vector<int> m){p=p_;mod=m;N=m.size()-1;Q=1;pw.assign(N+1,1);for(int i=0;i<N;i++){pw[i+1]=pw[i]*p;}Q=pw[N];}
  typedef ll E;
  E zero() const {return 0;} E one() const {return 1;}
  bool isz(E a) const {return a==0;}
  void dig(E a,int*d) const {for(int i=0;i<N;i++){d[i]=a%p;a/=p;}}
  E und(const int*d) const {E r=0;for(int i=N-1;i>=0;i--) r=r*p+d[i];return r;}
  E add(E a,E b) const {int x[24],y[24];dig(a,x);dig(b,y);for(int i=0;i<N;i++){x[i]+=y[i];if(x[i]>=p)x[i]-=p;}return und(x);}
  E neg(E a) const {int x[24];dig(a,x);for(int i=0;i<N;i++) x[i]=(p-x[i])%p;return und(x);}
  E sub(E a,E b) const {return add(a,neg(b));}
  E mul(E a,E b) const {if(!a||!b) return 0;int x[24],y[24];long r[48];dig(a,x);dig(b,y);for(int i=0;i<2*N;i++)r[i]=0;
    for(int i=0;i<N;i++) if(x[i]) for(int j=0;j<N;j++) r[i+j]+=x[i]*y[j];
    for(int i=0;i<2*N;i++) r[i]%=p;
    for(int d=2*N-2;d>=N;d--){long c=r[d]%p; if(c){for(int j=0;j<=N;j++){r[d-N+j]=((r[d-N+j]-c*mod[j])%p+p)%p;}} }
    int z[24];for(int i=0;i<N;i++) z[i]=r[i]%p;return und(z);}
  E pw_(E a,ll e) const {E r=1;while(e){if(e&1)r=mul(r,a);a=mul(a,a);e>>=1;}return r;}
  E inv(E a) const {return pw_(a,Q-2);}
  E fromint(ll c) const {c%=p;if(c<0)c+=p;return c;}
  E gen() const {return p;} // the class of x (index p means digit1=1) when N>1
};
struct GFlog {
  int p,N; ll Q; vector<int> ex; vector<int> lg; vector<int> zech; int negone; // log of -1
  typedef int E;
  void init(int p_, vector<int> m){p=p_;N=m.size()-1;GFdir d;d.init(p,m);Q=d.Q;ex.assign(Q-1,0);lg.assign(Q,-1);
    ll cur=1;ll g=(N>1)?p:0; if(N==1){ // find primitive root mod p
      for(ll c=2;c<p;c++){ll x=1;int ord=0;do{x=x*c%p;ord++;}while(x!=1);if(ord==p-1){g=c;break;}} }
    for(ll k=0;k<Q-1;k++){ex[k]=cur;if(lg[cur]!=-1){fprintf(stderr,"not primitive\n");exit(1);}lg[cur]=k;cur=d.mul(cur,g);}
    zech.assign(Q-1,-1);
    for(ll k=0;k<Q-1;k++){ll v=d.add(1,ex[k]);zech[k]=lg[v];}
    negone=lg[d.neg(1)];}
  E zero() const {return -1;} E one() const {return 0;}
  bool isz(E a) const {return a<0;}
  E mul(E a,E b) const {if(a<0||b<0) return -1; ll s=(ll)a+b; if(s>=Q-1) s-=Q-1; return s;}
  E add(E a,E b) const {if(a<0) return b; if(b<0) return a; ll k=(ll)b-a; if(k<0)k+=Q-1; int z=zech[k]; if(z<0) return -1; ll s=(ll)a+z; if(s>=Q-1)s-=Q-1; return s;}
  E neg(E a) const {return mul(a,negone);}
  E sub(E a,E b) const {return add(a,neg(b));}
  E pw_(E a,ll e) const {if(a<0) return e==0?0:-1; return (E)(( (__int128)a*e)%(Q-1));}
  E inv(E a) const {return a==0?0:(E)(Q-1-a);}
  E fromdig(ll idx) const {return lg[idx];}
  ll todig(E a) const {return a<0?0:ex[a];}
  E fromint(ll c) const {c%=p;if(c<0)c+=p;return lg[c];}
  E gen() const {return 1;}
};
// polynomial helpers (vector low->high), template on field
template<class F> using Poly = vector<typename F::E>;
template<class F> void trim(const F&f, Poly<F>&a){while(!a.empty()&&f.isz(a.back())) a.pop_back();}
template<class F> int deg(const F&f, const Poly<F>&a){int d=a.size()-1; while(d>=0&&f.isz(a[d])) d--; return d;}
template<class F> Poly<F> pmul(const F&f,const Poly<F>&a,const Poly<F>&b){ if(a.empty()||b.empty()) return {}; Poly<F> r(a.size()+b.size()-1,f.zero());
  for(size_t i=0;i<a.size();i++) if(!f.isz(a[i])) for(size_t j=0;j<b.size();j++) if(!f.isz(b[j])) r[i+j]=f.add(r[i+j],f.mul(a[i],b[j])); trim(f,r); return r;}
template<class F> Poly<F> psub(const F&f,const Poly<F>&a,const Poly<F>&b){Poly<F> r(max(a.size(),b.size()),f.zero());for(size_t i=0;i<a.size();i++)r[i]=a[i];for(size_t i=0;i<b.size();i++)r[i]=f.sub(r[i],b[i]);trim(f,r);return r;}
// multiply in place by (X - c)
template<class F> void mullin(const F&f, Poly<F>&a, typename F::E c){ a.push_back(f.zero()); for(int i=a.size()-1;i>=1;i--) a[i]=f.sub(a[i-1],f.mul(c,a[i])); a[0]=f.neg(f.mul(c,a[0])); }
template<class F> void pdivmod(const F&f, Poly<F> a, const Poly<F>&b, Poly<F>&q, Poly<F>&r){ int db=deg(f,b); auto ib=f.inv(b[db]); int da=deg(f,a); q.assign(max(0,da-db+1),f.zero());
  for(int d=da; d>=db; d--){ if(f.isz(a[d])) continue; auto c=f.mul(a[d],ib); q[d-db]=c; for(int j=0;j<=db;j++) a[d-db+j]=f.sub(a[d-db+j],f.mul(c,b[j])); }
  a.resize(max(0,db),f.zero()); trim(f,a); r=a; trim(f,q);}
template<class F> Poly<F> pgcd(const F&f, Poly<F> a, Poly<F> b){ trim(f,a); trim(f,b); while(!b.empty()){ Poly<F> q,r; pdivmod(f,a,b,q,r); a=b; b=r;} if(!a.empty()){auto c=f.inv(a.back()); for(auto&x:a) x=f.mul(x,c);} return a;}
template<class F> Poly<F> pderiv(const F&f,const Poly<F>&a){ Poly<F> r; for(size_t i=1;i<a.size();i++){ auto c=f.fromint(i); r.push_back(f.mul(c,a[i])); } trim(f,r); return r;}
