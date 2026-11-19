// Referee's own check of the VII.4 two-core computation (m = 8, first Walsh core, 7 primes of degree e).
// A four-set {p1..p4} of distinct irreducibles of degree e in F_2[t] is eligible iff deg (p1p2p3p4)' <= e/2 - 1.
// Mode "direct": test every 4-subset by forming the product and its derivative (e <= 12, product fits 64 bits).
// Mode "key": P = A^2 + t B^2 for pair products; eligible iff deg(AE+BC) <= b, b = floor((e/2-1)/2), which is
//             tested by bucketing on floor(t^{2e-b-1} B / A).  (Derived independently; see report.)
// Then count 2-dim subspaces (lines {a,b,a^b} of eligible sets) and 3-dim subspaces (Fano-complement systems).
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <array>
#include <wmmintrin.h>
#include <smmintrin.h>
using namespace std;
typedef uint64_t u64; typedef unsigned __int128 u128;
static inline u64 clmul(u64 x,u64 y){ __m128i a=_mm_set_epi64x(0,x), b=_mm_set_epi64x(0,y); __m128i r=_mm_clmulepi64_si128(a,b,0); return (u64)_mm_cvtsi128_si64(r); }
static inline u128 clmul128(u64 x,u64 y){ __m128i a=_mm_set_epi64x(0,x), b=_mm_set_epi64x(0,y); __m128i r=_mm_clmulepi64_si128(a,b,0); u64 lo=_mm_cvtsi128_si64(r); u64 hi=_mm_extract_epi64(r,1); return ((u128)hi<<64)|lo; }
static inline int dg(u64 x){ return x?63-__builtin_clzll(x):-1; }
static u64 pmod(u64 a,u64 f){ int df=dg(f); int d; while((d=dg(a))>=df) a^=f<<(d-df); return a; }
static u64 pdiv(u64 a,u64 f){ int df=dg(f); u64 q=0; int d; while((d=dg(a))>=df){ q|=1ULL<<(d-df); a^=f<<(d-df);} return q; }
static u64 mulmod(u64 a,u64 b,u64 f){ return pmod(clmul(a,b),f); } // deg f <= 16 ok
static u64 gcd2(u64 a,u64 b){ while(b){ u64 r=pmod(a,b); a=b; b=r;} return a; }
static bool irreducible(u64 f,int e){
  // t^(2^e) == t mod f, and gcd(t^(2^(e/q)) - t, f) = 1 for primes q | e
  u64 x=2; vector<u64> pw(e+1); pw[0]=x; for(int i=1;i<=e;i++){ x=mulmod(x,x,f); pw[i]=x; }
  if(pw[e]!=2) return false;
  for(int q=2;q<=e;q++){ if(e%q) continue; bool pr=true; for(int s=2;s*s<=q;s++) if(q%s==0) pr=false; if(!pr) continue;
    if(dg(gcd2(f, pw[e/q]^2))>0) return false; }
  return true;
}
static inline u64 deriv(u64 f){ return (f>>1)&0x5555555555555555ULL; } // coefficient of t^(k-1) for odd k
static inline int deriv_deg128(u128 f){ u64 lo=(u64)f, hi=(u64)(f>>64); u64 dlo=deriv(lo)|((hi&1)<<63); u64 dhi=deriv(hi>>0)>>0; // careful
  // derivative: bit k (k odd) -> bit k-1. For hi word bit j (k=64+j), odd iff j odd -> bit 63+j: j=1 -> lo bit 64?? handle simply
  (void)dlo;(void)dhi; return -2; }
static u64 sqrt_even(u64 x){ u64 r=0; for(int i=0;2*i<64;i++) if((x>>(2*i))&1) r|=1ULL<<i; return r; }

int main(int argc,char**argv){
  int e=atoi(argv[1]); const char*mode=argv[2];
  vector<u64> pr; for(u64 f=(1ULL<<e); f<(1ULL<<(e+1)); f++) if((f&1) && irreducible(f,e)) pr.push_back(f);
  int np=pr.size(); int b=(e/2-1)/2; int bound=e/2-1; if(argc>3) bound=atoi(argv[3]);
  printf("e=%d available primes=%d  derivative bound deg H' <= %d  (b=%d)\n",e,np,bound,b);
  vector<array<int,4>> elig;
  if(!strcmp(mode,"direct")){
    if(4*e>63){ printf("too big\n"); return 1; }
    vector<u64> pp((size_t)np*np); for(int i=0;i<np;i++) for(int j=0;j<np;j++) pp[(size_t)i*np+j]=clmul(pr[i],pr[j]);
    long long cnt=0;
    for(int i=0;i<np;i++) for(int j=i+1;j<np;j++){ u64 P=pp[(size_t)i*np+j];
      for(int k=j+1;k<np;k++){ u64 Pk=clmul(P,pr[k]); for(int l=k+1;l<np;l++){ u64 H=clmul(Pk,pr[l]); cnt++; if(dg(deriv(H))<=bound) elig.push_back({i,j,k,l}); } } }
    printf("direct: tested %lld four-sets\n",cnt);
  } else {
    // key bucketing
    int sh=2*e-b-1;
    vector<u64> packed; packed.reserve((size_t)np*(np-1)/2);
    for(int i=0;i<np;i++) for(int j=i+1;j<np;j++){ u64 P=clmul(pr[i],pr[j]); u64 ev=P&0x5555555555555555ULL, od=(P>>1)&0x5555555555555555ULL;
      u64 A=sqrt_even(ev), B=sqrt_even(od); // P = A^2 + t B^2
      if(dg(A)!=e){ printf("deg A error\n"); return 1; }
      u128 num=((u128)B)<<sh; // deg <= 3e - b - 2 < 64 for e<=16
      if(3*e-b-2>63){ printf("shift overflow\n"); return 1; }
      u64 key=pdiv((u64)num, A); if(dg(key)>=36-2*0 && 0){}
      if(key>>(64-2*15)){ printf("key too wide\n"); return 1; }
      packed.push_back((key<<30)|((u64)i<<15)|(u64)j); }
    sort(packed.begin(),packed.end());
    set<array<int,4>> S;
    for(size_t s=0;s<packed.size();){ size_t t=s; while(t<packed.size() && (packed[t]>>30)==(packed[s]>>30)) t++;
      for(size_t x=s;x<t;x++) for(size_t y=x+1;y<t;y++){
      int a=(packed[x]>>15)&0x7fff,c=packed[x]&0x7fff,d=(packed[y]>>15)&0x7fff,f=packed[y]&0x7fff; if(a==d||a==f||c==d||c==f) continue;
      array<int,4> q={a,c,d,f}; sort(q.begin(),q.end()); S.insert(q);} s=t; }
    { vector<u64>().swap(packed); }
    elig.assign(S.begin(),S.end());
    // sanity (e<=15): verify a few by direct derivative of the product (degree 4e may exceed 64 at e=16: use 128)
    int bad=0; for(auto&q:elig){ u128 H=clmul128(clmul(pr[q[0]],pr[q[1]]), clmul(pr[q[2]],pr[q[3]]));
      // derivative degree of 128-bit H: max odd k with bit set, minus 1
      int dd=-1; for(int k=127;k>=1;k-=1){ if((k&1) && ((H>>k)&1)){ dd=k-1; break; } } if(dd>bound) bad++; }
    printf("key: eligible sets re-verified by direct derivative; failures=%d\n",bad);
  }
  printf("admissible four-prime sets=%zu\n",elig.size());
  fflush(stdout);
  // fast search. pack a sorted 4-set as 4x15 bits
  auto pk=[&](const array<int,4>&q)->u64{ return ((u64)q[0]<<45)|((u64)q[1]<<30)|((u64)q[2]<<15)|(u64)q[3]; };
  unordered_map<u64,int> id; id.reserve(elig.size()*2); for(size_t i=0;i<elig.size();i++) id[pk(elig[i])]=i;
  auto sd=[&](const array<int,4>&x,const array<int,4>&y, array<int,4>&out)->bool{ int r[8]; int k=0,i=0,j=0;
    while(i<4||j<4){ if(j==4||(i<4&&x[i]<y[j])) r[k++]=x[i++]; else if(i==4||y[j]<x[i]) r[k++]=y[j++]; else {i++;j++;} if(k>4) return false; }
    if(k!=4) return false; for(int t=0;t<4;t++) out[t]=r[t]; return true; };
  unordered_map<u64, vector<int>> byPair; // pair of primes -> sets containing it
  for(size_t i=0;i<elig.size();i++){ auto&q=elig[i]; for(int x=0;x<4;x++) for(int y=x+1;y<4;y++) byPair[((u64)q[x]<<15)|q[y]].push_back(i); }
  set<array<int,3>> lines;
  for(auto&kv:byPair){ auto&v=kv.second; for(size_t x=0;x<v.size();x++) for(size_t y=x+1;y<v.size();y++){ array<int,4> s; if(!sd(elig[v[x]],elig[v[y]],s)) continue;
      auto it=id.find(pk(s)); if(it==id.end()) continue; array<int,3> L={v[x],v[y],it->second}; sort(L.begin(),L.end()); lines.insert(L);} }
  printf("two-dimensional subspaces (lines)=%zu\n",lines.size()); fflush(stdout);
  vector<vector<pair<int,int>>> linesAt(elig.size()); // element a -> (u, a^u)
  for(auto&L:lines){ linesAt[L[0]].push_back({L[1],L[2]}); linesAt[L[1]].push_back({L[0],L[2]}); linesAt[L[2]].push_back({L[0],L[1]}); }
  set<array<int,7>> planes;
  for(auto&L:lines){ int a=L[0],b=L[1],c=L[2];
    for(auto&pr2:linesAt[a]){ int u=pr2.first; if(u==b||u==c) continue; array<int,4> s1,s2;
      if(!sd(elig[u],elig[b],s1)) continue; auto i1=id.find(pk(s1)); if(i1==id.end()) continue;
      if(!sd(elig[u],elig[c],s2)) continue; auto i2=id.find(pk(s2)); if(i2==id.end()) continue;
      array<int,7> pl={a,b,c,u,pr2.second,i1->second,i2->second}; sort(pl.begin(),pl.end()); planes.insert(pl);} }
  printf("three-dimensional subspaces=%zu\n",planes.size());
  return 0;
}
