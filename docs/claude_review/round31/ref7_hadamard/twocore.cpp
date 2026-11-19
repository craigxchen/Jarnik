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
  int np=pr.size(); int b=(e/2-1)/2; int bound=e/2-1;
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
    unordered_map<u64, vector<pair<int,int>>> bucket;
    for(int i=0;i<np;i++) for(int j=i+1;j<np;j++){ u64 P=clmul(pr[i],pr[j]); u64 ev=P&0x5555555555555555ULL, od=(P>>1)&0x5555555555555555ULL;
      u64 A=sqrt_even(ev), B=sqrt_even(od); // P = A^2 + t B^2
      if(dg(A)!=e){ printf("deg A error\n"); return 1; }
      u128 num=((u128)B)<<sh; // deg <= 3e - b - 2 < 64 for e<=16
      if(3*e-b-2>63){ printf("shift overflow\n"); return 1; }
      u64 key=pdiv((u64)num, A); bucket[key].push_back({i,j}); }
    set<array<int,4>> S;
    for(auto&kv:bucket){ auto&v=kv.second; for(size_t x=0;x<v.size();x++) for(size_t y=x+1;y<v.size();y++){
      int a=v[x].first,c=v[x].second,d=v[y].first,f=v[y].second; if(a==d||a==f||c==d||c==f) continue;
      array<int,4> q={a,c,d,f}; sort(q.begin(),q.end()); S.insert(q);} }
    elig.assign(S.begin(),S.end());
    // sanity (e<=15): verify a few by direct derivative of the product (degree 4e may exceed 64 at e=16: use 128)
    int bad=0; for(auto&q:elig){ u128 H=clmul128(clmul(pr[q[0]],pr[q[1]]), clmul(pr[q[2]],pr[q[3]]));
      // derivative degree of 128-bit H: max odd k with bit set, minus 1
      int dd=-1; for(int k=127;k>=1;k-=1){ if((k&1) && ((H>>k)&1)){ dd=k-1; break; } } if(dd>bound) bad++; }
    printf("key: eligible sets re-verified by direct derivative; failures=%d\n",bad);
  }
  printf("admissible four-prime sets=%zu\n",elig.size());
  // set system as bitsets over np primes (np <= 4096) -> use sorted index arrays and a hash set of 4-sets
  set<array<int,4>> E(elig.begin(),elig.end());
  auto symdiff=[&](const array<int,4>&x,const array<int,4>&y, array<int,4>&out)->bool{
    vector<int> r; set_symmetric_difference(x.begin(),x.end(),y.begin(),y.end(),back_inserter(r)); if(r.size()!=4) return false; for(int i=0;i<4;i++) out[i]=r[i]; return true; };
  // index by prime for lines: a,b must share exactly 2 primes
  vector<vector<int>> byPrime(np); for(size_t i=0;i<elig.size();i++) for(int p:elig[i]) byPrime[p].push_back(i);
  set<array<int,3>> lines;
  for(size_t i=0;i<elig.size();i++){ set<int> cand; for(int p:elig[i]) for(int j:byPrime[p]) if((size_t)j>i) cand.insert(j);
    for(int j:cand){ array<int,4> s; if(!symdiff(elig[i],elig[j],s)) continue; if(E.count(s)){ int k=lower_bound(elig.begin(),elig.end(),s)-elig.begin(); array<int,3> L={(int)i,j,k}; sort(L.begin(),L.end()); lines.insert(L);} } }
  printf("two-dimensional subspaces (lines)=%zu\n",lines.size());
  set<array<int,7>> planes;
  for(auto&L:lines){ for(size_t u=0;u<elig.size();u++){ if((int)u==L[0]||(int)u==L[1]||(int)u==L[2]) continue;
      array<int,4> s1,s2,s3; if(!symdiff(elig[u],elig[L[0]],s1)||!E.count(s1)) continue; if(!symdiff(elig[u],elig[L[1]],s2)||!E.count(s2)) continue; if(!symdiff(elig[u],elig[L[2]],s3)||!E.count(s3)) continue;
      array<int,7> pl={L[0],L[1],L[2],(int)u,(int)(lower_bound(elig.begin(),elig.end(),s1)-elig.begin()),(int)(lower_bound(elig.begin(),elig.end(),s2)-elig.begin()),(int)(lower_bound(elig.begin(),elig.end(),s3)-elig.begin())}; sort(pl.begin(),pl.end()); planes.insert(pl);} }
  printf("three-dimensional subspaces=%zu\n",planes.size());
  return 0;
}
