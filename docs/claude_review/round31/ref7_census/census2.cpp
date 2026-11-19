// Exhaustive census of max clusters of monic z in F_9[t], deg z = D, sharing coefficients at degrees > floor(D/2),
// with equal norm z*conj(z) in F_3[t]. Reports max over all parents and over squarefree parents.
// Method: enumerate heads (coeffs at degrees floor(D/2)+1..D-1) up to the 12-element group
// (t->s t+a, conj), all low parts, group by exact norm via radix sort.
#include <bits/stdc++.h>
#include <thread>
#include <mutex>
using namespace std;
typedef uint64_t u64;
int MUL[9][9], ADD[9][9], CONJ[9], TR[9]; // TR(x)=x+conj(x) in F_3 (=2a)
int D, KL, NL, NH, NK; // NK = number of key trits = D+KL+1
const u64 ONES = 0x1249249249249249ULL; // bit 0 of each 3-bit field (21 fields)
static inline u64 tadd(u64 x, u64 y){ u64 s=x+y; u64 m=((s>>2)|((s>>1)&s))&ONES; return s-3*m; }
u64 enc(const vector<int>& t){ u64 r=0; for(size_t k=0;k<t.size();k++) r|=(u64)t[k]<<(3*k); return r; }
// poly over F_9 multiply
vector<int> pmul(const vector<int>&a,const vector<int>&b){ vector<int> r(a.size()+b.size()-1,0);
  for(size_t i=0;i<a.size();i++) if(a[i]) for(size_t j=0;j<b.size();j++) r[i+j]=ADD[r[i+j]][MUL[a[i]][b[j]]]; return r;}
// F_3 poly utilities (vector<int> low first, values 0..2)
static void trimv(vector<int>&p){ while(!p.empty()&&p.back()==0) p.pop_back(); }
static int inv3(int x){ return x; } // 1->1,2->2
vector<int> pmod3(vector<int> a, vector<int> b){ trimv(a); trimv(b); int il=inv3(b.back());
  while(a.size()>=b.size()){ int c=(a.back()*il)%3; size_t s=a.size()-b.size();
    for(size_t j=0;j<b.size();j++) a[s+j]=((a[s+j]-c*b[j])%3+3)%3; trimv(a);} return a;}
bool squarefree3(vector<int> N){ trimv(N); vector<int> d; for(size_t k=1;k<N.size();k++) d.push_back((int)(k%3)*N[k]%3); trimv(d);
  if(d.empty()) return false; vector<int> a=N,b=d; while(!b.empty()){ vector<int> r=pmod3(a,b); a=b; b=r;} return a.size()==1; }
// group transforms on full polys
vector<int> transform(const vector<int>& z, int a, int s, int cj){ // z(s t + a)*s^D, then conj
  int n=z.size(); vector<int> r(n,0); // Horner: r = sum z_k (s t + a)^k
  vector<int> lin={a, s}; // elements of F_3 subset of F_9: a in {0,1,2}, s in {1,2}
  vector<int> acc={0};
  for(int k=n-1;k>=0;k--){ acc=pmul(acc,lin); acc[0]=ADD[acc[0]][z[k]]; }
  acc.resize(n,0); int sd = (D%2==0)?1:s; for(int k=0;k<n;k++){ int v=MUL[acc[k]][sd]; r[k]=cj?CONJ[v]:v; } return r;}
struct Res { int maxAll=0, maxSF=0; map<int,u64> cntAll, cntSF; map<u64,int> parentsAll, parentsSF; vector<string> wit; };
mutex mtx; Res R; vector<u64> S; // S[L] = key of L*conj(L)
u64 codeN(const vector<int>& N){ u64 c=0; for(int k=(int)N.size()-1;k>=0;k--) c=c*3+N[k]; return c; }
void worker(vector<vector<int>>* heads, atomic<size_t>* nxt){
  size_t NLC=1; for(int k=0;k<NL;k++) NLC*=9;
  vector<u64> keys(NLC), tmp(NLC); vector<uint32_t> idx(NLC), idx2(NLC);
  int lo = min(NL,3); size_t NLO=1; for(int k=0;k<lo;k++) NLO*=9; size_t NHI=NLC/NLO;
  vector<u64> plo(NLO), phi(NHI);
  while(true){ size_t h=(*nxt)++; if(h>=heads->size()) break; const vector<int>& H=(*heads)[h];
    // T[k][v] = Tr(H * conj(v) t^k), low NK trits
    vector<array<u64,9>> T(NL);
    for(int k=0;k<NL;k++) for(int v=0;v<9;v++){ vector<int> t(NK,0);
      for(int j=KL+1;j<=D;j++){ int deg=j+k; if(deg<NK) t[deg]=TR[MUL[H[j]][CONJ[v]]]; } T[k][v]=enc(t);}
    for(size_t i=0;i<NLO;i++){ u64 s=0; size_t x=i; for(int k=0;k<lo;k++){ s=tadd(s,T[k][x%9]); x/=9;} plo[i]=s;}
    for(size_t i=0;i<NHI;i++){ u64 s=0; size_t x=i; for(int k=lo;k<NL;k++){ s=tadd(s,T[k][x%9]); x/=9;} phi[i]=s;}
    for(size_t L=0;L<NLC;L++){ keys[L]=tadd(tadd(plo[L%NLO],phi[L/NLO]),S[L]); idx[L]=L; }
    // LSD radix sort 4 passes of 16 bits (keys < 2^63)
    for(int pass=0;pass<4;pass++){ vector<uint32_t> cnt(65537,0); int sh=16*pass;
      for(size_t i=0;i<NLC;i++) cnt[((keys[i]>>sh)&0xFFFF)+1]++;
      for(int b=0;b<65536;b++) cnt[b+1]+=cnt[b];
      for(size_t i=0;i<NLC;i++){ size_t p=cnt[(keys[i]>>sh)&0xFFFF]++; tmp[p]=keys[i]; idx2[p]=idx[i]; }
      swap(keys,tmp); swap(idx,idx2); }
    // H*conj(H) full
    vector<int> Hc(H.size()); for(size_t j=0;j<H.size();j++) Hc[j]=CONJ[H[j]];
    vector<int> HH=pmul(H,Hc);
    Res loc; vector<pair<vector<int>,vector<uint32_t>>> big;
    for(size_t i=0;i<NLC;){ size_t j=i; while(j<NLC&&keys[j]==keys[i]) j++; int g=j-i;
      if(g>=(D<=5?2:3)){ vector<int> N(2*D+1); for(int k=0;k<=2*D;k++){ int v=HH[k]%3; // HH in F_3 (b=0)
          if(k<NK){ v=(v+(int)((keys[i]>>(3*k))&7))%3; } N[k]=v; }
        bool sf=squarefree3(N); u64 c=codeN(N);
        loc.cntAll[g]++; if(sf) loc.cntSF[g]++;
        if(g>=4||D<=5){ int&a=loc.parentsAll[c]; a=max(a,g); if(sf){int&b=loc.parentsSF[c]; b=max(b,g);} }
        if(g>loc.maxAll) loc.maxAll=g; if(sf&&g>loc.maxSF) loc.maxSF=g;
        if(g>=5){ vector<uint32_t> m(idx.begin()+i, idx.begin()+j); big.push_back({N,m}); } }
      i=j; }
    lock_guard<mutex> lk(mtx);
    R.maxAll=max(R.maxAll,loc.maxAll); R.maxSF=max(R.maxSF,loc.maxSF);
    for(auto&p:loc.cntAll) R.cntAll[p.first]+=p.second; for(auto&p:loc.cntSF) R.cntSF[p.first]+=p.second;
    // keep parents only for g>=4 groups; we record per-size later via cnt; store sets for counting
    for(auto&p:loc.parentsAll){int&a=R.parentsAll[p.first]; a=max(a,p.second);} for(auto&p:loc.parentsSF){int&b=R.parentsSF[p.first]; b=max(b,p.second);}
    for(auto& b:big){ string s="g="+to_string(b.second.size())+" sf="+to_string(squarefree3(b.first))+" N(high->low)=";
      for(int k=2*D;k>=0;k--) s+=char('0'+b.first[k]); s+=" head(low->high)="; for(int k=KL+1;k<=D;k++) s+=char('0'+H[k]);
      s+=" lows(base9 idx)="; for(auto x:b.second) s+=to_string(x)+","; R.wit.push_back(s);} }
}
int main(int argc,char**argv){ D=atoi(argv[1]); int nthr=argc>2?atoi(argv[2]):2;
  for(int x=0;x<9;x++) for(int y=0;y<9;y++){ int a=x%3,b=x/3,c=y%3,d=y/3; MUL[x][y]=((a*c-b*d)%3+3)%3+3*((a*d+b*c)%3); ADD[x][y]=(a+c)%3+3*((b+d)%3);}
  for(int x=0;x<9;x++){ int a=x%3,b=x/3; CONJ[x]=a+3*((3-b)%3); TR[x]=(2*a)%3; }
  KL=D/2; NL=KL+1; NH=D-1-KL; NK=D+KL+1; if(NK>21){fprintf(stderr,"too many trits\n");return 1;}
  size_t NLC=1; for(int k=0;k<NL;k++) NLC*=9; S.resize(NLC);
  for(size_t L=0;L<NLC;L++){ vector<int> l(NL); size_t x=L; for(int k=0;k<NL;k++){ l[k]=x%9; x/=9;}
    vector<int> lc(NL); for(int k=0;k<NL;k++) lc[k]=CONJ[l[k]]; vector<int> p=pmul(l,lc); vector<int> t(NK,0);
    for(size_t k=0;k<p.size();k++){ assert(p[k]<3); t[k]=p[k]; } S[L]=enc(t); }
  // heads
  size_t NHC=1; for(int k=0;k<NH;k++) NHC*=9; vector<vector<int>> heads;
  for(size_t h=0;h<NHC;h++){ vector<int> H(D+1,0); H[D]=1; size_t x=h; for(int k=KL+1;k<D;k++){ H[k]=x%9; x/=9; }
    auto code=[&](const vector<int>&z){ u64 c=0; for(int k=D-1;k>KL;k--) c=c*9+z[k]; return c; };
    u64 me=code(H); bool canon=true;
    for(int a=0;a<3&&canon;a++) for(int s=1;s<=2&&canon;s++) for(int cj=0;cj<2&&canon;cj++){
      vector<int> z=transform(H,a,s,cj); assert(z[D]==1); if(code(z)<me) canon=false; }
    if(canon) heads.push_back(H); }
  fprintf(stderr,"D=%d heads=%zu canonical=%zu lows=%zu\n",D,NHC,heads.size(),NLC);
  atomic<size_t> nxt(0); vector<thread> th; for(int i=0;i<nthr;i++) th.emplace_back(worker,&heads,&nxt); for(auto&t:th) t.join();
  // expand parent sets under affine maps t->s t+a on N (orbits)
  auto expand=[&](map<u64,int>& PP, int gmin){ set<u64> out; for(auto&pp:PP){ if(pp.second!=gmin) continue; u64 c=pp.first; vector<int> N(2*D+1); u64 x=c; for(int k=0;k<=2*D;k++){N[k]=x%3;x/=3;}
      for(int a=0;a<3;a++) for(int s=1;s<=2;s++){ int saveD=D; D=2*saveD; vector<int> M=transform(N,a,s,0); D=saveD; out.insert(codeN(M)); } } return out; };
  printf("D=%d maxAll=%d maxSF=%d\n",D,R.maxAll,R.maxSF);
  printf(" group-size counts (canonical heads, all parents):"); for(auto&p:R.cntAll) printf(" %d:%llu",p.first,(unsigned long long)p.second); printf("\n");
  printf(" group-size counts (canonical heads, squarefree):"); for(auto&p:R.cntSF) printf(" %d:%llu",p.first,(unsigned long long)p.second); printf("\n");
  // parents attaining max: need re-filter by size; store via witnesses for >=5, else via sets of size>=4
  { auto e=expand(R.parentsSF,R.maxSF); printf(" squarefree parents attaining maxSF=%d (orbit-expanded): %zu\n",R.maxSF,e.size()); }
  { auto e=expand(R.parentsAll,R.maxAll); printf(" all parents attaining maxAll=%d (orbit-expanded): %zu\n",R.maxAll,e.size()); }
  set<u64> p5all,p5sf; for(auto&w:R.wit){ printf(" WIT %s\n",w.c_str()); }
}
