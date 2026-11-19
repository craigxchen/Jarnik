// Referee's own brute force for VII.4 (binary model F_4[t]/F_2[t], theta^2+theta+1=0).
// Marks z_x = prod_{v in R} (chi_v(x)=+1 ? f_v : conj f_v), rows x in F_2^r, chi_v(x)=(-1)^{v.x}.
// Enumerate all monic nonreal f_v of degree 1..dmax; remove the full gcd of the marks; test the
// critical condition deg(z_x - z_y) <= D/2 for all pairs (D = reduced degree), distinctness.
// Report clusters, split by whether the base-prime supports of the columns are pairwise disjoint.
// Polynomials: p = a + theta*b, a,b in F_2[t] as bitmasks.
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <string>
using namespace std;
typedef uint64_t u64;
struct P { u64 a, b; };
static inline u64 clmul(u64 x, u64 y){ u64 r=0; while(y){ if(y&1) r^=x; x<<=1; y>>=1;} return r; }
static inline int dg(const P&p){ u64 m=p.a|p.b; return m? 63-__builtin_clzll(m) : -1; }
static inline P mul(const P&x,const P&y){ u64 ac=clmul(x.a,y.a), bd=clmul(x.b,y.b), ad=clmul(x.a,y.b), bc=clmul(x.b,y.a); return {ac^bd, ad^bc^bd}; }
static inline P conj(const P&p){ return {p.a^p.b, p.b}; }
static inline P add(const P&x,const P&y){ return {x.a^y.a, x.b^y.b}; }
// scalar e = (x + theta y), e in {1,theta,theta^2} encoded 1,2,3 as (x,y)
static inline P smul(int e,const P&p){ // e: 1 -> 1, 2 -> theta, 3 -> theta^2 = 1+theta
  if(e==1) return p;
  P th={p.b, p.a^p.b}; // theta*p
  if(e==2) return th;
  return {p.a^th.a, p.b^th.b}; // (1+theta)p
}
static inline int lead(const P&p,int d){ int x=(p.a>>d)&1, y=(p.b>>d)&1; return x+2*y; } // 1,2,3
static inline int inv(int e){ return e==1?1:(e==2?3:2); }
static inline int emul(int e,int f){ // in F_4 with codes 1=1,2=theta,3=theta^2
  int le = e==1?0:(e==2?1:2), lf = f==1?0:(f==2?1:2); int l=(le+lf)%3; return l==0?1:(l==1?2:3); }
static P pmod(P r,const P&s){ int ds=dg(s); int ls=inv(lead(s,ds)); int dr;
  while((dr=dg(r))>=ds){ int c=emul(lead(r,dr),ls); P t=smul(c,s); t.a<<=(dr-ds); t.b<<=(dr-ds); r=add(r,t);} return r; }
static P gcdp(P x,P y){ while(dg(y)>=0){ P r=pmod(x,y); x=y; y=r;} return x; }
static int degGcd(const P&x,const P&y){ return dg(gcdp(x,y)); }

int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage: r dmax col1 col2 ...\n"); return 1; }
  int r=atoi(argv[1]), dmax=atoi(argv[2]); int m=1<<r;
  vector<int> cols; for(int i=3;i<argc;i++) cols.push_back(atoi(argv[i]));
  int n=cols.size();
  vector<P> polys; for(int d=1; d<=dmax; d++) for(u64 lo=0; lo<(1ULL<<d); lo++) for(u64 b=1;b<(1ULL<<d);b++) polys.push_back({(1ULL<<d)|lo, b});
  int NP=polys.size();
  vector<vector<int>> sg(n, vector<int>(m));
  for(int j=0;j<n;j++) for(int x=0;x<m;x++) sg[j][x]=__builtin_parity(cols[j]&x)?-1:1;
  // rows distinguished?
  bool dist=true; for(int x=0;x<m;x++) for(int y=x+1;y<m;y++){ bool same=true; for(int j=0;j<n;j++) if(sg[j][x]!=sg[j][y]) same=false; if(same) dist=false; }
  printf("r=%d m=%d dmax=%d cols:",r,m,dmax); for(int c:cols) printf(" %d",c); printf("  rows distinguished by sign pattern: %s  #polys=%d\n", dist?"yes":"no", NP);
  vector<int> idx(n,0); long long tot=0, clDisj=0, clShared=0, printed=0;
  vector<P> z(m);
  while(true){
    tot++;
    for(int x=0;x<m;x++){ P acc={1,0}; for(int j=0;j<n;j++){ const P&f=polys[idx[j]]; acc=mul(acc, sg[j][x]>0? f: conj(f)); } z[x]=acc; }
    int D0=dg(z[0]);
    bool distinct=true; for(int x=0;x<m && distinct;x++) for(int y=x+1;y<m;y++) if(z[x].a==z[y].a && z[x].b==z[y].b){distinct=false;break;}
    if(distinct){
      int gub=degGcd(z[0],z[1]);
      bool ok=true;
      for(int x=0;x<m && ok;x++) for(int y=x+1;y<m;y++){ int dd=dg(add(z[x],z[y])); if(2*dd> D0+gub){ok=false;break;} }
      if(ok){
        P G=z[0]; for(int x=1;x<m;x++) G=gcdp(G,z[x]); int g=dg(G);
        for(int x=0;x<m && ok;x++) for(int y=x+1;y<m;y++){ int dd=dg(add(z[x],z[y])); if(2*dd> D0+g){ok=false;break;} }
        if(ok){
          // disjoint base-prime supports?  gcd of norms f fbar across columns
          bool disj=true;
          for(int j=0;j<n && disj;j++) for(int k=j+1;k<n;k++){ const P&f=polys[idx[j]]; const P&h=polys[idx[k]];
            if(degGcd(mul(f,conj(f)), mul(h,conj(h)))>0){disj=false;break;} }
          if(disj) clDisj++; else clShared++;
          if(printed<6){ printed++; printf("  CLUSTER (%s): D0=%d g=%d D=%d  f:", disj?"disjoint":"shared", D0, g, D0-g);
            for(int j=0;j<n;j++){ const P&f=polys[idx[j]]; printf(" [a=%llx b=%llx]",(unsigned long long)f.a,(unsigned long long)f.b);} printf("\n"); }
        }
      }
    }
    int j=0; while(j<n){ idx[j]++; if(idx[j]<NP) break; idx[j]=0; j++; } if(j==n) break;
  }
  printf("  assignments=%lld  critical clusters: disjoint-support=%lld  shared-support=%lld\n", tot, clDisj, clShared);
  return 0;
}
