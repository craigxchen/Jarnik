// Independent implementation of the "triple with common block" reduction for >=5 clusters on
// squarefree parents (constant ternary model). Phases are represented as truncated units
// z/conj(z) mod s^K (s=1/t) in the norm-one group, NOT as slopes; counting is done directly
// (no halving trick). Reports, for each D, the max cluster size found among clusters containing
// a triple with common block of degree >= ceil(D/10) (this captures every cluster of size >= 5).
#include <bits/stdc++.h>
using namespace std; typedef uint64_t u64;
int MUL[9][9],ADD[9][9],SUB[9][9],CONJ[9],INV[9];
int D,K,MAXB;
typedef vector<int> P; // low first, F_9 digits
u64 code(const P&p){u64 c=0;for(int k=p.size()-1;k>=0;k--)c=c*9+p[k];return c;}
P decode(u64 c){P p;while(c){p.push_back(c%9);c/=9;}return p;}
P mul(const P&a,const P&b){P r(a.size()+b.size()-1,0);for(size_t i=0;i<a.size();i++)if(a[i])for(size_t j=0;j<b.size();j++)r[i+j]=ADD[r[i+j]][MUL[a[i]][b[j]]];return r;}
bool pdivides(const P&f,const P&g,P&q){ // g | f, monic g
  P r=f; int n=f.size()-g.size(); if(n<0) return false; q.assign(n+1,0);
  for(int i=n;i>=0;i--){int c=r[i+g.size()-1]; q[i]=c; if(c) for(size_t j=0;j<g.size();j++) r[i+j]=SUB[r[i+j]][MUL[c][g[j]]];}
  for(size_t j=0;j+1<g.size();j++) if(r[j]) return false; return true;}
P conjP(const P&p){P r(p);for(auto&x:r)x=CONJ[x];return r;}
// phase group: truncated series in s, length K, constant 1; key = digits 1..K-1 base 9
int GS; unordered_map<u64,int> keyToIdx; vector<vector<int>> elems; vector<vector<uint16_t>> GM; vector<int> GINV;
vector<int> smul(const vector<int>&a,const vector<int>&b){vector<int> r(K,0);for(int i=0;i<K;i++)if(a[i])for(int j=0;i+j<K;j++)r[i+j]=ADD[r[i+j]][MUL[a[i]][b[j]]];return r;}
u64 skey(const vector<int>&a){u64 c=0;for(int i=K-1;i>=1;i--)c=c*9+a[i];return c;}
vector<int> phase(const P&p){ // p monic degree d: p*(s)/conj(p)*(s) mod s^K
  int d=p.size()-1; vector<int> a(K,0),b(K,0); for(int i=0;i<K&&i<=d;i++){a[i]=p[d-i];b[i]=CONJ[p[d-i]];}
  // a/b: b constant 1
  vector<int> q(K,0); for(int i=0;i<K;i++){int v=a[i]; for(int j=1;j<=i;j++) v=SUB[v][MUL[b[j]][q[i-j]]]; q[i]=v;} return q;}
int pidx(const vector<int>&s){ auto it=keyToIdx.find(skey(s)); int x= it==keyToIdx.end()? -1: it->second; if(x<0){fprintf(stderr,"phase not in group?\n");exit(1);} return x;}
struct Block{ u64 c; int deg; vector<u64> base; vector<int> plist; int ph; };
int main(int argc,char**argv){ D=atoi(argv[1]);
  for(int x=0;x<9;x++)for(int y=0;y<9;y++){int a=x%3,b=x/3,c=y%3,d=y/3;MUL[x][y]=((a*c-b*d)%3+3)%3+3*((a*d+b*c)%3);ADD[x][y]=(a+c)%3+3*((b+d)%3);SUB[x][y]=(a-c+3)%3+3*((b-d+3)%3);}
  for(int x=0;x<9;x++){CONJ[x]=x%3+3*((3-x/3)%3); for(int y=0;y<9;y++) if(MUL[x][y]==1) INV[x]=y;}
  K=(D+1)/2; int gmin=(D+9)/10; MAXB=D-gmin-K; int Gmax=D-(3*K+1)/2; int maxdeg=max(MAXB,Gmax);
  // build phase group by closure from phases of all monic linear and quadratic polys... safer: BFS closure of generators
  u64 KS=1; for(int i=1;i<K;i++) KS*=9; 
  // enumerate monic polys up to maxdeg, factor
  vector<vector<P>> irr(maxdeg+1); vector<Block> blocks[8];
  unordered_map<u64,vector<u64>> fac; fac.reserve(6000000); // code -> prime codes
  vector<vector<int>> gens;
  auto addElem=[&](const vector<int>&s){u64 k=skey(s); if(!keyToIdx.count(k)){keyToIdx[k]=elems.size(); elems.push_back(s);} };
  vector<int> one(K,0); one[0]=1; addElem(one);
  for(int d=1;d<=maxdeg;d++){ u64 lo=1; for(int i=0;i<d;i++) lo*=9;
    for(u64 c=lo;c<2*lo;c++){ P p=decode(c); vector<u64> pr; P q; bool found=false;
      for(int e=1;2*e<=d&&!found;e++) for(auto&f:irr[e]) if(pdivides(p,f,q)){ pr=fac[code(q)]; pr.push_back(code(f)); found=true; break; }
      if(!found){ irr[d].push_back(p); pr={c}; }
      fac[c]=pr; }
  }
  fac[1]={};
  // phases of irreducibles generate; closure
  vector<vector<int>> gen; for(int d=1;d<=maxdeg;d++) for(auto&f:irr[d]) { size_t before=elems.size(); addElem(phase(f)); if(elems.size()>before && gen.size()<40) gen.push_back(elems.back()); }
  for(size_t i=0;i<elems.size();i++) for(auto&g:gen) addElem(smul(elems[i],g));
  GS=elems.size(); fprintf(stderr,"D=%d K=%d group size %d (3^(K-1)=%d)\n",D,K,GS,(int)pow(3,K-1));
  GM.assign(GS,vector<uint16_t>(GS)); GINV.assign(GS,0);
  for(int i=0;i<GS;i++) for(int j=0;j<GS;j++){ GM[i][j]=pidx(smul(elems[i],elems[j])); if(GM[i][j]==0) GINV[i]=j; }
  // blocks
  for(int d=0;d<=maxdeg;d++){ u64 lo=1; for(int i=0;i<d;i++) lo*=9;
    for(u64 c=lo;c<2*lo;c++){ auto pr=fac[c]; bool ok=true; vector<u64> base;
      for(u64 f:pr){ u64 fb=code(conjP(decode(f))); if(fb==f){ok=false;break;} u64 k=min(f,fb); if(find(base.begin(),base.end(),k)!=base.end()){ok=false;break;} base.push_back(k);} 
      if(!ok) continue; sort(base.begin(),base.end());
      Block B; B.c=c; B.deg=d; B.base=base; B.ph=pidx(phase(decode(c)));
      vector<int> pl={0}; for(u64 f:pr){ int ph=pidx(phase(decode(f))), phi=GINV[ph]; vector<int> nl; for(int x:pl){nl.push_back(GM[x][ph]); nl.push_back(GM[x][phi]);} pl=nl; }
      B.plist=pl; blocks[d].push_back(B);} }
  auto disj=[&](const Block&a,const Block&b){ size_t i=0,j=0; while(i<a.base.size()&&j<b.base.size()){ if(a.base[i]==b.base[j]) return false; if(a.base[i]<b.base[j]) i++; else j++; } return true; };
  // bucket by phase of U (cluster condition on triple: phases of U_i equal)
  vector<vector<vector<int>>> bucket(maxdeg+1, vector<vector<int>>(GS));
  for(int d=0;d<=maxdeg;d++) for(size_t i=0;i<blocks[d].size();i++) bucket[d][blocks[d][i].ph].push_back(i);
  vector<int> hist(GS,0); vector<int> touched; int best=0; u64 ntrip=0,ncomp=0; map<int,u64> cnt; string bestw;
  for(int a=0;a<=MAXB;a++) for(int b=a;b<=MAXB;b++) for(int c=b;c<=MAXB;c++){ int g=D-a-b-c; if(a+b<K||g<gmin||g>maxdeg) continue;
    for(auto&U:blocks[a]) for(int vi:bucket[b][U.ph]){ const Block&V=blocks[b][vi]; if((a==b&&U.c>=V.c)||!disj(U,V)) continue;
      for(int ti:bucket[c][U.ph]){ const Block&T=blocks[c][ti]; if((b==c&&V.c>=T.c)||!disj(V,T)||!disj(U,T)) continue; ntrip++;
        for(int h:touched) hist[h]=0; touched.clear();
        for(int x:U.plist) for(int y:V.plist){ int xy=GM[x][y]; for(int z:T.plist){ int h=GM[xy][z]; if(hist[h]++==0) touched.push_back(h);} }
        int base=GM[GM[U.ph][GINV[V.ph]]][GINV[T.ph]]; // phase(U1 conj U2 conj U3)
        for(auto&G:blocks[g]){ if(!disj(G,U)||!disj(G,V)||!disj(G,T)) continue; ncomp++;
          int target=GM[base][G.ph]; int n=0; for(int y:G.plist) n+=hist[GM[target][GINV[y]]];
          cnt[n]++; if(n>best){best=n; char buf[200]; snprintf(buf,200,"U=%llu V=%llu T=%llu G=%llu",(unsigned long long)U.c,(unsigned long long)V.c,(unsigned long long)T.c,(unsigned long long)G.c); bestw=buf;} } } } }
  printf("D=%d best=%d triples=%llu completers=%llu  witness %s\n count-distribution:",D,best,(unsigned long long)ntrip,(unsigned long long)ncomp,bestw.c_str());
  for(auto&p:cnt) if(p.first>=4) printf(" %d:%llu",p.first,(unsigned long long)p.second); printf("\n");
}
