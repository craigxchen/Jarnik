// Referee BFS census of relation flats (independent of ii_min_flat.cpp).
// usage: ./flats [mid] [threshold]
// Every admissible flat of rank k+1 is the closure of an admissible rank-k flat plus one generating row
// (admissibility = no two formal points equal, no four formal points collinear; inherited by subflats).
// mod P = 10007; all ranks are of integer vectors (norm <= sqrt6, sum zero): minors <= 6^3.5 < P.
#include <bits/stdc++.h>
using namespace std;
const int P = 10007;
typedef array<int,8> Row;
int inv_[P];
vector<Row> GEN, ALL4, MIDS;
struct Basis { int r=0; int piv[8]; Row row[8]; };
Row red(const Basis&b, Row v){ for(int t=0;t<b.r;t++){int p=b.piv[t]; int c=v[p]; if(!c) continue; for(int k=0;k<8;k++) v[k]=(v[k]+ (long)(P-c)*b.row[t][k])%P;} return v;}
bool zero(const Row&v){for(int x:v) if(x) return false; return true;}
bool add(Basis&b, Row v){ v=red(b,v); int p=-1; for(int k=0;k<8;k++) if(v[k]){p=k;break;} if(p<0) return false;
  long s=inv_[v[p]]; for(int k=0;k<8;k++) v[k]=v[k]*s%P;
  for(int t=0;t<b.r;t++){int c=b.row[t][p]; if(c) for(int k=0;k<8;k++) b.row[t][k]=(b.row[t][k]+(long)(P-c)*v[k])%P;}
  b.piv[b.r]=p; b.row[b.r]=v; b.r++; return true;}
typedef bitset<128> Key;
Key key(const Basis&b){Key k; for(size_t i=0;i<GEN.size();i++) if(zero(red(b,GEN[i]))) k.set(i); return k;}
struct KH{size_t operator()(const Key&k)const{return hash<string>()(k.to_string());}};
bool admissible(const Basis&b){ Row n[8]; for(int i=0;i<8;i++){Row e{};e[i]=1;n[i]=red(b,e);}
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){ if(n[i]==n[j]) return false;}
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)for(int k=j+1;k<8;k++)for(int l=k+1;l<8;l++){
    Row D[3]; int js[3]={j,k,l}; for(int a=0;a<3;a++) for(int c=0;c<8;c++) D[a][c]=(n[js[a]][c]-n[i][c]+P)%P;
    bool r1=true; for(int a=0;a<3&&r1;a++)for(int bb=a+1;bb<3&&r1;bb++)for(int c=0;c<8&&r1;c++)for(int d=c+1;d<8;d++){
      long m=((long)D[a][c]*D[bb][d]-(long)D[a][d]*D[bb][c])%P; if(m){r1=false;break;}}
    if(r1) return false; }
  return true;}
int cnt(const Basis&b,const vector<Row>&R){int c=0; for(auto&v:R) if(zero(red(b,v))) c++; return c;}
Row mk(vector<pair<int,int>> d){Row v{}; for(auto [k,c]:d) v[k]=((v[k]+c)%P+P)%P; return v;}
int main(int argc,char**argv){ bool mid=argc>1&&string(argv[1])=="mid"; int thr=argc>2?atoi(argv[2]):11;
  inv_[1]=1; for(int i=2;i<P;i++) inv_[i]=(long)(P-P/i)*inv_[P%i]%P;
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)for(int k=j+1;k<8;k++)for(int l=k+1;l<8;l++){
    GEN.push_back(mk({{i,1},{l,1},{j,-1},{k,-1}}));
    ALL4.push_back(mk({{i,1},{l,1},{j,-1},{k,-1}}));ALL4.push_back(mk({{i,1},{j,1},{k,-1},{l,-1}}));ALL4.push_back(mk({{i,1},{k,1},{j,-1},{l,-1}}));}
  if(mid) for(int i=0;i<8;i++)for(int j=i+1;j<8;j++)for(int k=j+1;k<8;k++) GEN.push_back(mk({{i,1},{k,1},{j,-2}}));
  for(int j=0;j<8;j++)for(int i=0;i<8;i++)for(int k=i+1;k<8;k++) if(i!=j&&k!=j) MIDS.push_back(mk({{i,1},{k,1},{j,-2}}));
  vector<Basis> level{Basis()}; 
  for(int rank=1;rank<=5;rank++){
    unordered_set<Key,KH> seen; vector<Basis> next; map<tuple<int,int,int>,long> hist; long tot=0;
    for(auto&b:level){ Key k0=key(b);
      for(size_t r=0;r<GEN.size();r++){ if(k0[r]) continue; Basis c=b; add(c,GEN[r]); Key k=key(c);
        if(!seen.insert(k).second) continue; tot++;
        if(!admissible(c)) continue;
        int n70=k.count(), n4=cnt(c,ALL4), nm=cnt(c,MIDS);
        hist[{n70,n4,nm}]++;
        int score = mid? n4+nm : n4;
        if(score>=thr){ cout<<"REC rank "<<rank<<" n_gen "<<n70<<" n4 "<<n4<<" nmid "<<nm<<" rows";
          for(size_t i=0;i<GEN.size();i++) if(k[i]) cout<<' '<<i; cout<<"\n"; }
        if(rank<5) next.push_back(c); } }
    long adm=0; int mx4=0,mxall=0; for(auto&[t,n]:hist){adm+=n; mx4=max(mx4,get<1>(t)); mxall=max(mxall,get<1>(t)+get<2>(t));}
    cerr<<"rank "<<rank<<" distinct closures "<<tot<<" admissible "<<adm<<" max n4 "<<mx4<<" max n4+nmid "<<mxall<<"\n";
    cout<<"HIST rank "<<rank; for(auto&[t,n]:hist) cout<<" ("<<get<0>(t)<<","<<get<1>(t)<<","<<get<2>(t)<<"):"<<n; cout<<"\n";
    level.swap(next);} }
