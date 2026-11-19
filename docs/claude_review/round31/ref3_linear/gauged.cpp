// Gauged census: real affine maps act simply transitively on nonreal elements of F_{p^2}, so WLOG one
// edge root is i. Stars 0,1 = the two stars containing i. Count unlabelled solutions containing i.
// Moment order K (5 = real system).  Elements encoded x+p*y.
#include <bits/stdc++.h>
using namespace std;
int p,q,K; vector<int> conjv; vector<array<int,6>> imk;
int mulc(int a,int b){int x=a%p,y=a/p,X=b%p,Y=b/p;return ((x*X-y*Y)%p+p)%p + p*((x*Y+y*X)%p);}
typedef array<int,5> Sig;
Sig sg(int a){Sig s{};for(int k=0;k<K;k++)s[k]=imk[a][k+1];return s;}
Sig addS(Sig a,Sig b){for(int k=0;k<K;k++)a[k]=(a[k]+b[k])%p;return a;}
Sig subS(Sig a,Sig b){for(int k=0;k<K;k++)a[k]=(a[k]-b[k]+p)%p;return a;}
uint64_t hs(Sig s){uint64_t h=0;for(int k=0;k<K;k++)h=h*p+s[k];return h;}
int canon(int a){return min(a,conjv[a]);}
int main(int argc,char**argv){
  p=atoi(argv[1]);K=atoi(argv[2]);int NP=argc>3?atoi(argv[3]):1;q=p*p;conjv.resize(q);imk.resize(q);
  for(int a=0;a<q;a++){int x=a%p,y=a/p;conjv[a]=x+p*((p-y)%p);int v=1;for(int k=1;k<=5;k++){v=mulc(v,a);imk[a][k]=v/p;}}
  int I=p; // 0+1i
  vector<int> nr; for(int a=0;a<q;a++) if(a/p!=0) nr.push_back(a);
  // pair table: unordered pairs {x,y} of nonreal, non-conjugate, distinct; key = signature sum
  unordered_map<uint64_t,vector<pair<int,int>>> PT;
  for(size_t s=0;s<nr.size();s++)for(size_t t=s+1;t<nr.size();t++){int x=nr[s],y=nr[t]; if(canon(x)==canon(y))continue;
     PT[hs(addS(sg(x),sg(y)))].push_back({x,y});}
  // stars containing i: {i,a,b,c}
  unordered_map<uint64_t,vector<array<int,3>>> ST;
  vector<int> el; for(int a:nr) if(canon(a)!=canon(I)) el.push_back(a);
  long long nst=0; long long cnt=0, s01=0;
  for(int pass=0;pass<NP;pass++){ ST.clear(); nst=0;
  for(size_t s=0;s<el.size();s++)for(size_t t=s+1;t<el.size();t++)for(size_t u=t+1;u<el.size();u++){
    int a=el[s],b=el[t],c=el[u]; if(canon(a)==canon(b)||canon(a)==canon(c)||canon(b)==canon(c))continue;
    uint64_t hh=hs(addS(addS(addS(sg(I),sg(a)),sg(b)),sg(c))); if((int)(hh%NP)!=pass) continue; ST[hh].push_back({a,b,c}); nst++; }
  fprintf(stderr,"p=%d K=%d istars=%lld classes=%zu pairs=%zu\n",p,K,nst,ST.size(),PT.size());
  for(auto &kv:ST){ auto &L=kv.second; if(L.size()<2) continue;
    Sig sig; { // recompute signature from first
      auto &f=L[0]; sig=addS(addS(addS(sg(I),sg(f[0])),sg(f[1])),sg(f[2])); }
    for(size_t x=0;x<L.size();x++)for(size_t y=0;y<L.size();y++){ if(x==y)continue;
      auto A=L[x],D=L[y]; bool bad=false; set<int> cc={canon(I)};
      for(int z:A){ if(cc.count(canon(z)))bad=true; cc.insert(canon(z)); }
      for(int z:D){ if(cc.count(canon(z)))bad=true; cc.insert(canon(z)); }
      if(bad)continue; s01++;
      array<int,3> pm={0,1,2};
      do{ // star 2+j contains A[j], D[pm[j]]; inner edges r23,r24,r34
        vector<pair<int,int>> *Lst[3]; bool empty=false;
        for(int j=0;j<3;j++){ Sig need=subS(subS(sig,sg(A[j])),sg(D[pm[j]])); auto it=PT.find(hs(need)); if(it==PT.end()){empty=true;break;} Lst[j]=&it->second; }
        if(empty) continue;
        // star2 = {A0,D.,r23,r24}, star3={A1,D.,r23,r34}, star4={A2,D.,r24,r34}
        for(auto pr:*Lst[0]) for(int o=0;o<2;o++){
          int r23=o?pr.second:pr.first, r24=o?pr.first:pr.second;
          for(auto pr3:*Lst[1]) for(int o3=0;o3<2;o3++){
            int a3=o3?pr3.second:pr3.first, r34=o3?pr3.first:pr3.second; if(a3!=r23) continue;
            // star4 needs {r24,r34}
            Sig need=subS(subS(sig,sg(A[2])),sg(D[pm[2]]));
            Sig have=addS(sg(r24),sg(r34)); if(hs(need)!=hs(have)) continue;
            set<int> c2=cc; bool b2=false; for(int z:{r23,r24,r34}){ if(c2.count(canon(z)))b2=true; c2.insert(canon(z)); }
            if(b2) continue;
            cnt++;
            if(cnt<=2*3){ printf("sol: i"); for(int z:A)printf(" %d+%di",z%p,z/p); printf(" | i"); for(int z:D)printf(" %d+%di",z%p,z/p);
               printf(" | inner %d+%di %d+%di %d+%di\n",r23%p,r23/p,r24%p,r24/p,r34%p,r34/p);} }
        }
      }while(next_permutation(pm.begin(),pm.end()));
    }
  }
  }
  printf("p=%d K=%d ordered(star0,star1) pairs=%lld solutions_containing_i(unlabelled)=%lld (raw %lld /2)\n",p,K,s01,cnt/2,cnt);
}
