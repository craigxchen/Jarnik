// Split-form census over F_p: two independent root lists u_e, v_e in F_p, all 20 values distinct,
// sum_{e in star i}(u_e^k - v_e^k) independent of i for k=1..K. Gauge: affine maps act simply transitively on
// ordered pairs (u,v), u!=v, so WLOG one edge is (0,1). Same star-0/star-1 + pair-table search as gauged.cpp.
#include <bits/stdc++.h>
using namespace std;
int p,K;
typedef array<int,5> Sig;
vector<Sig> sgv; // sig of element e=(u,v): u^k - v^k
Sig addS(Sig a,Sig b){for(int k=0;k<K;k++)a[k]=(a[k]+b[k])%p;return a;}
Sig subS(Sig a,Sig b){for(int k=0;k<K;k++)a[k]=(a[k]-b[k]+p)%p;return a;}
uint64_t hs(Sig s){uint64_t h=0;for(int k=0;k<K;k++)h=h*p+s[k];return h;}
int U(int e){return e/p;} int V(int e){return e%p;}
bool disjointvals(initializer_list<int> es){ vector<int> vals; for(int e:es){vals.push_back(U(e));vals.push_back(V(e));}
  sort(vals.begin(),vals.end()); return adjacent_find(vals.begin(),vals.end())==vals.end(); }
int main(int argc,char**argv){
  p=atoi(argv[1]);K=atoi(argv[2]);int NP=argc>3?atoi(argv[3]):1;
  sgv.resize(p*p);
  for(int u=0;u<p;u++)for(int v=0;v<p;v++){Sig s{};long long a=1,b=1;for(int k=0;k<K;k++){a=a*u%p;b=b*v%p;s[k]=((a-b)%p+p)%p;}sgv[u*p+v]=s;}
  vector<int> el; for(int u=0;u<p;u++)for(int v=0;v<p;v++) if(u!=v) el.push_back(u*p+v);
  int I=0*p+1;
  unordered_map<uint64_t,vector<pair<int,int>>> PT;
  for(size_t s=0;s<el.size();s++)for(size_t t=s+1;t<el.size();t++){int x=el[s],y=el[t]; if(!disjointvals({x,y}))continue; PT[hs(addS(sgv[x],sgv[y]))].push_back({x,y});}
  vector<int> el2; for(int e:el) if(disjointvals({e,I})) el2.push_back(e);
  long long cnt=0,s01=0,nst=0;
  for(int pass=0;pass<NP;pass++){
    unordered_map<uint64_t,vector<array<int,3>>> ST; nst=0;
    for(size_t s=0;s<el2.size();s++)for(size_t t=s+1;t<el2.size();t++){ if(!disjointvals({el2[s],el2[t]}))continue;
      for(size_t w=t+1;w<el2.size();w++){int a=el2[s],b=el2[t],c=el2[w]; if(!disjointvals({a,b,c}))continue;
        uint64_t hh=hs(addS(addS(addS(sgv[I],sgv[a]),sgv[b]),sgv[c])); if((int)(hh%NP)!=pass)continue; ST[hh].push_back({a,b,c}); nst++;}}
    fprintf(stderr,"pass %d stars %lld\n",pass,nst);
    for(auto &kv:ST){auto &L=kv.second; if(L.size()<2)continue;
      Sig sig=addS(addS(addS(sgv[I],sgv[L[0][0]]),sgv[L[0][1]]),sgv[L[0][2]]);
      for(size_t x=0;x<L.size();x++)for(size_t y=0;y<L.size();y++){ if(x==y)continue; auto A=L[x],D=L[y];
        if(!disjointvals({I,A[0],A[1],A[2],D[0],D[1],D[2]}))continue; s01++;
        array<int,3> pm={0,1,2};
        do{ vector<pair<int,int>>* Ls[2]; bool emp=false;
          for(int j=0;j<2;j++){auto it=PT.find(hs(subS(subS(sig,sgv[A[j]]),sgv[D[pm[j]]]))); if(it==PT.end()){emp=true;break;} Ls[j]=&it->second;}
          if(emp)continue;
          uint64_t need4=hs(subS(subS(sig,sgv[A[2]]),sgv[D[pm[2]]]));
          for(auto pr:*Ls[0])for(int o=0;o<2;o++){int r23=o?pr.second:pr.first,r24=o?pr.first:pr.second;
            for(auto pr3:*Ls[1])for(int o3=0;o3<2;o3++){int a3=o3?pr3.second:pr3.first,r34=o3?pr3.first:pr3.second; if(a3!=r23)continue;
              if(hs(addS(sgv[r24],sgv[r34]))!=need4)continue;
              if(!disjointvals({I,A[0],A[1],A[2],D[0],D[1],D[2],r23,r24,r34}))continue;
              cnt++; if(cnt<=4){printf("sol: (0,1)"); for(int e:A)printf(" (%d,%d)",U(e),V(e)); printf(" |"); for(int e:D)printf(" (%d,%d)",U(e),V(e));
                 printf(" | inner (%d,%d) (%d,%d) (%d,%d)\n",U(r23),V(r23),U(r24),V(r24),U(r34),V(r34));}
            }}
        }while(next_permutation(pm.begin(),pm.end()));
      }}
  }
  printf("split p=%d K=%d ordered(star0,star1)=%lld solutions_containing_(0,1)=%lld (raw %lld /2)\n",p,K,s01,cnt/2,cnt);
}
