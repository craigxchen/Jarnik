// Independent census of the strict twisted linear ten-pair system over F_{p^2}, p = 3 mod 4.
// Twisted: barred root = Frobenius conjugate z^p. Strict: all 20 values r_e, r_e^p distinct (so 10 distinct
// conjugate pairs, none real). Moment condition: Im( sum_{e in star i} r_e^k ) equal for all five stars, k=1..K
// (K=5 is the real system; smaller K is used as a positive control).
// Method: enumerate all oriented 4-stars (4 distinct conjugate pairs, one orientation each), bucket by signature,
// then search 5 stars in a bucket that pairwise meet in exactly one element, no element in 3 stars,
// and 10 pairwise non-conjugate elements overall.  Independent of Astra's anchor/pivot elimination.
#include <bits/stdc++.h>
using namespace std;
int p,q,K;
struct E{int x,y;};
vector<int> conjv; vector<array<int,6>> imk; // imk[a][k]= Im(a^k) for k=1..5
int enc(int x,int y){return x+p*y;}
int mulc(int a,int b){int x=a%p,y=a/p,X=b%p,Y=b/p;return ((x*X-y*Y)%p+p)%p + p*((x*Y+y*X)%p);}
int main(int argc,char**argv){
  p=atoi(argv[1]); K=atoi(argv[2]); long long maxprint = argc>3?atoll(argv[3]):5; q=p*p;
  conjv.resize(q); imk.resize(q);
  for(int a=0;a<q;a++){int x=a%p,y=a/p; conjv[a]=x+p*((p-y)%p); int v=1; for(int k=1;k<=5;k++){v=mulc(v,a); imk[a][k]=v/p;}}
  // nonreal pair representatives: elements with y in 1..(p-1)/2 ; orientation gives a or conj(a)
  vector<int> reps; for(int a=0;a<q;a++){int y=a/p; if(y>=1 && y<=(p-1)/2) reps.push_back(a);}
  int P=reps.size(); fprintf(stderr,"p=%d pairs=%d\n",p,P);
  // star = 4 oriented elements. signature = K values mod p
  vector<pair<uint64_t,array<int,4>>> stars; stars.reserve((size_t)P*(P-1)/2*(P-2)/3*(P-3)/4*16+10);
  for(int a=0;a<P;a++)for(int b=a+1;b<P;b++)for(int c=b+1;c<P;c++)for(int d=c+1;d<P;d++)
   for(int m=0;m<16;m++){
     array<int,4> s={reps[a],reps[b],reps[c],reps[d]};
     for(int t=0;t<4;t++) if(m>>t&1) s[t]=conjv[s[t]];
     uint64_t sig=0; for(int k=1;k<=K;k++){int v=0; for(int t=0;t<4;t++) v+=imk[s[t]][k]; sig=sig*p+(v%p);} 
     sort(s.begin(),s.end()); stars.push_back({sig,s});
   }
  sort(stars.begin(),stars.end());
  size_t ns=stars.size(); size_t ncls=0,maxc=0; long long sols=0;
  auto canon=[&](int a){return min(a,conjv[a]);};
  for(size_t i=0;i<ns;){
    size_t j=i; while(j<ns&&stars[j].first==stars[i].first) j++;
    ncls++; maxc=max(maxc,j-i);
    int n=j-i; vector<array<int,4>> S(n); for(int t=0;t<n;t++) S[t]=stars[i+t].second;
    // compatibility: exactly one common element, and no element of one conjugate to an element of the other
    vector<vector<char>> ok(n,vector<char>(n,0)); vector<vector<int>> com(n,vector<int>(n,-1));
    for(int a=0;a<n;a++)for(int b=a+1;b<n;b++){
      int c=0,ce=-1; bool bad=false;
      for(int u:S[a])for(int v:S[b]){ if(u==v){c++;ce=u;} else if(canon(u)==canon(v)) bad=true; }
      if(c==1&&!bad){ok[a][b]=ok[b][a]=1;com[a][b]=com[b][a]=ce;}
    }
    // 5-cliques with distinct common elements (so no element in three stars)
    for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)if(ok[a][b])for(int c=b+1;c<n;c++)if(ok[a][c]&&ok[b][c])
     for(int d=c+1;d<n;d++)if(ok[a][d]&&ok[b][d]&&ok[c][d])for(int e=d+1;e<n;e++)if(ok[a][e]&&ok[b][e]&&ok[c][e]&&ok[d][e]){
       int id[5]={a,b,c,d,e}; set<int> el; for(int x=0;x<5;x++)for(int y=x+1;y<5;y++) el.insert(com[id[x]][id[y]]);
       if(el.size()!=10) continue;
       // all star elements must be shared ones (automatic when 10 distinct shared elements) ; verify
       bool good=true; for(int x=0;x<5;x++) for(int u:S[id[x]]) if(!el.count(u)) good=false;
       set<int> cc; for(int u:el) cc.insert(canon(u)); if(cc.size()!=10) good=false;
       if(!good) continue;
       sols++;
       if(sols<=maxprint){ printf("sol:"); for(int x=0;x<5;x++){printf(" {");for(int u:S[id[x]])printf("%d+%di ",u%p,u/p);printf("}");} printf("\n"); }
     }
    i=j;
  }
  printf("p=%d K=%d stars=%zu classes=%zu maxclass=%zu unlabelled_solutions=%lld\n",p,K,ns,ncls,maxc,sols);
}
