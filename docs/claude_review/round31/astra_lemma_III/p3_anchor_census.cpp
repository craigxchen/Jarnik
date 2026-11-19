// Affine-gauged strict rational-point census of the linear ten-pair moment system.
// Finite field element x+i*y encoded x+p*y, p=3 mod 4. No random sampling.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;
using Trip=array<int,3>;
int p,q;
vector<int> mul,plus_,minus_,conj_,inv_,can;
vector<array<int,5>> moments;
int mult(int a,int b){return mul[a*q+b];}
int add(int a,int b){return plus_[a*q+b];}
int sub(int a,int b){return minus_[a*q+b];}
int divide(int a,int b){assert(b);return mult(a,inv_[b]);}
int normval(int z,int w){return mult(sub(z,w),sub(z,conj_[w]));}
int key(int a,int b,int c){
  long long val=0;
  for(int j=0;j<5;j++){int v=moments[a][j]+moments[b][j]+moments[c][j]; if(v>=p)v-=p; if(v>=p)v-=p; val=val*p+v;}
  return int(val);
}
struct Cub{int c2,c1;};
Cub cubic(int a,int b,int c){return {sub(0,add(add(a,b),c)),add(add(mult(a,b),mult(a,c)),mult(b,c))};}
int pairroot(array<int,2> x,array<int,2> y){ // monic quadratics coefficients of t, 1
  int m=sub(x[0],y[0]); if(!m)return -1;
  int r=divide(sub(y[1],x[1]),m);
  if(add(mult(r,add(r,x[0])),x[1]))return -1;
  return r;
}
long long nt=0,nu=0,pairs_uv=0,one_completion=0,partial=0,full=0;
void process(const Trip &u,const Trip &v){
  for(int a:u)for(int b:v)if(can[a]==can[b])return;
  int h=mult(conj_[mult(mult(u[0],u[1]),u[2])],mult(mult(v[0],v[1]),v[2]))/p;
  assert(h!=0); // otherwise U and V have identical full signed root divisor
  ++pairs_uv;
  for(int a=p;a<q;++a){
    bool collision=false;for(int x:u)collision|=can[x]==can[a];for(int x:v)collision|=can[x]==can[a];if(collision)continue;
    int z=conj_[a],U=1,V=1,nv[3],nu_[3];
    for(int j=0;j<3;j++){U=mult(U,sub(z,u[j]));V=mult(V,sub(z,v[j]));nu_[j]=normval(z,u[j]);nv[j]=normval(z,v[j]);}
    int admiss[3]={0,0,0},A[3][3];
    bool skip=false;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){int L=mult(nv[j],U),R=mult(nu_[i],V); if(L==R)continue;
        int aij=divide(R,sub(R,L));if(aij<p){admiss[i]|=1<<j;A[i][j]=aij;}
      }
      if(!admiss[i]){skip=true;break;}
    }
    if(skip)continue;
    ++one_completion;
    Trip perm={0,1,2};
    do{
      int aa[3];bool ok=true;
      for(int j=0;j<3;j++){if(!(admiss[j]&(1<<perm[j]))){ok=false;break;} aa[j]=A[j][perm[j]];}
      if(!ok||aa[0]==aa[1]||aa[0]==aa[2]||aa[1]==aa[2])continue;
      ++partial;
      array<array<int,2>,3> C;
      for(int j=0;j<3;j++){
        int l=(j+1)%3,m=(j+2)%3;
        Cub L=cubic(conj_[v[perm[j]]],u[l],u[m]), R=cubic(conj_[u[j]],v[perm[l]],v[perm[m]]);
        int b=add(R.c2,mult(aa[j],sub(L.c2,R.c2))),c=add(R.c1,mult(aa[j],sub(L.c1,R.c1)));
        int c1=add(b,z);C[j]={c1,add(c,mult(z,c1))};
      }
      int r[3]; bool good=true;
      for(int j=0;j<3;j++){
         r[j]=pairroot(C[j],C[(j+1)%3]); if(r[j]<p||can[r[j]]==can[a]){good=false;break;}
         for(int t:u)if(can[t]==can[r[j]])good=false;
         for(int t:v)if(can[t]==can[r[j]])good=false;
         for(int n=0;n<j;n++)if(can[r[n]]==can[r[j]])good=false;
         if(!good)break;
      }
      if(!good)continue;
      ++full;
      if(full<=5){
        array<Trip,3> star={{u,v,Trip{r[0],r[1],r[2]}}};
        cout<<"sample a "<<a<<" u ";for(int x:u)cout<<x<<",";cout<<" v ";for(int x:perm)cout<<v[x]<<",";
        cout<<" inner ";for(int x:r)cout<<x<<",";cout<<"\n";
        // Check all five original moment orders in all five stars from scratch.
        for(int k=0;k<5;k++){int sig=(moments[a][k]+moments[u[0]][k]+moments[u[1]][k]+moments[u[2]][k])%p;
          assert(sig==(moments[a][k]+moments[v[0]][k]+moments[v[1]][k]+moments[v[2]][k])%p);
          for(int j=0;j<3;j++)assert(sig==(moments[u[j]][k]+moments[v[perm[j]]][k]+moments[r[j]][k]+moments[r[(j+2)%3]][k])%p);
        }
      }
    }while(next_permutation(perm.begin(),perm.end()));
  }
}
int main(int argc,char**argv){
  assert(argc==2);p=atoi(argv[1]);assert(p>3&&p%4==3&&p<=43);for(int d=2;d*d<=p;d++)assert(p%d);
  q=p*p; mul.resize(q*q);plus_.resize(q*q);minus_.resize(q*q);conj_.resize(q);inv_.resize(q);can.resize(q);moments.resize(q);
  for(int a=0;a<q;a++){int x=a%p,y=a/p;conj_[a]=x+p*((p-y)%p);can[a]=min(a,conj_[a]);
     for(int b=0;b<q;b++){int X=b%p,Y=b/p;int m=(x*X+p*p-y*Y)%p+p*((x*Y+y*X)%p);mul[a*q+b]=m;if(m==1)inv_[a]=b;
       plus_[a*q+b]=(x+X)%p+p*((y+Y)%p);minus_[a*q+b]=(x+p-X)%p+p*((y+p-Y)%p);
     }
  }
  for(int a=0;a<q;a++){int v=1;for(int j=0;j<5;j++){v=mult(v,a);moments[a][j]=v/p;}}
  unordered_map<int,vector<Trip>> U;
  int i=p;
  for(int a=p;a<q;a++)for(int b=p;b<a;b++)if(can[a]!=can[b]&&can[a]!=can[i]&&can[b]!=can[i]){U[key(i,a,b)].push_back({i,b,a});nu++;}
  for(int a=p;a<q;a++)for(int b=p;b<a;b++)if(can[b]!=can[a])for(int c=p;c<b;c++)if(can[c]!=can[a]&&can[c]!=can[b]){
    nt++;auto it=U.find(key(a,b,c));if(it!=U.end()){Trip v={c,b,a};for(const Trip&u:it->second)process(u,v);}
  }
  cout<<"p="<<p<<" anchored U="<<nu<<" V="<<nt<<" disjoint moment-pairs="<<pairs_uv
      <<" pivots passing all rows="<<one_completion<<" signed matchings="<<partial<<" strict completions="<<full<<"\n";
}
