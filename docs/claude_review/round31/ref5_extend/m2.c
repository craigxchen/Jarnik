// Exact M_2(D,c) for the anisotropic norm A^2+AB+B^2 of F_4(t)/F_2(t).
// Points: monic z=A+B*theta, deg A = D (lead 1), deg B < D.  Cluster: same N and same coeffs (of A and B) in degrees > K=floor(D/2+c).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned long long u64;
static u64 clmul(u64 a,u64 b){u64 r=0;while(b){if(b&1)r^=a;a<<=1;b>>=1;}return r;}
static u64 sq(u64 a){u64 r=0;int i=0;while(a){if(a&1)r|=1ULL<<(2*i);a>>=1;i++;}return r;}
static int cmp(const void*x,const void*y){u64 a=*(u64*)x,b=*(u64*)y;return a<b?-1:a>b;}
int main(int argc,char**argv){
  int D=atoi(argv[1]); int c2=argc>2?atoi(argv[2]):0; // c2 = 2c
  int K=(D+c2)/2; // floor(D/2+c)
  if(K>=D){printf("trivial\n");return 0;}
  int nb=K+1; u64 NB=1ULL<<nb; // bottom bits per coordinate
  int nt=D-K-1; u64 NT=1ULL<<nt;
  u64 *nbot=malloc(sizeof(u64)*NB*NB), *vals=malloc(sizeof(u64)*NB*NB);
  for(u64 a=0;a<NB;a++)for(u64 b=0;b<NB;b++) nbot[a*NB+b]=sq(a)^clmul(a,b)^sq(b);
  int best=0; u64 bestN=0,bestAt=0,bestBt=0; long long nbestclasses=0;
  long long hist[64]={0};
  for(u64 at=0;at<NT;at++)for(u64 bt=0;bt<NT;bt++){
    u64 At=((1ULL<<nt)|at)<<nb, Bt=bt<<nb;
    u64 ntop=sq(At)^clmul(At,Bt)^sq(Bt);
    // cross = At*Bb + Bt*Ab  (linear in bottom)
    u64 cA[64],cB[64];
    for(int k=0;k<nb;k++){cA[k]=clmul(Bt,1ULL<<k);cB[k]=clmul(At,1ULL<<k);}
    u64 n=0;
    for(u64 a=0;a<NB;a++){u64 xa=0;for(int k=0;k<nb;k++)if(a>>k&1)xa^=cA[k];
      for(u64 b=0;b<NB;b++){u64 xb=0;for(int k=0;k<nb;k++)if(b>>k&1)xb^=cB[k];
        vals[n++]=ntop^xa^xb^nbot[a*NB+b];}}
    qsort(vals,n,sizeof(u64),cmp);
    u64 i=0;while(i<n){u64 j=i;while(j<n&&vals[j]==vals[i])j++;int r=j-i;hist[r<63?r:63]++;
      if(r>best){best=r;bestN=vals[i];bestAt=At;bestBt=Bt;nbestclasses=0;} if(r==best)nbestclasses++; i=j;}
  }
  printf("D=%d c=%g K=%d M=%d #classes_attaining=%lld exampleN=%llx At=%llx Bt=%llx\n",D,c2/2.0,K,best,nbestclasses,bestN,bestAt,bestBt);
  printf("hist:");for(int r=1;r<64;r++)if(hist[r])printf(" %d:%lld",r,hist[r]);printf("\n");
  return 0;}
