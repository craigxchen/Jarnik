// Exact M_3(D,c): anisotropic norm A^2+B^2 of F_9(t)/F_3(t), z=A+iB monic (deg A=D lead 1, deg B<D).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u64;
static u64 ONE;
static u64 add3(u64 a,u64 b){u64 s=a+b;u64 m=((s+ONE)>>2)&ONE;return s-3*m;}
static u64 pack(int*d,int n){u64 r=0;for(int i=n-1;i>=0;i--)r=(r<<3)|(u64)d[i];return r;}
static void pmul(int*a,int na,int*b,int nb,int*r){memset(r,0,sizeof(int)*(na+nb));for(int i=0;i<na;i++)for(int j=0;j<nb;j++)r[i+j]=(r[i+j]+a[i]*b[j])%3;}
static int cmp(const void*x,const void*y){u64 a=*(u64*)x,b=*(u64*)y;return a<b?-1:a>b;}
int main(int argc,char**argv){
  int D=atoi(argv[1]); int c2=argc>2?atoi(argv[2]):0; int K=(D+c2)/2;
  int nd=2*D+1; ONE=0; for(int i=0;i<nd;i++)ONE|=((u64)1)<<(3*i);
  if(3*nd>127){printf("too big\n");return 1;}
  int nb=K+1, nt=D-K-1; int ndig=2*nb;
  u64 NBT=1; for(int i=0;i<ndig;i++)NBT*=3;
  u64 NTT=1; for(int i=0;i<2*nt;i++)NTT*=3;
  u64 *nbot=malloc(16*NBT),*vals=malloc(16*NBT);
  int dig[64],A[64],B[64],r1[128],r2[128],tmp[128];
  memset(dig,0,sizeof dig);
  for(u64 idx=0;idx<NBT;idx++){
    u64 x=idx; for(int j=0;j<ndig;j++){dig[j]=x%3;x/=3;}
    pmul(dig,nb,dig,nb,r1); pmul(dig+nb,nb,dig+nb,nb,r2);
    for(int i=0;i<2*nb-1;i++)tmp[i]=(r1[i]+r2[i])%3; for(int i=2*nb-1;i<nd;i++)tmp[i]=0;
    nbot[idx]=pack(tmp,nd);
  }
  int best=0; long long nbest=0; long long hist[64]={0}; u64 bestN=0; int bestTop[64];
  int td[64];
  for(u64 tix=0;tix<NTT;tix++){
    u64 x=tix; for(int j=0;j<2*nt;j++){td[j]=x%3;x/=3;}
    memset(A,0,sizeof A);memset(B,0,sizeof B);
    for(int k=0;k<nt;k++){A[nb+k]=td[k];B[nb+k]=td[nt+k];} A[D]=1;
    pmul(A,D+1,A,D+1,r1); pmul(B,D+1,B,D+1,r2);
    for(int i=0;i<nd;i++)tmp[i]=(r1[i]+r2[i])%3; u64 ntop=pack(tmp,nd);
    u64 v[64];
    for(int j=0;j<ndig;j++){ int *P=(j<nb)?A:B; int k=j%nb; memset(tmp,0,sizeof(int)*nd);
      for(int i=0;i<=D;i++)tmp[i+k]=(2*P[i])%3; v[j]=pack(tmp,nd);}
    // odometer
    memset(dig,0,sizeof dig); u64 cross=0;
    for(u64 idx=0;idx<NBT;idx++){
      vals[idx]=add3(add3(ntop,cross),nbot[idx]);
      int j=0; while(j<ndig){cross=add3(cross,v[j]); dig[j]++; if(dig[j]<3)break; dig[j]=0; j++;}
    }
    qsort(vals,NBT,16,cmp);
    u64 i=0;while(i<NBT){u64 j=i;while(j<NBT&&vals[j]==vals[i])j++;int r=j-i;hist[r<63?r:63]++;
      if(r>best){best=r;nbest=0;bestN=vals[i];memcpy(bestTop,td,sizeof td);} if(r==best)nbest++; i=j;}
  }
  printf("p=3 D=%d c=%g K=%d M=%d #classes=%lld\n",D,c2/2.0,K,best,nbest);
  printf("hist:");for(int r=1;r<64;r++)if(hist[r])printf(" %d:%lld",r,hist[r]);printf("\n");
}
