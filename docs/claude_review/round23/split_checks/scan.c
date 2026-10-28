/* Exhaustive scan of split A-points: primitive 5-tuples u (gcd 1) of distinct positive integers with lcm(u)=N<=NMAX,
   H=(max/min)*N, dist=(max-min)/min, gamma=log H/(-log dist) <= GMAX.  For each N, divisors sorted; for each i, the
   divisors in [d_i, d_i(1+N^{-1/GMAX})] (necessary since H>=N); all 5-subsets containing d_i as minimum; keep those with
   gcd 1 and lcm N.  Prints count and min gamma per dyadic range of N, and every tuple with gamma < GPRINT. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef unsigned long long u64;
static u64 gcdu(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
static int cmpu(const void*a,const void*b){u64 x=*(u64*)a,y=*(u64*)b;return x<y?-1:x>y;}
int main(int argc,char**argv){
  u64 NMAX=strtoull(argv[1],0,10); double GMAX=atof(argv[2]), GPRINT=atof(argv[3]);
  uint32_t*spf=calloc(NMAX+1,4);
  for(u64 i=2;i<=NMAX;i++) if(!spf[i]) for(u64 j=i;j<=NMAX;j+=i) if(!spf[j]) spf[j]=i;
  u64 dv[20000]; double minG[64]; long cnt[64]; for(int k=0;k<64;k++){minG[k]=1e9;cnt[k]=0;}
  for(u64 N=2;N<=NMAX;N++){
    int nd=1; dv[0]=1; u64 m=N;
    while(m>1){u64 p=spf[m]; int e=0; while(m%p==0){m/=p;e++;}
      int old=nd; u64 pk=1; for(int t=1;t<=e;t++){pk*=p; for(int q=0;q<old;q++) dv[nd++]=dv[q]*pk;}}
    if(nd<5) continue;
    qsort(dv,nd,8,cmpu);
    double w=pow((double)N,-1.0/GMAX);
    for(int i=0;i+4<nd;i++){
      int j=i+1; while(j<nd && (double)dv[j] <= dv[i]*(1+w)) j++;
      int c=j-i-1; if(c<4) continue;
      /* enumerate 4-subsets of indices i+1..j-1 */
      for(int a=i+1;a<j;a++)for(int b=a+1;b<j;b++)for(int cc=b+1;cc<j;cc++)for(int d=cc+1;d<j;d++){
        u64 t[5]={dv[i],dv[a],dv[b],dv[cc],dv[d]};
        u64 g=t[0]; for(int k=1;k<5;k++) g=gcdu(g,t[k]); if(g!=1) continue;
        u64 L=t[0]; for(int k=1;k<5;k++){ L=L/gcdu(L,t[k])*t[k]; } if(L!=N) continue;
        double dist=(double)(t[4]-t[0])/(double)t[0];
        double H=(double)N*(double)t[4]/(double)t[0];
        double G=log(H)/(-log(dist)); if(G>GMAX) continue;
        int k=(int)floor(log2((double)N)); cnt[k]++; if(G<minG[k]) minG[k]=G;
        if(G<GPRINT) printf("PT N=%llu u=(%llu,%llu,%llu,%llu,%llu) gamma=%.4f\n",N,t[0],t[1],t[2],t[3],t[4],G);
      }
    }
  }
  for(int k=0;k<64;k++) if(cnt[k]) printf("range 2^%d: count(gamma<=%.2f)=%ld  min gamma=%.4f\n",k,GMAX,cnt[k],minG[k]);
  return 0;
}
