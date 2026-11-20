// Exhaustive search: for every n in [2,N], every q >= n^{amin}, every class r mod q with gcd(r,q)=1,
// the number of positive divisors of n congruent to r mod q.  Records, for thresholds theta,
// the maximum count with log q / log n >= theta, and an example.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
static int *spf;
static int *cnt; static int *touched; static int ntouched;
static long long gcdll(long long a,long long b){while(b){long long t=a%b;a=b;b=t;}return a;}
#define NT 12
double th[NT]={0.20,0.21,0.22,0.23,0.24,0.245,0.25,0.26,0.28,0.30,0.3333333333,0.35};
int best[NT]; long long ex_n[NT], ex_q[NT], ex_r[NT];
int main(int argc,char**argv){
  int N=atoi(argv[1]); double amin=0.2;
  spf=calloc(N+1,sizeof(int)); cnt=calloc(N+1,sizeof(int)); touched=malloc((N+1)*sizeof(int));
  for(int i=2;i<=N;i++) if(!spf[i]) for(long long j=i;j<=N;j+=i) if(!spf[j]) spf[j]=i;
  int divs[20000]; int gd[20000];
  for(int n=2;n<=N;n++){
    // divisors of n
    int nd=1; divs[0]=1; int m=n;
    while(m>1){int p=spf[m],e=0; while(m%p==0){m/=p;e++;}
      int cur=nd; long long pp=1; for(int k=1;k<=e;k++){pp*=p; for(int t=0;t<cur;t++) divs[nd++]=divs[t]*pp;}}
    // sort
    for(int i=1;i<nd;i++){int v=divs[i],j=i-1; while(j>=0&&divs[j]>v){divs[j+1]=divs[j];j--;} divs[j+1]=v;}
    double ln=log((double)n); int qmin=(int)ceil(exp(amin*ln)); if(qmin<2)qmin=2;
    for(int a=0;a<nd;a++){
      int d0=divs[a]; ntouched=0;
      for(int b=a+1;b<nd;b++){
        int g=divs[b]-d0; if(g<qmin) continue;
        // divisors of g
        int ng=1; gd[0]=1; int mm=g;
        while(mm>1){int p=spf[mm],e=0; while(mm%p==0){mm/=p;e++;}
          int cur=ng; long long pp=1; for(int k=1;k<=e;k++){pp*=p; for(int t=0;t<cur;t++) gd[ng++]=gd[t]*pp;}}
        for(int t=0;t<ng;t++){int q=gd[t]; if(q<qmin) continue; if(gcdll(q,d0)!=1) continue;
          if(cnt[q]==0) touched[ntouched++]=q; cnt[q]++;}
      }
      for(int t=0;t<ntouched;t++){int q=touched[t]; int c=cnt[q]+1; cnt[q]=0;
        double al=log((double)q)/ln;
        for(int k=0;k<NT;k++) if(al>=th[k] && c>best[k]){best[k]=c;ex_n[k]=n;ex_q[k]=q;ex_r[k]=d0%q;}}
    }
  }
  for(int k=0;k<NT;k++) printf("theta=%.4f best=%d n=%lld q=%lld r=%lld alpha=%.4f\n",th[k],best[k],ex_n[k],ex_q[k],ex_r[k],
     log((double)ex_q[k])/log((double)ex_n[k]));
  return 0;
}
