/* Referee's independent exhaustive census: all 0=x0<x1<...<x_{k-1}<=B, gcd=1.
   Contacts by residue-class rule (both sides of a split), splits canonical = side without label 0.
   Output: count, coverage histogram, best eps among full coverage, number attaining max coverage + one example. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static int k,B; static int x[8]; static int primes[200],np=0;
static double ell[256]; static long long cnt=0, hist[130]; static int Vk;
static int maxcov=-1; static long long nmax=0; static int exmax[8];
static double besteps=2; static int exbest[8]; static double bestw=0;
static double bestw3eps=2; static int exb3[8];
static int gcdi(int a,int b){while(b){int t=a%b;a=b;b=t;}return a<0?-a:a;}
static void process(void){
  int full=(1<<k)-1;
  for(int s=0;s<(1<<k);s++) ell[s]=0;
  for(int pi=0;pi<np;pi++){int p=primes[pi]; long q=p; double lp=log((double)p);
    while(q<=B){ int seen=0; int any=0;
      for(int i=0;i<k;i++){ if(seen>>i&1) continue; int msk=0;
        for(int j=0;j<k;j++) if(((x[j]-x[i])%q+q)%q==0) msk|=1<<j;
        seen|=msk; int c=__builtin_popcount(msk);
        if(c>=2){ any=1; if(c<=k-2){ int S=(msk&1)? (full^msk):msk; ell[S]+=lp; } } }
      if(!any) break; q*=p; } }
  int cov=0; double lo=1e300,hi=0;
  for(int S=0;S<(1<<k);S++){ if(S&1) continue; int c=__builtin_popcount(S); if(c<2||c>k-2) continue;
    if(ell[S]>0) cov++; if(ell[S]<lo) lo=ell[S]; if(ell[S]>hi) hi=ell[S]; }
  hist[cov]++; cnt++;
  if(cov>maxcov){maxcov=cov;nmax=0;for(int i=0;i<k;i++)exmax[i]=x[i];}
  if(cov==maxcov) nmax++;
  if(cov==Vk){ double e=(hi-lo)/(hi+lo), w=(hi+lo)/2;
    if(w>=2 && e<besteps){besteps=e;bestw=w;for(int i=0;i<k;i++)exbest[i]=x[i];}
    if(w>=3 && e<bestw3eps){bestw3eps=e;for(int i=0;i<k;i++)exb3[i]=x[i];} }
}
static void rec(int i,int lo,int g){
  if(i==k){ if(g==1) process(); return; }
  for(int v=lo; v<=B-(k-1-i); v++){ x[i]=v; rec(i+1,v+1,gcdi(g,v)); }
}
int main(int argc,char**argv){k=atoi(argv[1]);B=atoi(argv[2]);
  for(int p=2;p<=B;p++){int ok=1;for(int d=2;d*d<=p;d++) if(p%d==0){ok=0;break;} if(ok) primes[np++]=p;}
  Vk=(1<<(k-1))-k-1; x[0]=0; rec(1,1,0);
  printf("k=%d B=%d configs=%lld V=%d\nhist:",k,B,cnt,Vk);
  for(int c=0;c<=Vk;c++) if(hist[c]) printf(" %d:%lld",c,hist[c]);
  printf("\nmaxcov=%d attained=%lld example:",maxcov,nmax); for(int i=0;i<k;i++) printf(" %d",exmax[i]);
  printf("\nfull=%lld besteps(w*>=2)=%.4f w*=%.3f ex:",hist[Vk],besteps,bestw); for(int i=0;i<k;i++) printf(" %d",exbest[i]);
  printf("\nbesteps(w*>=3)=%.4f ex:",bestw3eps); for(int i=0;i<k;i++) printf(" %d",exb3[i]); printf("\n");
}
