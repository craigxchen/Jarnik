// Referee's independent brute-force scan: for every n in [1, N], enumerate ALL lattice points of
// x^2+y^2=n directly (no factorization, no gcd reduction), take representatives x>0,y>=0 (one per
// unit orbit), angles mod pi/2, and compute for each k the least normalized span
//   C_k(n) = n^(1/4) * min_i (theta_{i+k-1} - theta_i)   (cyclic, wrapping by +pi/2 as often as needed).
// Since dividing a cluster by its Gaussian gcd lowers C and lands on a smaller n, min over n<=N of the
// raw C equals the min over primitive clusters with n<=N.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define KMAX 13
typedef long long ll;
static int cmpd(const void*a,const void*b){double x=*(double*)a,y=*(double*)b;return (x>y)-(x<y);}
int main(int argc,char**argv){
  ll N=atoll(argv[1]); ll S=1<<20;
  double thr[3]={0.5,1.0,1.4142135623730951};
  double bestC[KMAX+1]; ll bestn[KMAX+1];
  for(int k=0;k<=KMAX;k++){bestC[k]=1e9;bestn[k]=0;}
  int maxk_thr[3]={0,0,0}; ll first_thr[3][KMAX+1]; memset(first_thr,0,sizeof first_thr);
  ll nseg=(N+S)/S;
  #pragma omp parallel
  {
    double lbC[KMAX+1]; ll lbn[KMAX+1]; for(int k=0;k<=KMAX;k++){lbC[k]=1e9;lbn[k]=0;}
    ll lfirst[3][KMAX+1]; memset(lfirst,0,sizeof lfirst);
    int *cnt=malloc(sizeof(int)*(S+1)); int *off=malloc(sizeof(int)*(S+2));
    double *ang=NULL; size_t angcap=0;
    double tmp[4096];
    #pragma omp for schedule(dynamic,1)
    for(ll sg=0;sg<nseg;sg++){
      ll N0=sg*S+1, N1=N0+S; if(N1>N+1)N1=N+1; if(N0>N) continue;
      memset(cnt,0,sizeof(int)*(S+1));
      ll xmax=(ll)sqrt((double)N1)+2;
      // pass 1: count
      for(ll x=1;x<=xmax;x++){ ll x2=x*x; if(x2>=N1) break;
        ll lo=N0-x2; ll y=0; if(lo>0){ y=(ll)sqrt((double)lo); while(y*y<lo) y++; while(y>0&&(y-1)*(y-1)>=lo) y--; }
        for(;;y++){ ll n=x2+y*y; if(n>=N1) break; cnt[n-N0]++; } }
      ll tot=0; for(ll i=0;i<N1-N0;i++){off[i]=tot; tot+=cnt[i];} off[N1-N0]=tot;
      if((size_t)tot>angcap){angcap=tot*2; ang=realloc(ang,sizeof(double)*angcap);}
      memset(cnt,0,sizeof(int)*(S+1));
      for(ll x=1;x<=xmax;x++){ ll x2=x*x; if(x2>=N1) break;
        ll lo=N0-x2; ll y=0; if(lo>0){ y=(ll)sqrt((double)lo); while(y*y<lo) y++; while(y>0&&(y-1)*(y-1)>=lo) y--; }
        for(;;y++){ ll n=x2+y*y; if(n>=N1) break; ll i=n-N0; ang[off[i]+cnt[i]++]=atan2((double)y,(double)x); } }
      for(ll i=0;i<N1-N0;i++){
        int c=off[i+1]-off[i]; if(c<1) continue;
        ll n=N0+i; double n4=pow((double)n,0.25);
        // each k from 3..KMAX, also k up to beyond c with wrap
        memcpy(tmp,ang+off[i],sizeof(double)*c); qsort(tmp,c,sizeof(double),cmpd);
        for(int k=3;k<=KMAX;k++){
          double best=1e9;
          for(int s=0;s<c;s++){ int e=s+k-1; double th=tmp[e%c]+(e/c)*M_PI_2; double sp=th-tmp[s]; if(sp<best)best=sp; }
          double C=best*n4;
          if(C<lbC[k]){lbC[k]=C;lbn[k]=n;}
          for(int t=0;t<3;t++) if(C<=thr[t] && (lfirst[t][k]==0||n<lfirst[t][k])) lfirst[t][k]=n;
        }
      }
    }
    #pragma omp critical
    { for(int k=3;k<=KMAX;k++){ if(lbC[k]<bestC[k]||(lbC[k]==bestC[k]&&lbn[k]<bestn[k])){bestC[k]=lbC[k];bestn[k]=lbn[k];}
        for(int t=0;t<3;t++) if(lfirst[t][k] && (first_thr[t][k]==0||lfirst[t][k]<first_thr[t][k])) first_thr[t][k]=lfirst[t][k]; } }
    free(cnt);free(off);free(ang);
  }
  printf("N=%lld\n",N);
  for(int k=3;k<=KMAX;k++) printf("k=%2d bestC=%.6f at n=%lld\n",k,bestC[k],bestn[k]);
  for(int t=0;t<3;t++){ printf("C<=%.4f: first n per k:",thr[t]); for(int k=3;k<=KMAX;k++) if(first_thr[t][k]) printf(" k%d@%lld",k,first_thr[t][k]); printf("\n"); }
  return 0;
}
