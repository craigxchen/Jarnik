// Exhaustive two-class search: B integers 0=y0<y1<...<y_{B-1}=L.
// q(y) = lcm over all bipartitions sigma (incl. trivial) of g_sigma = gcd of same-side differences.
// Equivalently q = prod_p p^{f_p}, f_p = max f with y occupying <=2 classes mod p^f.
// Prints, for each L, the best log q / log L and the tuple.  (Floating logs only for ranking;
// the winners are re-verified exactly in python.)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static u64 gcd(u64 a, u64 b){ while(b){u64 t=a%b;a=b;b=t;} return a; }
int B, Lmax;
long long y[16];
double bestr; long long besty[16]; double bestlq;
static void evalq(int L){
  int nm = 1<<(B-1);
  u128 q=1;
  for(int m=0;m<nm;m++){
    long long r0=-1,r1=-1; u64 g=0;
    for(int i=0;i<B;i++){
      int side = (i==0)?0:((m>>(i-1))&1);
      if(side==0){ if(r0<0) r0=y[i]; else g=gcd(g,(u64)(y[i]-r0)); }
      else { if(r1<0) r1=y[i]; else g=gcd(g,(u64)(y[i]-r1)); }
      if(g==1) break;
    }
    if(g>1){ u64 h=gcd((u64)(q% (u128)g), g); q = q/(u128)h*(u128)g; }
  }
  double lq = log((double)q);
  double r = lq/log((double)L);
  if(r>bestr){ bestr=r; bestlq=lq; for(int i=0;i<B;i++) besty[i]=y[i]; }
}
static void rec(int k, int L){
  if(k==B-1){ y[k]=L; evalq(L); return; }
  for(long long v=y[k-1]+1; v<= L-(B-1-k); v++){ y[k]=v; rec(k+1,L); }
}
int main(int argc,char**argv){
  B=atoi(argv[1]); Lmax=atoi(argv[2]); int Lmin = argc>3?atoi(argv[3]):B;
  y[0]=0;
  double overall=0;
  for(int L=Lmin; L<=Lmax; L++){
    bestr=0; rec(1,L);
    printf("L=%d best_ratio=%.5f logq=%.4f y=",L,bestr,bestlq);
    for(int i=0;i<B;i++) printf("%lld%s",besty[i],i<B-1?",":"\n");
    fflush(stdout);
  }
  return 0;
}
