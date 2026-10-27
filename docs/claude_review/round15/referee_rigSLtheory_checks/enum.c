/* Referee's own enumerator. Input (stdin):
   r m mode cap
   par_1..par_r           (required parity of each coordinate)
   nv
   then nv lines: v_1..v_r mask   (mask bit s set => <v,c> == s*m/4 mod m is FORBIDDEN)
   ncls
   ncls codes (sorted ascending) : codes of class vectors mod n (base-n encoding, coordinate 0 most significant)
   mode 0: existence only (stop at first label), mode 1: count all + classify.
   Output: total, inside-class, outside-class counts; up to cap outside-class labels printed "L c1..cr". */
#include <stdio.h>
#include <stdlib.h>
int r,m,n,mode; long cap;
int par[8];
int nv; int (*V)[8]; int *mask; int *last;
int *byk[8]; int cnt[8];
long ncls; long long *cls;
int c[8];
long long tot=0,inside=0,outside=0,printed=0;
int found=0;
static int incls(void){ long long code=0; for(int j=0;j<r;j++){ code=code*n+ (( c[j]%n)+n)%n; }
  long lo=0,hi=ncls-1; while(lo<=hi){ long mid=(lo+hi)/2; if(cls[mid]==code) return 1; if(cls[mid]<code) lo=mid+1; else hi=mid-1;} return 0; }
static void dfs(int k){
  if(found && mode==0) return;
  if(k==r){ tot++; if(mode==0){found=1; printf("FOUND"); for(int j=0;j<r;j++)printf(" %d",c[j]); printf("\n"); return;}
    if(incls()) inside++; else { outside++; if(printed<cap){printed++; printf("L"); for(int j=0;j<r;j++)printf(" %d",c[j]); printf("\n");} }
    return; }
  for(int x=par[k]; x<m; x+=2){
    c[k]=x; int ok=1;
    for(int t=0;t<cnt[k] && ok;t++){ int i=byk[k][t]; long long s=0; for(int j=0;j<=k;j++) s+=(long long)V[i][j]*c[j];
      s%=m; if(s<0) s+=m; if(s%n==0){ int u=(int)(s/n); if(mask[i]&(1<<u)) ok=0; } }
    if(ok) dfs(k+1);
    if(found && mode==0) return;
  }
}
int main(){
  if(scanf("%d %d %d %ld",&r,&m,&mode,&cap)!=4) return 1; n=m/4;
  for(int j=0;j<r;j++) scanf("%d",&par[j]);
  scanf("%d",&nv); V=malloc(sizeof(int[8])*nv); mask=malloc(sizeof(int)*nv); last=malloc(sizeof(int)*nv);
  for(int j=0;j<r;j++){byk[j]=malloc(sizeof(int)*nv);cnt[j]=0;}
  for(int i=0;i<nv;i++){ int l=-1; for(int j=0;j<r;j++){scanf("%d",&V[i][j]); if(V[i][j]) l=j;} scanf("%d",&mask[i]); last[i]=l; if(l>=0 && mask[i]) byk[l][cnt[l]++]=i; }
  scanf("%ld",&ncls); cls=malloc(sizeof(long long)*(ncls+1)); for(long i=0;i<ncls;i++) scanf("%lld",&cls[i]);
  dfs(0);
  if(mode==0){ if(!found) printf("NONE\n"); }
  else printf("COUNT tot %lld inside %lld outside %lld\n",tot,inside,outside);
  return 0;
}
