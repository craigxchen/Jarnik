#include <stdlib.h>
#include <string.h>
/* free test: exists c in (Z/n)^r (parity: if n even and par[j]>=0 then c_j = par[j] mod 2) with <v,c> != 0 mod n for all v.
   vectors V (nv x r), grouped by last nonzero coordinate. c1 restricted to reps list. returns 1 if free (c filled), 0 if forced. */
static int R, N, NV; static int *Vv; static int *par; static int *last; static long long nodes;
static int c[16];
static int rec(int j){
  if (j==R) return 1;
  for (int x=0;x<N;x++){
    if (par[j]>=0 && (N%2==0) && (x%2)!=par[j]) continue;
    c[j]=x; int ok=1; nodes++;
    for (int k=0;k<NV;k++){ if (last[k]!=j) continue;
      long long s=0; for(int t=0;t<=j;t++) s+= (long long)Vv[k*R+t]*c[t];
      s%=N; if (s<0) s+=N; if (s==0){ok=0;break;} }
    if (ok && rec(j+1)) return 1;
  }
  return 0;
}
int isfree(int r,int n,int nv,int *V,int *parity,int nreps,int *reps,int *out){
  R=r;N=n;NV=nv;Vv=V;par=parity; last=malloc(sizeof(int)*nv);
  for(int k=0;k<nv;k++){int l=-1; for(int t=0;t<r;t++) if(V[k*r+t]) l=t; last[k]=l;}
  for(int i=0;i<nreps;i++){ c[0]=reps[i]; int ok=1;
    for(int k=0;k<nv;k++){ if(last[k]!=0) continue; long long s=(long long)V[k*r]*c[0]; s%=n; if(s<0)s+=n; if(s==0){ok=0;break;}}
    if(!ok) continue;
    if (rec(1)){ for(int t=0;t<r;t++) out[t]=c[t]; free(last); return 1;} }
  free(last); return 0;
}
