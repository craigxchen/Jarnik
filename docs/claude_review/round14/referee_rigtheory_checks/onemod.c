#include <stdlib.h>
/* enumerate labels l in (Z/m)^r with parity par[j]; V: nv x r (canonical +-); allowed[k*4+s]=1 if collision (v_k,s) affordable.
   genuine collision lattice: gcol[k]=1 if <v_k,L>=0 mod n. cls: bitmap over (Z/n)^r of genuine class (u*c, u unit mod n), index sum c_j n^j.
   counts: out[0]=survivors, out[1]=nonclass survivors, out[2]=nonclass F1 (some collision v with gcol=0), out[3]=nonclass F2. 
   also stores first nonclass survivor in wit. */
static int R,M,N,NV; static int *Vv,*par,*allowed,*last,*gcol; static unsigned char *cls; static long long out[4]; static int l[16]; static int *wit; static int haswit;
static void rec(int j){
  if (j==R){
    out[0]++;
    long long idx=0,p=1; for(int t=0;t<R;t++){idx+= (long long)(l[t]%N)*p; p*=N;}
    if (!cls[idx]){ out[1]++;
      int f1=0; for(int k=0;k<NV;k++){ long long s=0; for(int t=0;t<R;t++) s+=(long long)Vv[k*R+t]*l[t]; s%=M; if(s<0)s+=M;
         if (s%N==0 && !gcol[k]) {f1=1;break;} }
      if (f1) out[2]++; else out[3]++;
      if(!haswit){haswit=1; for(int t=0;t<R;t++) wit[t]=l[t];}
    }
    return;
  }
  for (int x=par[j]; x<M; x+=2){
    l[j]=x; int ok=1;
    for(int k=0;k<NV;k++){ if(last[k]!=j) continue;
      long long s=0; for(int t=0;t<=j;t++) s+=(long long)Vv[k*R+t]*l[t]; s%=M; if(s<0)s+=M;
      if (s%N==0){ int u=(int)(s/N); if(!allowed[k*4+u]){ok=0;break;} } }
    if (ok) rec(j+1);
  }
}
void run(int r,int m,int nv,int *V,int *parity,int *allow,int *gc,unsigned char *cl,long long *res,int *w){
  R=r;M=m;N=m/4;NV=nv;Vv=V;par=parity;allowed=allow;gcol=gc;cls=cl;wit=w;haswit=0;
  last=malloc(sizeof(int)*nv); for(int k=0;k<nv;k++){int L=-1; for(int t=0;t<r;t++) if(V[k*r+t]) L=t; last[k]=L;}
  for(int i=0;i<4;i++) out[i]=0; rec(0); for(int i=0;i<4;i++) res[i]=out[i]; free(last);
}
