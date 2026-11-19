// Exhaustive census of M(D,0) in the ramified order F_3[t,eta], eta^2 = Delta (default 2t^2+1).
// Point z = A + B eta, pole order D = max(deg A, deg B + deg(Delta)/2) (deg Delta even, nonsquare lead => inert inf).
// Cluster: distinct points, one norm N = A^2 - Delta B^2, pole order of each difference <= floor(D/2).
// Relation "difference has pole order <= H" is an equivalence (ultrametric): classes = fixed high coefficients.
// Hence M(D) = max over (class, N) of #points.  Sign flips A->-A, B->-B preserve everything, so the leading pair
// (a_D, b_{D-dd}) can be taken in {(1,0),(0,1),(1,1)} (dd = deg Delta / 2).
// Output: max over all N, and max over squarefree N, with a witness.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cstring>
using namespace std;
typedef uint64_t u64;
static const int MAXC = 21;           // up to degree 20 packed, 3 bits per coefficient
static u64 ONES;
static inline u64 addp(u64 a, u64 b){ u64 x=a+b; u64 m=((x+ONES)>>2)&ONES; return x-3*m; }
static inline u64 negp(u64 x){ u64 nz=(x|(x>>1))&ONES; return 3*nz-x; }
static inline int coef(u64 x,int i){ return (x>>(3*i))&7; }
static u64 pack(const int* c,int n){ u64 r=0; for(int i=0;i<n;i++){ int v=((c[i]%3)+3)%3; r|=(u64)v<<(3*i);} return r; }
static void polymul(const int*a,int na,const int*b,int nb,int*r){ for(int i=0;i<na+nb-1;i++) r[i]=0; for(int i=0;i<na;i++) if(a[i]) for(int j=0;j<nb;j++) r[i+j]=(r[i+j]+a[i]*b[j])%3; }
static int pdeg(const int*a,int n){ for(int i=n-1;i>=0;i--) if(a[i]%3) return i; return -1; }
static bool squarefree(u64 N){
    int a[MAXC+2], b[MAXC+2];
    for(int i=0;i<MAXC;i++) a[i]=coef(N,i);
    int da=pdeg(a,MAXC); if(da<=0) return true;
    for(int i=0;i<da;i++) b[i]=(a[i+1]*(i+1))%3; int db=pdeg(b,da);
    if(db<0) return false;  // N' = 0 -> N is a cube
    // gcd
    while(db>=0){
        // a = a mod b
        int inv = (b[db]==1)?1:2;
        while(da>=db){
            int c=(a[da]*inv)%3;
            if(c) for(int i=0;i<=db;i++) a[da-db+i]=((a[da-db+i]-c*b[i])%3+3)%3;
            da=pdeg(a,da+1);
            if(da<0) break;
        }
        // swap
        int tmp[MAXC+2]; memcpy(tmp,a,sizeof(int)*(MAXC+2)); memcpy(a,b,sizeof(int)*(MAXC+2)); memcpy(b,tmp,sizeof(int)*(MAXC+2));
        int t=da; da=db; db=t;
    }
    return da==0;
}
int main(int argc,char**argv){
    int D=atoi(argv[1]);
    // Delta from argv[2] as coefficient string low->high, default "102" = 1 + 0 t + 2 t^2
    const char* ds = argc>2? argv[2] : "102";
    int nd=strlen(ds); int Del[8]; for(int i=0;i<nd;i++) Del[i]=ds[i]-'0';
    int dD=nd-1; if(dD%2){fprintf(stderr,"odd deg\n");return 1;} int dd=dD/2;
    ONES=0; for(int i=0;i<MAXC;i++) ONES|=(u64)1<<(3*i);
    if(2*D>MAXC-1){fprintf(stderr,"D too big\n");return 1;}
    int H=D/2;
    // A: degrees 0..D ; B: degrees 0..D-dd. Low A: 0..H (H+1 coeffs); low B: degrees e with e+dd<=H: 0..H-dd.
    int nAl=H+1, nBl=H-dd+1; if(nBl<0) nBl=0;
    int nAh=D-H, nBh=(D-dd)-(H-dd+1)+1; // A: H+1..D ; B: max(0,H-dd+1)..D-dd
    int Bh0=H-dd+1; if(Bh0<0){Bh0=0;}
    nBh=(D-dd)-Bh0+1;
    long long nlo=1; for(int i=0;i<nAl+nBl;i++) nlo*=3;
    long long nA=1; for(int i=0;i<nAl;i++) nA*=3; long long nB=1; for(int i=0;i<nBl;i++) nB*=3;
    // Q table: Al^2 - Delta Bl^2
    vector<u64> Q(nlo);
    for(long long a=0;a<nA;a++){
        int Al[24]={0}; long long x=a; for(int i=0;i<nAl;i++){Al[i]=x%3;x/=3;}
        int A2[48]; polymul(Al,nAl,Al,nAl,A2);
        for(long long b=0;b<nB;b++){
            int Bl[24]={0}; long long y=b; for(int i=0;i<nBl;i++){Bl[i]=y%3;y/=3;}
            int r[48]={0};
            for(int i=0;i<2*nAl-1;i++) r[i]=A2[i];
            if(nBl>0){ int B2[48]; polymul(Bl,nBl,Bl,nBl,B2); int DB[64]; polymul(Del,nd,B2,2*nBl-1,DB); for(int i=0;i<nd+2*nBl-2;i++) r[i]=(r[i]-DB[i]+9)%3; }
            Q[a*nB+b]=pack(r,MAXC);
        }
    }
    int nfree = (nAh-1)+(nBh-1);
    long long nh=1; for(int i=0;i<nfree;i++) nh*=3;
    int leadA[3]={1,0,1}, leadB[3]={0,1,1};
    int best=0, bestsf=0; u64 wN=0,wNsf=0; long long hits=0, hitssf=0;
    vector<u64> keys(nlo); vector<u64> fA(nA), fB(nB);
    for(int ty=0;ty<3;ty++) for(long long h=0;h<nh;h++){
        int Ah[24]={0}, Bh[24]={0};   // full polys with only high coeffs
        long long x=h;
        for(int i=H+1;i<D;i++){Ah[i]=x%3;x/=3;}
        for(int i=Bh0;i<D-dd;i++){Bh[i]=x%3;x/=3;}
        Ah[D]=leadA[ty]; Bh[D-dd]=leadB[ty];
        // Nh = Ah^2 - Delta Bh^2
        int r[64]={0}, A2[64], B2[64], DB[64];
        polymul(Ah,D+1,Ah,D+1,A2); polymul(Bh,D-dd+1,Bh,D-dd+1,B2); polymul(Del,nd,B2,2*(D-dd)+1,DB);
        for(int i=0;i<=2*D;i++) r[i]=((A2[i]-DB[i])%3+3)%3;
        u64 Nh=pack(r,MAXC);
        // linear tables: fA[Al] = 2 Ah Al ; fB[Bl] = -2 Delta Bh Bl = Delta Bh Bl
        int DBh[64]; polymul(Del,nd,Bh,D-dd+1,DBh);
        vector<u64> ea(nAl), eb(nBl);
        for(int i=0;i<nAl;i++){ int c[64]={0}; for(int j=0;j<=D;j++) c[i+j]=(2*Ah[j])%3; ea[i]=pack(c,MAXC);}  // 2 Ah t^i
        for(int i=0;i<nBl;i++){ int c[64]={0}; for(int j=0;j<=D+dd;j++) c[i+j]=(DBh[j])%3; eb[i]=pack(c,MAXC);} // +Delta Bh t^i (=-2 Delta Bh t^i)
        // fill via base-3 increments
        fA[0]=0; for(long long a=1;a<nA;a++){ long long y=a; int i=0; while(y%3==0){y/=3;i++;} // lowest nonzero digit position i, digit d
            long long p3=1; for(int k=0;k<i;k++) p3*=3; long long prev=a-p3; fA[a]=addp(fA[prev],ea[i]); }
        fB[0]=0; for(long long b=1;b<nB;b++){ long long y=b; int i=0; while(y%3==0){y/=3;i++;}
            long long p3=1; for(int k=0;k<i;k++) p3*=3; long long prev=b-p3; fB[b]=addp(fB[prev],eb[i]); }
        for(long long a=0;a<nA;a++){ u64 fa=fA[a]; u64* kk=&keys[a*nB]; const u64* qq=&Q[a*nB];
            for(long long b=0;b<nB;b++) kk[b]=addp(addp(qq[b],fa),fB[b]); }
        sort(keys.begin(),keys.end());
        long long i=0;
        while(i<nlo){ long long j=i; while(j<nlo&&keys[j]==keys[i]) j++; int c=j-i;
            if(c>=best){ if(c>best){best=c;hits=0;wN=addp(Nh,keys[i]);} hits++; }
            if(c>=bestsf){ u64 N=addp(Nh,keys[i]); if(squarefree(N)){ if(c>bestsf){bestsf=c;hitssf=0;wNsf=N;} hitssf++; } }
            i=j; }
    }
    auto pr=[&](u64 N){ for(int i=2*D;i>=0;i--) printf("%d",coef(N,i)); };
    printf("D=%d Delta=%s  M_all=%d (classes attaining, up to sign flips: %lld) witnessN=",D,ds,best,hits); pr(wN);
    printf("  M_sf=%d (attaining %lld) witnessN=",bestsf,hitssf); pr(wNsf); printf("\n");
    return 0;
}
