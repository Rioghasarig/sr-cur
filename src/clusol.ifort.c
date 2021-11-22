#include "clusol.h"

// declarations for fortran function calls
void lusol_mp_lu1fac_(
  int64_t* m,
  int64_t* n,
  int64_t* nelem,
  int64_t* lena,
  int64_t* ap,
  int64_t* aq, 
  int64_t* frank,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* iploc,
  int64_t* iqloc,
  int64_t* ipinv,
  int64_t* iqinv,
  double* w,
  int64_t* inform);

void lusol_mp_lu6sol_(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void lu8rpc_(
  int64_t* mode1,
  int64_t* mode2,
  int64_t* m,
  int64_t* n,
  int64_t* jrep,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  double* lv,
  int64_t* li,
  int64_t* lj,
  int64_t* lenlv,
  int64_t* inform,
  double* diag,
  double* vnorm);

void lusol_mp_lu1or2_(
  int64_t* n, 
  int64_t* numa,
  int64_t* lena,
  double* a, 
  int64_t* inum,
  int64_t* jnum,
  int64_t* lenc,
  int64_t* locc);

void lu6mul_(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr);

void lu8adc_(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform,
  double* diag,
  double* vnorm);

void lu8adr_(
  int64_t* m,
  int64_t* n,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform,
  double* diag);

void lu8dlc_(
  int64_t* m,
  int64_t* n,
  int64_t* jdel,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);


void lu8mod_(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* beta,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform);

void lu8rpr_(
  int64_t* mode1,
  int64_t* mode2,
  int64_t* m,
  int64_t* n,
  int64_t* irep,
  double* v,
  double* w,
  double* wnew,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* lenlv,
  int64_t* li,
  int64_t* lj,
  double* lv,
  int64_t* inform);


void lu9clr_(
  int64_t* m,
  int64_t* n,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr, 
  double* c, 
  double* lmax,
  int64_t* inform); 



// c interface function definitions
void clu1fac(
  int64_t* m,
  int64_t* n,
  int64_t* nelem,
  int64_t* lena,
  int64_t* ap,
  int64_t* aq, 
  int64_t* frank,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* iploc,
  int64_t* iqloc,
  int64_t* ipinv,
  int64_t* iqinv,
  double* w,
  int64_t* inform) {
  lusol_mp_lu1fac_(m,n,nelem,lena,ap, aq, frank,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,iploc,iqloc,ipinv,iqinv,w,inform);
}

void clu6sol(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform) {
  lusol_mp_lu6sol_(mode,m,n,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

void clu8rpc(
  int64_t* mode1,
  int64_t* mode2,
  int64_t* m,
  int64_t* n,
  int64_t* jrep,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  double* lv,
  int64_t* li,
  int64_t* lj,
  int64_t* lenlv, 
  int64_t* inform,
  double* diag,
  double* vnorm) {
  lu8rpc_(mode1,mode2,m,n,jrep,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,lv,li,lj,lenlv,inform,diag,vnorm);
}

void clu1or2(
  int64_t* n, 
  int64_t* numa,
  int64_t* lena,
  double* a, 
  int64_t* inum,
  int64_t* jnum,
  int64_t* lenc,
  int64_t* locc) {
  lusol_mp_lu1or2_(n,numa,lena,a,inum,jnum,lenc,locc);
}

void clu6mul(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr) {
  lu6mul_(mode,m,n,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr);
}

void clu8adc(
  int64_t* mode,
  int64_t* m,
  int64_t* n,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform,
  double* diag,
  double* vnorm) {
  lu8adc_(mode,m,n,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform,diag,vnorm);
}

void clu8adr(
  int64_t* m,
  int64_t* n,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform,
  double* diag) {
  lu8adr_(m,n,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform,diag);
}

void clu8dlc(
  int64_t* m,
  int64_t* n,
  int64_t* jdel,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform) {
  lu8dlc_(m,n,jdel,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}


void clu9diagu(
  int64_t* m,
  int64_t* n, 
  int64_t* lena,
  int64_t* luparm,
  double* diagU,
  double* a, 
  int64_t* indc,
  int64_t* indr,
  int64_t* ip,
  int64_t* iq,
  int64_t* lenr,
  int64_t* locr,
  int64_t* inform) {
  lu9diagu_(m,n,lena,luparm,diagU,a,indc,indr,ip,iq,lenr,locr,inform);
}

void clu9mod(
  int64_t* m,
  int64_t* n,
  double* beta,
  double* v,
  double* w,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* inform) {
  lu9mod_(m,n,beta,v,w,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,inform);
}

void clu8rpr(
  int64_t* mode1,
  int64_t* mode2,
  int64_t* m,
  int64_t* n,
  int64_t* irep,
  double* v,
  double* w,
  double* wnew,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenc,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr,
  int64_t* lenlv,
  int64_t* li,
  int64_t* lj,
  double* lv,
  int64_t* inform) {
  lu8rpr_(mode1,mode2,m,n,irep,v,w,wnew,lena,luparm,parmlu,a,indc,indr,p,q,lenc,lenr,locc,locr,lenlv,li,lj,lv,inform);
}


void clu9clr(
  int64_t* m,
  int64_t* n,
  int64_t* lena,
  int64_t* luparm,
  double* parmlu,
  double* a,
  int64_t* indc,
  int64_t* indr,
  int64_t* p,
  int64_t* q,
  int64_t* lenr,
  int64_t* locc,
  int64_t* locr, 
  double* c,
  double* lmax,
  int64_t* inform) {
  lu9clr_(m, n, lena, luparm, parmlu, a, indc, indr, p, q, lenr, locc, locr, c, lmax, inform); 
}

