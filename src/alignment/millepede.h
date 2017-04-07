// $Id: millepede.h,v 1.4 2008/06/05 13:29:13 rgazda Exp $
#ifndef millepede_h
#define millepede_h

extern "C" {
  extern int initgl_(int*,int*,int*,int*);
  extern int parsig_(int*,float*);
  extern int constf_(float*,float*);
  extern int initun_(int*,float*);
  extern int zerloc_(float*,float*);
  extern int gener_(float*,float*,float*,float*,float*,float*,float*);
  extern int equloc_(float*,float*,float*,float*);
  extern int fitloc_();
  extern int fitglo_(float*);
  extern int parglo_(float*);
  extern int prtglo_(int*);
  extern float errpar_(int*);
}

// static float iPar; // moved to C_PARSIG in order to compile without warning "unused variable" (jj)
// static float sig; // moved to C_PARSIG & C_INITUN (jj)
// static float fRot[10000]; // defined but not used (jj)
// static int nGlobal, nLocal, nStdDev, printFlag;  // moved to C_INITGL (jj)
// static int lun; // moved to C_INITUN & C_PRTGLO (jj)
// static float dergb[10000],derlc[10000]; // defined but not used (jj)

#define C_PARSIG(A,B) {int iPar=A; float sig=B; parsig_(&iPar,&sig);}
#define C_CONSTF(A,B) {float fRot[10000]=A; float sig=B; constf_(fRot,&sig);}
#define C_INITUN(A,B) {int lun=A; float sig=B; initun_(&lun,&sig);}
#define C_ZERLOC(A,B) {dergb=a; derlc=B; zerloc_(dergb,derlc);}
#define C_PRTGLO(A)   {int lun=A; prtglo_(&lun);}

#define C_INITGL(A,B,C,D) {int nGlobal=A; int nLocal=B; int nStdDev=C; int printFlag=D; initgl_(&nGlobal,&nLocal,&nStdDev,&printFlag);}        

#endif
