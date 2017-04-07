// $Id: TLattice.cc 14080 2015-10-26 14:03:47Z lsilva $

/*!
  \file   TLattice.cc
  \brief  TLattice (==Tracks dictionary for QNewton fit in Traffic) utils
  \author Yann Bedfer
  */

/*
 * TLattice : Lattice constructor
 * get_coords   : Compute hit-coordinates from parameters using
 *               interpolation in the lattice
 * get_point    : Get grid point closest to input parameters
 * get_point_par: Get parameters of grid point closest to input parameters
 */

#include <stdio.h>
#include <string>
#include "CsInit.h"
#include "CsOpt.h"
#include "CsErrLog.h"
#include "CsGeom.h"
#include "CsEvent.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TLattice.h"

#if USE_RFIO
# include <shift.h>
#else
# include <sys/stat.h>
# include <unistd.h>
#endif

using namespace std;

TLattice* TLattice::address = 0;       // init static pointer
int TLattice::IplLast = 0;   // static data member

//  TLattice class constructor.
//  Lattice is loaded from file "TOpt::Dicofit"
TLattice::TLattice(int read) {

  if (address) return;                 // Protection against multiple instances
  address = this;

  //        ********** ARRAY of DETECTORS in DICO **********
  // (Note: All these are not simultaneously in any given dico. Their effective
  // presence is indicated for each dico by a bit pattern sieve, written in the
  // dico's header.)
  const char *dicoDets[] = {
    "SI03X1  ","SI03Y1  ","SI03U1  ","SI03V1  ", //M=0 H=0 
    "SI04X1  ","SI04Y1  ","SI04U1  ","SI04V1  ", //      4
    "SI05X1  ","SI05Y1  ","SI05U1  ","SI05V1  ", //      8
    "FI03X1  ","FI03Y1  ","FI03U1  ",            //  3
    "FI04X1  ","FI04Y1  ","FI04U1  ",            //  6
    "MM01X1  ","MM01Y1  ","MM01U1  ","MM01V1  ", // 10  12
    "MM02X1  ","MM02Y1  ","MM02U1  ","MM02V1  ", // 14  16
    "MM03X1  ","MM03Y1  ","MM03U1  ","MM03V1  ", // 18  20
    "DC00X1  ","DC00X2  ","DC00Y1  ","DC00Y2  ","DC00U1  ","DC00U2  ","DC00V1  ","DC00V2  ", // 26  28
    "DC01X1  ","DC01X2  ","DC01Y1  ","DC01Y2  ","DC01U1  ","DC01U2  ","DC01V1  ","DC01V2  ", // 34  36

    "GM01X1  ","GM01Y1  ","GM01U1  ","GM01V1  ", // 38  40
    "GM02X1  ","GM02Y1  ","GM02U1  ","GM02V1  ", // 42  44
    "GM03X1  ","GM03Y1  ","GM03U1  ","GM03V1  ", // 46  48
    "GM04X1  ","GM04Y1  ","GM04U1  ","GM04V1  ", // 50  52
    "GM05X1  ","GM05Y1  ","GM05U1  ","GM05V1  ", // 54  56
    "GM06X1  ","GM06Y1  ","GM06U1  ","GM06V1  ", // 58  60
    "DC02X1  ","DC02X2  ","DC02Y1  ","DC02Y2  ","DC02U1  ","DC02U2  ","DC02V1  ","DC02V2  ",
    "DC03X1  ","DC03X2  ","DC03Y1  ","DC03Y2  ","DC03U1  ","DC03U2  ","DC03V1  ","DC03V2  ",
    "DC04X1  ","DC04X2  ","DC04Y1  ","DC04Y2  ","DC04U1  ","DC04U2  ","DC04V1  ","DC04V2  ", // 66  68
    "DC05X1  ","DC05X2  ","DC05Y1  ","DC05Y2  ","DC05U1  ","DC05U2  ","DC05V1  ","DC05V2  ",
    "ST02X1ub","ST02X1db","ST02Y1ub","ST02Y1db","ST02U1ub","ST02U1db", // 72  74
    "ST02X2ub","ST02X2db","ST02Y2ub","ST02Y2db","ST02V1ub","ST02V1db", // 78  80
    "ST03X1ub","ST03X1db","ST03Y1ub","ST03Y1db","ST03U1ub","ST03U1db", // 84  86
    "ST03X2ub","ST03X2db","ST03Y2ub","ST03Y2db","ST03V1ub","ST03V1db", // 90  92
    "FI05X1  ","FI05Y1  ",                       // 92  94
    "FI55U1  ","FI55V1  ",                       // 94  94
    "FI06X1  ","FI06Y1  ","FI06V1  ",            // 97  97 (in DVCS)
    "PS01X1  ","PS01Y1  ","PS01U1  ","PS01V1  ", //101 101 
    "PA01X1  ","PA01U1  ","PA01V1  ",            //104 104
    "PA02X1  ","PA02U1  ","PA02V1  ",            //107 107

    "GM07X1  ","GM07Y1  ","GM07U1  ","GM07V1  ", //111 111
    "GM08X1  ","GM08Y1  ","GM08U1  ","GM08V1  ", //115 115
    "GM09X1  ","GM09Y1  ","GM09U1  ","GM09V1  ", //119 119
    "PA03X1  ","PA03U1  ","PA03V1  ",            //122 122
    "PA04X1  ","PA04U1  ","PA04V1  ",            //125 125
    "PA05X1  ","PA05U1  ","PA05V1  ",            //128 128
    "PA06X1  ",  // Only X!                        129 128
    "FI07X1  ","FI07Y1  ",                       //131 128
    "FI08X1  ","FI08Y1  ",                       //133 130

    // Subsidiary list
    "ST03X1u ","ST03X1d ","ST03Y1u ","ST03Y1d ","ST03U1u ","ST03U1d ",
    "ST03X2u ","ST03X2d ","ST03Y2u ","ST03Y2d ","ST03V1u ","ST03V1d ",

    // Note: For GPs, ordering matters: "X/U" first, "P" (for which one must
    // care not to increment the dico index) 2nd and "Y/V" last, so that the
    // 1st coordinate of "P" coincides w/ "X/U", and the 2nd one comes right
    // after (w/ the caveat that, for the 2nd coordinate of a "UV" pixel, one
    // has to take the opposite of "V").
    //  For MPs, "P" comes 2nd, and provides a slot for the 2nd coordinate.
    "GP01X1  ","GP01P2  ","GP01Y1  ",            //    132
    "GP02U1  ","GP02P1  ","GP02V1  ",            //    134
    "GP02X1  ","GP02P2  ","GP02Y1  ",            //    136
    "GP03U1  ","GP03P1  ","GP03V1  ",            //    138
    "GP03X1  ","GP03P2  ","GP03Y1  ",            //    140
    "MP00X1  ","MP00M1  "                        //    141
  };
  int nDicoDets = sizeof(dicoDets)/sizeof(char*);

  //          *************** GET DICO BOUNDS ***************

  //  Either from input file "dicofit.dat", when creating a new dico. Init the sieve.
  //  Or     from dico's header, when reading in an existing dico. In the latter
  // case, read also the bit pattern sieve.

  FILE *gridFile; int i, size = 0;
  if (!read) {   // ********** CASE of NEW DICO BEING CREATED **********

    char gridFileN[] = "dicofit.dat";
    if ((gridFile = fopen(gridFileN,(char*)"r"))==NULL) {
      printf("** TLattice::TLattice:\a Error Opening Input Grid File \"%s\"!\n",
	     gridFileN); fflush(stdout); perror((char*)" <=");
      exit(1);
    }
    parbound *bp;
    for (i = 0, bp = bound; i<5; i++, bp++) {
      char line[100];
      fscanf(gridFile,"%d%f%f",&bp->num,&bp->min,&bp->max);
      if (i<2) { bp->min /=10; bp->max /=10; } // mm -> cm
      fgets(line,100,gridFile);	
    }
    fclose(gridFile);
    for (i = 0; i<6; i++) sieve[i] = 0;
  }
  else {   // ********** CASE of EXISTING DICO BEING READ **********


    max_par[0] = 1000.; max_par[1] = 1000.; max_par[2] = 15.;
    max_par[3] = 10.; max_par[4] = 10.;
    
    //           ***** DETERMINE DICO FILENAME *****
    // - Filename to open derives from the knowledge of:
    //   - path specified in options file: can be
    //     - either file: obvious case.
    //     - or directory: have to guess most appropriate file therein.
    //   + target magnet info: sign of solenoid field, dipole on/off.
    // - First: retrieve target magnet info.
    CsMagInfo *mag = CsGeom::Instance()->getCsField()->getMagInfo();
    // - Polarity: If options file name already contains reference to a target
    //  setting (extension "plus/minus/transv"): keep it as is and just emit a
    //  non-fatal error if it diagrees w/ target magnet info.
    // - Otherwise search for filename w/ appropriate extension.
    // - Note: the derivation is very similar to that done in "CsInit" ctor for
    //  the best appropriate detectors.dat. A priori, one would like that the
    //  dico file here determined corresponds to that detectors.dat (i.e. was
    //  compiled  based on it). This is not a strong requirement though: dicos
    //  are expressed in detector coordinates corrected for offsets, they depend
    //  on alignment only as far as angle, picth and Z abscissae are concerned,
    //  the latter parameter being stored in the dico file and X-checked against
    //  coral's abscissae (themselves read from det.dat).
    //  => Therefore, let's not be too finicky...

    int runNum = // Run# is actual one, at variance w/ what's done in CsInit...
      // ...ctor, where it's the one extracted from the input file name.
      CsEvent::Instance()->getRunNumber();
    int solSign, dipole, smSign;
    if (mag[0].flg1>=3 /* I.e. SMC or OD dipole */) {
      dipole = mag[0].fsc ? 1 : 0; solSign = 0;
    }
    else {
      dipole = 0; solSign = mag[0].fsc ? (mag[0].fsc>0 ? +1 : -1) : 0;
    }
    smSign =  mag[1].fsc ? (mag[1].fsc>0 ? +1 : -1) : 0;

    string GetDicoFile(string path, int run, int jobType,
		       int solSign = 0, int dipole = 0, int smSign = 0);// Prototype
    string Dicofit; CsInit *init = CsInit::Instance();
    if      (init->IsAHadronJob())
      // Hadron job: ignore target setting, in line w/ what's done in CsInit.
      Dicofit = GetDicoFile(TOpt::Dicofit,runNum,1);
    else if (init->IsABCSJob())
      Dicofit = GetDicoFile(TOpt::Dicofit,runNum,2,solSign,dipole,smSign);
    else
      Dicofit = GetDicoFile(TOpt::Dicofit,runNum,0,solSign,dipole);
    
    if (!init->IsAMonteCarloJob()) {             // ***** CHECK FILE NAME SYNTAX
      // MC's setup files (det.dat and dico) do not usually follow the same
      // syntactic rules as RD's: they need not to. => Bypass them.
      if      (solSign==+1 && (int)Dicofit.find(".plus")<0 ||
	       solSign==-1 && (int)Dicofit.find(".minus")<0)
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "File retrieved owing to option \"TraF Dicofit\" = \"%s\" while "
	  "expected polarity is \"%s\"!",
		      Dicofit.c_str(),solSign>0?"plus":"minus");
      else if (dipole      && (int)Dicofit.find(".transv")<0)
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "File retrieved owing to option \"TraF Dicofit\" = \"%s\" while "
	  "file name w/ \"transv\" extension is expected!",Dicofit.c_str());
    }

    //           ***** OPEN DICO FILE *****
    if ((gridFile = fopen((char*)Dicofit.c_str(),(char*)"r"))==NULL) {
      cerr<<"** TLattice::TLattice:\a Error Opening Input Dico File \""
	  <<Dicofit<<"\"!\n";
      perror((char*)" <="); exit(1);
    }

    //         ***** READ IN BIT PATTERN SIEVE *****
    fread(sieve,sizeof(sieve),1,gridFile);

    //       ***** READ IN BOUNDS => LATTICE SIZE *****
    fread(bound,sizeof(bound),1,gridFile);
    for (i = 4, size = 1; i>=0; i--) { jump[i] = size; size *= bound[i].num;}

    // ***** COMPUTE THE MULTI-JUMP TABLE FOR DIRECT JUMPING *****
    for (i = 0; i<32; i++) {
      int jp, j;
      for (jp=j = 0; j<5; j++) if (i & (1<<j)) jp += jump[j];
      multi_jump[i] = jp;
    }

    //   ***** COMPUTE THE INVERSE LATTICE SPACING'S *****
    il_diff.x    = (bound[0].num-1)/(bound[0].max-bound[0].min); 
    il_diff.y    = (bound[1].num-1)/(bound[1].max-bound[1].min); 
    il_diff.pxzi = (bound[2].num-1)/(bound[2].max-bound[2].min); 
    il_diff.tx   = (bound[3].num-1)/(bound[3].max-bound[3].min); 
    il_diff.ty   = (bound[4].num-1)/(bound[4].max-bound[4].min); 

    //                     ***** ALLOCATE *****
    if ((grid = (Hits *)malloc(size*sizeof(Hits)))==NULL) {
      cerr<<"** TLattice::TLattice:\a Cannot Allocate "
	  <<size*sizeof(Hits)<<" Bytes in Memory!\n";
      exit(1);
    }
  }

  //          *************** DICO INDEX ***************

  // It's <0 for all dets not belonging to the dico.
  // Otherwise it's index in the array "dicoDets" of all detectors liable to
  // enter a dico, sieved by:
  //  - "detectors.dat"      in the case of a new dico being created,
  //  - the bit patter sieve in the case of an existing dico being read.

  const TSetup& setup = TSetup::Ref();
  for (int ipl = 0; ipl<(int)setup.vDetect().size(); ipl++)  // ***** INITIALISE
    dico_idx.push_back(-1);

  bool muonSetup =  // Pinned down approximatively by requiring !=0 target field
    setup.MagScale[0];
  int idx, jdet; char TBprv[8];
  for (jdet = 0, idx = -1, TBprv[0] = '\0'; jdet<nDicoDets;
       jdet++) {
    const char *dicoDet = dicoDets[jdet]; char coord = dicoDet[4];
    // Special dico entries:
    //  I) Pixel (as opposed to stripped) piece of a CsPixelGEM: Shares the same
    // dico slot as its "X" or "U" counterpart, which must precede it in the
    // "dicoDets" array, and gets its second coordinate from the next dico slot,
    // which must be "Y" or "V" in "dicoDets".
    bool isGPP = strncmp(dicoDet,"GP",2)==0 && coord=='P';
    // II) M-pixel piece of CsPixelMM: needs only one dico slot, which we
    // allocate as supra, viz.: "M" shares the same dico slot as its stripped
    // counterpart. Earlier CsPMs had square P-pixels, providing precise 2D
    // info: they are not considered here.
    bool isMPM = strncmp(dicoDet,"MP",2)==0 && coord=='M';
    // III) MMs: can be replaced by corresponding MP
    bool isMM = strncmp(dicoDet,"MM",2)==0;
    if (read) {                              // ***** CASE READING EXISTING DICO
      if (!(sieve[jdet/32]&1<<jdet%32)) continue;
      // The index is determined from the bit pattern sieve read from the dico's
      // header, independent on whether the detector actually is present in the
      // detector table or not. So that the dico can be used on a partial setup.
      if (!isGPP && !isMPM) idx++;
    }
    for (int ipl = 0; ipl<(int)setup.vDetect().size(); ipl++) {
      const TPlane &p = setup.vPlane()[ipl];
      const TDetect &d = setup.vDetect()[ipl];
      const char *TBname = d.Name.c_str();
      bool isDicoDet = strcmp(TBname,dicoDet)==0;
      bool dicoOuterST = false; if (!isDicoDet) { // Not in "dicoDets", yet...
	// ...may have to be in the dico. It's the case of outer ("a/c" slices)
	// STs. They may be added to the dico of the muon setup, if space left
	// in dico allows.
	dicoOuterST = !strncmp(TBname,dicoDet,7) &&// 1st 7 chars as "a/c" slice
	  dicoDet[7]==' ' && (TBname[7]=='a' || TBname[7]=='c') && muonSetup;
      }
      // Several other special cases:
      // I) SI03: can be in spectro (2004), hence in "dicoDets", or out (2006)
      isDicoDet &= d.X(0)>setup.TargetCenter[0];
      // II) PA06: allowed in muon setup
      isDicoDet &= strcmp(TBname,"PA06X1  ") || muonSetup;
      // III) MP stripped piece: it fills the same slot as corresponding MM
      isDicoDet |= isMM && strncmp(TBname,"MP",2)==0 &&
	strcmp(TBname+2,dicoDet+2)==0;
      if (isDicoDet || dicoOuterST) {
	if (!read) {                             // ***** CASE CREATING NEW DICO
	  // Index is the rank in the sub-set of "dicoDets" present in the
	  // detector table. Several detector planes can share a common index:
	  // - OuterST "a/c" slices: skip 2nd encountered (which follows suit).
	  // - GP/MP "P" and either "X" or "U"
	  if ((!dicoOuterST || strncmp(TBname,TBprv,7)) && !isGPP && !isMPM)
	    idx++;
	  strncpy(TBprv,TBname,7); sieve[jdet/32] |= 1<<jdet%32;
	}
	dico_idx[ipl] = idx;
	if (ipl>IplLast) IplLast = ipl;
	if (TOpt::Print[7]) printf("%s: %d\n",TBname,idx);
      }
    }
    // Forget about last detectors in "dicoDets", if dico is already full (
    // they should, in principle, fall in the ``subsidiary list'').
    if (idx>=(NROWS-1)) break;
  }
  for (int ipl = (int)setup.vDetect().size()-1; ipl>=0; ipl--) {
    if (dico_idx[ipl]>=0) break;
    dico_idx[ipl] = -2;  // Fill last part of "dico_idx" w/ "-2"
  }
  IplLast++;  // "IplLast" is in fact last det in dico plus one

  if (read) {

    // ********** CASE of EXISTING DICO BEING READ: EXTRA STEPS **********

    //              ***** CHECK GEOMETRY *****
    Hits dico_abscs;   // As many abscissa as there are hits
    if (fread(dico_abscs,sizeof(Hits),1,gridFile)!=1) {
      CsErrLog::mes(elFatal,"Inconsistent Lattice file!\n");
    }
    double xpl; int ipl; for (ipl = 0, xpl = setup.TargetCenter[0];
			      ipl<(int)setup.vDetect().size(); ipl++) {
      // loop over planes
      double xprv = xpl; xpl = setup.vDetect()[ipl].X(0);
      if (xpl==xprv) continue;         // ...consecutive det's @ same abscissa
      int idx = dico_idx[ipl];
      if (idx==-1) continue;   // ...outside dico
      if (idx==-2) break;      // ...beyond dico
      if (fabs(dico_abscs[idx]-xpl)>.01) { // .01 cm x 200 mrd = 20um
	// Non fatal, to allow looking @ Dico response on a wrong geometry
	CsErrLog::msg(elError,__FILE__,__LINE__,
	  "Geometry: Dico's (=%f) != current (=%f) @ index == %d %s",
	  dico_abscs[idx],xpl,idx,setup.vDetect()[ipl].Name.c_str());
      }
    }

    //                        ***** READ *****
    Hits *grid_p; 
    printf("Reading dico file..."); fflush(stdout);
    for (i = 0, grid_p = grid; i<size; i++)
      if (fread(grid_p++,sizeof(Hits),1,gridFile)!=1) {
	CsErrLog::mes(elFatal,"Inconsistent Lattice file!\n");
      }
    fclose (gridFile);
    printf("\n");
  }

  // ***** DERIVE STARTING Z from SETUP
  if (setup.MagScale[0]==0)
    dZStart = 0;
  else
    // If the target field is active, we avoid incorporating it into the dico, 
    // cf. "TLattice::GenInit". We still set the origin at the center of the
    // target, we move in straight line from this origin to a plane situated
    // "dZStart" away from it and only then start extrapolating in the magnetic
    // field.
    dZStart = 70;

}

//! Accessor to the object
TLattice *TLattice::Ptr()
{
  return address;
}

unsigned int TLattice:: getSize()
{
  unsigned int i, size;
  for (i = 0, size = 1; i<5; i++)
    size *= bound[i].num;
  return size;
}

void TLattice::get_coords(Track *tp)
{
  int    i, grid_point;
  Hits  *grid_p, hits;
  float *pp;     /* pointer to parameter section in structure *tp */

  pp = &(tp->par.x);

  if (!grid) {
    printf("** get_coords_from_grid:\a No lattice loaded !\n"); exit(1);
  }

  /* now there is a grid pointed to by `grid` and described by 'bound' */
   
  /* for the actual parameters, find the nearest gridpoint */
  /* assure minimum distance to boundary of valid phase-space of 1 */

  grid_point = global_grid_point;

  grid_p = grid + grid_point;

  /* CHECK PARAMETERS AGAINST ALLOWED VALUES */

  for (i = 0; i<5; i++) if (fabs(pp[i])>max_par[i]) return;

  for (i = 0; i<5; i++) {
    float x;
    int idx;

    x=(pp[i] - bound[i].min)/(bound[i].max-bound[i].min)*(bound[i].num -1);

    idx = grid_point / jump[i];
    grid_point      %= jump[i];
 
    x -= idx; al_v[i] = 1.-x; ar_v[i] = x;
  }

  //********** INTERPOLATE **********

#define INTERPOLATE1 interpolate_lin0

  INTERPOLATE1 (hits,grid_p);

  for(i = 0; i<NROWS; i++) 
    if (fabs(hits[i])>1e9) hits[i] = 1e9;  /* avoid overflows */

  for (i = 0; i<NROWS; i++) tp->hits.h[i] = hits[i];
}


void TLattice::INTERPOLATE1(Hits ret_p,Hits* grid_p)
{
#undef DIM
#undef INTERPOLATE1

#define DIM 0
#define INTERPOLATE1 interpolate_lin1

  Hits ql, qr;
  int i;

  INTERPOLATE1(ql,grid_p          );
  INTERPOLATE1(qr,grid_p+jump[DIM]);

  for (i = 0; i<NROWS; i++)
    ret_p[i] = al_v[DIM]*ql[i] + ar_v[DIM]*qr[i];
}

void TLattice::INTERPOLATE1(Hits ret_p,Hits* grid_p)
{
#undef DIM
#undef INTERPOLATE1

#define DIM 1
#define INTERPOLATE1 interpolate_lin2

  Hits ql, qr;
  int i;

  INTERPOLATE1(ql,grid_p          );
  INTERPOLATE1(qr,grid_p+jump[DIM]);

  for (i = 0; i<NROWS; i++)
    ret_p[i] = al_v[DIM]*ql[i] + ar_v[DIM]*qr[i];
}

void TLattice::INTERPOLATE1(Hits ret_p,Hits* grid_p)
{
#undef DIM
#undef INTERPOLATE1

#define DIM 2
#define INTERPOLATE1 interpolate_lin3

  Hits ql, qr;
  int i;

  INTERPOLATE1(ql,grid_p          );
  INTERPOLATE1(qr,grid_p+jump[DIM]);

  for (i=0; i<NROWS; i++)
     ret_p[i] = al_v[DIM]*ql[i] + ar_v[DIM]*qr[i];
}

void TLattice::INTERPOLATE1(Hits ret_p,Hits* grid_p)
{
#undef DIM
#undef INTERPOLATE1

#define DIM 3
#define INTERPOLATE1 interpolate_lin4

  static Hits ql, qr;
  int i;

  INTERPOLATE1(ql,grid_p          );
  INTERPOLATE1(qr,grid_p+jump[DIM]);

  for (i=0; i<NROWS; i++)
    ret_p[i] = al_v[DIM]*ql[i] + ar_v[DIM]*qr[i];
}

/* this one here is optimized for speed, cause
   thats where ALL the cpu-time goes !!! */

void TLattice::INTERPOLATE1(Hits ret_p,Hits* grid_p)
{
#undef DIM
#define DIM 4

  register float *flp,*frp,*rp;

  rp = &ret_p[0];
  flp = &(*grid_p)[0];
  frp = &(*(grid_p+1))[0];

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);

  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp++);
  *rp++ = al_v[DIM]*(*flp++) + ar_v[DIM]*(*frp  );
}

int TLattice::get_point(Parameters *p)
{
  int   grid_point = 0, idx, i;
  float *pp;   /* pointer to parameter section in structure *tp */

  pp = &(p->x);

  for(i = 0; i<5; i++) { 
    float x = 
      (pp[i]-bound[i].min)/(bound[i].max-bound[i].min)*(bound[i].num-1);

    idx = (int)x;
    if (idx<0             ) idx = 0;
    if (idx>bound[i].num-2) idx = bound[i].num-2;

    grid_point += idx*jump[i];

  }
  return grid_point;
}

int TLattice::get_point_par(Parameters *p,Parameters *q)
{
  int grid_point = 0, idx, i;
  float *pp = &(p->x);
  float *qq = &(q->x);

  for(i = 0; i<5;i ++) {
    float x = 
      (pp[i]-bound[i].min)/(bound[i].max-bound[i].min)*(bound[i].num -1);

    idx = (int)x;
    if (idx<0             ) idx = 0;
    if (idx>bound[i].num-2) idx = bound[i].num-2;
    grid_point += idx*jump[i];

    qq[i] = x - idx;
  }
  return grid_point;
}

//  **************************************************************************
//  *****  FIND MOST RELEVANT DICO FILE FOR GIVEN RUN and JOB TYPE  ******
// - Among all detectors.dat in given path, w/ <run#> and <type> extensions.
// - MOST RELEVANT is CLOSEST EARLIER <run#> 
//                 w/ MAGNETS SETTING RELEVANT for <type> 
// - JOB TYPE is:
//   0: spin asymmetry:   relevant depends on target (solenoid,dipole) field
//                   - (+/-1,0) <-> <type> = "plus" or "minus" 
//                   - (  0 ,1) <-> <type> = "transv"
//   1: hadron                      <type> = "hadron"
//   2: charge asymmetry: relevant depends on SM1/2 polarity
//                   - +/-1     <-> <type> = "mu+" or "mu-"
//  **********************************************************************
#include <dirent.h>
string GetDicoFile(string path, // Path to either file or directory
		   int run,
		   int jobType,
		   int solSign, int dipole, int smSign)
{
  cout<<endl;
  //                                          ***** IS "path" DIRECTORY OR FILE?
  if (path[path.length()-1]!='/') { // Path not ending w/ "/"
    cout<<"GetDicoFile ==> \""<<path<<"\" will be used\n";
    return path;                           // ...=> FILE RETURN INPUT "path" ARG
  }
  else {                                   // ...ELSE DIRECTORY
    if (run<0) {
      cout<<"GetDicoFile ==> Run number can't be extracted from data.\n";
      cout<<"Please specify 'dico' explicitly in \"TraF Dicofit\" option.\n";
      exit(1);
    }
  }
  
  printf  ("GetDicoFile ==> Run# to be processed "
	   "(as extracted from input data file) is %d\n\n",run);
  if      (jobType==0) {
    if (solSign)
      printf("GetDicoFile ==> This run is assumed to have solenoid field: \"%c\"\n",
	     solSign==+1 ? '+' : '-');
    else if (dipole)
      printf("GetDicoFile ==> This run is assumed to have target dipole ON\n");
    else
      printf("GetDicoFile ==> This run is assumed to have target field OFF\n");
  }
  else if (jobType==2)
    printf("GetDetectorsDat ==> This run is assumed to have SM1/2 polarity %c0\n",
	   smSign>0?'>':'<');

  string targetedType;
  if      (jobType==0) {
    if      (solSign==+1) targetedType = "plus";
    else if (solSign==-1) targetedType = "minus";
    else if (dipole)      targetedType = "transv";
    else                  targetedType = "zero";
  }
  else if (jobType==1)    targetedType = "hadron";
  else if (jobType==2) {
    if      (smSign==+1)  targetedType = "mu+";
    else if (smSign==-1)  targetedType = "mu-";
  }

  DIR *dp;                                            // ***** OPENING DIRECTORY
  if (!(dp = opendir(const_cast<char *>(path.c_str()))))
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
		  "Can't open directory \"%s\"",path.c_str());
  struct dirent *ep; map<int,string> mRun2File;      // ***** PARSING FILE NAMES
  while ((ep = readdir(dp))) { // Loop over file names
    string fnam(ep->d_name); 
    //                           ***** REQUIRE "dico[.<anything>].<run#>.<type>"
    if (fnam.find("dico.")!=0) continue;
    int dot1 = fnam.rfind('.');
    int dot2 = fnam.rfind('.', dot1-1); if (dot2<0) continue;  

    string type(fnam,dot1+1,fnam.length());  // ***** FILTER TYPE w/ TARGET INFO
    if (type!=targetedType) continue;

    string runNumField(fnam,dot2+1,dot1-dot2-1);           // ***** EXTRACT RUN#
    char *end, **endptr = &end;     // Check run# field houses a numerical value
    int num = strtol(runNumField.c_str(),endptr,0);
    if (**endptr!='\0') continue;

    mRun2File[num] = fnam;                       // ***** MAP RUN# <-> FILE NAME
  } // End of loop over file names
  closedir (dp);
  if (mRun2File.empty())
    CsErrLog::msg(elFatal,__FILE__,__LINE__,
      "No files w/ syntax \"dico[.<anything>].<run#>.%s\""
		  "in directory \"%s\"",targetedType.c_str(),path.c_str());

  // Files w/ well formed syntax are expected to be ordered by increasing run#. 
  int bestRun; map<int,string>::reverse_iterator irun;
  for (irun = mRun2File.rbegin(), bestRun = -1; irun!=mRun2File.rend();
       irun++) {
    if ((*irun).first>run) continue;
    bestRun = (*irun).first; break;
  }

  if (bestRun!=-1)
    printf("GetDicoFile ==> \"%s%s\" will be used to process it.\n",
	   path.c_str(),mRun2File[bestRun].c_str());
  else {
    printf("GetDicoFile ==> Files w/ syntax "
	   "\"dico[.<anything>].<run#>.%s\" in \"%s\" are:\n",
	   targetedType.c_str(),path.c_str());
    for (irun = mRun2File.rbegin(); irun!=mRun2File.rend(); irun++)
      printf("\"%s\" ",(*irun).second.c_str());
    cout<<endl;
    CsErrLog::mes(elFatal,"No appropriate \"dico\" could be found.\n"
      "=> Check your options file. In last resort, "
      "you can always specify an explicit \"TraF Dicofit\" file name.");
  }
  cout<<endl;
  return path+mRun2File[bestRun];
}
