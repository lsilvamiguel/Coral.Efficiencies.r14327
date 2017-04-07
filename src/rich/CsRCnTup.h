// $Id: CsRCnTup.h,v 1.9 2010/06/18 10:44:21 tnagel Exp $

/*!
   \file    CsRCnTup.h
   \brief   Output event info to an HBOOK ntuple.
   \author  Dmitri Pechekhonov, Benigno Gobbo, Alexandre Korzevev
   \date    $Date: 2010/06/18 10:44:21 $
   \version $Revision: 1.9 $
   
   \par History:
   2000/05/19 Added Root Tree <br>
Giuia per root ntuple 070702
*/

//--- recontructed by Paolo 01/10/19.


#ifndef CsRCnTup_h
#define CsRCnTup_h

#include "coral_config.h"

#ifdef COMPASS_USE_OSPACE
# include <ospace/std/string>
#else
# include <string>
#endif

#include "TTree.h" 
#include "TVector3.h"

class CsRCnTup {

 public:

  //! singleton
  static CsRCnTup* Instance();

  void CsRCnTupNm();

  //! fill hbook cw ntuple or root tree
  void fillOutput();

 protected:

  CsRCnTup();          //!< The Protected Singleton Constructor

 private:

  static CsRCnTup* instance_; //!< The singleton static attribute 

  bool block_[5];    //!< true - fills, false - does not fill.


  // This is the HBOOK part

  // arrays size
  static const int nRmax =  30;
  static const int nPmax = 300;

  // ntuple path name
  const char * pathname;
  int ntupNum;

  struct { 
    int align; 
    float Run;
    float Ev;
    float Index;
  } Header;   //!< Event and Run numbers + Refr. index
  
  struct { 
    int align;
    int nTracks;
    float xIn[nRmax]; float yIn[nRmax]; float zIn[nRmax];
    float tgX[nRmax]; float tgY[nRmax];
    float cmom[nRmax];
    float zHx0[nRmax]; float zHx1[nRmax];
  } Track;   //!< Event Tracks

  struct { 
    int align;
    int nRings;
    float theta[nRmax];
    int cPaDe[nRmax]; float xPaDe[nRmax]; float yPaDe[nRmax];
    int nPhoRing[nRmax]; int nPhoBkx[nRmax];
  } Ring;   //!< Event Rings

  struct { 
    int align;
    int nPhotons;
    float the[nPmax]; float phi[nPmax]; float PH[nPmax];
    int cath[nPmax]; float xc[nPmax]; float yc[nPmax];
    int nPads[nPmax];
    float paPH[nPmax];
  } Photon;   //!< Ring Photons

 struct { 
    int align;
    int nLike;
    float likeBkg[nRmax];
    float likePi[nRmax];
    float likeKa[nRmax];
    float likePr[nRmax];
    float dLiDnPi[nRmax];
    float dLiDnKa[nRmax];
    float dLiDnPr[nRmax];
    float theLike[nRmax];
    float theRec[nRmax];
    float nPhoLk[nRmax];
    float theFit[nRmax];
    float chiRing[nRmax];
    float chiPi[nRmax];
    float chiKa[nRmax];
    float chiPr[nRmax];
    float likeEl[nRmax];
    float likeMu[nRmax];
    float dLiDnEl[nRmax];
    float dLiDnMu[nRmax];
    float chiEl[nRmax];
    float chiMu[nRmax];
  } Like;   //!< Ring Likelihoods

 struct { 
    int align;
    int kEvGood;
    int kMuPart;
    float zVertex;
    float phiMass;
  } Phys;   //!< Physics

  //! fill ntuple
  void fillOutNtuple();


  // And this is the ROOT one
  TTree*    MyTree_;
  //! fill tree
  void fillTree();

  void print();

};
#endif //CsRCnTup_h
