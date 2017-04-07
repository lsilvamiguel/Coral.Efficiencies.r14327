// $Id: CsRPDetector.h,v 1.6 2010/01/24 16:10:41 suhl Exp $

/*!
   \file   CsRPDetector.h
   \brief   Compass Recoil Proton detector Class.
   \author  S. Gerassimov
   \version $Revision: 1.6 $
   \date    $Date: 2010/01/24 16:10:41 $
*/

#ifndef CsRPDetector_h
#define CsRPDetector_h

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <iostream>
#include <vector>

#include "CsDetector.h"
#include "CsHistograms.h"

class CsZone;

class TString;
class TLorentzVector;
#include "DaqDataDecoding/ChipSADC.h"
#include "DaqDataDecoding/ChipADC.h"
#include "DaqDataDecoding/ChipF1.h"


/*! \class CsRPDetector 
    \brief Compass Recoil proton detector Class.
*/

// copied from straw detector class

class CsRPDetector : public CsDetector {

 public:

   /*! \fn CsRPDetector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate,
    const double vel, const double t0, const double thRes, 
    const double spSli, const double tiSli );
    \brief Constructor for tracking detector types.
    \param row   Detector file raw number
    \param id    Detector identification number
    \param name  Detector name (see Comgeant)
    \param unit  Detector number in station
    \param type  Detector type (see Comgeant)
    \param rdLen Radiation lenght
    \param xsiz  Detector X size (DRS) (mm)
    \param ysiz  Detector Y size (DRS) (mm)
    \param zsiz  Detector Z size (DRS) (mm)
    \param xcm   Detector centre X position (MRS) (mm)
    \param ycm   Detector centre Y position (MRS) (mm)
    \param zcm   Detector centre Z position (MRS) (mm)
    \param rotDRS Rotation matrix of DRS w.r.t MRS
    \param rotWRS Rotation matrix of WRS w.r.t MRS
    \param wirD  1st wire offset (WRS) (mm)
    \param ang   Rotation angle of WRS w.r.t MRS
    \param nWir  Number of wires
    \param wirP  Wires pitch  (mm)
    \param eff   Detector efficiency
    \param bkg   Detector background
    \param tGate Detector gate  (ns)
  */
  CsRPDetector( const int    row,
        const int    id,    const char* name,  const char *TBname,
        const int    unit,  const int    type,
        const double rdLen, const double xsiz,  
        const double ysiz,  const double zsiz,
        const double xcm,   const double ycm,   
        const double zcm,   const CLHEP::HepMatrix rotDRS,
        const CLHEP::HepMatrix rotWRS,
        const double wirD,  const double ang,   
        const int    nWir,  const double wirP, 
        const double eff,   const double bkg,
	const double tGate);

  ~CsRPDetector();

  virtual void AddSubDetector( const int    row,
                               const int    id,    const char* name,
                               const char *TBname,
                               const int    unit,  const int    type,
                               const double rdLen, const double xsiz,
                               const double ysiz,  const double zsiz,
                               const double xcm,   const double ycm,
                               const double zcm,   const CLHEP::HepMatrix rotDRS,
                               const CLHEP::HepMatrix rotWRS,
                               const double wirD,  const double ang,
                               const int    nWir,  const double wirP,
                               const double eff,   const double bkg,
                               const double tGate );

  public:

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( true ); }

  inline bool hasDrift() const { return( false ); }

  //! Returns the TDC resolution (where applicable)
  inline double    getTDCResol() const { return(tRes_); }

  //! Decode the MC data for this detector
  void makeMCDecoding();

  //! Clusterize the digits
  void clusterize();

  /// Decode raw data
  void          DecodeChipDigits                (const CS::Chip::Digits &digits);
  void          DecodeChipDigit                 (const CS::Chip::Digit &digit);

  /// Read calibration data.
  virtual void readCalibration(time_t timePoint);

  /// Book histograms
  void BookHistograms();

  /// Clear
  void Clear( void );

  private:
  TString     ChannelToPMT(const char *RpdPlanes, int chId, int &thisPMT);
  void       InitRPD(void);
  void       ClearRPD(void);
  float      treatADC(const CS::ChipSADC::Digit *ds);
  float       correctEnergyLoss(float pRpd, float sinTheta, float phi, float xVertex, float yVertex, int debug);
  float       correction_energy(float p, int mat, float thickness, int debug);
  float       terpol3(float xx1,float xx2,float xx3,float yy1,float yy2,float yy3,float x);
  void        Search(void);


 private:
  
  double tRes_;                    //!<          TDC resolution
  bool           decodeCard_;      //!<         Decoding set by cards
  bool           ishorizontal_;    //!<         Has horizontal elements
  bool           isup_;            //!<         Is the upper part 
  float  hittMin_, hittMax_;       //!< Cut on Hit Time (TDC ticks)

  std::map<std::string, CsHist1D*>   mH1;                     //! 1D histogram pointers
  std::map<std::string, CsHist2D*>   mH2;                     //! 2D histogram pointers

private:
int fBestProtonTrack;
std::vector <TLorentzVector> fRPDTrack; // RPD Tracks are stored here
double fMomentum, fMom_x, fMom_y, fMom_z, fTheta, fCosTheta, fSinTheta, fPhi;
std::vector <double> fzT; // primary vertex seen by the RPD for every track
std::vector <double> dEA; // energy loss in A for every track
std::vector <double> dEB;// energy loss in B for every track
std::vector <double> zA_vec; //z Position of Ring A hit
std::vector <double> zB_vec; //z Position of Ring B hit
std::vector <bool> fCalibratedEvent; // Helper to see whether beta is between 0 and 1 (for every track)
std::vector <std::pair <int,int> > hits; // pair of A and B hits for every track
bool something; // RPD saw something

float mP;
int iAup, iBup, iAdo, iBdo, nTracks;
int fDataSentAupt[12], fDataSentAdot[12], fDataSentBupt[24], fDataSentBdot[24];
bool fDataSentAupq[12], fDataSentAdoq[12], fDataSentBupq[24], fDataSentBdoq[24];
float fAdcAup[12], fAdcAdo[120], fAdcBup[24], fAdcBdo[24];
float fTdcAup[12][5], fTdcAdo[12][5], fTdcBup[24][5], fTdcBdo[24][5];
float positionOffsetForA[12], positionOffsetForB[24];
float tofOffset[12][24];
float lightSpeedInA[12], lightSpeedInB[24];
float dEBcusp[24], dEAcusp[12], zVertexOffset[12][24];

float zA, tAoff, dtA, tA, zB, tBoff, dtB, tB, beta, dof, fdEA, fdEB;
float rA, rB, zTargetCenter;
float deltaR;
float lightSpeed, dEAcuspMeV, dEBcuspMeV, tofGlobalOffset;

bool initialized;

};

#endif //CsRPDetector_h
// /*
// RPD Helper for Phast
// 
// comments & questions to 
// Etienne Burtin <Etienne.Burtin@cea.fr>
// Johannes Bernhard <Johannes.Bernhard@cern.ch>
// 
// Enjoy!
// 
// 09-02-12 18:50 jbernhar : Changed defintion of phi according to root definition (i.e. -pi < phi <= pi)
// 09-02-18 00:18 jbernhar : Implemented Etienne's corrections (energy loss in the target, calibration constants, vertex correction) 
// 09-02-21 18:30 jbernhar : Introduced multiple RPD Tracks
// 09-02-23 18:50 jbernhar : guarded several variables against division by 0
// 09-02-26 23:53 jbernhar : changed class structure to be according to Phast
// 
// */
// 
// #ifndef Phast_RPD_Helper_h
// #define Phast_RPD_Helper_h
// 
// #include <iostream>
// #include <sstream>
// #include <cmath>
// 
// 
// #include "TLorentzVector.h"
// #include "TH1.h"
// #include "TH2.h"
// /*
// #include "TTree.h"
// #include "TLorentzVector.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "Phast.h"
// #include "PaSetup.h"
// #include "PaEvent.h"
// #include "G3part.h"
// #include "PaAlgo.h"
// #include "PaParticle.h"
// #include "PaTPar.h"
// #include "PaTrack.h"
// #include "PaVertex.h"
// */
// #include <typeinfo>
// 
// 
// //using namespace std;
// 
// class RPD : public CsMiscDetector
// {
// 
// public:
// 
// RPD::RPD(void);
// 
// RPD::~RPD() {}
// 
// void RPD::Clear(void);
// 
// // void RPD::DecodePhast( const PaEvent& e );
// void RPD::DecodePhast( void );
// 
// const std::vector <TLorentzVector> &RPD::vTracks() const {
// return fRPDTrack;
// }
// 
// //void RPD::Search(const PaEvent& e, bool correct_momentum = false, double PV_x = 0., double PV_y = 0., bool correct_overlaps = false);
// //void RPD::Search(const PaEvent& e);
// void RPD::Search(void);
// 
// int fBestProtonTrack;
// std::vector <TLorentzVector> fRPDTrack; // RPD Tracks are stored here
// double fMomentum, fMom_x, fMom_y, fMom_z, fTheta, fCosTheta, fSinTheta, fPhi;
// std::vector <double> fzT; // primary vertex seen by the RPD for every track
// std::vector <double> dEA; // energy loss in A for every track
// std::vector <double> dEB;// energy loss in B for every track
// std::vector <double> zA_vec; //z Position of Ring A hit
// std::vector <double> zB_vec; //z Position of Ring B hit
// std::vector <bool> fCalibratedEvent; // Helper to see whether beta is between 0 and 1 (for every track)
// std::vector <std::pair <int,int> > hits; // pair of A and B hits for every track
// bool something; // RPD saw something
// TLorentzVector fRPD4Vector; // one RPD track
// 
// std::vector <bool> RPD::CheckCalibration() const {
// return fCalibratedEvent;
// } // checks whether tracks have good calibrations
// 
// int RPD::nTrack() const {
// return nTracks;
// } // returns the number of tracks
// 
// std::vector <double> RPD::PV_Z() const {
// return fzT;
// } // returns the z component of the primary vertex seen by the RPD
// 
// std::vector <std::pair <int,int> > RPD::Hits() const {
// return hits;
// } // returns pairs of A and B hits
// 
// bool RPD::HasTracks() const {
// return something;
// } // returns whether the RPD saw something
// 
// int RPD::iBestProtonTrack() const {
// return fBestProtonTrack;
// } // returns a guess for the best track
// 
// std::vector <double> RPD::dE_A() const {
// return dEA; 
// } // energy loss in A for every track
// 
// std::vector <double> RPD::dE_B() const {
// return dEB;
// } // energy loss in B for every track
// 
// std::vector <double> RPD::ZA() const {
// return zA_vec;
// }
// 
// std::vector <double> RPD::ZB() const {
// return zB_vec;
// }
// 
// void RPD::Dump() const {
// //std::cout << "\n--------------------------------------------------------------------------------\n   ::RPD EVENT DUMP:: \n--------------------------------------------------------------------------------\n\nRPD has an event: " << something << std::endl << std:flush;
// std::cout << "\n--------------------------------------------------------------------------------\n   ::RPD EVENT DUMP:: \n--------------------------------------------------------------------------------\n\nRPD has an event: " << something << std::endl;
// if (something) {
// std::cout << "\nTracks found: " << fRPDTrack.size() << "\n\nTrack parameters:\n----------------\n\n" << std::flush;
// for (unsigned int i=0; i<fRPDTrack.size(); i++) {
// 	std::cout << "Track " << i << "	Momentum: " << fRPDTrack[i].Vect().Mag() << "\n      	Px: " << fRPDTrack[i].Px() << "\n      	Py: " << fRPDTrack[i].Py() << "\n      	Pz: " << fRPDTrack[i].Pz() << "\n      	E: " << fRPDTrack[i].E() << "\n      	Pt: " << fRPDTrack[i].Pt() << "\n      	Theta: " << fRPDTrack[i].Theta() << "\n      	Phi: " << fRPDTrack[i].Phi() << "\n      	Beta: " << fRPDTrack[i].Beta() << "\n      	PV_Z: " << fzT[i] << "\n      	Hit A in Element: " << hits[i].first << "\n      	Hit B in 
// 	}
// }
// }
// 
// private:
// 
// TString RPD::ChannelToPMT(const char *RpdPlanes, int chId, int &thisPMT);
// 
// float RPD::correctEnergyLoss(float pRpd, float sinTheta, float phi, float xVertex, float yVertex, int debug);
// 
// float RPD::correction_energy(float p, int mat, float thickness, int debug);
// 
// float RPD::terpol3(float xx1,float xx2,float xx3,float yy1,float yy2,float yy3,float x);
// 
// //float RPD::treatADC(const PaDigit& d);
// float RPD::treatADC(void);
// 
// private:
// 
// float mP;
// int iAup, iBup, iAdo, iBdo, nTracks;
// int fDataSentAupt[12], fDataSentAdot[12], fDataSentBupt[24], fDataSentBdot[24];
// bool fDataSentAupq[12], fDataSentAdoq[12], fDataSentBupq[24], fDataSentBdoq[24];
// float fAdcAup[12], fAdcAdo[120], fAdcBup[24], fAdcBdo[24];
// float fTdcAup[12][5], fTdcAdo[12][5], fTdcBup[24][5], fTdcBdo[24][5];
// float positionOffsetForA[12], positionOffsetForB[24];
// float tofOffset[12][24];
// float lightSpeedInA[12], lightSpeedInB[24];
// float dEBcusp[24], dEAcusp[12], zVertexOffset[12][24];
// 
// float zA, tAoff, dtA, tA, zB, tBoff, dtB, tB, beta, dof, fdEA, fdEB;
// float rA, rB, zTargetCenter;
// float deltaR;
// float lightSpeed, dEAcuspMeV, dEBcuspMeV, tofGlobalOffset;
// 
// bool initialized;
// 
// };
// 
// #endif
