#include "CsPixelMumegaDetector.h"

//-----------------------------------------------------------------------------
#include <math.h>
#include <string.h>
#include <strings.h>
#include <functional>
#include "CsZone.h"
#include "CsDigit.h"
#include "CsOpt.h"
#include "CsEvent.h"
#include "CsMCDigit.h"
#include "CsComgNtCommons.h"
#include "CsGeom.h"
#include "CsInit.h"
#include <cstdlib>
#include "CsMCTrkHit.h"
#include "CsRandom.h"
#include "CsGeant3.h"
#include "CsMCTrack.h"
#include "CDB.h"

#include "DaqDataDecoding/ChipAPV.h"

using namespace std;
using namespace CLHEP;
using CS::DetID;

extern QhitType Ghit;

//=============================================================================

CsPixelMumegaDetector::CsPixelMumegaDetector( const int    row,    const int    id,
		const char*  name,   const char *TBname,
		const int    unit,   const int    type,
		const double rdLen,  const double xsiz,
		const double ysiz,   const double zsiz,
		const double xcm,    const double ycm,
		const double zcm,    const HepMatrix rotDRS,
		const HepMatrix rotWRS,
		const double wirD,   const double ang,
		const int    nWir,   const double wirP,
		const double eff,    const double bkg,
		const double tGate,  const double spSig,
		const double eGain,  const double eGSig,
		const double sWidth, const double tRes):
		CsDetector(row, id, name, TBname, unit, type, rdLen, xsiz, ysiz, zsiz,
				xcm, ycm, zcm, rotDRS, rotWRS, wirD, ang,
				nWir, wirP, eff, bkg, tGate),
				isMaster_(false), associateDet_(NULL),
				spSig_(spSig),
				eGain_(eGain),
				eGSig_(eGSig),
				sWidth_(sWidth),
				tRes_(tRes),
				stripplane(NULL), pixelplane(NULL), rectpixelplane(NULL),
                                pixmmobject(NULL),
				apv_chan_cals_mapped(false), compl_det_searched(false)
{


	nWirV_ = 0;
	if(tRes_ == 0.0){
		tRes_ = 12;
	}

	this->sparse = true;
	this->do_clustering = 1; // 0 : primitive clustering, 1 : Full clustering
	// may be superseded by option

	string TB = GetTBName().substr(0,2);


	// === Get the pixelMM design version (1 by default, i.e. 2010-2012 prototypes)  ===

	string key = "PixelMM_design_version";
        int pixelsVersion = 0;
	if (CsOpt::Instance()->getOpt(TBname, key, pixelsVersion) ||
	    CsOpt::Instance()->getOpt(TB, key, pixelsVersion)) {
          pixmmobject = pixmm::Instance(pixelsVersion);
          pixmmobject->initconndata();
	  CsErrLog::msg(elInfo, __FILE__, __LINE__, "%s : PixelMM_design_version is no more useful for real data, you can delete it in this case",
			TBname);
	}

	// === Check if correction in U should be applied  ===

	corr_U_bool = false;
	if (CsOpt::Instance()->getOpt(TBname, "corr_U") ||
			CsOpt::Instance()->getOpt(TB, "corr_U")) {
		corr_U_bool = true;
		corr_U_tab = NULL;
		CsErrLog::msg(elWarning, __FILE__, __LINE__,
				"%s : correction in U will be applied", TBname);
	}
	else {
		CsErrLog::msg(elWarning, __FILE__, __LINE__,
				"%s: correction in U will NOT be applied!", GetTBName().c_str());
	}

	decodeCard_ = false;
	string tag = "";
	key = "make decoding";
	bool status = CsOpt::Instance()->getOpt(tag, key);
	if (status) {
		list<string> options;
		list<string>::iterator Is;
		status = CsOpt::Instance()->getOpt(tag, key, options);
		if (status) {
			for ( Is = options.begin(); Is != options.end(); Is++ ) {
				if ( *Is == TB || *Is == "PixelMumega" || *Is == "all" ) {
					decodeCard_ = true;
				}
			}
		}
		else {
			decodeCard_ = true;
		}
	}

	key = "Threshold";
	fThresholdHit_ = 3.5;
	fThresholdClu_ = 5.;
	vector<double> v;
	if (CsOpt::Instance()->getOpt(TBname, key, v) ||
			CsOpt::Instance()->getOpt(TB, key, v)) {
		if (v.size() != 2) {
			CsErrLog::msg(elFatal, __FILE__, __LINE__,
					"%s : Syntax error in specifying thresholds !", TBname);
		}
		else {
			fThresholdHit_ = v[0]; fThresholdClu_ = v[1];
		}
	}
	CsErrLog::msg(elInfo, __FILE__, __LINE__,
			"%s : Strip amplitude >= %.3f sigma, Cluster amplitude >= %.3f sigma",
			TBname, fThresholdHit_, fThresholdClu_);
	/*
	key = "AmplitudeRatio";
	if (CsOpt::Instance()->getOpt(TBname, key, v) ||
			CsOpt::Instance()->getOpt(TB, key, v)) {
		if ( v.size()%2 != 0 ) {
			CsErrLog::msg(elFatal, __FILE__, __LINE__,
					"%s : Syntax error in specifying cut on amplitude ratio !",
					TBname);
		}
		else {
			int np = v.size() / 2;
			string mes = "%s : Polygon for amplitude ratio cut \n";
			fAmpRatio23_.clear();
			fAmpRatio13_.clear();
			for ( int i = 0; i < np; i++ ) {
				fAmpRatio23_.push_back(v[2*i]);   // x coordinates
				fAmpRatio13_.push_back(v[2*i+1]); // y coordinates
				char coord[100];
				sprintf(coord, "(%.3f, %.3f)\n", v[2*i], v[2*i+1]);
				mes += coord;
			}
			// Close polygon
			fAmpRatio23_.push_back(v[0]);
			fAmpRatio13_.push_back(v[1]);
			CsErrLog::msg(elInfo, __FILE__, __LINE__,
					mes.c_str(), TBname);
		}
	}
	*/

	key = "AmplitudeRatioCuts";
	int valuearc;
	do_amplituderatiocuts = 1;
	if (CsOpt::Instance()->getOpt(TBname, key, valuearc) ||
			CsOpt::Instance()->getOpt(TB, key, valuearc)) {
		do_amplituderatiocuts = valuearc;
	}
	CsErrLog::msg(elInfo, __FILE__, __LINE__,
			"%s : AmplitudeRatioCuts %d",
			TBname, do_clustering);


	key = "Multiplicity";
	fUppMult_ = 1000.;
	fLowMult_ = 0.;
	if ( CsOpt::Instance()->getOpt(TBname, key, v) ||
			CsOpt::Instance()->getOpt(TB, key, v)) {
		if ( v.size() != 2 ) {
			CsErrLog::msg(elFatal, __FILE__, __LINE__,
					"%s : Syntax error in specifying multiplicity cut !", TBname);
		}
		else {
			fLowMult_ = v[0];
			fUppMult_ = v[1];
		}
	}
	CsErrLog::msg(elInfo, __FILE__, __LINE__,
			"%s : %.3f <= Multiplicity <= %.3f",
			TBname, fLowMult_, fUppMult_);

	key = "Clustering";
	int valuei;
	do_clustering = 1;
	if (CsOpt::Instance()->getOpt(TBname, key, valuei) ||
			CsOpt::Instance()->getOpt(TB, key, valuei)) {
		do_clustering = valuei;
	}
	CsErrLog::msg(elInfo, __FILE__, __LINE__,
			"%s : Clustering %d",
			TBname, do_clustering);

	key = "ClusterConfigMask";
	if (CsOpt::Instance()->getOpt(TBname, key, valuei) ||
			CsOpt::Instance()->getOpt(TB, key, valuei)) {
		fClusterConfigMask = valuei;
		CsErrLog::msg(elInfo, __FILE__, __LINE__,
				"%s : ClusterConfigMask %d",
				TBname, fClusterConfigMask);
	}
	else {
		fClusterConfigMask = 0xb3;
		CsErrLog::msg(elInfo, __FILE__, __LINE__,
				"%s : No ClusterConfigMask specified : using standard mask %d",
				TBname, fClusterConfigMask);
	}

	// ***** MONTE CARLO : CUT ON Hit Time *****

	fMCHitTMin = -tGate_/2; fMCHitTMax = tGate_/2;
	key = "MCHitTime";
	if (CsOpt::Instance()->getOpt(TBname, key, v) ||
			CsOpt::Instance()->getOpt(TB, key, v)) {
		if (v.size() != 2)
			CsErrLog::msg(elFatal, __FILE__, __LINE__,
					"%s : Syntax error in specifying MC hit time cut !", TBname);
		fMCHitTMin = v[0]; fMCHitTMax = v[1];
		CsErrLog::msg(elInfo, __FILE__, __LINE__, "%s : %f <= MC hit time <= %f",
				TBname, fMCHitTMin, fMCHitTMax, tGate_);
	}

	fDebugLevel = 0;
	key = "debug level";
	if (CsOpt::Instance()->getOpt(TBname, key, fDebugLevel) ||
			CsOpt::Instance()->getOpt(TB, key, fDebugLevel)) {
	}

	// === Check if amplitudes simulation in MC decoding should be used ===

	amplitudeMCDecoding_ = false;
	if (CsOpt::Instance()->getOpt(TBname, "amplitudeMCDecoding") ||
			CsOpt::Instance()->getOpt(TB, "amplitudeMCDecoding")) {
		amplitudeMCDecoding_ = true;
		CsErrLog::msg(elWarning, __FILE__, __LINE__,
				"%s : amplitudes simulation in MC decoding will be used", TBname);
	}

	// ***** MONTE CARLO : Amplitudes simulation parameters *****
	//
	// spSig_ detector space resolution (mm)
	// eGain_ effective gain - should be tuned to reproduce cluster amplitudes
	// eGsig_ gain sigma (a.u.) for amplitude correlation, for example.
	// sWidth_ signal width (mm) (should be tuned to have correct number of strips/cluster)
	// tRes_ time resolution (ns)
	// These parameters can be made dependent upon detector's TB name, cf. supra:
	// MP ampParsMC [0-4] 0.05 2000. 20. 0.30 9.
	// MP00X1__ ampParsMC [0-4] 0.05 2000. 20. 0.30 9.
	key = "ampParsMC";
	if (CsOpt::Instance()->getOpt(TBname, key, v) ||
			CsOpt::Instance()->getOpt(TB, key, v)) {

		if (v.size() != 5)
			CsErrLog::msg(elFatal, __FILE__, __LINE__,
					"%s : Syntax error in specifying MC amplitude parameters !", TBname);
		spSig_ = v[0]; eGain_ = v[1]; eGSig_ = v[2]; sWidth_ = v[3]; tRes_ = v[4];
		CsErrLog::msg(elInfo, __FILE__, __LINE__, "%s : space res. %f mm, gain %f, "
				"sigma gain %f, signal width %f mm, time res. %f ns",
				TBname, spSig_, eGain_, eGSig_, sWidth_, tRes_);
	}

	if (amplitudeMCDecoding_) {
		if (TBname[4]!='P' && TBname[4]!='M') {
			if ( sWidth_ <= 0. || eGain_ <= 0. ) {
				CsErrLog::msg(elFatal, __FILE__, __LINE__,
						"%s : amplitude simulation impossible, wrong "
						"width = %f or eGain = %f  ",
						TBname, sWidth_, eGain_);
			}
		}
	}

	//=== Check if amplitudes correlation in MC should be used ===

	string vs;
	ampCorrelationMC_ = false;
	isMaster_ = false;
	if (CsOpt::Instance()->getOpt(TBname, "ampCorrelationMC") ||
			CsOpt::Instance()->getOpt(TB, "ampCorrelationMC")) {
		ampCorrelationMC_ = true;
		CsErrLog::msg(elWarning, __FILE__, __LINE__,
				"%s : amplitudes correlation in MC decoding will be used", TBname);
		key = "Master";
		if (CsOpt::Instance()->getOpt(TB, key, vs) ) {
			if (vs.size() == 2) {
				if (GetTBName()[4] == vs[0] || GetTBName()[4] == vs[1] ) {
					isMaster_ = true;
				}
			}
		}
		key = "Master";
		if (CsOpt::Instance()->getOpt(TBname, key)) {
			isMaster_ = true;
		}
		key = "Slave";
		if ( CsOpt::Instance()->getOpt(TBname, key) ) {
			isMaster_ = false;
		}
		CsErrLog::msg(elInfo, __FILE__, __LINE__,
				"%s : detector is master if non-zero : %i ", TBname, int(isMaster_) );
	}

	fAmpCorr[0] = 0; fAmpCorr[1] = 1; fAmpCorr[2] = 0;
	if (CsOpt::Instance()->getOpt(TBname, "AmplitudeCorr", v)) {
		if (v.size() != 3)
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "%s : Syntax error in specifying AmplitudeCorr !", TBname);
		// ***** Correction to the amplitude in order to match that of counterpart
		fAmpCorr[0] = v[0]; fAmpCorr[1] = v[1]; fAmpCorr[2] = v[2];
		CsErrLog::msg(elInfo, __FILE__, __LINE__, "%s : Amplitude Corr %f %f %f",
				TBname, fAmpCorr[0], fAmpCorr[1], fAmpCorr[2]);
	}
	/*
  // read parameters for cross talk suppression algorithm
  fCrossTalkParams[0] = 0.30; fCrossTalkParams[1] = 0.55;
  if (CsOpt::Instance()->getOpt(TBname, "CrodssTalkParams", v) ||
      CsOpt::Instance()->getOpt(TB, "CrossTalkParams", v)) {
    if (v.size() != 2)
      CsErrLog::msg(elFatal, __FILE__, __LINE__,
		    "%s : Syntax error in specifying CrossTalkParams !", TBname);
    fCrossTalkParams[0] = v[0]; fCrossTalkParams[1] = v[1];
    CsErrLog::msg(elInfo, __FILE__, __LINE__, "%s : CrossTalkParams %f %f",
		  TBname, fCrossTalkParams[0], fCrossTalkParams[1]);
  }
  //  double valued;
	 */

	// activation and setting of merging of clusters
	float valuef;
	key = "ClusterMergingAmpratioThr";
	if (CsOpt::Instance()->getOpt(TBname, key, valuef) ||
			CsOpt::Instance()->getOpt(TB, key, valuef)) {
		fClusterMergingAmpratioThr = valuef;
		CsErrLog::msg(elInfo, __FILE__, __LINE__,
				"%s : ClusterMergingAmpratioThr %.3f",
				TBname, fClusterMergingAmpratioThr);
	}
	else {
		fClusterMergingAmpratioThr = 0.35;
		CsErrLog::msg(elWarning, __FILE__, __LINE__,
				"%s : No ClusterMergingAmpratioThr specified : using standard value %.3f",
				TBname, fClusterMergingAmpratioThr);
	}

	
	key = "TimeCuts";
	std::vector<float> vtc;
	if (CsOpt::Instance()->getOpt(TBname, key, vtc) ||
	    CsOpt::Instance()->getOpt(TB, key, vtc)) {
	  if (vtc.size() != 2){
	    CsErrLog::msg(elError, __FILE__, __LINE__, "%s : Syntax error in specifying Additional Time Cuts ! (Set 2 parameters)", TBname);
	    fTimeCuts.resize(2);
	    fTimeCuts[0] = -1e9;
	    fTimeCuts[1] = 1e9;
	  }
	  else {
	    fTimeCuts = vtc;
	    CsErrLog::msg(elInfo, __FILE__, __LINE__,
			  "%s : TimeCuts %.3f ns, %.3f ns",
			  TBname, fTimeCuts[0], fTimeCuts[1]);}
	}
	else {
	  fTimeCuts.resize(2);
	  fTimeCuts[0] = -1e9;
	  fTimeCuts[1] = 1e9;
	  CsErrLog::msg(elInfo, __FILE__, __LINE__,
			"%s : No Additional Time Cuts specified",
			TBname);
	}
	

	// Create CsMumegaPlane (to be used in CsPM::clusterize)
	if (TBname[4] == 'P') {  // PGEM like pixels (2009 proto)
		ChipsPerPlane = ChipsPerPixelPlane;
		pixelplane = new CsPixelMumegaPlane(GetTBName().c_str()); // pixel part
		pixelplane->GetPar()->SetClusConfigMask(fClusterConfigMask);
		pixelplane->GetPar()->SetClusMeth(do_clustering);
		pixelplane->GetPar()->SetThrClus(fThresholdClu_);
		pixelplane->GetPar()->SetThrHit(fThresholdHit_);
		pixelplane->GetPar()->SetSample(2); // sample used for clustering

		//    if (!pixelplane->GetPar()->SetPosCorrs(fPosCorrs))
		//  CsErrLog::msg(elError, __FILE__, __LINE__, "%s : Wrong format of position corrections", TBname);
	} else if (TBname[4] == 'M') {  // PMM like pixels (since 2010)
		ChipsPerPlane = ChipsPerRectPixelPlane;
		rectpixelplane = new CsRectPixelMumegaPlane(GetTBName().c_str()); // pixel part
		rectpixelplane->GetPar()->SetClusConfigMask(fClusterConfigMask);
		rectpixelplane->GetPar()->SetClusMeth(do_clustering);
		rectpixelplane->GetPar()->SetThrClus(fThresholdClu_);
		rectpixelplane->GetPar()->SetThrHit(fThresholdHit_);
		rectpixelplane->GetPar()->SetMergingAmpratioThreshold(fClusterMergingAmpratioThr);
		rectpixelplane->GetPar()->SetTimeCuts(fTimeCuts);
	} else {  // strips
		ChipsPerPlane = ChipsPerStripPlane;
		stripplane = new CsMumegaPlane(GetTBName().c_str());
		stripplane->GetPar()->SetClusMeth(do_clustering);
		stripplane->GetPar()->SetThrClus(fThresholdClu_);
		stripplane->GetPar()->SetThrHit(fThresholdHit_);
		stripplane->GetPar()->SetSample(2); // sample used for clustering
		stripplane->GetPar()->SetMergingAmpratioThreshold(fClusterMergingAmpratioThr);
		stripplane->GetPar()->SetTimeCuts(fTimeCuts);
	}

	// ***** MONTE CARLO : Clusters simulation parameters *****
	// MP clusParsMC [0-13] pix_res strip_res1 strip_res2 strip_res3 clus_size_mean clus_size_mean_const clus_size_mean_slope clus_size_sig elec_spread_mean elec_spread_sigma elec_spread_sigmabg amp2_mean amp2_sig amp2_fact
	// MP clusParsMC [0-13]
	key = "clusParsMC";
	bool def = true;
	if (CsOpt::Instance()->getOpt(TBname, key, v) ||
			CsOpt::Instance()->getOpt(TB, key, v)) {

		if (v.size() != 14){
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s : Syntax error in specifying MC clusters parameters ! Default one will be taken", TBname);
		}
		else{
			/// detector resolution in U for strips /////
			pix_res = v[0];// detector resolution for pixel zone
			strip_res1 = v[1];// detector resolution for short strip zone
			strip_res2 = v[2];// detector resolution for long strip zone (left : U < 0.0)
			strip_res3 = v[3];// detector resolution for long strip zone (right : U > 0.0)
			///////////

			//// cluster size's parameters of the detector ////
			//cluster size mean
			clus_size_mean = v[4];
			clus_size_mean_const = v[5];
			clus_size_mean_slope = v[6];
			////
			clus_size_sig = v[7]; //cluster size sigma
			///////////

			//// coefficient for electron shower size and sigma ///
			elec_spread_mean = v[8]; //sigma mean in general
			elec_spread_sigma = v[9]; //sigma dependant of sqrt(amplitude of chanel)
			elec_spread_sigmabg = v[10]; //sigma constant due to electronic noise
			///////////

			//// amplitude of the third sample ////
			amp2_mean = v[11]; //mean amplitude 2 for 2.2kev deposit
			amp2_sig = v[12]; //sigma amplitude 2
			amp2_fact = v[13];// factor between amplitude 2 and deposited energy
			///////////

			CsErrLog::msg(elInfo, __FILE__, __LINE__, "%s : pix_res %f, strip_res1 %f, strip_res2 %f, strip_res3 %f,"
					" clus_size_mean %f, clus_size_mean_const %f, clus_size_mean_slope %f, clus_size_sig %f,"
					" elec_spread_mean %f, elec_spread_sigma %f, elec_spread_sigmabg %f, amp2_mean %f, amp2_sig %f, amp2_fact %f",
					TBname, pix_res, strip_res1, strip_res2, strip_res3, clus_size_mean,clus_size_mean_const,clus_size_mean_slope,clus_size_sig,
					elec_spread_mean, elec_spread_sigma, elec_spread_sigmabg, amp2_mean, amp2_sig, amp2_fact);
			def = false;
		}
	}
	if (def) {
		if (amplitudeMCDecoding_){
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s : Default MC clusters parameters is taken", TBname);
		}
		/// detector resolution in U for strips /////
		pix_res = 0.092;// detector resolution for pixel zone the clustering add 38um
		strip_res1 = 0.048;// detector resolution for short strip zone the clustering add 57um
		strip_res2 = 0.050;// detector resolution for long strip zone (left : U < 0.0)
		strip_res3 = 0.050;// detector resolution for long strip zone (right : U > 0.0)
		///////////

		//// cluster size's parameters of the detector ////
		//cluster size mean
		clus_size_mean = 5.28;
		clus_size_mean_const = 1.35;
		clus_size_mean_slope = -0.0029;
		////
		clus_size_sig = 0.20; //cluster size sigma
		///////////

		//// coefficient for electron shower size and sigma ///
		elec_spread_mean = 0.300; //sigma mean in general
		elec_spread_sigma = 0.15; //sigma dependant of sqrt(amplitude of chanel)
		elec_spread_sigmabg = 5.; //sigma constant due to electronic noise
		///////////

		//// amplitude of the second sample ////
		amp2_mean = 480.0; //amplitude 2 mean for 2.2kev deposit
		amp2_sig = 200.0; //amplitude 2 sigma
		amp2_fact = 3.0; // factor between amplitude 2 and deposited energy
		///////////
	}
	///////**********************************************///////
}


// int CsPixelMumegaDetector::pix_address [CsPixelMumegaDetector::nbpixX][CsPixelMumegaDetector::nbpixY];

void CsPixelMumegaDetector::AddSubDetector(const int    row,   const int    id,
		const char*  name,  const char   *TBname,
		const int    unit,  const int    type,
		const double rdLen, const double xsiz,
		const double ysiz,  const double zsiz,
		const double xcm,   const double ycm,
		const double zcm,   const HepMatrix rotDRS,
		const HepMatrix rotWRS,
		const double wirD,  const double ang,
		const int    nWir,  const double wirP,
		const double eff,   const double bkg,
		const double tGate)
{
	// ***** 2ND PART OF THE CONSTRUCTION OF A CsPM OF THE PIXELISED KIND *****

	if (TBname[4]!='P' && TBname[4]!='M') {// Stripped, as opposed to pixelised, piece of pixelMM...
		// ...fall back to the base class method.
		CsDetector::AddSubDetector(row, id, name, TBname, unit,  type, rdLen,
				xsiz,  ysiz, zsiz,
				xcm,   ycm,  zcm,
				rotDRS, rotWRS, wirD, ang, nWir, wirP,
				eff,   bkg, tGate);
		return;
	}

	CsErrLog::msg(elInfo, __FILE__, __LINE__,
			"Adding subdetector to %s", GetTBName().c_str());

	if (TBname[4]=='P') {  // GEM like pixelMM
		// The following parameters must be the same in all subdetectors.

		if (strcmp(TBname, GetTBName().c_str()))
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s += %s : Dissimilar TBnames", GetTBName().c_str(), TBname);
		if ( type != getType() )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s (%d) += %s (%d) : Dissimilar types", GetTBName().c_str(), getType(), TBname, type);
		if ( rdLen - rdLen_ )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s, #%d(%f) += #%d(%f) : Dissimilar rad. lengths", GetTBName().c_str(), GetID().GetNumber(), rdLen_, id, rdLen);
		if ( !(rotDRS == rotDRS_) )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s, #%d += #%d : Added piece's MRS->WRS != this CsPixelMumega's + pi/2",
					GetTBName().c_str(), GetID().GetNumber(), id);
		HepMatrix x2y(3,3,0); x2y[1][0] = -1; x2y[0][1] = 1; x2y[2][2] = 1;
		HepMatrix M = rotWRS*x2y; M -= rotWRS_;
		double diff; int i, j; for (i = 0, diff = 0; i < 3; i++) for (j = 0; j < 3; j++)
			if ( fabs(M[i][j]) > diff ) diff = fabs(M[i][j]);
		if ( diff > 1e-6 )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s, #%d +=#%d : Added piece's MRS->WRS != this CsPixelMumega's + pi/2",
					GetTBName().c_str(), GetID().GetNumber(), id);

		// Setting parameters specific to 2nd dimension

		wirDV_ = wirD;                 // Wires distance   in 2nd dimension
		angV_  = ang;                  // Wires angle      in 2nd dimension
		sinAngV_ = sin(angV_*M_PI/180); cosAngV_ = cos(angV_*M_PI/180);
		nWirV_ = nWir;                 // Number of wires  in 2nd dimension
		wirPV_ = wirP;                 // Wire pitch       in 2nd dimension
	}

	if (TBname[4]=='M') {  // pixelMM with rectangular pixels
		// The following parameters must be the same in all subdetectors.

		if (strcmp(TBname, GetTBName().c_str()))
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s += %s : Dissimilar TBnames", GetTBName().c_str(), TBname);
		if ( type != getType() )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s (%d) += %s (%d) : Dissimilar types", GetTBName().c_str(), getType(), TBname, type);
		if ( rdLen - rdLen_ )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s, #%d(%f) += #%d(%f) : Dissimilar rad. lengths", GetTBName().c_str(), GetID().GetNumber(), rdLen_, id, rdLen);
		if ( !(rotDRS == rotDRS_) )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s, #%d += #%d : Added piece's MRS->WRS != this CsPixelMumega's + pi/2",
					GetTBName().c_str(), GetID().GetNumber(), id);
		HepMatrix x2y(3,3,0); x2y[1][0] = -1; x2y[0][1] = 1; x2y[2][2] = 1;
		HepMatrix M = rotWRS*x2y; M -= rotWRS_;
		double diff; int i, j; for (i = 0, diff = 0; i < 3; i++) for (j = 0; j < 3; j++)
			if ( fabs(M[i][j]) > diff ) diff = fabs(M[i][j]);
		if ( diff > 1e-6 )
			CsErrLog::msg(elFatal, __FILE__, __LINE__, "AddSubDetector : %s, #%d +=#%d : Added piece's MRS->WRS != this CsPixelMumega's + pi/2",
					GetTBName().c_str(), GetID().GetNumber(), id);

		// Setting parameters specific to 2nd dimension

		wirDV_ = wirD;                 // Wires distance   in 2nd dimension
		angV_  = ang;                  // Wires angle      in 2nd dimension
		sinAngV_ = sin(angV_*M_PI/180); cosAngV_ = cos(angV_*M_PI/180);
		nWirV_ = nWir;                 // Number of wires  in 2nd dimension
		wirPV_ = wirP;                 // Wire pitch       in 2nd dimension

	}
	// std::cout<<"AddSubDetector: TBname "<<TBname<<" wirDV_ "<<wirDV_<<" angV_ "<<angV_<<" sinAngV_ "<<sinAngV_<<" nWirV_ "<<nWirV_<<" wirPV_ "<<wirPV_<<std::endl;
}


void CsPixelMumegaDetector::BookHistograms( void ) {

	// Check if histograms are to be booked
	CsDetector::ReadHistLevel();
	if ( hLevel_ == None ) return;

	string tbn = GetTBName();
	CsHistograms::SetCurrentPath("/PixelMumegas");

	// Normal level
	if ( hLevel_ >= Normal ) {
		mH1[tbn+"_cA0"] = new CsHist1D(tbn+"_cA0", tbn+" cluster amplitude 0 ", 256, -2., 1022.);
		mH1[tbn+"_cA1"] = new CsHist1D(tbn+"_cA1", tbn+" cluster amplitude 1 ", 256, -2., 1022.);
		mH1[tbn+"_cA2"] = new CsHist1D(tbn+"_cA2", tbn+" cluster amplitude 2 ", 256, -2., 1022.);
		mH1[tbn+"_cSize"] = new CsHist1D(tbn+"_cSize", tbn+" cluster size ", 32, -0.5, 31.5);
		mH1[tbn+"_nClus"] = new CsHist1D(tbn+"_nClus", tbn+" number of clusters ", 32, -0.5, 31.5);
		mH1[tbn+"_cTime"] = new CsHist1D(tbn+"_cTime", tbn+" cluster time ", 100, -45.5, 50.5);
		mH2[tbn+"_cAmpRatio"] = new CsHist2D(tbn+"_cAmpRatio", tbn+" cluster amplitude ratio ", 25, -0.1, 2.4, 25, -0.1, 2.4);
		//mH1[tbn+"_cCTProbC"] = new CsHist1D(tbn+"_cCTProbC", tbn+" cluster cross talk probability after cut", 100, 0., 1.001);

		// Pixel region, designated by a "P" in the projection field of TBname
		if (tbn.find("P", 2, 1) != string::npos || tbn.find("M", 2, 1) != string::npos) {
			mH2[tbn+"_cMap"] = new CsHist2D(tbn+"_cMap", tbn+" cluster map ", 200, -30, 30, 200, -30, 30);
		} else { // strip region
			mH1[tbn+"_cPos"] = new CsHist1D(tbn+"_cPos", tbn+" cluster position ", nWir_, 0., nWir_-1);
			//	mH1[tbn+"_hPos"] = new CsHist1D(tbn+"_hPos", tbn+" hit position ", 2*nWir_+1, -0.5-nWir_, 0.5+nWir_);
		}

	}

	// High level
	if ( hLevel_ >= High ) {
	  /*	mH1[tbn+"_cA0"] = new CsHist1D(tbn+"_cA1", tbn+" cluster amplitude 0 all", 256, -2., 1022.);
		mH1[tbn+"_cA1"] = new CsHist1D(tbn+"_cA2", tbn+" cluster amplitude 1 all", 256, -2., 1022.);
		mH1[tbn+"_cA2"] = new CsHist1D(tbn+"_cA3", tbn+" cluster amplitude 2 all", 256, -2., 1022.);
		mH1[tbn+"_nHit"] = new CsHist1D(tbn+"_nHit", tbn+" number of hit strips all", 128, -0.5, 127.5);
		mH1[tbn+"_nClus"] = new CsHist1D(tbn+"_nClus", tbn+" number of clusters all", 32, -0.5, 31.5);
		mH2[tbn+"_cAmpRatio"] = new CsHist2D(tbn+"_cAmpRatio", tbn+" cluster amplitude ratio all", 25, -0.1, 2.4, 25, -0.1, 2.4);
		mH1[tbn+"_cCTProbT"] = new CsHist1D(tbn+"_cCTProbT", tbn+" cluster cross talk probability time cut", 100, 0., 1.001);
	  */
		mH1[tbn+"_ha0"] = new CsHist1D(tbn+"_ha0", tbn+" hit amplitude 0 ", 256, -2., 1022.);
		mH1[tbn+"_ha1"] = new CsHist1D(tbn+"_ha1", tbn+" hit amplitude 1 ", 256, -2., 1022.);
		mH1[tbn+"_ha2"] = new CsHist1D(tbn+"_ha2", tbn+" hit amplitude 2 ", 256, -2., 1022.);
		mH1[tbn+"_nHit"] = new CsHist1D(tbn+"_nHit", tbn+" number of hit strips ", 128, -0.5, 127.5);

		mH1[tbn+"_haiovaip1"] = new CsHist1D(tbn+"_haiovaip1", tbn+"ratio ai / ai+1 of consecutive hits in multiplexed signal ratio ", 100, 0., 1.);
		mH1[tbn+"_haip1ovai"] = new CsHist1D(tbn+"_haip1ovai", tbn+"ratio ai+1 / ai of consecutive hits in multiplexed signal ratio ", 100, 0., 1.);

		// Pixel region, designated by a "P" in the projection field of TBname
		if (tbn.find("P", 2, 1) != string::npos || tbn.find("M", 2, 1) != string::npos) {
			mH2[tbn+"_hMap"] = new CsHist2D(tbn+"_hMap", tbn+" hit map ", 200, -30, 30, 200, -30, 30);
			//	mH2[tbn+"_cMap"] = new CsHist2D(tbn+"_cMap", tbn+" cluster map all", nWirU_, -0.5, nWirU_-0.5, nWirV_, -0.5, nWirV_-0.5);
		} else { // strip region
		  //	mH1[tbn+"_cPos"] = new CsHist1D(tbn+"_cPos", tbn+" cluster position all", 2*nWir_+1, -0.5-nWir_, nWir_+0.5);
		  	mH1[tbn+"_hPos"] = new CsHist1D(tbn+"_hPos", tbn+" hit position ", nWir_, 0., nWir_-1);
		}

		if( CsInit::Instance()->IsAMonteCarloJob() )  {
		  mH1[tbn+"_MCdA1"] = new CsHist1D(tbn+"_MCdA1", tbn+" MC digit amplitude 1 ", 256, -2., 1022.);
		  mH1[tbn+"_MCdA2"] = new CsHist1D(tbn+"_MCdA2", tbn+" MC digit amplitude 2 ", 256, -2., 1022.);
		  mH1[tbn+"_MCdA3"] = new CsHist1D(tbn+"_MCdA3", tbn+" MC digit amplitude 3 ", 256, -2., 1022.);
		  mH1[tbn+"_MCdnHit"] = new CsHist1D(tbn+"_MCdnHit", tbn+" number of MC hits per digit", 32, -0.5, 31.5);
		  mH2[tbn+"_MCdAmpRatio"] = new CsHist2D(tbn+"_MCdAmpRatio", tbn+" MC digit amplitude ratio ", 25, -0.1, 2.4, 25, -0.1, 2.4);
		  mH2[tbn+"_MCnDvsT"] = new CsHist2D(tbn+"_MCnDvsT", tbn+" number of MC digits per hit vs time", 40, -220, 100, 32, -0.5, 31.5);
		  mH2[tbn+"_MChMap"] = new CsHist2D(tbn+"_MChMap", tbn+" MC hit map ", 2000, -220, 220, 2000, -220, 220);
		  mH1[tbn+"_MCHitTime"] = new CsHist1D(tbn+"_MCHitTime", tbn+" MC Hit Time ", 801, -400, 400);
		}
	}

	CsHistograms::SetCurrentPath("/");
}


void CsPixelMumegaDetector::DecodeChipDigit(const CS::Chip::Digit &digit) {
	string tbn = GetTBName();
	char mess[500];

	if (fDebugLevel > 9) {
		std::cout << "CsPixelMumegaDetector::DecodeChipDigit : " << tbn << std::endl;
	}

	// Apply mapping to all channel calibrations if not done yet
	// cannot be done in CsMumegaDetector::readCalibration because
	// CS::Chip::Maps are not yet available

	if ( !apv_chan_cals_mapped ) apv_chan_cals_mapped = mapChannelCal();

	// Check digit
	const CS::ChipAPV::Digit *d = dynamic_cast<const CS::ChipAPV::Digit *>(&digit);
	if ( d == NULL )
		throw CS::Exception("CsPixelMumegaDetector::DecodechipDigit() : Wrong digit type");

	// Check if pixel digit
	const CS::ChipAPV::DigitPixel *dp = dynamic_cast<const CS::ChipAPV::DigitPixel*>(d);
	const CS::ChipAPV::DigitPixelMM *dpmm = dynamic_cast<const CS::ChipAPV::DigitPixelMM*>(d);

	// Strip digit
	if ( dp == NULL && dpmm == NULL ) {

		if ( d->GetChannel() >= int(getNWir()) ) {
			d->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Strip) ==> Unexpected strip number %d for PixelMumega %d %s",
					d->GetChannel(), GetID().GetNumber(), GetTBName().c_str());
		}

		if ( d->GetChipChannel() >= ChipChannels ) {
			d->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Strip) : bad channel number! chip=%d, chan=%d(max=%d)",
					d->GetChip(), d->GetChipChannel(), ChipChannels-1);
		}

		if( d->IsSingleFrame() ) {
			d->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Strip): single frame is not supported yet!");
		}

		if( d->IsSparsifed() )
			sparse=true;
		else {
			d->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Strip): found digits without sparsification. Latch all events not yet supported!");
		}

		int chan_position = d->GetChanPos();
		if ( chan_position > 1 || chan_position < -1 )
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Strip) ==> Unexpected channel position %d for PixelMumega %d %s",
					d->GetChanPos(), GetID().GetNumber(), GetTBName().c_str());

		int iwir = d->GetChannel();

		// Get data
		vector<double> data;
		data.push_back(d->GetAmplitude()[0]);
		data.push_back(d->GetAmplitude()[1]);
		data.push_back(d->GetAmplitude()[2]);
		data.push_back(d->GetChanPos());
		data.push_back(d->GetChipChannel());
		data.push_back(d->GetStripConnNb());
		data.push_back(d->GetChip());

// cerr<<"strip "<<tbn<<" GetStripConnNb "<<d->GetStripConnNb()<<" GetChip "<<d->GetChip()<<endl;

		// Fill histograms
		if ( hLevel_ >= High ) {
			if ( mH1[tbn+"_hPos"]!=NULL)  mH1[tbn+"_hPos"]->Fill(iwir); // Hit position
			if(mH1[tbn+"_ha0"]!=NULL)  mH1[tbn+"_ha0"]->Fill(data[0]); // amp sample 2
			if(mH1[tbn+"_ha1"]!=NULL)  mH1[tbn+"_ha1"]->Fill(data[1]); // amp sample 2
			if(mH1[tbn+"_ha2"]!=NULL)  mH1[tbn+"_ha2"]->Fill(data[2]); // amp sample 2
		}

		//
		// Create CORAL digit
		//
		myDigits_.push_back( new CsDigit(*this, iwir, &data[0], data.size()) );

		// Log
		if (fDebugLevel>9) {
			sprintf(mess,
					"CsPixelMumegaDetector::DecodeChipDigit(Strip) : %s \n"
					"Created digit at %d, hem %d, amplitudes %d, %d, %d",
					tbn.c_str(),
					iwir,
					chan_position,
					d->GetAmplitude()[0],
					d->GetAmplitude()[1],
					d->GetAmplitude()[2]);
			std::cout << mess << std::endl;
		}
	}

	// Pixel digit
	if ( dp ) {

		if (dp->GetPixelX()>=nWirU_ || dp->GetPixelY()>=nWirV_ ||
				dp->GetPixelX()<0 || dp->GetPixelY()<0) {
			dp->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Pixel) "
					"==> Unexpected pixel number %d, %d for PixelMumega %d %s",
					dp->GetPixelX(), dp->GetPixelY(),
					GetID().GetNumber(), GetTBName().c_str());
		}

		if (dp->GetChannel()>=int(getNWir())*int(getNWirV())) {
			dp->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Pixel) "
					"==> Unexpected channel number %d for PixelMumega %d %s",
					dp->GetChannel(), GetID().GetNumber(), GetTBName().c_str());
		}

		if( dp->GetChipChannel()>=ChipChannels ){
			dp->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Pixel): "
					"bad channel number! chip=%d, chan=%d(max=%d)",
					dp->GetChip(), dp->GetChipChannel(), ChipChannels-1);
		}

		if( dp->IsSingleFrame() ) {
			dp->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Pixel): "
					"single frame is not supported yet!");
		}

		if( dp->IsSparsifed() )
			sparse=true;
		else {
			dp->Print(cerr, "BAD DIGIT:  ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Pixel): "
					"found digits without sparsification. "
					"Latch all events not yet supported!");
		}

		int iwir = d->GetChannel();
		int pixX = dp->GetPixelX();
		int pixY = dp->GetPixelY();

		// Get data
		vector<double> data;
		data.push_back((double)dp->GetAmplitude()[0]);
		data.push_back((double)dp->GetAmplitude()[1]);
		data.push_back((double)dp->GetAmplitude()[2]);
		data.push_back((double)(pixX+(pixY<<10)));
		data.push_back((double)dp->GetChipChannel());
		data.push_back((double)dp->GetChip());

		// Fill histograms
		if (hLevel_>=High) {
			if(mH2[tbn+"_hMap"]!=NULL)  mH2[tbn+"_hMap"]->Fill(pixX, pixY); // Hit map
			if(mH1[tbn+"_ha0"]!=NULL)  mH1[tbn+"_ha0"]->Fill(data[0]); // amp sample 2
			if(mH1[tbn+"_ha1"]!=NULL)  mH1[tbn+"_ha1"]->Fill(data[1]); // amp sample 2
			if(mH1[tbn+"_ha2"]!=NULL)  mH1[tbn+"_ha2"]->Fill(data[2]); // amp sample 2
		}

		//
		// Create CORAL digit, encode pixel address
		//
		myDigits_.push_back( new CsDigit(*this, iwir, &data[0], data.size()) );

		// Log
		if (fDebugLevel>9) {
			sprintf(mess,
					"CsPixelMumegaDetector::DecodeChipDigit(Pixel) : %s \n"
					"Created digit at %d, %d, %d, amplitudes %d, %d, %d",
					tbn.c_str(), iwir, pixX, pixY,
					dp->GetAmplitude()[0],
					dp->GetAmplitude()[1],
					dp->GetAmplitude()[2]);
			cout << mess << endl;
		}
	}

	// Rectangular pixel digit of pixelMM
	if ( dpmm ) {

		if( dpmm->IsSingleFrame() ) {
			dpmm->Print(cerr, "BAD DIGIT : ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Pixel): "
					"single frame is not supported yet!");
		}
		if( dpmm->IsSparsifed() )
			sparse=true;
		else {
			dpmm->Print(cerr, "BAD DIGIT:  ");
			throw CS::Exception("CsPixelMumegaDetector::DecodeChipDigit(Pixel): "
					"found digits without sparsification. "
					"Latch all events not yet supported!");
		}

                if (!pixmmobject) pixmmobject = pixmm::Instance(dpmm->GetPixelMMversion());

		int ich = d->GetChannel() + ChipChannels*dpmm->GetConnNb();
		double pixX = dpmm->GetPixelX();
		double pixY = dpmm->GetPixelY();

		// Get data
		vector<double> data;
		data.push_back((double)dpmm->GetAmplitude()[0]);
		data.push_back((double)dpmm->GetAmplitude()[1]);
		data.push_back((double)dpmm->GetAmplitude()[2]);
		data.push_back(pixX);
		data.push_back(pixY);
		data.push_back((double)dpmm->GetPixelNb());
		data.push_back((double)dpmm->GetChipChannel());
// 		data.push_back((double)dpmm->GetConnNb());
		data.push_back((double)pixmmobject->GetNConn(dpmm->GetConnNb())); // connectors 0 to 9
		data.push_back((double)dpmm->GetChip());

		// Fill histograms
		if (hLevel_>=High) {
			if(mH2[tbn+"_hMap"]!=NULL)  mH2[tbn+"_hMap"]->Fill(pixX, pixY); // Hit map
			if(mH1[tbn+"_ha0"]!=NULL)  mH1[tbn+"_ha0"]->Fill(data[0]); // amp sample 2
			if(mH1[tbn+"_ha1"]!=NULL)  mH1[tbn+"_ha1"]->Fill(data[1]); // amp sample 2
			if(mH1[tbn+"_ha2"]!=NULL)  mH1[tbn+"_ha2"]->Fill(data[2]); // amp sample 2
		}

		//
		// Create CORAL digit
		//
		myDigits_.push_back( new CsDigit(*this, ich, &data[0], data.size()) );

		// Log
		if (fDebugLevel>9) {
			sprintf(mess,
					"CsPixelMumegaDetector::DecodeChipDigit(Pixel) : %s \n"
					"Created MM pixel digit at ich %d, pixnb %d, chipch %d, connnb %d, pX %f, pY %f, amplitudes %d, %d, %d",
					tbn.c_str(), ich, dpmm->GetPixelNb(),dpmm->GetChipChannel(), dpmm->GetConnNb(), pixX, pixY,
					dpmm->GetAmplitude()[0],
					dpmm->GetAmplitude()[1],
					dpmm->GetAmplitude()[2]);
			cout << mess << endl;
		}

	}
	return;
}


void CsPixelMumegaDetector::clusterize() {

          //Get Event properties
  const CsEvent *event = CsEvent::Instance();
  int runnb = event->getRunNumber();
  int evnb = event->getEventNumberInRun();
  int date = 0;
  double tcsphase = event->getTCSPhaseTime();

	// Clear list of CORAL clusters
	clearClusterList();

	if (pixelplane) { // Pixel region
		clusterizePixels();
	} else if (rectpixelplane) { // Pixel region of new pixelMM
		rectpixelplane->GetHeader()->Set(evnb, runnb, date, tcsphase);
		clusterizeRectPixels();
	} else if (stripplane) { // Strip region
		stripplane->GetHeader()->Set(evnb, runnb, date, tcsphase);
		clusterizeStrips();
	}
}


void CsPixelMumegaDetector::clusterizeRectPixels() {

	string tbn = GetTBName();

	// CORAL digits
	list<CsDigit*>::iterator Id;
	list<CsDigit*> digits = getMyDigits();
	vector<list<CsDigit*>::iterator> iterators;
	iterators.clear();

	// Protection
	if ( digits.empty() ) return;
	HepMatrix iRotM = getRotWRSInv();

	double wireDCorrU = wirDU_ + iRotM(1,1) * _deltax + iRotM(1,2) * _deltay;
	double wireDCorrV = wirDV_ + iRotM(2,1) * _deltax + iRotM(2,2) * _deltay;

	// Map to store CORAl digits
	map<int, CsDigit*> digitmap;

	// Loop over digits and fill them as hits into the plane
	for ( Id = digits.begin(); Id != digits.end(); Id++ ) {

		// Make sure that digit belongs to this detector
		assert( (*Id)->getDet() == this );

		// Skip if already clusterized
		if ( (*Id)->isClusterized() ) continue;

		// Get amplitudes
		const float a0 = (float)(*Id)->getData()[0];
		const float a1 = (float)(*Id)->getData()[1];
		const float a2 = (float)(*Id)->getData()[2];

		// Detector channel
		const int channel = (int)(*Id)->getAddress();
		const int pixnb = (int)(*Id)->getData()[5];

		// Connector Number
		// const int connnb = (int)(*Id)->GetConnNb()

		// Decode pixel address
		const float pixX = (*Id)->getData()[3];
		const float pixY = (*Id)->getData()[4];

		int debug(1);

		if ( (pixX < -60.) || (pixX > 60.) || (pixY <  -60.) || (pixY > 60.) ) {
			if ( debug )
				cout << "----------------------------------------------------------\n"
				<< "PixelMumega : Hit with pixX = " << pixX << ", pixY = " << pixY << endl
				<< "Channel = " << (int)(*Id)->getAddress() << endl
				<< "Chip channel = " << (int)(*Id)->getData()[4] << endl
				<< "Hit ignored !\n"
				<< "--------------------------------------------------------- \n";
			continue;
		}

		// Add hit to plane
		rectpixelplane->AddHit(channel, pixnb, pixX, pixY, a0, a1, a2);

		// Add reference to CsDigit
		digitmap[channel] = *Id;

	} // end of loop over digits

	// Fill histogram with number of hits
	if ( hLevel_ >= High ) {
		if ( mH1[tbn+"_nHit"] != NULL )  mH1[GetTBName()+"_nHit"]->Fill(rectpixelplane->GetNhit());

	}
	// Clustering
	rectpixelplane->Clusterize();

	//Fill histograms with multiplexed signal amp. ratios
	if ( hLevel_ >= High ) {
		if ((mH1[tbn+"_haiovaip1"] != NULL) && (mH1[tbn+"_haip1ovai"] != NULL)){
			std::list<CsMumegaHit*> hits = rectpixelplane->GetHits();
			std::list<CsMumegaHit*>::iterator ithit;
			for (ithit = hits.begin(); ithit != hits.end(); ithit++){
				if ((*ithit)->GetMultiplexAmpRatio() != -1){
					mH1[tbn+"_haiovaip1"]->Fill((*ithit)->GetMultiplexAmpRatio());
					mH1[tbn+"_haip1ovai"]->Fill(1./((*ithit)->GetMultiplexAmpRatio()));
				}
			}
		}
	}
	// Display and print clusters
	if ( fDebugLevel > 9 ) {
		rectpixelplane->PrintClusters();
		rectpixelplane->Display();
	}

	// Interface to CORAL clusters
	int nclus = rectpixelplane->GetNcluster();
	//int ngoodclus = 0;

	static double triggerOffset; // MC : Time offset of the trigger w.r.t. the event
	const CsEvent *event = CsEvent::Instance();
	bool isMC = event->isAMonteCarloEvent();
	if (isMC) triggerOffset = event->getTriggerMCOffset();

	// Get found clusters
	list<CsRectPixelMumegaCluster*> clusters = rectpixelplane->GetClusters();

	// Loop over clusters to fill CORAL clusters
	list<CsRectPixelMumegaCluster*>::iterator itclus;
	for (itclus = clusters.begin(); itclus != clusters.end(); itclus++) {

		float pixU = (*itclus)->GetPositionX();
		float pixV = (*itclus)->GetPositionY();

		// Get cluster amplitudes
		std::vector<Float_t> Amp = (*itclus)->GetAmp();
		float Amp13 = (Amp[2] == 0) ? 0 : (Amp[0]/Amp[2]);
		float Amp23 = (Amp[2] == 0) ? 0 : (Amp[1]/Amp[2]);

		// Fill histograms //before time cut
		if (hLevel_>=Normal){
			if (mH1[tbn+"_cA0"] != NULL)  mH1[tbn+"_cA0"]->Fill(Amp[0]);
			if (mH1[tbn+"_cA1"] != NULL)  mH1[tbn+"_cA1"]->Fill(Amp[1]);
			if (mH1[tbn+"_cA2"] != NULL)  mH1[tbn+"_cA2"]->Fill(Amp[2]);
			if (mH2[tbn+"_cMap"] != NULL)  mH2[tbn+"_cMap"]->Fill(pixU, pixV);
			if (mH2[tbn+"_cAmpRatio"] != NULL)  mH2[tbn+"_cAmpRatio"]->Fill(Amp23, Amp13);
		}

		/*	// Cut on amplitude ratios of clusters
		if (fAmpRatio13_.size() > 0 && fAmpRatio23_.size() > 0) {
			if ( !CsPixelMumega::IsInside(Amp23, Amp13, fAmpRatio23_, fAmpRatio13_) ) continue;
			}*/

		/*	// Fill histograms after time cut
		if (hLevel_>=Normal){
			if (mH1[tbn+"_cA0T"] != NULL)  mH1[tbn+"_cA0T"]->Fill(Amp[0]);
			if (mH1[tbn+"_cA1T"] != NULL)  mH1[tbn+"_cA1T"]->Fill(Amp[1]);
			if (mH1[tbn+"_cA2T"] != NULL)  mH1[tbn+"_cA2T"]->Fill(Amp[2]);
			if (mH2[tbn+"_cMapT"] != NULL)  mH2[tbn+"_cMapT"]->Fill(pixU, pixV);
			if (mH1[tbn+"_cSizeT"] != NULL) mH1[tbn+"_cSizeT"]->Fill((*itclus)->GetSize());
			if (mH2[tbn+"_cAmpRatioT"] != NULL)  mH2[tbn+"_cAmpRatioT"]->Fill(Amp23, Amp13);
			}*/

		//Print cross talk probability
		// std::cout<<"pix xtalk prob ="<<(*itclus)->GetXTalk()<<std::endl;
		// Cut on cross talk probability
		/*	if ( (*itclus)->GetXTalk() > fCrossTalkParams[1] )
		  continue;*/

		// fill histogram after cut on crosstalk probability
		/*	if ( hLevel_>=Normal )
			if (mH1[tbn+"_cCTProbC"] != NULL) mH1[tbn+"_cCTProbC"]->Fill( (*itclus)->GetXTalk() );
		*/
		//	ngoodclus++;

		// Set coordinates in detector reference system
		//     const double u(pixU*wirPU_+wireDCorrU); // column
		//     const double v(pixV*wirPV_+wireDCorrV); // rows
		const double u(pixU*wirPU_/pixmmobject->GetPitchU()+wireDCorrU-pixmmobject->GetFirstWireU()); // column, pixU was already in mm with center of det at 0,0
		const double v(pixV*wirPV_/pixmmobject->GetPitchV()+wireDCorrV-pixmmobject->GetFirstWireV()); // rows, pixU was already in mm with center of det at 0,0
		double w = zcm_;
		// std:cout<<"clusterizeRectPixels: pixU "<<pixU<<" wirPU_ "<<wirPU_<<" pixV "<<pixV<<" wirPV_ "<<wirPV_<<" w "<<w
		//         <<" wireDCorrU "<<wireDCorrU<<" wireDCorrV "<<wireDCorrV<<" _deltax "<<_deltax<<" _deltay "<<_deltay
		//         <<" wirDU_ "<<wirDU_<<" wirDV_ "<<wirDV_<<std::endl;

		// Set errors
		HepMatrix cov(3,3,0); // Zero matrix
		if (GetTBName() == "MP00M1__"){
			cov(1,1) = pow( 0.160, 2 ); //awful patch for worse resolution of resistive PMM Prototype
			cov(2,2) = pow( wirPV_ / sqrt(12.0), 2 );
		} else{
			if ( (*itclus)->GetPositionErr() < 0 ) { // if cluster error is unknown
				cov(1,1) = pow( wirPU_ / sqrt(12.0), 2 );
				cov(2,2) = pow( wirPV_ / sqrt(12.0), 2 );
			} else {
				cov(1,1) = pow( (*itclus)->GetPositionXErr()*wirPU_, 2 );
				cov(2,2) = pow( (*itclus)->GetPositionYErr()*wirPV_, 2 );
			}
		}
		cov(3,3) = 1.; // assume 1 mm resolution in Z

		// Create CsCluster
		CsCluster* cluster = new CsCluster( u, v, w, cov );

		// Save CsDigits the CsCluster is made of
		std::list<CsMumegaHit*> hits = (*itclus)->GetHits();
		std::list<CsMumegaHit*>::iterator ithit;
		for ( ithit = hits.begin(); ithit != hits.end(); ithit++ ) {
			cluster->addDigit( *(digitmap[(*ithit)->GetChan()->GetId()->GetDetectorChannel()]) );
		}

		// Store three amplitudes
		cluster->addAnalogData( (*itclus)->GetAmp()[0], (*itclus)->GetNoise() );
		cluster->addAnalogData( (*itclus)->GetAmp()[1], (*itclus)->GetNoise() );
		cluster->addAnalogData( (*itclus)->GetAmp()[2], (*itclus)->GetNoise() );

		// Store x and y pos
		cluster->addAnalogData( pixU );
		cluster->addAnalogData( pixV );

		// Store crosstalk probability
		cluster->addAnalogData( (*itclus)->GetXTalk() );

		// Set detector
		cluster->addDet( *this );

		// Store cluster time (TCS phase corrected)
		double time, etime;
		if ( (*itclus)->GetTime(time, etime) ) {
			time = time - event->getTCSPhaseTime();
			if (isMC) time -= triggerOffset; // MC : Correct for trigger offset
			cluster->setTime(time, etime);
		} else
			cluster->setTime(0., 1.e9);

		//if(rectpixelplane->GetName() == "MP01MX__") std::cout<<rectpixelplane->GetName()<<" Time : "<<time<<" tcs "<<event->getTCSPhaseTime()<<std::endl;

		// Fill histograms 
		if ( hLevel_ >= Normal ) {
			if ( mH1[tbn+"_cTime"] != NULL ) mH1[tbn+"_cTime"]->Fill(time);
		}

		// Add cluster
		addCluster( *cluster );

		// std::cout<<"clusterizeRectPixels: getU "<<cluster->getU()<<" getV "<<cluster->getV()<<" cov "<<" cov "<<cluster->getCov()<<std::endl;
	}

	// fill histogram for cluster multiplicities before and after time cut
	/*	if ( hLevel_ >= Normal ) {
		if ( mH1[tbn+"_nClusT"] != NULL )  mH1[tbn+"_nClusT"]->Fill(ngoodclus);
		}*/
	if ( hLevel_>=Normal ){
		if ( mH1[tbn+"_nClus"] != NULL )  mH1[tbn+"_nClus"]->Fill(nclus);
	}

	sortClusters();
	setClusteringDone();

	// Clean up
	//  delete plane
	rectpixelplane->Clear();

}


void CsPixelMumegaDetector::clusterizePixels() {

	string tbn = GetTBName();

	// CORAL digits
	list<CsDigit*>::iterator Id;
	list<CsDigit*> digits = getMyDigits();
	vector<list<CsDigit*>::iterator> iterators;
	iterators.clear();

	// Protection
	if ( digits.empty() ) return;
	HepMatrix iRotM = getRotWRSInv();
	double wireDCorrU = wirDU_ + iRotM(1,1) * _deltax + iRotM(1,2) * _deltay;
	double wireDCorrV = wirDV_ + iRotM(2,1) * _deltax + iRotM(2,2) * _deltay;

	// Map to store CORAl digits
	map<int, CsDigit*> digitmap;

	// Loop over digits and fill them as hits into the plane
	for ( Id = digits.begin(); Id != digits.end(); Id++ ) {

		// Make sure that digit belongs to this detector
		assert( (*Id)->getDet() == this );

		// Skip if already clusterized
		if ( (*Id)->isClusterized() ) continue;

		// Get amplitudes
		const float a0 = (float)(*Id)->getData()[0];
		const float a1 = (float)(*Id)->getData()[1];
		const float a2 = (float)(*Id)->getData()[2];

		// Detector channel
		const int channel = (int)(*Id)->getAddress();

		// Decode pixel address
		int pixU = ((int)(*Id)->getData()[3]) & 0x3ff;
		int pixV = (((int)(*Id)->getData()[3])>>10) & 0x3ff;

		int debug(0);

		if ( pixU < 0 || pixU > 31 || pixV <  0 || pixV > 31 ) {
			if ( debug )
				cout << "----------------------------------------------------------\n"
				<< "PixelMumega : Hit with pixu = " << pixU << ", pixV = " << pixV << endl
				<< "Channel = " << (int)(*Id)->getAddress() << endl
				<< "Chip channel = " << (int)(*Id)->getData()[4] << endl
				<< "Hit ignored !\n"
				<< "--------------------------------------------------------- \n";
			continue;
		}

		// Add hit to plane
		pixelplane->AddHit(channel, pixU, pixV, a0, a1, a2);

		// Add reference to CsDigit
		digitmap[channel] = *Id;

	} // end of loop over digits

	// Fill histogram with number of hits
	if ( hLevel_ >= High ) {
		//  if ( mH1[tbn+"_nHit"] != NULL )  mH1[GetTBName()+"_nHit"]->Fill(pixelplane->GetNhit());

	}

	// Clustering
	pixelplane->Clusterize();

	// Display and print clusters
	if ( fDebugLevel > 9 ) {
		pixelplane->PrintClusters();
		pixelplane->Display();
	}

	// Interface to CORAL clusters
	int nclus = pixelplane->GetNcluster();
	//int ngoodclus = 0;

	static double triggerOffset; // MC : Time offset of the trigger w.r.t. the event
	const CsEvent *event = CsEvent::Instance();
	bool isMC = event->isAMonteCarloEvent();
	if (isMC) triggerOffset = event->getTriggerMCOffset();

	// Get found clusters
	list<CsPixelMumegaCluster*> clusters = pixelplane->GetClusters();

	// Loop over clusters to fill CORAL clusters
	list<CsPixelMumegaCluster*>::iterator itclus;
	for (itclus = clusters.begin(); itclus != clusters.end(); itclus++) {

		float pixU = (*itclus)->GetPositionX();
		float pixV = (*itclus)->GetPositionY();

		// Get cluster amplitudes
		std::vector<Float_t> Amp = (*itclus)->GetAmp();
		float Amp13 = (Amp[2] == 0) ? 0 : (Amp[0]/Amp[2]);
		float Amp23 = (Amp[2] == 0) ? 0 : (Amp[1]/Amp[2]);

		// Fill histograms 
		if (hLevel_>=Normal){
			if (mH1[tbn+"_cA0"] != NULL)  mH1[tbn+"_cA0"]->Fill(Amp[0]);
			if (mH1[tbn+"_cA1"] != NULL)  mH1[tbn+"_cA1"]->Fill(Amp[1]);
			if (mH1[tbn+"_cA2"] != NULL)  mH1[tbn+"_cA2"]->Fill(Amp[2]);
			if (mH2[tbn+"_cMap"] != NULL)  mH2[tbn+"_cMap"]->Fill(pixU, pixV);
			if (mH2[tbn+"_cAmpRatio"] != NULL)  mH2[tbn+"_cAmpRatio"]->Fill(Amp23, Amp13);
		}

		/*	// Cut on amplitude ratios of clusters
		if (fAmpRatio13_.size() > 0 && fAmpRatio23_.size() > 0) {
			if ( !CsPixelMumega::IsInside(Amp23, Amp13, fAmpRatio23_, fAmpRatio13_) ) continue;
		}

		// Fill histograms after time cut
		if (hLevel_>=Normal){
			if (mH1[tbn+"_cA0T"] != NULL)  mH1[tbn+"_cA0T"]->Fill(Amp[0]);
			if (mH1[tbn+"_cA1T"] != NULL)  mH1[tbn+"_cA1T"]->Fill(Amp[1]);
			if (mH1[tbn+"_cA2T"] != NULL)  mH1[tbn+"_cA2T"]->Fill(Amp[2]);
			if (mH2[tbn+"_cMapT"] != NULL)  mH2[tbn+"_cMapT"]->Fill(pixU, pixV);
			if (mH1[tbn+"_cSizeT"] != NULL) mH1[tbn+"_cSizeT"]->Fill((*itclus)->GetSize());
			if (mH2[tbn+"_cAmpRatioT"] != NULL)  mH2[tbn+"_cAmpRatioT"]->Fill(Amp23, Amp13);
		}

		// fill histogram before cut on crosstalk probability
		if ( hLevel_>=High )
			if (mH1[tbn+"_cCTProbT"] != NULL) mH1[tbn+"_cCTProbT"]->Fill( (*itclus)->GetXTalk() );

		// Cut on cross talk probability
		if ( (*itclus)->GetXTalk() > fCrossTalkParams[1] )
			continue;

		// fill histogram after cut on crosstalk probability
		if ( hLevel_>=Normal )
			if (mH1[tbn+"_cCTProbC"] != NULL) mH1[tbn+"_cCTProbC"]->Fill( (*itclus)->GetXTalk() );
		*/
		//	ngoodclus++;

		// Set coordinates in detector reference system
		const double u(pixU*wirPU_+wireDCorrU); // column
		const double v(pixV*wirPV_+wireDCorrV); // rows
		double w = zcm_;

		// Set errors
		HepMatrix cov(3,3,0); // Zero natrix
		if ( (*itclus)->GetPositionErr() < 0 ) { // if cluster error is unknown
			cov(1,1) = pow( wirPU_ / sqrt(12.0), 2 );
			cov(2,2) = pow( wirPV_ / sqrt(12.0), 2 );
		} else {
			cov(1,1) = pow( (*itclus)->GetPositionXErr()*wirPU_, 2 );
			cov(2,2) = pow( (*itclus)->GetPositionYErr()*wirPV_, 2 );
		}
		cov(3,3) = 1.; // assume 1 mm resolution in Z

		// Create CsCluster
		CsCluster* cluster = new CsCluster( u, v, w, cov );

		// Save CsDigits the CsCluster is made of
		std::list<CsMumegaHit*> hits = (*itclus)->GetHits();
		std::list<CsMumegaHit*>::iterator ithit;
		for ( ithit = hits.begin(); ithit != hits.end(); ithit++ ) {
			cluster->addDigit( *(digitmap[(*ithit)->GetChan()->GetId()->GetDetectorChannel()]) );
		}

		// Store three amplitudes
		cluster->addAnalogData( (*itclus)->GetAmp()[0], (*itclus)->GetNoise() );
		cluster->addAnalogData( (*itclus)->GetAmp()[1], (*itclus)->GetNoise() );
		cluster->addAnalogData( (*itclus)->GetAmp()[2], (*itclus)->GetNoise() );

		// Store x and y pos
		cluster->addAnalogData( pixU );
		cluster->addAnalogData( pixV );

		// Store crosstalk probability
		cluster->addAnalogData( (*itclus)->GetXTalk() );

		// Set detector
		cluster->addDet( *this );

		// Store cluster time (TCS phase sorrected)
		double time, etime;
		if ( (*itclus)->GetTime(time, etime) ) {
		  time = time - event->getTCSPhaseTime();
			if (isMC) time -= triggerOffset; // MC : Correct for trigger offset
			cluster->setTime(time, etime);
		} else
			cluster->setTime(0., 1.e9);


		// Fill histograms after time cut
		if ( hLevel_ >= Normal ) {
			if ( mH1[tbn+"_cTime"] != NULL ) mH1[tbn+"_cTime"]->Fill(time);
		}

		// Add cluster
		addCluster( *cluster );
	}

	// fill histogram for cluster multiplicities before and after time cut
	/*	if ( hLevel_ >= Normal ) {
		if ( mH1[tbn+"_nClusT"] != NULL )  mH1[tbn+"_nClusT"]->Fill(ngoodclus);
		}*/
	if ( hLevel_>=Normal ){
		if ( mH1[tbn+"_nClus"] != NULL )  mH1[tbn+"_nClus"]->Fill(nclus);
	}

	sortClusters();
	setClusteringDone();

	// Clean up
	//  delete plane
	pixelplane->Clear();
}


void CsPixelMumegaDetector::clusterizeStrips() {
	string tbn = GetTBName();

	// CORAL digits
	list<CsDigit*>::iterator Id;
	list<CsDigit*> digits = getMyDigits();
	vector<list<CsDigit*>::iterator> iterators;
	iterators.clear();

	// Protection
	if ( digits.empty() ) return;

	int err;
	HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
	double wireDCorr = wirD_ + iRotM(1,1) * _deltax + iRotM(1,2) * _deltay;

	// Map to store CORAL digits
	map<pair<int, int>, CsDigit*> digitmap;

	// Loop over digits and fill them as hits into the plane
	for ( Id = digits.begin(); Id != digits.end(); Id++ ) {

		// Make sure that digit belongs to this detector
		assert( (*Id)->getDet() == this );

		// Skip if already clusterized
		if ( (*Id)->isClusterized() ) continue;
		(*Id)->setClusterized();

		// Get amplitudes
		const float a0 = (float)(*Id)->getData()[0];
		const float a1 = (float)(*Id)->getData()[1];
		const float a2 = (float)(*Id)->getData()[2];

		// Detector channel
		const int channel = (int)(*Id)->getAddress();

		// Hemisphere
		int hem = (int)(*Id)->getData()[3];

		// Add hit to plane
		stripplane->AddHit(channel, hem, a0, a1, a2);

		// Add reference to csDigit
		pair<int, int> id(channel, hem);
		digitmap[id] = *Id;

	} // end of interface

	// Fill histogram with number of hits
	if ( hLevel_ >= High ) {
		if ( mH1[tbn+"_nHit"] != NULL )  mH1[tbn+"_nHit"]->Fill(stripplane->GetNhit());
	}

	// Clustering
	stripplane->Clusterize();

	//Fill histograms with multiplexed signal amp. ratios
	if ( hLevel_ >= High ) {
		if ((mH1[tbn+"_haiovaip1"] != NULL) && (mH1[tbn+"_haip1ovai"] != NULL)){
			std::list<CsMumegaHit*> hits = stripplane->GetHits();
			std::list<CsMumegaHit*>::iterator ithit;
			for (ithit = hits.begin(); ithit != hits.end(); ithit++){
				if ((*ithit)->GetMultiplexAmpRatio() != -1){
					mH1[tbn+"_haiovaip1"]->Fill((*ithit)->GetMultiplexAmpRatio());
					mH1[tbn+"_haip1ovai"]->Fill(1./((*ithit)->GetMultiplexAmpRatio()));
				}
			}
		}
	}

	// Display and dump
	if ( fDebugLevel > 9 ) {
		stripplane->PrintClusters();
		stripplane->Display();
	}

	// Interface to CORAL clusters
	int nclus = stripplane->GetNcluster();
	//	int ngoodclus = 0;

	static double triggerOffset; // MC : Timeoffset of the trigger w.r.t. the event
	const CsEvent *event = CsEvent::Instance();
	bool isMC = event->isAMonteCarloEvent();
	if (isMC) triggerOffset = event->getTriggerMCOffset();

	// Get found clusters
	list<CsMumegaCluster*> clusters = stripplane->GetClusters();

	// Loop over clusters to fill CORAL clusters
	list<CsMumegaCluster*>::iterator itclus;
	for ( itclus = clusters.begin(); itclus != clusters.end(); itclus++ ) {

		// Get cluster amplitudes
		std::vector<Float_t> Amp = (*itclus)->GetAmp();
		float Amp13 = (Amp[2] == 0) ? 0 : (Amp[0]/Amp[2]);
		float Amp23 = (Amp[2] == 0) ? 0 : (Amp[1]/Amp[2]);

		// fill histograms 
		if ( hLevel_ >= Normal ) {
			if ( mH1[tbn+"_cA0"] != NULL )  mH1[tbn+"_cA0"]->Fill(Amp[0]);
			if ( mH1[tbn+"_cA1"] != NULL )  mH1[tbn+"_cA1"]->Fill(Amp[1]);
			if ( mH1[tbn+"_cA2"] != NULL )  mH1[tbn+"_cA2"]->Fill(Amp[2]);
			if ( mH1[tbn+"_cPos"] != NULL )  mH1[tbn+"_cPos"]->Fill((*itclus)->GetPosition());
			if ( mH2[tbn+"_cAmpRatio"] != NULL )  mH2[tbn+"_cAmpRatio"]->Fill(Amp23, Amp13);
		}

		/*	// Cut on amplitude ratios of clusters
		if ( fAmpRatio13_.size() > 0 && fAmpRatio23_.size() > 0 ) {
			if ( !CsPixelMumega::IsInside(Amp23, Amp13, fAmpRatio23_, fAmpRatio13_) ) continue;
		}

		// Fill histograms after time cut
		if ( hLevel_ >= Normal ) {
			if ( mH1[tbn+"_cA0T"] != NULL )  mH1[tbn+"_cA0T"]->Fill(Amp[0]);
			if ( mH1[tbn+"_cA1T"] != NULL )  mH1[tbn+"_cA1T"]->Fill(Amp[1]);
			if ( mH1[tbn+"_cA2T"] != NULL )  mH1[tbn+"_cA2T"]->Fill(Amp[2]);
			if ( mH1[tbn+"_cPosT"] != NULL )  mH1[tbn+"_cPosT"]->Fill( (*itclus)->GetPosition() );
			if ( mH1[tbn+"_cSizeT"] != NULL ) mH1[tbn+"_cSizeT"]->Fill( (*itclus)->GetSize() );
			if ( mH2[tbn+"_cAmpRatioT"] != NULL )  mH2[tbn+"_cAmpRatioT"]->Fill(Amp23, Amp13);
		}

		// fill histogram before cut on crosstalk probability
		if ( hLevel_ >= High )
			if ( mH1[tbn+"_cCTProbT"] != NULL ) mH1[tbn+"_cCTProbT"]->Fill( (*itclus)->GetXTalk() );

		//Print cross talk probability
		// std::cout<<" xtalk prob ="<<(*itclus)->GetXTalk()<<std::endl;
		// Cut on cross talk probability
		if ( (*itclus)->GetXTalk() > fCrossTalkParams[1] )
			continue;

		// fill histogram after cut on crosstalk probability
		if ( hLevel_ >= Normal )
			if ( mH1[tbn+"_cCTProbC"] != NULL ) mH1[tbn+"_cCTProbC"]->Fill( (*itclus)->GetXTalk() );

			ngoodclus++;*/

		//; Set coordinates in detector reference system
		double u;
		if ( IsVarPitch() )
			u = wireDCorr + Wire2Pos( (*itclus)->GetPosition() );  // perp. to wire
		else
			u = wireDCorr + (*itclus)->GetPosition() * wirP_;  // perp. to wire
		double v = 0;                                      // parallel to wire
		double w = zcm_;

		// Set errors
		HepMatrix cov(3,3,0);  // Zero matrix
		if ( (*itclus)->GetPositionErr() < 0 ) { // if cluster error is unknown
			if ( IsVarPitch() ) cov(1,1) = pow( Pitch((*itclus)->GetPosition()) / sqrt(12.0), 2 );
			else                cov(1,1) = pow( wirP_ / sqrt(12.0), 2 );
		} else {
		  if ( IsVarPitch() ) cov(1,1) = pow( (*itclus)->GetPositionErr()*Pitch((*itclus)->GetPosition()), 2 );
		  else                cov(1,1) = pow( (*itclus)->GetPositionErr()*wirP_, 2 );
		}
		cov(2,2) = pow( max(getXsiz(), getYsiz()), 2 ); // just very big error
		cov(3,3) = 10.;

		// Create CsCluster
		CsCluster* cluster = new CsCluster( u, v, w, cov );

		// Save CsDigits the CsCluster is made of
		list<CsMumegaHit*> hits = (*itclus)->GetHits();
		list<CsMumegaHit*>::iterator ithit;
		for ( ithit=hits.begin(); ithit!=hits.end(); ithit++ ) {
			pair<int,int> id( (*ithit)->GetChan()->GetId()->GetDetectorChannel(),
					(*ithit)->GetChan()->GetId()->GetHemisphere() );
			cluster->addDigit( *(digitmap[id]) );
		}

		// Store three amplitudes
		cluster->addAnalogData( (*itclus)->GetAmp()[0], (*itclus)->GetNoise() );
		cluster->addAnalogData( (*itclus)->GetAmp()[1], (*itclus)->GetNoise() );
		cluster->addAnalogData( (*itclus)->GetAmp()[2], (*itclus)->GetNoise() );

		// Store wire position
		cluster->addAnalogData( (*itclus)->GetPosition(), (*itclus)->GetHemisphere() );

		// Store crosstalk probability
		cluster->addAnalogData( (*itclus)->GetXTalk() );

		cluster->addDet( *this );

		// Store cluster time (TCS phase corrected)
		double time, etime;
		if ( (*itclus)->GetTime(time, etime) ) {
			time = time - event->getTCSPhaseTime();
			if (isMC) time -= triggerOffset;   // MC : Correct for trigger offset
			cluster->setTime(time, etime);
		} else
			cluster->setTime(0., 1.e9);

		// Fill histograms 
		if ( hLevel_ >= Normal ) {
			if ( mH1[tbn+"_cTime"] != NULL ) mH1[tbn+"_cTime"]->Fill(time);
		}

		// Add cluster
		addCluster( *cluster );

	} // end of loop over found clusters

	// Fill histogram for cluster mutliplicities before and after time cut
	if ( hLevel_ >= Normal ) {
		if ( mH1[tbn+"_nClus"] != NULL )  mH1[tbn+"_nClus"]->Fill(nclus);
	}
	/*	if ( hLevel_ >= High ) {
		if ( mH1[tbn+"_nClusT"] != NULL )  mH1[tbn+"_nClusT"]->Fill(ngoodclus);
		}*/

	sortClusters();
	setClusteringDone();

	// Clean up
	stripplane->Clear();
}

//
// Local class for detector response simulation
// (Used in CsPixelMumegaDetector::makeMCDecoding())
//

// temporary digits
class pGMdig
{
public:
	int     wire; // wire #
	int     hemi; // hemisphere (-1 for lower half, +1 for upper half)
	double  amp1; // amp1
	double  amp2; // amp2
	double  amp3; // amp3
	CsMCHit* ref; // reference to MCHit
	bool operator < (const pGMdig& gd) const
	{
		return ( wire == gd.wire ? hemi < gd.hemi : wire < gd.wire );
	};
};

void getAmp(const double& t, double amps[]);

void CsPixelMumegaDetector::getMCAmp(const double& t, double energy, double amps[] ){
	/// inspired from PixelGEM /////

	//energy in keV
	// timing calibrations , default ones.
	const double t0_1  =  -40.;
	const double sl_1  =  22.;
	const double tc_1  = 1.65;

	const double t0_2  = 0.0;
	const double sl_2  = 27.;
	const double tc_2  = 1.23;

	amps[2]= (amp2_mean + amp2_sig*CsRandom::gauss())*(energy/2.2)*amp2_fact;


	// inverted ratios from getGEMtime()
	double r13=tc_1/(exp((t-t0_1)/sl_1)+1.);
	double r23=tc_2/(exp((t-t0_2)/sl_2)+1.);

	amps[0]=amps[2]*r13;
	amps[1]=amps[2]*r23;

	return;

}

void CsPixelMumegaDetector::makeMCDecoding() {
	// Apply mapping to all channel calibrations if not done yet
	if ( !apv_chan_cals_mapped ) apv_chan_cals_mapped = mapMCChannelCal();

	//Find the complementary detector, if it exist, for sharing the MChits between the detectors
	if(!compl_det_searched){
		list<CsDetector*>::iterator id;
		list<CsDetector*> det = CsGeom::Instance()->getDetectors();
		const char* TBname_comp;
		const char* TBname_det = this->GetTBName().c_str();
		for( id=det.begin(); id!=det.end(); id++ ) {
			if(strncmp( (*id)->GetTBName().c_str(), TBname_det, 4 ) == 0){
				TBname_comp = (*id)->GetTBName().c_str();

				if((TBname_det[4] == TBname_comp[5])||(TBname_det[5] == TBname_comp[4])){
					CsPixelMumegaDetector *comp_det = dynamic_cast<CsPixelMumegaDetector *>(*id);
					if ( comp_det == NULL ){
						CsErrLog::msg(elError, __FILE__, __LINE__,
								"%s : Wrong complementary detector found : %s", this->GetTBName().c_str(), TBname_comp);
					}
					else{
						comp_det->setComplementaryDet(*this);
						comp_det->setComplDetSearched(true);
						complementary_det = comp_det;

						CsErrLog::msg(elError, __FILE__, __LINE__,
								"complementary detectors found are : %s and %s", this->GetTBName().c_str(), this->getComplementaryDet()->GetTBName().c_str());

						compl_det_searched = true;
						break;
					}
				}
			}
		}

		if(!compl_det_searched){
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s : complementary detectors not found !!", this->GetTBName().c_str());
			complementary_det = NULL;
			compl_det_searched = true;
		}
	}
	////////////

	if ( !decode_ && !decodeCard_ ) return;
	if ( decodingDone_ ) return;
	myDigits_.clear();

	if ( CsGeant3::Instance()->isAnNtFile() ) {
		// COMGeant NTuples (an obsolete type of output used in earlier COMPASS
		// times) are not supported by the CsPG class
		CsErrLog::msg(elFatal, __FILE__, __LINE__,
				"%s : No MC decoding COMGeant NTuple", GetTBName().c_str());
	}

	if ( nWirV_ ) { //                 Pixelised CsPixelMumegaDetector
		if(complementary_det != NULL) MChits_Xdet_pix(); //crossing the MChit between the two complementaries detectors from strip to pixel
		if ( !amplitudeMCDecoding() ){ // Use simplistics decoding
			makeMCPixelsSimple();
		}
		else {//                          Use amplitude decoding
			makeMCPixelsAmplitude();
		}
	} else { //                        csPixelMumegaDetector of strip kind
		if(complementary_det != NULL) MChits_Xdet_strip(); //crossing the MChit between the two complementaries detectors from pixel to strip
		if ( !amplitudeMCDecoding() )
			makeMCStripsSimple();
		else
			makeMCStripsAmplitude();
	}

	decodingDone_ = true;
}


void CsPixelMumegaDetector::makeMCPixelsSimple() {
	list<CsMCHit*>::iterator Ih;
	for ( Ih = myMCHits_.begin(); Ih != myMCHits_.end(); Ih++ ) { // loop on hits
		if ( ( ( (*Ih)->getMCTrack() )->getParticle() )->getCharge() == 0 && (*Ih)->getOrigin() == 0 )
			continue;

		// ***** LOOP ON HITS from CHARGED PARTICLES or CHARGED PRODUCTS *****

		CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
		if ( thit == 0 ) return;

		double t  = thit->getDTime(); // Delay time (ns)
		double ui = thit->getUin();   // Hit in point (DRS)
		double vi = thit->getVin();
		double wi = thit->getWin();
		double uo = thit->getUout();  // hit out point (DRS)
		double vo = thit->getVout();
		double wo = thit->getWout();

		// ***** TIMING : SMEARING, TIME CUT
		double tdc = t+tRes_*CsRandom::gauss();	
		string tbn = GetTBName();
		if ( mH1[tbn+"_MCHitTime"]!=NULL)  mH1[tbn+"_MCHitTime"]->Fill(tdc); // MCHit time
		if ( fMCHitTMin > tdc || tdc > fMCHitTMax )
			continue;

		HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
		double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+ xcm_ - _xshift;
		double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+ ycm_ - _yshift;
		double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+ zcm_;
		double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+ xcm_ - _xshift;
		double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+ ycm_ - _yshift;
		double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+ zcm_;
		double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
		double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
		double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
		double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
		double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
		double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;

		if ( mH2[tbn+"_MChMap"]!=NULL)  mH2[tbn+"_MChMap"]->Fill((Ui+Uo)/2.,(Vi+Vo)/2.); // Hit map

		int iu, iv, pixUF, pixVF; // First X/Y raw for this hit
		if ( (Ui - wirDU_)/wirPU_ < 0 ) pixUF = int( (Ui - wirDU_)/wirPU_ - 0.5 );
		else                            pixUF = int( (Ui - wirDU_)/wirPU_ + 0.5 );
		if ( (Vi - wirDV_)/wirPV_ < 0 ) pixVF = int( (Vi - wirDV_)/wirPV_ - 0.5 );
		else                            pixVF = int( (Vi - wirDV_)/wirPV_ + 0.5 );
		int pixUL, pixVL; // Last X/Y raw for this hit
		if ( (Uo - wirDU_)/wirPU_ < 0 ) pixUL = int( (Uo - wirDU_)/wirPU_ - 0.5 );
		else                            pixUL = int( (Uo - wirDU_)/wirPU_ + 0.5 );
		if ( (Vo - wirDV_)/wirPV_ < 0 ) pixVL = int( (Vo - wirDV_)/wirPV_ - 0.5 );
		else                            pixVL = int( (Vo - wirDV_)/wirPV_ + 0.5 );
		if ( pixUL < pixUF ) {
			int tmp = pixUL; pixUL = pixUF; pixUF = tmp;
		}
		if ( pixVL < pixVF ) {
			int tmp = pixVL; pixVL = pixVF; pixVF = tmp;
		}

		for ( iu = pixUF; iu <= pixUL; iu++ ) for ( iv = pixVF; iv <= pixVL; iv++ ) {

			// ****** DOES A CsDigit ALREADY EXIST ON THIS PIXEL ? ******

			list<CsDigit*>::iterator Id; bool found;
			for ( Id = myDigits_.begin(), found = false; Id != myDigits_.end(); Id++ ) {
				int raw = (*Id)->getAddress(), col = raw%nWirU_; raw /= nWirU_;
				if ( iu == col && iv == raw ) { // Here it is...
					found = true;                 // ... ADD THIS HIT TO EXISTING CsDigit...
					dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
					double info[3];
					getAmp(tdc, info);
					for ( unsigned int infoidx = 0; infoidx < 3; infoidx++ )
						(*Id)->getData()[infoidx] += info[infoidx]*20000.;
					break;
				}
			}
			if ( !found ) { // no digits found with these wire and detector
				if ( iu < 0 || nWirU_ <= iu || iv < 0 || nWirV_ <= iv ) {
					CsErrLog::msg(elAnomaly, __FILE__, __LINE__,
							"%s : Pixel (= %d, %d) outside range [0, %d[x[0, %d[",
							GetTBName().c_str(), iu, iv, nWirU_, nWirV_);
					continue;
				} else {   // ***** NEW CsDigit...

					// Lump the 2 pixel coord's into a single CsDigit's address
					double info[4];
					getAmp(tdc, info);
					for ( unsigned int infoidx = 0; infoidx < 3; infoidx++ )
						info[infoidx] *= 20000.;
					info[3] = (iv<<10)+iu;
					CsDigit *digit = new CsMCDigit(*this, iu+iv*nWirU_, info, 4);
					dynamic_cast<CsMCDigit*>(digit)->addHit(*(*Ih));
					myDigits_.push_back( digit );  // ***** ADD IT TO LIST
				}
			}
		} // End of loop on hits associated to current MC hit
	} // End of loop on MC hits

	// Fill histograms
	if ( hLevel_ >= High ) {
		list<CsDigit*>::iterator Id;
		for (Id = myDigits_.begin(); Id != myDigits_.end(); Id++) {
			double amp1 = (*Id)->getData()[0];
			double amp2 = (*Id)->getData()[1];
			double amp3 = (*Id)->getData()[2];
			string tbn  = GetTBName();
			if ( mH1[tbn+"_MCdA1"] != NULL )  mH1[tbn+"_MCdA1"]->Fill(amp1);
			if ( mH1[tbn+"_MCdA2"] != NULL )  mH1[tbn+"_MCdA2"]->Fill(amp2);
			if ( mH1[tbn+"_MCdA3"] != NULL )  mH1[tbn+"_MCdA3"]->Fill(amp3);
			if ( mH2[tbn+"_MCdAmpRation"] != NULL )  mH2[tbn+"_MCdAmpRatio"]->Fill(amp2/amp3,  amp1/amp3);
		}
	}
}


void CsPixelMumegaDetector::makeMCPixelsAmplitude() {

	list<CsMCHit*>::iterator Ih;
	for ( Ih = myMCHits_.begin(); Ih != myMCHits_.end(); Ih++ ) { // loop on hits

		if ( ( ( (*Ih)->getMCTrack() )->getParticle() )->getCharge() == 0 && (*Ih)->getOrigin() == 0 )
			continue;

		// ***** LOOP ON HITS from CHARGED PARTICLES or CHARGED PRODUCTS *****

		CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
		if ( thit == 0 ) return;

		double t  = thit->getDTime(); // Delay time (ns)
		double ui = thit->getUin();   // Hit in point (DRS)
		double vi = thit->getVin();
		double wi = thit->getWin();
		double uo = thit->getUout();  // hit out point (DRS)
		double vo = thit->getVout();
		double wo = thit->getWout();


		// ***** TIMING : SMEARING, TIME CUT
		double tdc = t+tRes_*CsRandom::gauss();
		string tbn = GetTBName();
		if ( mH1[tbn+"_MCHitTime"]!=NULL)  mH1[tbn+"_MCHitTime"]->Fill(tdc); // MCHit time
		if ( fMCHitTMin > tdc || tdc > fMCHitTMax )
			continue;

		HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
		double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+ xcm_ - _xshift;
		double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+ ycm_ - _yshift;
		double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+ zcm_;
		double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+ xcm_ - _xshift;
		double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+ ycm_ - _yshift;
		double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+ zcm_;
		double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
		double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
		double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
		double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
		double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
		double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;
		double uShift = iRotM(1,1)*xcm_+iRotM(1,2)*ycm_;  // shift of plane wrt WRS, for cuts based on detector geometry
		double vShift = iRotM(2,1)*xcm_+iRotM(2,2)*ycm_;

		double pathX = fabs((Ui - Uo)/2.);

		if ( mH2[tbn+"_MChMap"]!=NULL)  mH2[tbn+"_MChMap"]->Fill((Ui+Uo)/2.,(Vi+Vo)/2.); // Hit map

		///// computation of ADC amplitude ///
		double amplitude[3];
		double Eloss = thit->getELos()*1.E6;
		getMCAmp(tdc, Eloss, amplitude);
		///////////////////////////

		double U_final = ((Ui+Uo)/2.) + pix_res*CsRandom::gauss(); // WRS
		double V_final = ((Vi+Vo)/2.) + pix_res*CsRandom::gauss();
		//	double clus_size = clus_size_mean + clus_size_sig*CsRandom::gauss(); //(mm)
		double mean = (clus_size_mean - exp(clus_size_mean_const + (clus_size_mean_slope*amplitude[2])));
		double clus_size = max((mean*0.4 + clus_size_sig*CsRandom::gauss()),0.); //(mm)
		double clus_rad = clus_size/2.; //cluster radius (cluster is supposed circular) (mm)

		double U_clust_max = (U_final + clus_rad);
		double U_clust_min = (U_final - clus_rad);
		double V_clust_max = (V_final + clus_rad);
		double V_clust_min = (V_final - clus_rad);

		int iu, iv, pixUT, pixVT; // max top-left pixel of the "cluster"
		if ( (U_clust_min - wirDU_)/wirPU_ < 0 ) pixUT = int( (U_clust_min - wirDU_)/wirPU_ - 0.5 );
		else                            pixUT = int( (U_clust_min - wirDU_)/wirPU_ + 0.5 );
		if ( (V_clust_max - wirDV_)/wirPV_ < 0 ) pixVT = int( (V_clust_max - wirDV_)/wirPV_ - 0.5 );
		else                            pixVT = int( (V_clust_max - wirDV_)/wirPV_ + 0.5 );
		int pixUB, pixVB; // max bottom-right of the "cluster"
		if ( (U_clust_max - wirDU_)/wirPU_ < 0 ) pixUB = int( (U_clust_max - wirDU_)/wirPU_ - 0.5 );
		else                            pixUB = int( (U_clust_max - wirDU_)/wirPU_ + 0.5 );
		if ( (V_clust_min - wirDV_)/wirPV_ < 0 ) pixVB = int( (V_clust_min - wirDV_)/wirPV_ - 0.5 );
		else                            pixVB = int( (V_clust_min - wirDV_)/wirPV_ + 0.5 );

		/// border pixels area TEST (under-range and over-range) ////
		if(( pixUB < 0)&&(pixUT < 0)) continue; //all cluster under-range in U
		if(( pixUB >= pixmmobject->GetNbPixU())&&(pixUT >= pixmmobject->GetNbPixU())) continue;//all cluster over-range in U
		if(( pixVB < 0)&&(pixVT < 0)) continue;//all cluster under-range in V
		if(( pixVB >= pixmmobject->GetNbPixV())&&(pixVT >= pixmmobject->GetNbPixV())) continue;//all cluster over-range in V

		double out_U_fact = 0.; //per cent of cluster charge outside the detector in U
		if(U_clust_min < -(getXsiz()*0.5)+uShift){
			out_U_fact += getErfFact(U_clust_min, U_final, -(getXsiz()*0.5)+uShift, clus_size);
		}
		if(U_clust_max > (getXsiz()*0.5)+uShift){
			out_U_fact += getErfFact((getXsiz()*0.5)+uShift, U_final, U_clust_max, clus_size);
		}

		double out_V_fact = 0.; //per cent of cluster charge outside the detector in V
		if(V_clust_min < -(getYsiz()*0.5)+vShift){
			out_V_fact += getErfFact(V_clust_min, V_final, -(getYsiz()*0.5)+vShift, clus_size);
		}
		if(V_clust_max > (getYsiz()*0.5)+vShift){
			out_V_fact += getErfFact((getYsiz()*0.5)+vShift, V_final, V_clust_max, clus_size);
		}
		double out_fact = (out_V_fact + out_U_fact - (out_U_fact*out_V_fact));
		////////////////

		////border limitation for the table pix
		pixUB = max(pixUB,0);
		pixUT = max(pixUT,0);
		pixVT = max(pixVT,0);
		pixVB = max(pixVB,0);

		pixUB = min(pixUB,(pixmmobject->GetNbPixU() - 1));
		pixUT = min(pixUT,(pixmmobject->GetNbPixU() - 1));
		pixVT = min(pixVT,(pixmmobject->GetNbPixV() - 1));
		pixVB = min(pixVB,(pixmmobject->GetNbPixV() - 1));
		/////////////////////////////////////

		if ( pixUB < pixUT ) { //not really needed but in case of error (top pixel must be higher than bottom one).
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s: MC pixels decoding wrong in U : pix bottom : %d , pix top : %d",
					GetTBName().c_str(), pixUB, pixUT);
			int tmp = pixUB; pixUB = pixUT; pixUT = tmp;
		}
		if ( pixVT < pixVB) {
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s: MC pixels decoding wrong in V : pix bottom : %d , pix top : %d",
					GetTBName().c_str(), pixVB, pixVT);
			int tmp = pixVB; pixVB = pixVT; pixVT = tmp;
		}

		// table of weight for amplitude repartition (-1 = digit not selected)
		int nb_dig_V = pixVT - pixVB + 1;
		int nb_dig_U = pixUB - pixUT + 1;
		int nb_digit = nb_dig_V*nb_dig_U;
		double W_tab[nb_digit];
		double W_total = 0.0;
		double W_tmp = 0.0;
		for (int i = 0; i<nb_digit; i++){
			W_tab[i] = -1.0;
		}

		//check if the pixel is in the cluster area and attribute it a weight
		for ( iv = pixVB; iv <= pixVT; iv++ ) for ( iu = pixUT; iu <= pixUB; iu++ ) {
			if (dist_pixtopoint(iu,iv,U_final,V_final) < clus_rad){
				W_tmp = weight_for_pix(iu,iv,U_final,V_final,clus_size, pathX);
				W_total += W_tmp;
				W_tab[(iu - pixUT) + (pixVT - iv)*nb_dig_U] = W_tmp;
			}
		}
		W_total = (W_total/(1. - out_fact)); //take into account the percent of charge outside the detector
		///////

		map<int, double> digitmap; //first is Address of pixel and second is the weight for energy loss
		map<int, double>::iterator it_dig;
		int address;

		//fusion and weight afectation for each real pixel///
		for ( iv = pixVB; iv <= pixVT; iv++ ) for ( iu = pixUT; iu <= pixUB; iu++ ) {
			W_tmp = W_tab[(iu - pixUT) + (pixVT - iv)*nb_dig_U];
			if((W_tmp > 0.)&& (iu < pixmmobject->GetNbPixU())&&(iv < pixmmobject->GetNbPixV())){ // if the pixel is hited and belonging to det
// 				address = CsPixelMumegaDetector::pix_address[iu][iv];
// 				address = pixmmobject->PixelAddresses()[(iu,iv)];
				address = pixmmobject->PixelAddr(iu,iv);
				if(address != -1){
					if (digitmap.count(address)> 0){ // we fusion the half pixel as real pixel (long and short)
						digitmap[address] += (W_tmp/W_total);
					}
					else{
						digitmap[address] = (W_tmp/W_total);
					}
				}
			}
		}
		///////

		// ****** DOES A CsDigit ALREADY EXIST ON THIS PIXEL ? ******
		double weight = 0.;
		list<CsDigit*>::iterator Id;
		bool found;
		int address_found;

		for (it_dig = digitmap.begin(); it_dig != digitmap.end(); it_dig++){
			weight = (*it_dig).second;
			found = false;
			for ( Id = myDigits_.begin(); Id != myDigits_.end(); Id++ ) {
				address_found = (*Id)->getAddress();
				if ((*it_dig).first == address_found) { // Here it is...
					found = true;                 // ... ADD THIS HIT TO EXISTING CsDigit...
					dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
					for ( unsigned int i = 0; i < 3; i++ ){
						(*Id)->getData()[i] += max((amplitude[i]*weight + elec_spread_sigmabg*CsRandom::gauss()),0.);
					}
					break;
				}
			}
			if ( !found ) { // no digits found with these wire and detector
				// ***** NEW CsDigit...
				int ich = (*it_dig).first;
				int pix_nb = pixmmobject->GetPixNb(ich);
				double pixX = pixmmobject->GetXPix(pix_nb);
				double pixY = pixmmobject->GetYPix(pix_nb);

				// Get data
				vector<double> data;
				data.push_back(max(amplitude[0]*weight + elec_spread_sigmabg*CsRandom::gauss(),0.));
				data.push_back(max(amplitude[1]*weight + elec_spread_sigmabg*CsRandom::gauss(),0.));
				data.push_back(max(amplitude[2]*weight + elec_spread_sigmabg*CsRandom::gauss(),0.));
				data.push_back(pixX);
				data.push_back(pixY);
				data.push_back((double)pix_nb);

				CsDigit *digit = new CsMCDigit(*this, ich, &data[0], data.size());
				dynamic_cast<CsMCDigit*>(digit)->addHit(*(*Ih));
				myDigits_.push_back( digit );  // ***** ADD IT TO LIST
			}
		} // End of loop on hits associated to current MC hit
	} // End of loop on MC hits


	// Fill histograms
	if ( hLevel_ >= High ) {
		list<CsDigit*>::iterator Id;
		for (Id = myDigits_.begin(); Id != myDigits_.end(); Id++) {
			double amp1 = (*Id)->getData()[0];
			double amp2 = (*Id)->getData()[1];
			double amp3 = (*Id)->getData()[2];
			string tbn  = GetTBName();
			if ( mH1[tbn+"_MCdA1"] != NULL )  mH1[tbn+"_MCdA1"]->Fill(amp1);
			if ( mH1[tbn+"_MCdA2"] != NULL )  mH1[tbn+"_MCdA2"]->Fill(amp2);
			if ( mH1[tbn+"_MCdA3"] != NULL )  mH1[tbn+"_MCdA3"]->Fill(amp3);
			if ( mH2[tbn+"_MCdAmpRation"] != NULL )  mH2[tbn+"_MCdAmpRatio"]->Fill(amp2/amp3,  amp1/amp3);
		}
	}
}


void CsPixelMumegaDetector::makeMCStripsSimple() {
	list<CsMCHit*>::iterator Ih;
	for ( Ih = myMCHits_.begin(); Ih != myMCHits_.end(); Ih++ ) { // loop on hits

		// Only charged particles or charged products
		if ( ( ( (*Ih)->getMCTrack() )->getParticle() )->getCharge() || (*Ih)->getOrigin() ) {

			CsMCTrkHit* thit = dynamic_cast<CsMCTrkHit*>(*Ih);
			if ( thit == 0 ) return;

			double t  = thit->getDTime(); // Delay time (ns)
			double ui = thit->getUin();   // Hit in point (DRS)
			double vi = thit->getVin();
			double wi = thit->getWin();
			double uo = thit->getUout();  // Hit out point (DRS)
			double vo = thit->getVout();
			double wo = thit->getWout();
			int    id = thit->GetID();

			// ***** TIME SMEARING :
			double tdc = t + tRes_*CsRandom::gauss();
			string tbn = GetTBName();
			if ( mH1[tbn+"_MCHitTime"]!=NULL)  mH1[tbn+"_MCHitTime"]->Fill(tdc); // MCHit time
			if ( fMCHitTMin > tdc || tdc > fMCHitTMax ) continue;
		
			// Find center of this hit's subdetector
			double xcm=xcm_, ycm=ycm_, zcm=zcm_;
			map<int, double>::iterator i;
			if((i=xcms_.find(id))!=xcms_.end())
			  xcm = i->second;
			//	else
			//	  cout<<"Warning : CsPixelMumegaDetector::makeMCDecoding : No xcm for subdetector with id "<<id<<endl;
			if((i=ycms_.find(id))!=ycms_.end())
			  ycm = i->second;
			//	else
			//	  cout<<"Warning : CsPixelMumegaDetector::makeMCDecoding : No ycm for subdetector with id "<<id<<endl;
			if((i=zcms_.find(id))!=zcms_.end())
			  zcm = i->second;
			//	else
			//	  cout<<"Warning : CsPixelMumegaDetector::makeMCDecoding : No zcm for subdetector with id "<<id<<endl;


			int err;
			HepMatrix iRotM(3,3); iRotM = rotWRS_.inverse( err );
			double xi = rotDRS_(1,1)*ui + rotDRS_(1,2)*vi + rotDRS_(1,3)*wi + xcm - _xshift;
			double yi = rotDRS_(2,1)*ui + rotDRS_(2,2)*vi + rotDRS_(2,3)*wi + ycm - _yshift;
			double zi = rotDRS_(3,1)*ui + rotDRS_(3,2)*vi + rotDRS_(3,3)*wi + zcm;
			double xo = rotDRS_(1,1)*uo + rotDRS_(1,2)*vo + rotDRS_(1,3)*wo + xcm - _xshift;
			double yo = rotDRS_(2,1)*uo + rotDRS_(2,2)*vo + rotDRS_(2,3)*wo + ycm - _yshift;
			double zo = rotDRS_(3,1)*uo + rotDRS_(3,2)*vo + rotDRS_(3,3)*wo + zcm;
			double Ui = iRotM(1,1)*xi + iRotM(1,2)*yi + iRotM(1,3)*zi; // WRS
			double Vi = iRotM(2,1)*xi + iRotM(2,2)*yi + iRotM(2,3)*zi;
			double Wi = iRotM(3,1)*xi + iRotM(3,2)*yi + iRotM(3,3)*zi;
			double Uo = iRotM(1,1)*xo + iRotM(1,2)*yo + iRotM(1,3)*zo; // WRS
			double Vo = iRotM(2,1)*xo + iRotM(2,2)*yo + iRotM(2,3)*zo;
			double Wo = iRotM(3,1)*xo + iRotM(3,2)*yo + iRotM(3,3)*zo;

			if ( mH2[tbn+"_MChMap"]!=NULL)  mH2[tbn+"_MChMap"]->Fill((Ui+Uo)/2.,(Vi+Vo)/2.); // Hit map

			int wireF; // First wire for this hit
			if ( ( Ui-wirD_ )/wirP_ < 0 ) {
				wireF = int( ( Ui-wirD_ )/wirP_ - 0.5 );
			} else {
				wireF = int( ( Ui-wirD_ )/wirP_ + 0.5 );
			}
			int wireL; // Last  wire for this hit
			if( ( Uo-wirD_ )/wirP_ < 0 ) {
				wireL = int( ( Uo-wirD_ )/wirP_ - 0.5 );
			} else {
				wireL = int( ( Uo-wirD_ )/wirP_ + 0.5 );
			}
			if( wireL < wireF ) {
				int tmp = wireL;
				wireL = wireF;
				wireF = tmp;
			}

			for( int i = wireF; i <= wireL; i++ ) {

				list<CsDigit*>::iterator Id;
				bool found = false;
				for ( Id = myDigits_.begin(); ( Id != myDigits_.end() ) && ( !found ); Id++ ) {

					if ( i == (*Id)->getAddress() ) {
						found = true;
						dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
					}
				}
				if ( !found ) { // No digits found with these wire and detector
					if ( i < 0 || i >= nWir_ ) {
						ostringstream ost;
						ost << "Unreliable wire number : " << i << " (0, " << nWir_ << "), "
								<< "detector : " << GetID() << " " << unit_ << " " << type_ << ".";
						CsErrLog::Instance()->mes( elAnomaly, ost.str() );
					} else {
						CsDigit* digit = new CsMCDigit( *this, i );
						dynamic_cast<CsMCDigit*>(digit)->addHit(*(*Ih));
						// Add this digit to my list
						myDigits_.push_back( digit );
					}
				}
			}
		}
	} // End of MCHits loop
}


void CsPixelMumegaDetector::makeMCStripsAmplitude() {
	list<CsMCHit*>::iterator Ih;
	for ( Ih = myMCHits_.begin(); Ih != myMCHits_.end(); Ih++ ) { // loop on hits

		if ( ( ( (*Ih)->getMCTrack() )->getParticle() )->getCharge() == 0 && (*Ih)->getOrigin() == 0 )
			continue;

		// ***** LOOP ON HITS from CHARGED PARTICLES or CHARGED PRODUCTS *****

		CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
		if ( thit == 0 ) return;

		double t  = thit->getDTime(); // Delay time (ns)
		double ui = thit->getUin();   // Hit in point (DRS)
		double vi = thit->getVin();
		double wi = thit->getWin();
		double uo = thit->getUout();  // hit out point (DRS)
		double vo = thit->getVout();
		double wo = thit->getWout();
		int    id = thit->GetID();


		// ***** TIMING : SMEARING, TIME CUT
		double tdc = t+tRes_*CsRandom::gauss();
		string tbn = GetTBName();
		if ( mH1[tbn+"_MCHitTime"]!=NULL)  mH1[tbn+"_MCHitTime"]->Fill(tdc); // MCHit time
		if ( fMCHitTMin > tdc || tdc > fMCHitTMax )
		  continue;
		
			// Find center of this hit's subdetector
			double xcm=xcm_, ycm=ycm_, zcm=zcm_;
			map<int, double>::iterator i;
			if((i=xcms_.find(id))!=xcms_.end())
			  xcm = i->second;
			//	else
			//	  cout<<"Warning : CsPixelMumegaDetector::makeMCDecoding : No xcm for subdetector with id "<<id<<endl;
			if((i=ycms_.find(id))!=ycms_.end())
			  ycm = i->second;
			//	else
			//	  cout<<"Warning : CsPixelMumegaDetector::makeMCDecoding : No ycm for subdetector with id "<<id<<endl;
			if((i=zcms_.find(id))!=zcms_.end())
			  zcm = i->second;
			//	else
			//	  cout<<"Warning : CsPixelMumegaDetector::makeMCDecoding : No zcm for subdetector with id "<<id<<endl;


		HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
		double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+ xcm - _xshift;
		double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+ ycm - _yshift;
		double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+ zcm;
		double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+ xcm - _xshift;
		double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+ ycm - _yshift;
		double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+ zcm;
		double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
		double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
		double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
		double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
		double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
		double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;
		double uShift = iRotM(1,1)*xcm_+iRotM(1,2)*ycm_;  // shift of plane wrt WRS, for cuts based on detector geometry
		double vShift = iRotM(2,1)*xcm_+iRotM(2,2)*ycm_;

		if ( mH2[tbn+"_MChMap"]!=NULL)  mH2[tbn+"_MChMap"]->Fill((Ui+Uo)/2.,(Vi+Vo)/2.); // Hit map

		double pathX = fabs((Ui - Uo)/2.);

		/// find the resolution associated with the MChit position (calculate the center of the cluster)
		double U_final = ((Ui+Uo)/2.);
		double V_final = ((Vi+Vo)/2.);
		bool error;
		int wire_final = DistToWire(U_final, error);
		if((wire_final < 0)||(wire_final > 895)){
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s: Wrong MC strip decoding : bad initial strip wire nb or not found (U_clust: %f, wire nb : %d)",
					GetTBName().c_str(),U_final, wire_final);
			continue;
		}
		if(wire_final < 256){
			U_final = U_final + strip_res2*CsRandom::gauss(); // WRS
			V_final = V_final + strip_res2*CsRandom::gauss();
		}
		else if(wire_final > 639){
			U_final = U_final + strip_res3*CsRandom::gauss(); // WRS
			V_final = V_final + strip_res3*CsRandom::gauss();
		}
		else{
			U_final = U_final + strip_res1*CsRandom::gauss(); // WRS
			V_final = V_final + strip_res1*CsRandom::gauss();
		}
		//////////////

		///// computation of ADC amplitude ///
		double amplitude[3];
		double Eloss = thit->getELos()*1.E6;
		getMCAmp(tdc, Eloss, amplitude);
		///////////////////////////

		//double clus_size = clus_size_mean + clus_size_sig*CsRandom::gauss(); //(mm)
		double mean = (clus_size_mean - exp(clus_size_mean_const + (clus_size_mean_slope*amplitude[2])));
		const double clus_size = max((mean*0.4 + clus_size_sig*CsRandom::gauss()),0.); //(mm)
		const double clus_rad = clus_size/2.; //cluster radius (cluster is supposed circular) (mm)

		const double U_clust_max = (U_final + clus_rad);
		const double U_clust_min = (U_final - clus_rad);
		const double V_clust_max = (V_final + clus_rad);
		const double V_clust_min = (V_final - clus_rad);

		/// border detector area TEST (under-range and over-range) ////
		if((U_clust_max <= -(getXsiz()*0.5) + uShift)|| (U_clust_min >= getXsiz()*0.5 + uShift)) continue;
		if((V_clust_max <= -(getYsiz()*0.5) + vShift)|| (V_clust_min >= getYsiz()*0.5 + vShift)) continue;

		double out_U_fact = 0.; //per cent of cluster charge outside the detector in U
		if(U_clust_min < -(getXsiz()*0.5) + uShift){
			out_U_fact += getErfFact(U_clust_min, U_final, -(getXsiz()*0.5) + uShift, clus_size);
		}
		if(U_clust_max > (getXsiz()*0.5) + uShift){
			out_U_fact += getErfFact((getXsiz()*0.5) + uShift, U_final, U_clust_max, clus_size);
		}
		/////////////

		/////// definition of wire number for last and first wire
		bool max_found;
		bool min_found;
		int wireF = DistToWire(U_clust_min, min_found);  // first wire for this hit
		int wireL = DistToWire(U_clust_max, max_found); // last  wire for this hit
		max_found = !max_found;//logic for found is inverse than for error
		min_found = !min_found;//logic for found is inverse than for error

		if((wireF < 0)||(wireF > 895)){ //check the number of the strip found
			min_found = false;
		}
		if((wireL < 0)||(wireL > 895)){ //check the number of the strip found
			max_found = false;
		}

		if((!max_found) && (!min_found)){
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s: Wrong MC strip decoding : strip wires nb not found or incorrect (U_min : %f, U_max : %f)",
					GetTBName().c_str(), U_clust_min, U_clust_max);
			continue;
		}
		if((!max_found) && min_found) wireL = wireF; //avoid some bug due to over (under) -range of the cluster position
		if((!min_found) && max_found) wireF = wireL; //avoid some bug due to over (under) -range of the cluster position
		///////////////////////////

		if ( wireL < wireF ) {//not really needed but in case of error (last wire must be higher than first one).
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s: MC strip decoding wrong strip value : First wire : %d , Last  wire : %d",
					GetTBName().c_str(), wireF, wireL);
			int tmp = wireL;
			wireL = wireF;
			wireF = tmp;
		}

		// find the hemisphere(s) hited
		unsigned int hemneg(0), hempos(0);
		if ( V_clust_max > 0. ) hempos=1;
		if ( V_clust_min < 0. ) hemneg=1;
		////////////////////////////


		// table of weight for amplitude repartition (-1 = strip not selected)
		const int nb_hem = (hempos + hemneg);
		const int nb_wire = (wireL - wireF + 1);
		const int nb_digit = nb_wire*nb_hem;
		double W_tab[nb_digit];
		double W_total = 0.0;
		for (int i = 0; i<nb_digit; i++){
			W_tab[i] = -1.0;
		}
		//////////////

		//check if the strip is in the cluster area and attribute it a weight it take into account
		//the over and under range in V and pixel area (and in the end the over and under range in U)
		for (int iu = wireF; iu <= wireL; iu++ ){
			//calculate the weight for all this strip position (for taking account the charge outside the detector in V
			//and charge in pixel area)
			W_total += weight_for_strip(iu,U_final,V_final,uShift,vShift,clus_size,pathX,3);
			if (hempos){
				if (dist_striptopoint(iu,U_final,V_final,uShift,vShift,1) < clus_rad){
					W_tab[(iu - wireF)] = weight_for_strip(iu,U_final,V_final,uShift,vShift,clus_size,pathX,1);
				}
			}
			if (hemneg){
				if (dist_striptopoint(iu,U_final-uShift,V_final-vShift,uShift,vShift,-1) < clus_rad){
					W_tab[(iu - wireF) + (nb_hem - 1)*nb_wire] = weight_for_strip(iu,U_final,V_final,uShift,vShift,clus_size,pathX,-1);
				}
			}
		}
		W_total = (W_total/(1. - out_U_fact)); //take into account the percent of charge outside the detector in U
		///////

		map<pair<int,int>, double> digitmap; //first is nb and hemisphere of strip, and second is the weight for energy loss
		map<pair<int,int>, double>::iterator it_dig;
		pair<int,int> ref_strip;
		int hem;
		int orientation = 1;
		double W_tmp;
		/// orientation convention
		if ( GetTBName()[4] == 'X' ) {
			orientation = 1;
		}
		else if ( GetTBName()[4] == 'Y' ) { //for keeping the hem -1 on saleve and hem 1 on jura
			orientation = -1;
		}
		else if ( GetTBName()[4] == 'U' ) { //for keeping the hem -1 on saleve and hem 1 on jura
			orientation = -1;
		}
		else if ( GetTBName()[4] == 'V' ) {
			orientation = 1;
		}
		else {
			CsErrLog::msg(elFatal, __FILE__, __LINE__,
					"%s: MC strip decoding detector orientation not implemented !!!",
					GetTBName().c_str());
		}
		/////

		//fusion the hemisphere for hemisphere 0 and weight affectation for each strip///
		for (int iu = wireF; iu <= wireL; iu++ ){
			//find if the wire is hemisphere 0 or not
			if ((iu < 256) || (iu > 639)){ //long strips
				hem = 0;
			}
			else {//short strips
				hem = 1*orientation;
			}
			//////////
			if (hempos){
				W_tmp = W_tab[(iu - wireF)];
				if(W_tmp > 0.){
					ref_strip = make_pair(iu, hem);
					if (digitmap.count(ref_strip) > 0){ // we fusion the hemisphere 0
						digitmap[ref_strip] += (W_tmp/W_total);
					}
					else{
						digitmap[ref_strip] = (W_tmp/W_total);
					}
				}
			}
			if (hemneg){
				W_tmp = W_tab[(iu - wireF) + (nb_hem - 1)*nb_wire];
				if(W_tmp > 0.){
					ref_strip = make_pair(iu, -hem);
					if (digitmap.count(ref_strip)> 0){ // we fusion the hemisphere 0
						digitmap[ref_strip] += (W_tmp/W_total);
					}
					else{
						digitmap[ref_strip] = (W_tmp/W_total);
					}
				}
			}
		}
		///////

		if(digitmap.size() == 0) continue; //check if we have at least one digit

		// ****** DOES A CsDigit ALREADY EXIST ON THIS STRIP ? ******
		pair<int,int> ref_found;
		list<CsDigit*>::iterator Id;
		bool found;
		double weight = 0.;

		for (it_dig = digitmap.begin(); it_dig != digitmap.end(); it_dig++){
			weight = (*it_dig).second;
			found = false;
			for ( Id = myDigits_.begin(); Id != myDigits_.end(); Id++ ) {
				ref_found = make_pair((*Id)->getAddress(),int((*Id)->getData()[3]));
				if ((*it_dig).first == ref_found) { // Here it is...
					found = true;                 // ... ADD THIS HIT TO EXISTING CsDigit...
					dynamic_cast<CsMCDigit*>(*Id)->addHit(*(*Ih));
					for ( unsigned int i = 0; i < 3; i++ ){
						(*Id)->getData()[i] += max((amplitude[i]*weight + elec_spread_sigmabg*CsRandom::gauss()),0.);
					}
					break;
				}
			}
			if ( !found ) { // no digits found with these wire and detector
				// ***** NEW CsDigit...
				int iwire = ((*it_dig).first).first;
				hem = ((*it_dig).first).second;

				// Get data
				vector<double> data;
				data.push_back(max(amplitude[0]*weight + elec_spread_sigmabg*CsRandom::gauss(),0.));
				data.push_back(max(amplitude[1]*weight + elec_spread_sigmabg*CsRandom::gauss(),0.));
				data.push_back(max(amplitude[2]*weight + elec_spread_sigmabg*CsRandom::gauss(),0.));
				data.push_back(hem);

				CsDigit *digit = new CsMCDigit(*this, iwire, &data[0], data.size());
				dynamic_cast<CsMCDigit*>(digit)->addHit(*(*Ih));
				myDigits_.push_back( digit );  // ***** ADD IT TO LIST
			}
		} // End of loop on hits associated to current MC hit
	} // End of loop on MC hits


	// Fill histograms
	if ( hLevel_ >= High ) {
		list<CsDigit*>::iterator Id;
		for (Id = myDigits_.begin(); Id != myDigits_.end(); Id++) {
			double amp1 = (*Id)->getData()[0];
			double amp2 = (*Id)->getData()[1];
			double amp3 = (*Id)->getData()[2];
			string tbn  = GetTBName();
			if ( mH1[tbn+"_MCdA1"] != NULL )  mH1[tbn+"_MCdA1"]->Fill(amp1);
			if ( mH1[tbn+"_MCdA2"] != NULL )  mH1[tbn+"_MCdA2"]->Fill(amp2);
			if ( mH1[tbn+"_MCdA3"] != NULL )  mH1[tbn+"_MCdA3"]->Fill(amp3);
			if ( mH2[tbn+"_MCdAmpRation"] != NULL )  mH2[tbn+"_MCdAmpRatio"]->Fill(amp2/amp3,  amp1/amp3);
			//if ( mH1[tbn+"_MChPos"]!=NULL)  mH1[tbn+"_MChPos"]->Fill(channel); // Hit position
		}
	}
}


istream & operator >> (istream &in, vector<CsPixelMumegaDetector::APVCal> &c) {

	c.clear();

	try {

		unsigned chip = 9999999;
		while(1) {

			char s[111];
			in.getline(s, sizeof(s)); // read line with information about chips
			if ( !in )
				break; // end of data stream

			if ( s[0] == '#' ) {

				char s2[sizeof(s)];
				int a, src_id, adc_id, chip_id;
				if ( 5!=sscanf(s, "%s %d %d %d %d", s2, &a, &src_id, &adc_id, &chip_id) )
					throw CS::Exception("CsPixelMumegaDetector::APVCal::operator>>: bad string \"%s\".", s);
				chip--;

				CsPixelMumegaDetector::APVCal cal;
				cal.src_id = src_id;
				cal.adc_id = adc_id;
				cal.chip_id = chip_id;

				for ( int i=0; i<CsPixelMumegaDetector::ChipChannels; i++ )
				{
					in.getline(s, sizeof(s)); // read line with information about chips
					if ( !in )
						throw CS::Exception("CsPixelMumegaDetector::APVCal::operator>>: Error in reading data.");
					cal.channels.push_back(CsPixelMumegaDetector::APVCal::Channel(s));
				}

				c.push_back(cal);
			}
			else
				throw CS::Exception("CsPixelMumegaDetector::APVCal::operator>> : format error.");
		}
	}
	catch(...)
	{
		cerr << "CsPixelMumegaDetector::APVCal::operator>> error!\n";
		throw;
	}

	in.clear(); // Why should I do this? I don't know...   Neither do I...

	return in;
}


void CsPixelMumegaDetector::readCalibration(time_t timePoint){
	// Clear all existing calibrations
	if ( stripplane ) stripplane->ClearChannels();
	else if ( pixelplane ) pixelplane->ClearChannels();
	else if ( rectpixelplane ) rectpixelplane->ClearChannels();
	apv_chan_cals.clear();
	time_cals.Clear();
	apv_chan_cals_mapped = false;

	CDB::Time tp(timePoint, 0);
	tm *t = localtime( &tp.first );

	// Read calibration data for each strip
	string strdata("");
	cdb_->read(GetTBName(), strdata, tp);
	istringstream istrdata(strdata);
	istrdata >> apv_chan_cals;

	// Check if calibrations for correct number of chips was found
	bool apv_chan_cals_ok = true;
	if ( (int)apv_chan_cals.size() == ChipsPerPlane ){
		CsErrLog::msg(elInfo, __FILE__, __LINE__,
				"%s calibration data found in CDB for local time %s",
				GetTBName().c_str(), asctime(t));

		// Check if correct number of channels for each chip
		vector<CsPixelMumegaDetector::APVCal>::iterator itcal;
		for ( itcal = apv_chan_cals.begin(); itcal != apv_chan_cals.end(); itcal++ ) {
			vector<CsPixelMumegaDetector::APVCal::Channel>& vchn = (*itcal).channels;
			if ( vchn.size() != ChipChannels ) {
				CsErrLog::msg(elError, __FILE__, __LINE__,
						"%s wrong calibration data size (# channels) in CDB for local time %s", GetTBName().c_str(), asctime(t));
				apv_chan_cals_ok = false;
			}
		}

		// Debugging output
		if ( getenv("COND_DB_DEBUG") != 0 || getenv("COND_DB_DEBUG_MUMEGA") != 0 ) {
			cout << GetTBName() <<endl;
			cout << apv_chan_cals.size() << endl;
			vector<CsPixelMumegaDetector::APVCal>::iterator itcal;
			vector<CsPixelMumegaDetector::APVCal::Channel>::iterator itchn;
			for ( itcal = apv_chan_cals.begin(); itcal != apv_chan_cals.end(); itcal++ ) {
				cout<<"next APVCal structure:\n";
				vector<CsPixelMumegaDetector::APVCal::Channel>& vchn = (*itcal).channels;
				for ( itchn = vchn.begin(); itchn != vchn.end(); itchn++ ) {
					(*itchn).Print();
				}
			}
			cout<<endl;
		}
	}

	// Calibration data not found for all chips
	else {
		CsErrLog::msg(elError, __FILE__, __LINE__,
				"%s wrong calibration data size (# chips) in CDB for local time %s",
				GetTBName().c_str(), asctime(t));
		apv_chan_cals_ok = false;
	}

	// Calibrations not found for every strip of plane
	if ( !apv_chan_cals_ok ) apv_chan_cals.clear();


	// Calibration for time, cluster position correction and cross talk in time
	Bool_t time_cal_ok = false, pos_corr_cal_ok = false, txt_cal_ok = false, amp_cut_cal_ok = false;
	unsigned int n_params_time_cal = ChipsPerPlane, n_params_pos_corr_cal = 16, n_params_txt_cal = ChipsPerPlane;


	// Read calibration data for conversion of amplitude->time for each plane
	string strdata2("");
	cdb_->read(GetTBName(), strdata2, tp, "timing_per_chip");
	istringstream istrdata2(strdata2);
	std::vector<Float_t> f_time_cal_offset, f_time_cal_slope, f_time_cal_time_offset ;
	if (strdata2 == "") {
	  time_cal_ok = false;
	  CsErrLog::msg(elError, __FILE__, __LINE__,
	                  "%s no timing_per_chip calibration found in CDB for local time %s",
	                  GetTBName().c_str(), asctime(t));
	} else {
	  string tmp_os = "", tmp_sl = "", tmp_tos = "";
	  while (!istrdata2.eof()) {
	    istrdata2 >> tmp_os >> tmp_sl >> tmp_tos;
	    f_time_cal_offset.push_back(atof(tmp_os.c_str()));
	    f_time_cal_slope.push_back(atof(tmp_sl.c_str()));
	    f_time_cal_time_offset.push_back(atof(tmp_tos.c_str()));
	  }
	
	  // std::cout<<"sizes : "<<s_time_cal_offset.size()<<" "<<s_time_cal_slope.size()<<" "<<s_time_cal_time_offset.size()<<std::endl;

	  if ( (f_time_cal_offset.size() != n_params_time_cal) || (f_time_cal_slope.size() != n_params_time_cal) || (f_time_cal_time_offset.size() != n_params_time_cal) ){  //test on string objects to detect double dec errors in file (ex : 0.21.4)
	    time_cal_ok = false;
	    CsErrLog::msg(elError, __FILE__, __LINE__,
	                    "%s wrong number of parameters for time calibration",
	                    GetTBName().c_str());
	  } else time_cal_ok = true;
	}

	// Read calibration data for cluster position correction
	string strdata3("");
	cdb_->read(GetTBName(), strdata3, tp, "pos_corr");
	istringstream istrdata3(strdata3);
	std::vector<string> s_pos_corr;
	std::vector<Float_t> f_pos_corr;

	if(strdata3 == ""){
		pos_corr_cal_ok = false;
		CsErrLog::msg(elInfo, __FILE__, __LINE__,
				"%s no calibration for correction of clusters position found in CDB for local time %s",
				GetTBName().c_str(), asctime(t));
	}
	else{
		string tmp = "";
		while (!istrdata3.eof()){
			istrdata3 >> tmp;
			s_pos_corr.push_back(tmp);
			f_pos_corr.push_back(atof(tmp.c_str()));
		}
		tmp = "";
		if (s_pos_corr.size() != n_params_pos_corr_cal){  //test on string objects to detect double dec errors in file (ex : 0.21.4)
			pos_corr_cal_ok = false;
			CsErrLog::msg(elInfo, __FILE__, __LINE__,
					"%s wrong number of calibration parameters for correction of clusters position",
					GetTBName().c_str());
		}
		else pos_corr_cal_ok = true;
	}

	// Read calibration data for crosstalk in time parameters for all planes
	string strdata4("");
	cdb_->read(GetTBName(), strdata4, tp, "crosstalk_in_time");
	istringstream istrdata4(strdata4);
	std::vector<string> s_param_txt_i, s_param_txt_ip1;
	std::vector<Float_t> f_param_txt_i, f_param_txt_ip1;
	if(strdata4 == "") {
		txt_cal_ok = false;
		CsErrLog::msg(elError, __FILE__, __LINE__,
				"%s no calibration for correction of cross-talk in time found in CDB for local time %s",
				GetTBName().c_str(), asctime(t));
	}
	else{
		string tmpi = "", tmpip1 = "";
		while (!istrdata4.eof()){
			istrdata4 >> tmpi >> tmpip1;
			s_param_txt_i.push_back(tmpi);
			s_param_txt_ip1.push_back(tmpip1);
			f_param_txt_i.push_back(atof(tmpi.c_str()));
			f_param_txt_ip1.push_back(atof(tmpip1.c_str()));
		}
		tmpi = "";
		tmpip1 = "";
		if ((s_param_txt_i.size() != n_params_txt_cal) || (s_param_txt_ip1.size() != n_params_txt_cal) ){  //test on string objects to detect double dec errors in file (ex : 0.21.4)
			txt_cal_ok = false;
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s wrong number of calibration parameters for correction of cross talk in time",
					GetTBName().c_str());
		}
		else txt_cal_ok = true;
	}

	
	// Read calibration data for Amplitude Ratio cuts for all planes
	string strdata5("");
	cdb_->read(GetTBName(), strdata5, tp, "amplitude_ratio_cuts");
	istringstream istrdata5(strdata5);
	std::vector<string> s_a1a2, s_a0a2, s_a1a2tcs, s_tcs;
	std::vector<Float_t> f_a1a2, f_a0a2, f_a1a2tcs, f_tcs ;
	if(strdata5 == "") {
	  amp_cut_cal_ok = false;
	  CsErrLog::msg(elError, __FILE__, __LINE__,
			"%s no calibration for cuts on cluster amplitude ratios found in CDB for local time %s",
			GetTBName().c_str(), asctime(t));
	}
	else{
	  string tmpsa0a2a1a2 = "", tmpsa1a2tcs = "", tmpa1a2 = "", tmpa0a2 = "", tmpa1a2tcs = "", tmptcs = "";
	  istrdata5 >> tmpsa0a2a1a2 >> tmpsa1a2tcs;
	  int size_a0a2a1a2 = atoi(tmpsa0a2a1a2.c_str());
	  int size_a1a2tcs = atoi(tmpsa1a2tcs.c_str());
	  tmpsa0a2a1a2 = ""; tmpsa1a2tcs = "";

	  //while (!istrdata5.eof()){
	  for (int i = 0;i<size_a0a2a1a2;i++){
	    istrdata5 >> tmpa1a2 >> tmpa0a2;
	    s_a1a2.push_back(tmpa1a2);
	    s_a0a2.push_back(tmpa0a2);
	    f_a1a2.push_back(atof(tmpa1a2.c_str()));
	    f_a0a2.push_back(atof(tmpa0a2.c_str()));
	    }
	  tmpa1a2 = "";
	  tmpa0a2 = "";

	  for (int i = 0;i<size_a1a2tcs;i++){
	    istrdata5 >> tmptcs >> tmpa1a2tcs;
	    s_a1a2tcs.push_back(tmpa1a2tcs);
	    s_tcs.push_back(tmptcs);
	    f_a1a2tcs.push_back(atof(tmpa1a2tcs.c_str()));
	    f_tcs.push_back(atof(tmptcs.c_str()));
	  }
	  tmpa1a2tcs = "";
	  tmptcs = "";

	  if ( (s_a1a2.size() != s_a0a2.size()) || (s_a1a2tcs.size() != s_tcs.size()) ){  //test on string objects to detect double dec errors in file (ex : 0.21.4)
	    amp_cut_cal_ok = false;
	    CsErrLog::msg(elError, __FILE__, __LINE__,
			  "%s wrong number of calibration parameters for cluster amplitude ratios cuts",
			  GetTBName().c_str());
	  }
	  else amp_cut_cal_ok = true;
	}
	

	if (stripplane) {
	  if (time_cal_ok) {
	    stripplane->GetPar()->EnableTimeCal();
	    stripplane->GetPar()->SetTimeCal(f_time_cal_offset, f_time_cal_slope, f_time_cal_time_offset);
	  } else stripplane->GetPar()->TimeCalDefault();
	  if (txt_cal_ok) {
	    stripplane->GetPar()->EnableTxtCal();
	    stripplane->GetPar()->SetTimeCrossTalkHitRatioi(f_param_txt_i);
	    stripplane->GetPar()->SetTimeCrossTalkHitRatioip1(f_param_txt_ip1);
	  } else stripplane->GetPar()->TxtCalDefault();
	  if (amp_cut_cal_ok && do_amplituderatiocuts) {
	    stripplane->GetPar()->EnableAmpCuts();
	    stripplane->GetPar()->SetAmpRatioCutsParams(f_a1a2, f_a0a2, f_tcs, f_a1a2tcs);//std::cout<<"OK cal amp strips"<<std::endl;
          } else stripplane->GetPar()->DisableAmpCuts();
	}
	else if (rectpixelplane) {
	  if (time_cal_ok) {
	    rectpixelplane->GetPar()->EnableTimeCal();
	    rectpixelplane->GetPar()->SetTimeCal(f_time_cal_offset, f_time_cal_slope, f_time_cal_time_offset);
	  } else rectpixelplane->GetPar()->TimeCalDefault();
	  if (pos_corr_cal_ok) {
	    rectpixelplane->GetPar()->EnablePosCorrCal();
	    rectpixelplane->GetPar()->SetPosCorrs(f_pos_corr);
	  } else rectpixelplane->GetPar()->PosCorrCalDefault();
	  if (txt_cal_ok) {
	    rectpixelplane->GetPar()->EnableTxtCal();
	    rectpixelplane->GetPar()->SetTimeCrossTalkHitRatioi(f_param_txt_i);
	    rectpixelplane->GetPar()->SetTimeCrossTalkHitRatioip1(f_param_txt_ip1);
	  } else rectpixelplane->GetPar()->TxtCalDefault();
	  if (amp_cut_cal_ok && do_amplituderatiocuts) {
	    rectpixelplane->GetPar()->EnableAmpCuts();
	    rectpixelplane->GetPar()->SetAmpRatioCutsParams(f_a1a2, f_a0a2, f_tcs, f_a1a2tcs);//std::cout<<"OK cal amp pixels"<<std::endl;
	  } else rectpixelplane->GetPar()->DisableAmpCuts();
	}
	
	// Just debug printouts
	if ( getenv("COND_DB_DEBUG") != 0 || getenv("COND_DB_DEBUG_MUMEGA") != 0) {
		cout << endl << "Timing calibrations read in for " << GetTBName() << endl;
		if ( stripplane ) {
			cout<<"const : "<<stripplane->GetPar()->GetTimeCalOffset()[0]<<" slope : "<<stripplane->GetPar()->GetTimeCalSlope()[0]<<endl;
		}
		/*    else if ( pixelplane ) {
      cout<<"const : "<<pixelplane->GetPar()->GetTimeCal()[0]<<" slope : "<<pixelplane->GetPar()->GetTimeCal()[1]<<endl;
      }*/
		else if ( rectpixelplane ) {
			cout<<"const : "<<rectpixelplane->GetPar()->GetTimeCalOffset()[0]<<" slope : "<<rectpixelplane->GetPar()->GetTimeCalSlope()[0]<<endl;
		}
	}

	// Read calibration data for correction in U
	if (corr_U_bool){
		string str_corrU("");
		int binX,binY;
		int tab_size = 0;
		double value;
		int indice;
		size_t ptr, ptr2;
		std::string line;
		int cpt = 1;
		cdb_->read(GetTBName(), str_corrU, tp, "corr_U");
		istringstream istr_corrU(str_corrU);
		if(!istr_corrU.fail()){

			while(getline(istr_corrU,line)){
				if ((line[0] == '#') || (line.length() == 0)) continue;

				std::istringstream is(line);

				if(cpt == 1){
					is >> corr_U_shift_u;
					is >> corr_U_shift_v;
					is >> corr_U_phi;
					cpt++;
					continue;
				}
				if(cpt == 2){
					is >> corr_U_size_u;
					is >> corr_U_size_v;
					cpt++;
					continue;
				}
				if(cpt == 3){
					is >> corr_U_nb_u;
					is >> corr_U_nb_v;
					tab_size = corr_U_nb_u * corr_U_nb_v;

					if(corr_U_tab != NULL) delete [] corr_U_tab;
					corr_U_tab = new double[tab_size];
					cpt++;
					continue;
				}

				if (cpt > 3){
					is >> binX;
					is >> binY;
					is >> value;

					indice = ((binX-1) + (binY-1)*corr_U_nb_u);
					if ((indice < tab_size) && (indice >= 0)){
						corr_U_tab[indice] = value;
						value = 0.;
					}
					else {
						corr_U_bool = false;
						CsErrLog::msg(elError, __FILE__, __LINE__,
								"%s wrong calibration data in CDB for cluster position correction (out of range) %s : correction in U will NOT be applied!", GetTBName().c_str(), asctime(t));
					}
				}
			}
			if(cpt < 3){
				corr_U_bool = false;
				CsErrLog::msg(elError, __FILE__, __LINE__,
						"%s wrong calibration data in CDB for cluster position correction (bad file) %s, correction in U will NOT be applied!", GetTBName().c_str(), asctime(t));
			}
			else{

				corr_U_angle = (ang_ * (M_PI) / 180.) +  corr_U_phi;
				corr_U_cos_angle = cos(corr_U_angle);
				corr_U_sin_angle = sin(corr_U_angle);
				corr_U_center_u = (getXcm()*corr_U_cos_angle + getYcm()*corr_U_sin_angle)*0.1 + corr_U_shift_u;
				corr_U_center_v = (- getXcm()*corr_U_sin_angle + getYcm()*corr_U_cos_angle)*0.1 + corr_U_shift_v;
			}

		}
		else {
			corr_U_bool = false;
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s wrong calibration data in CDB for cluster position correction (read failled) %s, correction in U will NOT be applied!", GetTBName().c_str(), asctime(t));
		}

		// Just debug printouts
		if ( getenv("COND_DB_DEBUG") != 0 || getenv("COND_DB_DEBUG_MUMEGA") != 0 ) {
			cout << endl << "cluster's position correction read in for " << GetTBName() << endl;
			cout << "tab_size : "<< tab_size << endl;
			for ( int i=0; i<corr_U_nb_v; i++) {
				for ( int j=0; j<corr_U_nb_u; j++) {
					cout << (j+1) << "\t" << (i+1) << "\t" << corr_U_tab[i + j*corr_U_nb_u] << endl;
				}
			}
		}
	}
}


bool CsPixelMumegaDetector::mapChannelCal() {
	// All decoding maps
	const CS::Chip::Maps &daq_maps =
			CsInit::Instance()->getDaqMaps();

	// Create a short name for iterator through map
	typedef CS::Chip::Maps::const_iterator m_it;

	if ( daq_maps.size() != 0 ) {

		// Apply mapping to calibrations for this plane
		unsigned  mumega_chan_cals_ok = 0;
		if ( apv_chan_cals.size() != 0 ) {

			// Loop over all chips and get calibration data for each chip
			for ( unsigned chip = 0; chip < apv_chan_cals.size(); chip++ ) {
				APVCal &c = apv_chan_cals[chip];

				// Get Source ID and ADC ID and Chip ID of chip
				int src_id = c.src_id, adc_id = c.adc_id, chip_id = c.chip_id;
				int wire = -9999; int hem = 0;
				int pixX = -1; int pixY = -1;
				float pixXMM = -1; float pixYMM = -1;
				int conn = -1, conn09 = -1, pixnb = -1, chandet = -1, stripconn = -1;
				float amplfactor = 1;

				// Loop over all channels of chip
				for ( unsigned chip_chan = 0; chip_chan < c.channels.size(); chip_chan++ ) {

					// All maps with given data ID
					CS::ChipAPV::DataID data_id(src_id, adc_id, chip_id, chip_chan);

					const pair<m_it,m_it> m_range =
							daq_maps.equal_range(data_id);

					// Read map
					bool mapped = false;
					for ( m_it cc=m_range.first; cc!=m_range.second; cc++ ) {
						const CS::ChipAPV::Digit *digit1 =
								dynamic_cast<CS::ChipAPV::Digit*>(cc->second);

						if ( digit1==NULL ){
							CsErrLog::msg(elWarning, __FILE__, __LINE__,
									"%s: ChipAPV wrong map!", GetTBName().c_str());
							continue;
						}

						if ( GetTBName()!=digit1->GetDetID().GetName() ){
							CsErrLog::msg(elWarning, __FILE__, __LINE__,
									"%s: Inconsistency between mapping file and calibration file!", GetTBName().c_str());
							CsErrLog::msg(elWarning, __FILE__, __LINE__,
									"%s from calibration file does not match TBName",
									digit1->GetDetID().GetName().c_str());
							continue;
						}

						// Map found

						// Check if pixel or rect pixel digit
						const CS::ChipAPV::DigitPixel *digit1p =
								dynamic_cast<const CS::ChipAPV::DigitPixel*>(digit1);
						const CS::ChipAPV::DigitPixelMM *digit1pmm =
								dynamic_cast<const CS::ChipAPV::DigitPixelMM*>(digit1);

						// Strip digit
						if ( digit1p == NULL && digit1pmm == NULL ) {

							// Apply mapping
							try {
								wire = digit1->GetChannel();
								hem = digit1->GetChanPos();
							}
							catch(...) {
								CsErrLog::msg(elError, __FILE__, __LINE__,
										"%s: Could not apply mapping ! Can't access to GetChannel() or GetChanPos()",
										GetTBName().c_str());
								continue;
							}
							
							amplfactor = digit1->GetAmplFactor();
							stripconn = digit1->GetStripConnNb();

							mapped = true;
							break;
						}

						// Pixel digit
						else if ( digit1p ) {

							// Apply mapping
							try {
								wire = digit1p->GetChannel();
							}
							catch(...) {
								CsErrLog::msg(elError, __FILE__, __LINE__,
										"%s: Could not apply mapping ! Can't access to GetChannel()",
										GetTBName().c_str());
								continue;
							}

							// Get pixel coordinates
							pair<int,int> xy = pixgem::detch2xy(wire);
							pixX = xy.first;
							pixY = xy.second;

							// Pixel coordinates invalid
							if ( pixX == -1 && pixY == -1 )
							{
								CsErrLog::msg(elError, __FILE__, __LINE__,
										"%s: Could not apply mapping ! Pixel coordinates invalid (-1)",
										GetTBName().c_str());
								continue;
							}

							// Invert X pads if detector orientation negative (det facing upstream)
							if( digit1p->GetDetOrientation() < 0 )
							{
								pixX = 31 - pixX;
							}
							mapped = true;
							break;
						}

						// Rectangular pixel digit
						else if ( digit1pmm ) {

							// Apply mapping
							try {
								wire = digit1pmm->GetChannel();
							}
							catch(...) {
								CsErrLog::msg(elError, __FILE__, __LINE__,
										"%s: Could not apply mapping ! Can't access to GetChannel()",
										GetTBName().c_str());
								continue;
							}

							// Get pixel infos
							conn = digit1pmm->GetConnNb();
                                                        register int mapPixelVersion = digit1pmm->GetPixelMMversion();
                                                        if (!pixmmobject) pixmmobject = pixmm::Instance(mapPixelVersion);
							pixnb = pixmmobject->GetPixNb(conn,wire);
                                                        conn09 = pixmmobject->GetNConn(conn); // convert from 0-19 to 0-9 numbering

							// Pixel coordinates invalid
							if ( pixnb == -1 ) {
								CsErrLog::msg(elError, __FILE__, __LINE__,
										"%s: Could not apply mapping ! Pixel number invalid (-1)",
										GetTBName().c_str());
								continue;
							}
							pixXMM = pixmmobject->GetXPix(pixnb);
							pixYMM = pixmmobject->GetYPix(pixnb);
							chandet = wire + conn*ChipChannels;
							if( digit1pmm->GetDetOrientation() < 0 ) {
								pixXMM = -pixXMM;
							}

							amplfactor = digit1pmm->GetAmplFactor();

							mapped = true;
							break;
						}

					} // Loop over maps

					// Calibrations mapped
					if ( mapped ){

						mumega_chan_cals_ok++;

						// Create strip channel with calibrations
						if ( stripplane ) {
							stripplane->AddChan(wire, hem, c.channels[chip_chan].flag,
									c.channels[chip_chan].pedestal_mean,
									c.channels[chip_chan].pedestal_sigma*amplfactor,
									c.chip_id, chip_chan, stripconn);
// cerr<<"mapChannelCal "<<GetTBName()<<" stripconn "<<stripconn<<" chip "<<chip<<" wire "<<wire<<endl;
						}

						// Create pixel channel with calibrations
						else if ( pixelplane ) {
							pixelplane->AddChan(wire, pixX, pixY, c.channels[chip_chan].flag,
									c.channels[chip_chan].pedestal_mean,
									c.channels[chip_chan].pedestal_sigma,
									c.chip_id, chip_chan);
						}

						// Create rectangular pixel channel with calibrations
						else if ( rectpixelplane ) {
							rectpixelplane->AddChan(chandet, pixnb, pixXMM, pixYMM, c.channels[chip_chan].flag,
									c.channels[chip_chan].pedestal_mean,
									c.channels[chip_chan].pedestal_sigma*amplfactor,
										c.chip_id, chip_chan, conn09);
						}

					}
					else {
						CsErrLog::msg(elWarning, __FILE__, __LINE__,
								"%s: Mapping failed for ChipId: %i, ChipChannel: %i!",
								GetTBName().c_str(), chip_id, chip_chan);
					}
				} // Loop over chip channels
			} // Loop over chips
		} // Mapping of APV calibrations data

		// Apply default calibrations
		if ( mumega_chan_cals_ok != daq_maps.GetWires( DetID(GetTBName()) ) ) {

		  //testest
		  //std::cout<<GetTBName()<<" "<<mumega_chan_cals_ok<<" "<<daq_maps.GetWires( DetID(GetTBName()) )<<std::endl;
		  
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s: Using default calibrations!", GetTBName().c_str());

			// Delete any existing calibrations
			if ( stripplane ) {
				if ( stripplane->GetNchannel() != 0 ) {
					stripplane->ClearChannels();
				}
			}
			// Delete any existing calibrations
			if ( pixelplane ) {
				if ( pixelplane->GetNchannel() != 0 ) {
					pixelplane->ClearChannels();
				}
			}
			// Delete any existing calibrations
			if ( rectpixelplane ) {
				if ( rectpixelplane->GetNchannel() != 0 ) {
					rectpixelplane->ClearChannels();
				}
			}

			set<uint16> srcIDs;
			daq_maps.GetSrcIDs( DetID( GetTBName() ), srcIDs );
			if ( srcIDs.size() != 1 ) {
				CsErrLog::msg(elFatal, __FILE__, __LINE__,
						"%s: One plane must be connected to only one source ID!",
						GetTBName().c_str() );
			}
			uint16 srcID = *srcIDs.begin();

			// All maps with given data ID
			CS::ChipAPV::DataID start_id (srcID,   0, 0, 0);
			CS::ChipAPV::DataID finish_id(srcID+1, 0, 0, 0);

			const pair<m_it, m_it> m_start  = daq_maps.equal_range(start_id);
			const pair<m_it, m_it> m_finish = daq_maps.equal_range(finish_id);

			// Fill calibrations for all wires of plane
			for ( m_it cc = m_start.first; cc != m_finish.second; cc++ ) {
				const CS::ChipAPV::Digit *digit1 =
						dynamic_cast<CS::ChipAPV::Digit*>(cc->second);

				if ( digit1 == NULL ) {
					CsErrLog::msg(elWarning, __FILE__, __LINE__,
							"%s: ChipAPV wrong map!", GetTBName().c_str() );
					continue;
				}

				if ( GetTBName() != digit1->GetDetID().GetName() )
					continue;

				int flag = 1;
				float ped = 750.;
				float sigma = 4./*0.0*/; //test 6/01/12
				int wire, hem, PX, PY;
				float PXMM = -1; float PYMM = -1;
				int conn = -1, conn09 = -1, pixnb = -1, chandet = -1, chipid = -1, chipchan = -1;
				if ( stripplane ) {
					try {
						wire = digit1->GetChannel();
						hem  = digit1->GetChanPos();
						chipid = digit1->GetChip();
						chipchan = digit1->GetChipChannel();
					}
					catch (...) {
						CsErrLog::msg(elError, __FILE__, __LINE__,
								"%s: Could not apply mapping!", GetTBName().c_str());
						continue;
					}
					stripplane->AddChan(wire, hem, flag, ped, sigma, chipid, chipchan);
				}
				else if ( pixelplane ) {
					const CS::ChipAPV::DigitPixel *digit1p =
							dynamic_cast<const CS::ChipAPV::DigitPixel*>(digit1);

					if ( digit1p == NULL ) {
						CsErrLog::msg(elWarning, __FILE__, __LINE__,
								"%s: ChipAPV wrong map!", GetTBName().c_str());
						continue;
					}

					// Apply mapping
					try {
						wire = digit1p->GetChannel();
						chipid = digit1p->GetChip();
						chipchan = digit1p->GetChipChannel();
					}
					catch (...) {
						CsErrLog::msg(elError, __FILE__, __LINE__,
								"%s: Could not apply mapping!", GetTBName().c_str());
						continue;
					}

					// Get pixel coordinates
					pair<int,int> xy = pixgem::detch2xy(wire);
					PX = xy.first;
					PY = xy.second;

					// Pixel coordinates invalid
					if ( PX == -1 && PY == -1 ) {
						CsErrLog::msg(elError, __FILE__, __LINE__,
								"%s: Could not apply mapping!",
								GetTBName().c_str());
						continue;
					}

					// Invert X pads if detector orientation negative (det facing upstream)
					if ( digit1p->GetDetOrientation() < 0 ) {
						PX = 31 - PX;
					}

					pixelplane->AddChan(wire, PX, PY, flag, ped, sigma, chipid, chipchan);
				}
				else if ( rectpixelplane ) {
					const CS::ChipAPV::DigitPixelMM *digit1pmm =
							dynamic_cast<const CS::ChipAPV::DigitPixelMM*>(digit1);

					if ( digit1pmm == NULL ) {
						CsErrLog::msg(elWarning, __FILE__, __LINE__,
								"%s: ChipAPV wrong map!", GetTBName().c_str());
						continue;
					}

					// Apply mapping
					try {
						wire = digit1pmm->GetChannel();
						chipid = digit1pmm->GetChip();
						chipchan = digit1pmm->GetChipChannel();
					}
					catch (...) {
						CsErrLog::msg(elError, __FILE__, __LINE__,
								"%s: Could not apply mapping!", GetTBName().c_str());
						continue;
					}

					// Get pixel coordinates
					conn = digit1pmm->GetConnNb();
                                        register int mapPixelVersion = digit1pmm->GetPixelMMversion();
                                        if (!pixmmobject) pixmmobject = pixmm::Instance(mapPixelVersion);
					pixnb = pixmmobject->GetPixNb(conn,wire);
                                        conn09 = pixmmobject->GetNConn(conn); // convert from 0-19 to 0-9 numbering

					// Pixel coordinates invalid
					if ( pixnb == -1 ) {
						CsErrLog::msg(elError, __FILE__, __LINE__,
								"%s: Could not apply mapping!",
								GetTBName().c_str());
						continue;
					}
					PXMM = pixmmobject->GetXPix(pixnb);
					PYMM = pixmmobject->GetYPix(pixnb);
					chandet = wire + conn*ChipChannels;
					// Invert X pads if detector orientation negative (det facing upstream)
					if ( digit1pmm->GetDetOrientation() < 0 ) {
						PXMM = -PXMM;
					}

					rectpixelplane->AddChan(chandet, pixnb, PXMM, PYMM, flag, ped, sigma, chipid, chipchan, conn09);
				}
			}
		}
	} // Daq maps found

	// No daq maps found
	else {
		CsErrLog::msg(elFatal, __FILE__, __LINE__,
				"%s: Daq maps not found!", GetTBName().c_str());
		return false;
	}

	// Debugging output
	if ( getenv("COND_DB_DEBUG") != 0 || getenv("COND_DB_DEBUG_MUMEGA") != 0 ) {
		if ( stripplane ) stripplane->PrintChannels();
		else if ( pixelplane ) pixelplane->PrintChannels();
		else if ( rectpixelplane ) rectpixelplane->PrintChannels();
	}
	if (fDebugLevel > 9) {
		if ( rectpixelplane ) rectpixelplane->PrintChannelsSimple();
	}



	return true;
}


bool CsPixelMumegaDetector::mapMCChannelCal() {
	if ( !stripplane && !pixelplane && !rectpixelplane ) {
		CsErrLog::msg(elFatal, __FILE__, __LINE__,
				"%s: mapMCChannelCal() was called, but plane objects do not exist!",
				GetTBName().c_str());

		return false;
	}

	CsErrLog::msg(elError, __FILE__, __LINE__,
			"%s: Using default calibrations!", GetTBName().c_str());

	// Delete any existing calibrations
	if ( stripplane ) {
		if ( stripplane->GetNchannel() != 0 ) {
			stripplane->ClearChannels();
		}
	}
	// Delete any existing calibrations
	if ( pixelplane ) {
		if ( pixelplane->GetNchannel() != 0 ) {
			pixelplane->ClearChannels();
		}
	}
	// Delete any existing calibrations
	if ( rectpixelplane ) {
		if ( rectpixelplane->GetNchannel() != 0 ) {
			rectpixelplane->ClearChannels();
		}
	}

	// Fill calibrations for all wires of plane
	int flag = 1;
	float ped = 750.;
	float sigma = 4.; 
	int   hem   =   1;
	int wire, PX, PY;
	float PXMM = -1; float PYMM = -1;
	int pixnb = -1, chandet = -1;
	int shortpix_start, shortpix_end;
	if ( stripplane ) {
		if(this->hasVarP_){
			if (wirPs_.size() == 3){
				map<int, double>::iterator it_pitch = wirPs_.begin();
				it_pitch++;
				shortpix_start = (*it_pitch).first; it_pitch++;
				shortpix_end = (*it_pitch).first;
			}
			else
			{
				CsErrLog::msg(elError, __FILE__, __LINE__,
						"%s: Using default strips geometry!", GetTBName().c_str());
				shortpix_start = 256;
				shortpix_end = 640;
			}
			for ( int wire = 0; wire < nWir_; wire++ ) { //wire = address
				if((wire>=shortpix_start)&&(wire<shortpix_end)){
					stripplane->AddChan(wire, hem, flag, ped, sigma);
					stripplane->AddChan(wire, -hem, flag, ped, sigma);
				}
				else{
					stripplane->AddChan(wire, 0, flag, ped, sigma);
				}
			}
		}
		else{ //no variable pitch
			for ( int wire = 0; wire < nWir_; wire++ ) {
				stripplane->AddChan(wire,  hem, flag, ped, sigma);
				stripplane->AddChan(wire, -hem, flag, ped, sigma);
			}
		}

		// set the standard time calibration, this time calibration is used
		// in MC to simulate the amplitude of the three samples, though either
		// change all occurences or no!!!
		CsMumegaTimeCals time_cals; time_cals.Clear();
		CsMumegaTimeCalOld *cal0 = new CsMumegaTimeCalOld(0.2,1.1,-40.,22.,1.65);
		CsMumegaTimeCalOld *cal1 = new CsMumegaTimeCalOld(0.3,1.0,0.,27.,1.23);
		time_cals.Set(cal0,cal1);
		stripplane->SetTimeCals(time_cals);
	}
	else if ( pixelplane ) {
		for ( int pixU = 0; pixU < nWir_; pixU++ )
			for ( int pixV = 0; pixV < nWirV_; pixV++ )
				pixelplane->AddChan(pixU + pixV*nWir_, pixU, pixV, flag, ped, sigma);

		// set the standard time calibration, this time calibration is used
		// in MC to simulate the amplitude of the three samples, though either
		// change all occurences or no!!!
		CsMumegaTimeCals time_cals; time_cals.Clear();
		CsMumegaTimeCalOld *cal0 = new CsMumegaTimeCalOld(0.2,1.1,-40.,22.,1.65);
		CsMumegaTimeCalOld *cal1 = new CsMumegaTimeCalOld(0.3,1.0,0.,27.,1.23);
		time_cals.Set(cal0, cal1);
		pixelplane->SetTimeCals(time_cals);
	}
	else if ( rectpixelplane ) {
		for ( int i = 0; i < pixmmobject->GetNbPix(); i++ ){
			// Get pixel coordinates
			chandet = i;
			pixnb = pixmmobject->GetPixNb(chandet);
			PXMM = pixmmobject->GetXPix(pixnb);
			PYMM = pixmmobject->GetYPix(pixnb);

			rectpixelplane->AddChan(chandet, pixnb, PXMM, PYMM, flag, ped, sigma);
		}
	}

	return true;
}

//____________________________________________________________________________
double CsPixelMumegaDetector::getCorrU(const CsCluster *c,
		const double x, const double y, // MRS, mm
		const double tt,  // Offset w.r.t. *event* time
		bool &error)
{
	if(corr_U_bool == true){
		// *** This is a fast method to correct the alignment of the detector (cluster's position)
		// due to geometry detector.

		static bool first=true;
		if(first){
			cout<<"CsPixelMumegaDetector::getCorrU : Correction of MP cluster position is being applied!"<<endl;
			first=false;
		}

		double angle, cos_angle, sin_angle;
		double center_u, center_v;
		double utrk, vtrk;

		/// get the position in the detector of the track associated to the cluster c
		utrk = (x*corr_U_cos_angle + y*corr_U_sin_angle)*0.1;
		vtrk = (-x*corr_U_sin_angle + y*corr_U_cos_angle)*0.1;

		int indice; // Indice of the table correction
		double u_in_det = ((utrk - corr_U_center_u)+ (getXsiz()*0.05));
		double v_in_det = ((vtrk - corr_U_center_v)+ (getYsiz()*0.05));
		int indice_u = (int) floor(u_in_det/corr_U_size_u);
		int indice_v = (int) floor(v_in_det/corr_U_size_v);

		indice = (indice_u + (indice_v*corr_U_nb_u));

		/// applied correction to the postion of the cluster c
		if ((indice < (corr_U_nb_u * corr_U_nb_v)) && (indice >= 0) ){ //check the indice
			double corrU = c->getU() + corr_U_tab[indice]*10;
			error=false;
			return corrU ; //in mm
		}
		else {
			CsErrLog::msg(elError, __FILE__, __LINE__,
					"%s: wrong indice for correction U table! indice (U,V): %d (%d,%d)",
					GetTBName().c_str(),indice,indice_u,indice_v);
			error = true;
			return c->getU();
		}
	}
	else {
		error = true;
		return c->getU();
	}
}


/// return the closest distant between a pixel with coordinate (pixU, pixV) and a point with coordinate (pointU, pointV)
double CsPixelMumegaDetector::dist_pixtopoint(double pixU, double pixV, double pointU, double pointV){
	int pixUf, pixVf; // pixel belonging to the point
	if ( (pointU - wirDU_)/wirPU_ < 0 ) pixUf = int( (pointU - wirDU_)/wirPU_ - 0.5 );
	else                            pixUf = int( (pointU - wirDU_)/wirPU_ + 0.5 );
	if ( (pointV - wirDV_)/wirPV_ < 0 ) pixVf = int( (pointV - wirDV_)/wirPV_ - 0.5 );
	else                            pixVf = int( (pointV - wirDV_)/wirPV_ + 0.5 );

	double posU_pix, posV_pix;
	double distU = 0.0;
	double distV = 0.0;

	posU_pix = pixU*wirPU_ + wirDU_;
	posV_pix = pixV*wirPV_ + wirDV_;

	if (posV_pix < pointV){
		distV = (pointV - (posV_pix + (wirPV_/2)));
	}
	else{
		distV = ((posV_pix - (wirPV_/2)) - pointV);
	}

	if(pixUf == pixU){ //same coordinate in U : minimal distance is distV
		return distV;
	}

	if (posU_pix < pointU){
		distU = (pointU - (posU_pix + (wirPU_/2)));
	}
	else{
		distU = ((posU_pix - (wirPU_/2)) - pointU);
	}

	if(pixVf == pixV){ //same coordinate in V : minimal distance is distU
		return distU;
	}

	return sqrt(distU*distU + distV*distV);
}

/// return the closest distance between a strip with number wire and hemisphere hem and a point with coordinate (pointU, pointV)
/// change by DN 21/9/2015: added uShift,vShift args to be able to center to the detector !
double CsPixelMumegaDetector::dist_striptopoint(const int wire, const double pointU, const double pointV,
                                                const double uShift, const double vShift, const int hem) {
	int wiref, hemf;  // Strip belonging to the point
	double posU, posV;
	bool error;
	double distU = 0.0;
	double distV = 0.0;

	if (pointV > 0.){
		hemf = 1;
	}
	else {
		hemf = -1;
	}

	posU = WireToDist(wire, error);
	if (fabs(pointV-vShift) <= (getYsiz()*0.5)) { //in detector (in V coordinate)
		if ((wire < 384) || (wire > 511)){ //long strips
			posV = 0.0;
		}
		//short strips
		else if ((wire < 392) || (wire > 503)){
			posV = 12.5*hem;
		}
		else if ((wire < 408) || (wire > 487)){
			posV = 18.75*hem;
		}
		else{
			posV = 25.*hem;
		}
	}
	else if (pointV > 0.) { //over-range
		posV = (getYsiz()*0.5)+vShift;
		hemf = 3;
	}
	else { //under-range
		posV = -(getYsiz()*0.5)-vShift;
		hemf = -3;
	}
	distV = fabs(pointV - posV);
	if (fabs(pointU-uShift) <= (getXsiz()*0.5)) { //in detector in U coordinate
		wiref = DistToWire(pointU, error);
		if((wiref == wire) && (hemf != hem)){ //same wire (coordinate in U) : minimal distance is distV
			return distV;
		}
	}

	distU = (fabs(pointU - posU) - (getWirePitch(wire)/2));

	if((hemf == hem) && (fabs(posV) <= fabs(pointV))){ //minimal distance is distU
		return distU;
	}

	return sqrt(distU*distU + distV*distV);
}


double CsPixelMumegaDetector::weight_for_pix(const double pixU, const double pixV, const double pointU, const double pointV, const double cluster_size, const double pathX){
	double posU_pix, posV_pix;
	double amp_Bot, amp_Cen, amp_Top; // the three amplitudes for the average amplitude of the pixel
	double distU_2, distV_2, W_sup, W_inf, distU;


	double sigma = elec_spread_mean;

	posU_pix = pixU*wirPU_ + wirDU_;
	posV_pix = pixV*wirPV_ + wirDV_;

	distU = max((fabs(posU_pix - pointU) - pathX),0.);
	distU_2 = distU*distU;

	distV_2 = (posV_pix - (wirPV_/2) - pointV)*(posV_pix - (wirPV_/2) - pointV);
	amp_Bot = exp(-(distU_2+distV_2)/(2*sigma*sigma));

	distV_2 = (posV_pix - pointV)*(posV_pix - pointV);
	amp_Cen = exp(-(distU_2+distV_2)/(2*sigma*sigma));

	distV_2 = (posV_pix + (wirPV_/2) - pointV)*(posV_pix + (wirPV_/2) - pointV);
	amp_Top = exp(-(distU_2+distV_2)/(2*sigma*sigma));

	W_sup = ((amp_Top + amp_Cen)/2);
	W_inf = ((amp_Bot + amp_Cen)/2);

	return (max(((W_sup + W_inf)+(sqrt(W_sup + W_inf)*elec_spread_sigma*CsRandom::gauss())),0.));
}


/// change by DN 21/9/2015: added uShift,vShift args to be able to center to the detector !
double CsPixelMumegaDetector::weight_for_strip(const int wire, const double pointU, const double pointV,
                                               const double uShift, const double vShift,
                                               const double cluster_size, const double pathX, const int hem) {
	double posU, lim_sup, lim_inf;
	bool error;
	double distU_2, weight, distU;
	double distUmin_2;

	double clust_dim = (cluster_size / 2.); //cluster dimension in the strip
	double radius_2 = clust_dim*clust_dim; //global radius square
	double sigma = elec_spread_mean;

	posU = WireToDist(wire, error);

	distU = max((fabs(posU - pointU) - pathX),0.);
	distU_2 = distU*distU;


	if(hem == 3){ //if hem == 3 we return the total weight for this strip
		return exp(-(distU_2)/(2*sigma*sigma))*getWirePitch(wire);
	}
	else{
		weight = exp(-(distU_2)/(2*sigma*sigma))*getWirePitch(wire);
	}


	////// find the cluster dimension in the strip
	distUmin_2 = (distU - (getWirePitch(wire)/2.)) * (distU - (getWirePitch(wire)/2.));
	clust_dim = sqrt(radius_2 - distUmin_2);//Pythagore equation
	////////////////:

	/////////////find the strip limits
	if ((wire < 384) || (wire > 511)){ //long strips
		lim_inf = 0.0;
	}
	//short strips
	else if ((wire < 392) || (wire > 503)){
		lim_inf = 12.5*hem;
	}
	else if ((wire < 408) || (wire > 487)){
		lim_inf = 18.75*hem;
	}
	else{
		lim_inf = 25.*hem;
	}

	lim_sup = (getYsiz()*0.5)*hem;

// 	if (lim_inf != 0.0)//long strips
// 		lim_inf += vShift;
	lim_inf += vShift;
	lim_sup += vShift;
	if ( lim_sup < lim_inf ) {//lim_sup should be higher than lim_inf
		double tmp = lim_sup;
		lim_sup = lim_inf;
		lim_inf = tmp;
	}
	///////////////////////

	return (max((weight + (sqrt(weight)*elec_spread_sigma*CsRandom::gauss())),0.)*getErfFact(lim_inf, pointV, lim_sup, (2*clust_dim)));
}

// useless method, as initilization done in pixmm object of the decoding library
void CsPixelMumegaDetector::initialise_pix(void){
// 	int conn_nb, pixU, pixV;
// 	double pos_u = 0.0;
// 	double pos_v = 0.0;
// 	const double wirDU_DRS = -((wirPU_ * (nWirU_-1))/2);
// 	const double wirDV_DRS = -((wirPV_ * (nWirV_-1))/2);
// 
// 	for(int iu=0; iu<CsPixelMumegaDetector::nbpixX; iu++) for(int iv=0; iv<CsPixelMumegaDetector::nbpixY; iv++){
// 		CsPixelMumegaDetector::pix_address[iu][iv]= -1;
// 	}
// 
// 	for(int i=0; i<pixmm::Getpix_size(); i++){
// 		conn_nb = pixmm::GetPix(i);
// 		pos_u = pixmm::GetXPix(conn_nb);
// 		pos_v = pixmm::GetYPix(conn_nb);
// 		if ( (pos_u - wirDU_DRS)/wirPU_ < 0 ) pixU = int( (pos_u - wirDU_DRS)/wirPU_ - 0.5 );
// 		else                            pixU = int( (pos_u - wirDU_DRS)/wirPU_ + 0.5 );
// 		if ( (pos_v - wirDV_DRS)/wirPV_ < 0 ) pixV = int( (pos_v - wirDV_DRS)/wirPV_ - 0.5 );
// 		else                            pixV = int( (pos_v - wirDV_DRS)/wirPV_ + 0.5 );
// 
// 		if((pixU < 128)&&(pixU >=0)&&(pixV < 38)&&(pixV >= 2)){
// 			if((pixU>31)&&(pixU<96)&&(pixV>9)&&(pixV<30)){
// 				CsPixelMumegaDetector::pix_address[pixU][pixV]= i;
// 				if (pixV%2 == 0){
// 					CsPixelMumegaDetector::pix_address[pixU][pixV+1]= i;
// 				}
// 				else
// 					CsPixelMumegaDetector::pix_address[pixU][pixV-1]= i;
// 			}
// 			else{
// 				CsPixelMumegaDetector::pix_address[pixU][pixV-2]= i;
// 				CsPixelMumegaDetector::pix_address[pixU][pixV-1]= i;
// 				CsPixelMumegaDetector::pix_address[pixU][pixV]= i;
// 				CsPixelMumegaDetector::pix_address[pixU][pixV+1]= i;
// 				CsPixelMumegaDetector::pix_address[pixU][pixV+2]= i;
// 			}
// 		}
// 		else {
// 			CsErrLog::msg(elError, __FILE__, __LINE__,
// 					"%s wrong initialization of pixel address for MC :\nPixU : %d \nPixV : %d ", GetTBName().c_str(),pixU,pixV);
// 		}
// 	}

}

void CsPixelMumegaDetector::MChits_Xdet_pix(void){
	list<CsMCHit*>::const_iterator Ih;
	const list<CsMCHit*> & lMChits = complementary_det->getMyMCHits();

	for ( Ih = lMChits.begin(); Ih != lMChits.end(); Ih++ ) { // loop on hits

		if(strcmp((*Ih)->getDet()->GetTBName().c_str(),this->GetTBName().c_str()) == 0) continue; //avoid adding twice the same MChit in the list

		const CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
		if ( thit == 0 ) return;

		double t  = thit->getDTime(); // Delay time (ns)
		double ui = thit->getUin();   // Hit in point (DRS)
		double vi = thit->getVin();
		double wi = thit->getWin();
		double uo = thit->getUout();  // hit out point (DRS)
		double vo = thit->getVout();
		double wo = thit->getWout();

		HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
		double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+ xcm_ - _xshift;
		double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+ ycm_ - _yshift;
		double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+ zcm_;
		double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+ xcm_ - _xshift;
		double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+ ycm_ - _yshift;
		double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+ zcm_;
		double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
		double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
		double Wi = iRotM(3,1)*xi+iRotM(3,2)*yi+iRotM(3,3)*zi;
		double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
		double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
		double Wo = iRotM(3,1)*xo+iRotM(3,2)*yo+iRotM(3,3)*zo;
		double uShift = iRotM(1,1)*xcm_+iRotM(1,2)*ycm_;  // shift of plane wrt WRS, for cuts based on detector geometry
		double vShift = iRotM(2,1)*xcm_+iRotM(2,2)*ycm_;
		if(((fabs(Ui-uShift) < (getXsiz()*0.5 + 3.))||(fabs(Uo-uShift) < (getXsiz()*0.5 + 3.)))&&((fabs(Vi-vShift) < (getYsiz()*0.5 + 3.))||(fabs(Vo-vShift) < (getYsiz()*0.5 + 3.)))){
			addMCHit(*(*Ih)); // no changes on the MChit, in this way it keep it's original CsDet
			continue;
		}
	}
}

void CsPixelMumegaDetector::MChits_Xdet_strip(void){
	list<CsMCHit*>::const_iterator Ih;
	const list<CsMCHit*> & lMChits = complementary_det->getMyMCHits();

	for ( Ih = lMChits.begin(); Ih != lMChits.end(); Ih++ ) { // loop on hits

		if(strcmp((*Ih)->getDet()->GetTBName().c_str(),this->GetTBName().c_str()) == 0) continue; //avoid adding twice the same MChit in the list

		const CsMCTrkHit *thit = dynamic_cast<CsMCTrkHit*>(*Ih);
		if ( thit == 0 ) return;

		double t  = thit->getDTime(); // Delay time (ns)
		double ui = thit->getUin();   // Hit in point (DRS)
		double vi = thit->getVin();
		double wi = thit->getWin();
		double uo = thit->getUout();  // hit out point (DRS)
		double vo = thit->getVout();
		double wo = thit->getWout();

		HepMatrix iRotM(3,3); iRotM = getRotWRSInv();
		double xi = rotDRS_(1,1)*ui+rotDRS_(1,2)*vi+rotDRS_(1,3)*wi+ xcm_ - _xshift;
		double yi = rotDRS_(2,1)*ui+rotDRS_(2,2)*vi+rotDRS_(2,3)*wi+ ycm_ - _yshift;
		double zi = rotDRS_(3,1)*ui+rotDRS_(3,2)*vi+rotDRS_(3,3)*wi+ zcm_;
		double xo = rotDRS_(1,1)*uo+rotDRS_(1,2)*vo+rotDRS_(1,3)*wo+ xcm_ - _xshift;
		double yo = rotDRS_(2,1)*uo+rotDRS_(2,2)*vo+rotDRS_(2,3)*wo+ ycm_ - _yshift;
		double zo = rotDRS_(3,1)*uo+rotDRS_(3,2)*vo+rotDRS_(3,3)*wo+ zcm_;
		double Ui = iRotM(1,1)*xi+iRotM(1,2)*yi+iRotM(1,3)*zi; // WRS
		double Vi = iRotM(2,1)*xi+iRotM(2,2)*yi+iRotM(2,3)*zi;
		double Uo = iRotM(1,1)*xo+iRotM(1,2)*yo+iRotM(1,3)*zo; // WRS
		double Vo = iRotM(2,1)*xo+iRotM(2,2)*yo+iRotM(2,3)*zo;
		double uShift = iRotM(1,1)*xcm_+iRotM(1,2)*ycm_;  // shift of plane wrt WRS, for cuts based on detector geometry
		double vShift = iRotM(2,1)*xcm_+iRotM(2,2)*ycm_;
		if ((fabs(Ui-uShift)> (25.6 - 3))||(fabs(Vi-vShift)> (25. - 3))||(fabs(Uo-uShift)> (25.6 - 3))||(fabs(Vo-vShift)> (25. - 3))){
			addMCHit(*(*Ih)); // no changes on the MChit, in this way it keep it's original CsDet
			continue;
		}
		if (((fabs(Ui-uShift)> (22.4 - 3))&&(fabs(Vi-vShift)> (12.5 - 3)))||((fabs(Uo-uShift)> (22.4 - 3))&&(fabs(Vo-vShift)> (12.5 - 3)))){
			addMCHit(*(*Ih)); // no changes on the MChit, in this way it keep it's original CsDet
			continue;
		}
		if (((fabs(Ui-uShift)> (16. - 3))&&(fabs(Vi-vShift)> (18.75 - 3)))||((fabs(Uo-uShift)> (16. - 3))&&(fabs(Vo-vShift)> (18.75 - 3)))){
			addMCHit(*(*Ih)); // no changes on the MChit, in this way it keep it's original CsDet
			continue;
		}
	}
}


bool CsPixelMumegaDetector::in_pix(const double u, const double v){

	double wirDU_pix;
	double wirDV_pix;
	double wirPU_pix;
	double wirPV_pix;
	if(complementary_det != NULL){
// 		wirDU_pix = complementary_det->getWirD();
// 		wirDV_pix = complementary_det->getWirP();
// 		wirPU_pix = complementary_det->getWirDV();
// 		wirPV_pix = complementary_det->getWirPV();
		wirDU_pix = complementary_det->getWirD();
		wirPU_pix = complementary_det->getWirP();
		wirDV_pix = complementary_det->getWirDV();
		wirPV_pix = complementary_det->getWirPV();

	}
	else {
		wirDU_pix = pixmmobject->GetFirstWireU();
		wirDV_pix = pixmmobject->GetFirstWireV();
		wirPU_pix = pixmmobject->GetPitchU();
		wirPV_pix = pixmmobject->GetPitchV();
	}

	if ((fabs(u)> (25.6))||(fabs(v)> (25.))) return false;

	int pixU, pixV; // center of the "cluster"
	if ( (u - wirDU_pix)/wirPU_pix < 0 ) pixU = int( (u - wirDU_pix)/wirPU_pix - 0.5 );
	else                            pixU = int( (u - wirDU_pix)/wirPU_pix + 0.5 );
	if ( (v - wirDV_pix)/wirPV_pix < 0 ) pixV = int( (v - wirDV_pix)/wirPV_pix - 0.5 );
	else                            pixV = int( (v - wirDV_pix)/wirPV_pix + 0.5 );

// 	if(CsPixelMumegaDetector::pix_address[pixU][pixV] == -1) return false;
// 	if (pixmmobject->PixelAddresses()[(pixU,pixV)] == -1) return false;
	if (pixmmobject->PixelAddr(pixU,pixV) == -1) return false;
	else return true;
}

int CsPixelMumegaDetector::DistToWire(const double U_pos, bool &error){
	error = false;

	if(!hasVarP_){
		if (( U_pos - wirD_ )/wirP_ < 0 ) {
			return int(( U_pos - wirD_ )/wirP_ - 0.5 );
		}
		else {
			return int(( U_pos - wirD_ )/wirP_ + 0.5 );
		}
	}
	else{
		map<int,double>::const_iterator i;
		double dim = (wirD_ - (wirPs_[0]/2));  // Hold dimension covered (WRS)
		double dim_tmp = (wirD_ - (wirPs_[0]/2)); // better than (-getXsiz()/2)
		int lastw = 0; double lastp = 0;
		for (i = wirPs_.begin(); i!=wirPs_.end(); i++) {
			dim_tmp += ((i->first-lastw)*lastp);
			if((dim_tmp >= U_pos)&&(lastp != 0.)){
				if (( U_pos - dim )/lastp < 0 ) {
					return int (lastw + (( U_pos - dim )/lastp - 0.5 ));
				}
				else {
					return int (lastw + (( U_pos - dim )/lastp + 0.5 ));
				}
			}
			dim += (((i->first-lastw)*lastp - (lastp/2)) + (i->second/2));
			lastw = i->first; lastp = i->second;
		}
		dim_tmp += ((nWir_ - lastw)*lastp);
		if(dim_tmp >= U_pos){
			// first wire for this hit
			if (( U_pos - dim )/lastp < 0 ) {
				return int (lastw + (( U_pos - dim )/lastp - 0.5 ));
			}
			else {
				return int (lastw + (( U_pos - dim )/lastp + 0.5 ));
			}
		}
	}
	error = true;
	return -1;
}


double CsPixelMumegaDetector::WireToDist(const int wire, bool &error ){
	error = false;

	if((wire < 0) || (wire >= nWir_)){
		error = true;
		return -1;
	}
	if(!hasVarP_){
		return double(wire*wirP_ + wirD_);
	}
	else{
		map<int,double>::const_iterator i;
		double dim = (wirD_ - (wirPs_[0]/2));  // Hold dimension covered (WRS)
		int lastw = 0; double lastp = 0;
		for (i = wirPs_.begin(); i!=wirPs_.end(); i++) {
			if(i->first > wire){
				return double (dim + (wire - lastw)*lastp);
			}
			dim += (((i->first-lastw)*lastp - (lastp/2)) + (i->second/2));
			lastw = i->first; lastp = i->second;
		}
		if(nWir_ > wire){
			return double (dim + (wire - lastw)*lastp);
		}
	}
	error = true;
	return -1;
}

double CsPixelMumegaDetector::getWirePitch(const int wire){
	double pitch = 0; map<int,double>::iterator ipitch;
	for (ipitch = wirPs_.begin(); ipitch!=wirPs_.end(); ipitch++ ){
		if (wire >= ipitch->first) pitch = ipitch->second;
	}
	return pitch;
}

double CsPixelMumegaDetector::getErfFact(const double min, const double center, const double max, const double clust_size){
	double radius = (clust_size / 2.);
	double sigma = elec_spread_mean;
	double fact1 = 0.5;
	double fact2 = 0.5;
	double result = 0.;
	double dist1 = fabs(center - min);
	double dist2 = fabs(max - center);
	double range = fabs(min - max);

	if(dist1 < radius){
		fact1 = (0.5*(1 + erf(dist1/(sigma*sqrt(2.)))) - fact1);
	}
	if(dist2 < radius){
		fact2 = (0.5*(1 + erf(dist2/(sigma*sqrt(2.)))) - fact2);
	}

	if((dist1 > range) || (dist2 > range)){
		result = fabs(fact1 - fact2);
		if ((result < 0.)||(result > 0.5)){//beware of the double limits
			return 0.;
		}
		else return result;
	}
	else{
		return (fact1 + fact2);
	}
}
/*
//-----------------------------------------------------------------------------
// Function which returns true if point xp,yp lies inside the
// polygon defined by the points in vectors x and y, false otherwise
// NOTE that the polygon must be a closed polygon (1st and last point
// must be identical) (taken from ROOT)
// Put to separate namespace CsPixelMumega
//-----------------------------------------------------------------------------
bool CsPixelMumega::IsInside(const float xp, const float yp, const std::vector<float>& x, const std::vector<float>& y) {
	register double xint;
	register int i;
	register int inter = 0;
	register int np = x.size();
	for (i=0;i<np-1;i++) {
		if (y[i] == y[i+1]) continue;
		if (yp <= y[i] && yp <= y[i+1]) continue;
		if (y[i] < yp && y[i+1] < yp) continue;
		xint = x[i] + (yp-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]);
		if (xp < xint) inter++;
	}
	if (inter%2) return true;
	return false;
}
*/
