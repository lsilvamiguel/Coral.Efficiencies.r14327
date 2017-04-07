#ifndef CsPixelMumegaDetector_h
#define CsPixelMumegaDetector_h

//---------------------------------------------------------------------------

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include "CsDetector.h"
#include "CsHistograms.h"
#include "TCutG.h"
#include "DaqDataDecoding/Chip.h"
#include "DaqDataDecoding/ChipAPV.h"
#include "DaqDataDecoding/PixelMMDecoding.h"
#include "CsMumegaPlane.h"
#include "CsPixelMumegaPlane.h"
#include "CsRectPixelMumegaPlane.h"

class CsZone;

class CsPixelMumegaDetector : public CsDetector {

public:

	int ChipsPerPlane;
	enum {ChipsPerStripPlane=10, ChipsPerPixelPlane=8, ChipsPerRectPixelPlane=10, ChipChannels=128};//modif 12/11/12 (4, 8, 10)

	class APVCal
	{
	public:

		class Channel
		{
		public:
			Channel(int f, double pm, double ps, double cm, double cs) :
				flag(f),
				pedestal_mean(pm), pedestal_sigma(ps),
				calibration_mean(cm), calibration_sigma(cs)
			{}

			Channel(const char *s)
			{
				if ( 5 != sscanf(s, "%d%lf%lf%lf%lf", &flag, &pedestal_mean, &pedestal_sigma, &calibration_mean, &calibration_sigma) )
					throw CS::Exception("CsPixelMumegaDetector::Calib::Channel::Channel() : bad line \"%s\"",s);
			}

			void Print(std::ostream &o=std::cout, const std::string &prefix="") const
			{
				o << prefix;
				o << "flag = " << flag
						<< "PED = ("
						<< std::setw(8) << std::setprecision(3) << pedestal_mean
						<< ", "
						<< std::setw(8) << std::setprecision(3) << pedestal_sigma
						<< ") CALIB = ("
						<< std::setw(8) << std::setprecision(3) << calibration_mean
						<< ", "
						<< std::setw(8) << std::setprecision(3) << calibration_sigma
						<< ")" << std::endl;
			}

			int flag;
			double pedestal_mean;
			double pedestal_sigma;
			double calibration_mean;
			double calibration_sigma;
		};

		int src_id;
		int adc_id;
		int chip_id;

		std::vector<Channel> channels;

		friend std::istream & operator >> (std::istream &, std::vector<APVCal> &c);
	};

public:

	CsPixelMumegaDetector(const int    row,
			const int    id,     const char*  name, const char *TBname,
			const int    unit,   const int    type,
			const double rdLen,  const double xsiz,
			const double ysiz,   const double zsiz,
			const double xcm,    const double ycm,
			const double zcm,    const CLHEP::HepMatrix rotDRS,
			const CLHEP::HepMatrix rotWRS,
			const double wirD,   const double ang,
			const int    nWir,   const double wirP,
			const double eff,    const double bkg,
			const double tGate,  const double spSig,
			const double eGain,  const double eGSig,
			const double sWidth, const double tRes);

	void AddSubDetector(const int    row,
			const int    id,    const char*  name,   const char *TBname,
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

	//! Returns the first pixel position along V (= WRS + pi/2) (mm)
	double getWirDV() const { return(wirDV_); }
	//! Returns the pixel pitch along V (mm)
	double getWirPV() const { return(wirPV_); }
	//! Returns number of pixels along V
	int    getNWirV() const { return(nWirV_); }

private:

	bool operator==( const CsPixelMumegaDetector& ) const; //!< "equal to" operator

	bool operator<( const CsPixelMumegaDetector& ) const; //!< "less than" operator

public:

	//! book all histograms (they should not be created in the constructor)
	void BookHistograms();

	//! Returns true if TCD informations are available
	inline bool hasTDC() const { return false; }


	//! Returns true if Drift informations are available
	inline bool hasDrift() const { return false; }

	//! Decode the MC data for this detector
	void makeMCDecoding();

private: // decode the MC data for the two kind
	// of detectors (pixel/strip) and the
	// two different options in hit
	// creation (simplistic/with amplitude)
	void makeMCPixelsSimple();
	void makeMCPixelsAmplitude();
	void makeMCStripsSimple();
	void makeMCStripsAmplitude();

	/*! \fn double dist_pixtopoint(double pixU, double pixV, double pointU, double pointV)
	  		\brief return the distance between the closest corner of a pixel and the point (in WRS, in mm)
			\param pixU, pixV 	pixel numbers with reference 0,0 in bottom-left corner (-25.4mm,-21.875mm : WRS)
			\param pointU, pointV 	coordinates of the point of interest (in WRS, in mm)
	 */
	double dist_pixtopoint(double pixU, double pixV, double pointU, double pointV);

	/*! \fn double dist_striptopoint(const int wire, const double pointU, const double pointV,
                                        const double uShift, const double vShift,  const int hem)
		  		\brief return the distance between the closest corner of a strip and the point (in WRS, in mm)
				\param wire		strip numbers with reference 0 in the strip the most negative in U (in WRS)
				\param pointU, pointV	coordinates of the point of interest (in WRS, in mm)
				\param uShift,vShift	shift to center the detector (in mm)
				\param hem		hemisphere of the strip (-1, 0, 1)
	 */
	double dist_striptopoint(const int wire, const double pointU, const double pointV,
                                 const double uShift, const double vShift,  const int hem);

	/*! \fn double weight_for_pix(double pixU, double pixV, double pointU, double pointV, double cluster_size)
		  	\brief return the integral of the electronic shower for a pixel
			\param pixU, pixV 	pixel numbers with reference 0,0 in bottom-left corner (-25.4mm,-21.875mm : WRS)
			\param pointU, pointV 	coordinates of the center of the cluster (shower) (in WRS, in mm)
			\param cluster_size 	size of the cluster diameter in (mm)
	 */
	double weight_for_pix(const double pixU, const double pixV, const double pointU, const double pointV,
                              const double cluster_size, const double pathX);

	/*! \fn double weight_for_strip(const int wire, const double pointU, const double pointV,
                                        const double uShift, const double vShift, const double cluster_size, const int hem)
			  	\brief return the integral of the electronic shower for a strip in accordance with its pitch
				\param wire 	strip numbers with reference 0 in the strip the most negative in U (in WRS)
				\param pointU, pointV 	coordinates of the center of the cluster (shower) (in WRS, in mm)
				\param uShift,vShift	shift to center the detector (in mm)
				\param cluster_size 	size of the cluster diameter in (mm)
				\param hem		hemisphere of the strip (-1, 0, 1) if hem == 3 then the function return the total weight for this strip (size infinite : no cut)
	 */
	double weight_for_strip(const int wire, const double pointU, const double pointV,
                                const double uShift, const double vShift,
                                const double cluster_size, const double pathX, const int hem);

	/* \fn static void initialise_pix(void)
	 	 \brief initialise the pix_adress table.
	 */
	void initialise_pix(void);

	//! Returns the amplitude from an deposited energy
	void getMCAmp(const double& t, double energy, double amps[] );

	/*! store MChit from the complementary detector in the MChit vector
	 * 	for pixel detector (MChits from strip detector to pixel detector)
	 * 	In this way we escape dependency on detector'limits of the MChit (simplest way
	 * 	found to reproduced the pixel shape)
	 */
	void MChits_Xdet_pix(void);

	/*! store MChit from the complementary detector in the MChit vector
	 * 	for strip detector (MChits from pixel detector to strip detector)
	 * 	In this way we escape dependency on detector'limits of the MChit (simplest way
	 * 	found to reproduced the pixel shape)
	 */
	void MChits_Xdet_strip(void);

public:

	//! Amplitude correlation is simulated in MC
	inline bool doAmpCorrelationMC() const { return (ampCorrelationMC_); }

	//! "Master" detector for amp. correlation in MC
	inline bool isMaster() const { return (isMaster_); }

	//! Returns "associate" detector (e.g. for MC amplitude correlation)
	inline CsDetector* getAssociateDet() const { return (associateDet_); }

	//! Sets "associate" detector (e.g. for MC amplitude correlation)
	inline void setAssociateDet(CsDetector &det) { associateDet_ = &det; }

	//! Returns true if amplitude simulation in MC decoding is used
	inline bool amplitudeMCDecoding() const { return (amplitudeMCDecoding_); }

	//! Clusterize the digits
	void clusterize();
	void clusterizePixels();
	void clusterizeRectPixels();
	void clusterizeStrips();

	/// Decode raw data
	virtual void DecodeChipDigit(const CS::Chip::Digit &digit);

	//! Read calibration data
	virtual void readCalibration(time_t timePoint);

	//! Map channel calibrations
	bool mapChannelCal();

	//! Map channel calibrations for MC
	bool mapMCChannelCal();

	/*! \fn double getCorrU(const CsCluster &c, const double x, const double y, const double t, bool &error)
	      \brief Returns U corrected for propagation time, given MRS (x,y), and track's latency tt (fast method).
	      \param c    Input CsCluster
	      \param x,y  MRS, in mm
	      \param tt   Time offset
	 */
	double getCorrU(const CsCluster *c,
			const double x, const double y, const double tt,
			bool &error);

	//! return the complementary detector (the detector corresponding to the strips for pixel one and conversely)
	CsPixelMumegaDetector* getComplementaryDet() const
	{
		return complementary_det;
	}

	//! set the complementary detector
	void setComplementaryDet(CsPixelMumegaDetector &complementaryDet)
	{
		complementary_det = &complementaryDet;
	}

	//! return the complementary detector flag
	bool isComplDetSearched() const {
		return compl_det_searched;
	}

	//! set the complementary detector flag
	void setComplDetSearched(bool complDetSearched) {
		compl_det_searched = complDetSearched;
	}

	//! returns time resolution
	inline double getTRes() const { return (tRes_); }

	//! Correction to the amplitude in order to match counterpart (X vs. Y or V vs. U)
	inline const float *getAmpCorr() const { return fAmpCorr; }


	/*! Check if the point with the coordinates (u;v) is in pix area
	 * 	if complementary detector is found this function will be
	 * 	based on the pix_address table for the pix area shape
	 * 	else it will be based on default area and send an Error message.
	 	\brief return true if the point is in pix area, else return false
	 	\param u, v the coordinates of the point to be tested
	 */
	bool in_pix(const double u, const double v);

	/*! calculate the wire number for a coordinate in U in WRS (mm)
		 	\brief return the wire number if it found, -1 else
		 	\param t the coordinate of the point
	 */
	int DistToWire(const double t, bool &error );

	/*! calculate the coordinate in U in WRS (mm) from wire number
			 	\brief return the coordinate in U in WRS (mm)
			 	\param t the wire number
	 */
	double WireToDist(const int wire, bool &error );

	//! \returns the pitch associated to a given wire
	double getWirePitch(const int wire);

	//! \returns the percent of charge inside [min;max] for the cluster
	double getErfFact(const double min, const double center, const double max, const double clust_size);

private:
	bool          decodeCard_; //!< Decoding set by cards
	bool amplitudeMCDecoding_; //!< Amplitude simulation in MC
	bool    ampCorrelationMC_; //!< Amplitude correlation in MC
	bool            isMaster_; //!< True if "master" detector
	CsDetector* associateDet_; //!< Pointer to the coupled plane

// 	int       fPixelsVersion_; //!< Version of pixel scheme (1 normal, 2 simplified)
	double             spSig_; //!< space resolution, mm
	double             eGain_; //!< effective gain, a.u.
	double             eGSig_; //!< gain sigma, a.u.
	double            sWidth_; //!< sigma width, mm
	double              tRes_; //!< time resolution, ns

	std::vector<APVCal>    apv_chan_cals; //!< pedestals and sigmas for each APV chip
	CsMumegaTimeCals           time_cals; // timing calibrations

	std::map<std::string, CsHist1D*> mH1; //! 1D histogram pointers
	std::map<std::string, CsHist2D*> mH2; //! 2D histogram pointers
	bool                          sparse; //! true, if zero suppression is ON
	int                    do_clustering; //! clustering procedure switch
	int               fClusterConfigMask; //! clustering configuration
	int            do_amplituderatiocuts; //! amplitude ratio cuts switch

	float                 fThresholdHit_; //! Threshold for hits in units of sigma
	float                 fThresholdClu_; //! Threshold for clusters in units of sigma
	float           fLowMult_, fUppMult_; //! Cuts on multiplicity
	std::vector<float> fAmpRatio13_, fAmpRatio23_; //! Polygon for cut on amplitude ratio
	std::vector<float> fRollOverCor13_1, fRollOverCor23_1; //! Polygon for roll over correction
	std::vector<float> fRollOverCor13_2, fRollOverCor23_2; //! Polygon for roll over correction
	float         fMCHitTMin, fMCHitTMax; //! Cuts on hit time (for Monte Carlo)
	float                    fAmpCorr[3]; //! Correction to the amplitude in order to match counterpart (X vs. Y or V vs. U)
	float            fCrossTalkParams[2]; //! Cross talk suppression algorithm parameters
	// float            fTimeCrossTalkParam; //! Cross talk in multiplexed amalog signal suppression alogorithm parameter
	std::vector<float> fTimeCrossTalkParamsThrPix; //! Cross talk in multiplexed amalog signal suppression alogorithm parameters for pixel plane chips (Threshold)
	std::vector<float> fTimeCrossTalkParamsThrStrips; //! Cross talk in multiplexed amalog signal suppression alogorithm parameters for strip plane chips (Threshold)
	std::vector<float> fTimeCrossTalkParamsPeakPix; //! Cross talk in multiplexed amalog signal suppression alogorithm parameters for pixel plane chips (Peak)
	std::vector<float> fTimeCrossTalkParamsPeakStrips; //! Cross talk in multiplexed amalog signal suppression alogorithm parameters for strip plane chips (Peak)
	Float_t fClusterMergingAmpratioThr; //! Threshold on amplitude ratio for cluster merging
	std::vector<float> fTimeCuts;

	// const int fNchips;
	// float fTimeCrossTalkParamsPix[fNchips];
	// float fTimeCrossTalkParamsStrips[fNchips];

	CsMumegaPlane            *stripplane; //! Contains hits and clusters
	CsPixelMumegaPlane       *pixelplane; //! Contains hits and clusters
	CsRectPixelMumegaPlane  *rectpixelplane; //! Contains hits and clusters
	pixmm                   *pixmmobject; //! pixmm object in DaqDataDecoding/PixelMMDecoding.h for pixel positions decoding
	int                      fDebugLevel; //! Debug level

	bool            apv_chan_cals_mapped; //! Flag to indicate mapping of channel calibrations

#define wirDU_   wirD_   //!< Wires distance = diff 1st_wire-center
#define angU_    ang_    //!< Wires angle
#define sinAngU_ sinAng_ //!< sin ( wires angle )
#define cosAngU_ cosAng_ //!< cos ( wires angle )
#define nWirU_   nWir_   //!< Number of wires
#define wirPU_   wirP_   //!< Wire pitch

	double                        wirDV_; //!< Wires distance      in 2nd dimension
	double                         angV_; //!< Wires angle         in 2nd dimension
	double                      sinAngV_; //!< sin ( wires angle ) in 2nd dimension
	double                      cosAngV_; //!< cos ( wires angle ) in 2nd dimension
	int                           nWirV_; //!< Number of wires     in 2nd dimension
	double                        wirPV_; //!< Wire pitch          in 2nd dimension

	double *corr_U_tab; //!< correction table for the cluster's position
	int corr_U_nb_u; //!< number of bin in u for the correction
	int corr_U_nb_v; //!< number of bin in v for the correction
	double corr_U_size_u;//!< size of the bin in u for the correction
	double corr_U_size_v;//!< size of the bin in u for the correction
	double corr_U_shift_u; //!< additional correction for the detector position in u
	double corr_U_shift_v;//!< additional correction for the detector position in v
	double corr_U_phi; //!< Additional correction for the detector angle
	double corr_U_center_u; //!< detector centre U position (DRS) (cm)
	double corr_U_center_v; //!< detector centre V position (DRS) (cm)
	double corr_U_angle; //!< the Wires rotation angle of WRS respect to MRS after correction
	double corr_U_cos_angle;//!< cos ( wires angle ) after correction
	double corr_U_sin_angle; //!< sin ( wires angle ) after correction
	bool corr_U_bool; //!< boolean for enable correction

	/// default value for MC reconstruction /// (mm)

	/// detector resolution in U /////
	double pix_res;
	double strip_res1;// detector resolution for short strip zone (top : V > 0.0)
	double strip_res2;// detector resolution for short strip zone (bottom : V < 0.0)
	double strip_res3;// detector resolution for long strip zone (left : U < 0.0)
	double strip_res4;// detector resolution for long strip zone (right : U > 0.0)
	///////////

	//// cluster size of the detector ////
	double clus_size_mean; //cluster size mean (mm)
	double clus_size_mean_const;
	double clus_size_mean_slope;
	double clus_size_sig; //cluster size sigma (mm)
	///////////

	//// coefficient for relation between sigma of electron shower and cluster size ////
	double elec_spread_mean;
	double elec_spread_sigma;
	double elec_spread_sigmabg;
	///////////

	//// table for half pixel to pixel adress ////
// 	static const int nbpixY=40, nbpixX=128;
// 	static int pix_address [nbpixX][nbpixY];
// 	static bool pix_initialized;
	///////////

	//// amplitude of the third sample ////
	double amp2_mean; //amplitude 2 mean for 2.2kev deposit (for a MIP)
	double amp2_sig; //amplitude 2 sigma
	double amp2_fact; // factor between amplitude 2 and deposited energy (ref = 1.0 for 300V 600V 900V 1200V)
	///////////

	CsPixelMumegaDetector* complementary_det;
	bool compl_det_searched;
};

namespace CsPixelMumega {
  bool IsInside(const float xp, const float yp, const std::vector<float>& x, const std::vector<float>& y);
} // CsPixelMumega

#endif
