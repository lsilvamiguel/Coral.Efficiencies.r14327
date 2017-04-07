// $Id: CsSiTrackerDetector.h 13197 2012-01-13 16:31:22Z suhl $

/*!
   \file    CsSiTrackerDetector.h
   \brief   Compass Silicon Tracker like detector Class.
   \author  Benigno Gobbo
   \version $Revision: 13197 $
   \date    $Date: 2012-01-13 17:31:22 +0100 (Fri, 13 Jan 2012) $
*/

#ifndef CsSiTrackerDetector_h
#define CsSiTrackerDetector_h

#include "coral_config.h"

#include <CLHEP/Matrix/Matrix.h>
#include <list>
#include "CsDetector.h"
#include "CsHist.h"
#include "TTree.h"

#include "silicon_timing.h"

class CsZone;
class SiDigit;
struct cluster_timing;

/*! \class CsSiTrackerDetector 
    \brief Compass Silicon Tracker like detector Class.
*/

class CsSiTrackerDetector : public CsDetector {

 public:

  enum {ChipChannels=128};    // channels per chip

  class Calib
    {
    public:
      
      class Channel
        {
	public:
	  Channel(int f,double pm,double ps,double cm,double cs) :
	    flag(f),
	    pedestal_mean(pm), pedestal_sigma(ps),
	    calibration_mean(cm), calibration_sigma(cs)
	    {}
	  
	  Channel(const char *s)
            {
              if( 5!=sscanf(s,"%d%lf%lf%lf%lf",&flag,&pedestal_mean,&pedestal_sigma,
			    &calibration_mean,&calibration_sigma) )
                throw CS::Exception("CsSiTrackerDetector::Calib::Channel::Channel(): bad line \"%s\"",s);
            }
	  
	  void Print(std::ostream &o=std::cout,const std::string &prefix="") const
            {
              o<<prefix;
	      o << "flag=" << flag
		<< " PED=(" 
		<< std::setw(8) << std::setprecision(3) << pedestal_mean
		<< ","
		<< std::setw(8) << std::setprecision(3) << pedestal_sigma
		<< "  CALIB=("
		<< std::setw(8) << std::setprecision(3) << calibration_mean
		<< ","
		<< std::setw(8) << std::setprecision(3) << calibration_sigma
		<< ")" << std::endl;
            }

	  int    flag;
	  double pedestal_mean;
	  double pedestal_sigma;
	  double calibration_mean;
	  double calibration_sigma;
        };


      // constructor
      Calib() : initialised(false) {};


      // accessors
      bool IsInitialised() const { return initialised; };
      int GetSrcId() const {
        if (!initialised) throw CS::Exception("uninitialised use of CsSiTrackerDetector::Calib()");
        return src_id;
      };
      int GetAdcId() const {
        if (!initialised) throw CS::Exception("uninitialised use of CsSiTrackerDetector::Calib()");
        return adc_id;
      };
      const std::vector<Channel>& GetChannels() const {
        if (!initialised) throw CS::Exception("uninitialised use of CsSiTrackerDetector::Calib()");
        return channels;
      }

    private:
      bool initialised;     // true when object has been fully initialised
      int src_id;
      int adc_id;
      
      std::vector<Channel> channels;
      
      friend std::istream & operator >> (std::istream &,std::vector<Calib> &c);
    };


  /*! \fn CsSiTrackerDetector( const int row, 
    const int id, const char* name, const int unit,  
    const int type, const double rdLen, const double xsiz,  
    const double ysiz,  const double zsiz, const double xcm,   
    const double ycm, const double zcm, const HepMatrix rotDRS,
    const HepMatrix rotWRS, const double wirD, const double ang, 
    const int nWir, const double wirP, const double eff, const double bkg,
    const double tGate );
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
    \param wirP  Wires pitch
    \param eff   Detector efficiency
    \param bkg   Detector background
    \param tGate Detector gate
  */
  CsSiTrackerDetector( const int    row,
	      const int    id,    const char* name,    const char *TBname,
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


 private:
  bool operator==( const CsSiTrackerDetector& ) const; //!< "equal to" operator
  
  bool operator<( const CsSiTrackerDetector& ) const; //!< "less than" operator

  //! Returns true if TCD informations are available
  inline bool hasTDC() const { return( false ); }

  //! Returns true if Drift informations are available
  inline bool hasDrift() const { return( false ); }

  //! Book histograms if needed
  void BookHistograms();

  //! Decode the MC data for this detector
  void makeMCDecoding();


  //! Returns true if amplitude simulation in MC decoding is used
  inline bool amplitudeMCDecoding() const { return( amplitudeMCDecoding_ ); }


  bool  amplitudeMCDecoding_;      //!< Amplitude simulation in MC 

  double spSig_; //!< space resolution, mm, for simulation
  double eGain_; //!< deffective gain, a.u., for simulation
  double eGSig_; //!< gain sigma, a.u., for simulation
  double sWidth_; //!< signal width, mm, for simulation
  double tRes_; //!< time resolution, ns, for simulation

  double eta_[8]; // charge sharing distribution parameters, for simulation
  double cls_[7]; // clustersize distribution parameters, for simulation
  double cls3_[7]; // clustersize distribution parameters, for simulation

  double cls( double *, double *); // used to determine clustersize in simulation
  double eta( double *, double *); // used to determine charge sharing in simulation
  double cls3( double *, double *); // used to determine clustersize in simulation



  //! Clusterize the digits
  void clusterize();
  void clusterizeCind();
  void clusterizeRatio();
  
  //! Calculate time from 3 ADC values of cluster and store it into the same cluster 
  double decodeTime(CsCluster* cluster);

  //! Decode raw data
  virtual void DecodeChipDigit(const CS::Chip::Digit &digit);

  //! Read calibration data from CDB
  virtual void readCalibration(time_t timePoint);

  void SetupAlignmentPatch();
  bool align_patch;               // true when Sergei's alignment patch is enabled
  unsigned align_patch_run;       // run number for which alignment patch has been setup

  void SetupCinderellaClustering();
  bool cind_clus_setup;           // true when Cinderella clustering has been setup
  double params_noise[1300];       // input for Cinderella clustering
  double params_pedestal[1300];
  double pedestal_mean_chip[10];
  
  // Computes convariance matrix
  void comp_cov(CLHEP::HepMatrix&, CLHEP::HepMatrix&, double &);
  
  // checks if time jumps were found for the processed event
  int find_jump();
  
  // Timing Calibration and coefficients of the eta correction are written
  // in the according structures
  void get_calib(int);
  
  // for MC case: depending on the double and detector type, charge
  // sharing factor is calculated
  double get_charge_sharing(double);

			
   // options if mate clutering is used
   bool cluster_mate; 	
   bool cluster_mate_opt;		

   double shft0;			 // timing shift depending on cluster size
   double shft1;			 // timing shift depending on cluster size

   // structure used in cinderella clustering
   time_reconstruction_option_t params_time; 
 
	std::list<SiDigit> sort_digits(const std::list<CsDigit*>& digits);
	
	// Digits are read and filled in array of strips,
	// needed for cinderella clustering
	inline void fill_strips(const std::list<SiDigit>& siDigits, strip_t*);
	
	// apply timing, eta correction and postion to clusters. 
	// output is a vector of clusters which can be used in mate_clustering	
	void finalize_clusters(cluster_timing&, CLHEP::HepMatrix&, double, std::list<CsCluster*>&, const std::list<SiDigit>&);
	
	// depending on detector type and clustersize, the according covariance matrix is returned
	inline void get_spatialerror(CLHEP::HepMatrix&);
	
	// correcting for residual against eta distribution
	inline void apply_etacorr(double& u, double r);
	
	// mate clustering
	void mate_clustering(std::list<CsCluster*> & , std::list<CsCluster*> & ,  double);
	
	// cluster is splitted according to ration of given ints
	void split_cluster(std::list<CsCluster*> & , std::list<CsCluster*>::iterator , int , int , double );
	
	// if no appropriate "splitting partners" are found, 3 strip clusters are splitting with factor 0.5
	void split_cluster(std::list<CsCluster*> & , std::list<CsCluster*>::iterator , double );
	
	// if no appropriate "splitting partners" are found, and clustersize is 4 or lager, two new clusters 
	// are created using the two first and the two last strips of the old cluster
	void split_delta(std::list<CsCluster*> &, std::list<CsCluster*>::iterator, double);
	
	// pointer to mate detector, which is the plane on the same waver. this is used for the amplitude
	// correalation. In constructor, this initialized to NULL. The pointer is set inside clusterize()
	CsSiTrackerDetector* mateDet;
        
    // for a given cluster, one looks for clusters which match the amplitude sum. returns false
    // if match is poor    
	bool find_best_candidates(std::list<CsCluster*> & ,std::list<std::list<CsCluster*>::iterator>::iterator , int &, int& );
	
	double wireDCorr_;

  //Tree for debugung (hist_level high)
  TTree *T_SIdeb_Ccl;

  double amp0; // amplitude of SI digit 1. sample
  double amp1; // amplitude of SI digit 2. sample
  double amp2; // amplitude of SI digit 3. sample
  
  double r0; //ratio of amplitudes amp0/amp2 
  double r1; //ratio of amplitudes amp1/amp2
  
  double hit_time_0; double hit_time_err_0; // time of hit in ratio clustering of ratio 0
  double hit_time_1; double hit_time_err_1; // time of hit in ratio clustering of ratio 1
  double hit_time;   double hit_time_err; // time of hit in ratio clustering of hit
  
  double tcs_cor;
  double tcs_cor_c;
  double time_deb; // time of cluster 
  double time_deb_t0;
  double time_deb_t1;

  double time_clus;
  double time_clus_t0;
  double time_clus_t1;
  double pos_out;
  double time_out;
  double clustertime_err_out;
  double amp2_out; //amplidude 2 of cluster
  double amp1_out; //amplidude 1 of cluster
  double amp0_out; //amplidude 0 of cluster
  double r1_out; //ratio of cluster amlplitude  amp1_out/amp2_out
  double r0_out; //ratio of cluster amlplitude  amp0_out/amp2_out

  unsigned int first_wire_out; // wire number of first (most little) strip in one cluster after CINDERELLA clustering 
  unsigned int cluster_size_out; // cluster size after CINDERELLA clustering 
  
  int cSize; // cluster size after ratio clustering
  
  param_t p_ratio1;
  param_t p_ratio2;
  param_t *param[2];
  
  eta_t p_eta;
  

  double time_thiemo;
  bool cluster_cind_opt; // value set in the option file, to swich on and of cinderella clustering
  bool time_calib_munich;
  bool time_calib_db;
  bool time_12nsjumps;      // set in option file using "SI Time12nsJumps on"
  bool eta_correction;
  bool eta_calib;
  std::string calib_dir;

  unsigned int params_strip_count; // number to be processed for this event

  int i_cl; //cluster counter in cindarella clustering


 
  bool           decodeCard_;      //!<         Decoding set by cards

  double         fMultCut;         //!< For detector checks:  clean event

  bool           fTimeCut;          //!< switch for ratio cut (default: on)
  bool           fRatioCut;         //!< switch for expl. time cut (def: off)
  double         fLCut,fRCut;       //!< left, right limits for time cut (ns)
  
  unsigned int ChipsPerPlane;       ///number of chips calculated for each plan
  uint32 fEventNr[10];               //!< memory for Common Mode Noise
  bool sparse;                       //!< switch set by DaqDataDecoding

  bool                       has_calib_data;
  std::vector<Calib>         calib_data;    //!< channel-wise calibs (pedestals, sigmas etc), grouped by chips,
                                            // ordered by consecutive chip #
  std::vector<double>        timing_calib;  //!< time measurement related constants
  std::vector<int>           timing_spill_jumps;     // Correction per spill for planes with "jumping" timing
  std::vector<double> eta_corr; 

  class Histograms {

    std::map<std::string,CsHist1D*> hist1;
    std::map<std::string,CsHist2D*> hist2;

  public:

    CsHist1D* operator()(std::string,int,double,double);
    CsHist2D* operator()(std::string,int,double,double,int,double,double);
    void operator()(std::string,double);
    void operator()(std::string,double,double);

    ~Histograms() {std::cout<<"CsSiTrackerDetector::Histograms::~Histograms"<<std::endl;}

  } h;
  
  // boolinas used in clustering method and different spatial resolution 
  bool isX;
  bool isU;
  bool isY;
  bool isV;

  double rescor;   // Sergei's correction

};





#endif // CsSiTrackerDetector_h
