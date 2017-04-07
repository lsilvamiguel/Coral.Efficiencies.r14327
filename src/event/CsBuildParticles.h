// $Id: CsBuildParticles.h,v 1.12 2010/11/17 16:42:12 tnagel Exp $

/*!
   \file    BuildParticles.h
   \brief   Track and Calorimeter Cluster Association (Singleton).
   \version $Revision: 1.12 $
   \date    $Date: 2010/11/17 16:42:12 $
*/

#ifndef CsBuildParticles_h
#define CsBuildParticles_h

#include <vector>

#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "CsHist.h"

class CsCalorimeter;
class CsHelix;
class CsParticle;
class CsRecoEvent;
class CsTrack;
namespace Reco {
  class CalorimeterParticle;
}

/// Track and Calorimeter Cluster Association (Singleton)
class CsBuildParticles {
 public:
  static CsBuildParticles* Instance();

  /// Associate tracks with calorimeter clusters.
  void Build( const CsRecoEvent& event, std::vector<CsParticle*>& parts );

  /// Calculate calorimeter cluster Z-positions
  void SetClusterZ( std::vector<CsParticle*>& parts ) const;

 private:

  /// constructor is private (singleton), use CsBuildParticles::Instance() to
  /// access the instance
  CsBuildParticles();

  static CsBuildParticles *_instance;

  class BPCalo;
  std::vector<BPCalo> _calos;


  /// Everything that belongs to a calorimeter: options, histograms, pointer
  /// to CsCalorimeter
  class BPCalo {
  public:
    /// tiny class to bundle a calorimeter cluster together with the information
    /// whether it has been associated already
    class BPCluster {
    public:
      BPCluster(Reco::CalorimeterParticle* _p) : clus(_p), associated(false) {};
      
      Reco::CalorimeterParticle *clus;
      bool associated;
    };  // end class BPCluster


    //
    // member functions
    //
    
    BPCalo( CsCalorimeter *c, int levelBPHistos_, int levelCellTimeHistos_ );

    void FillClustersHists( void );
    void FillClusterTrackHists( const Reco::CalorimeterParticle* clus, 
				double exclus, double eyclus,
				double xtrk, double ytrk,
				double momentum );
    void FillClusterTrackAssocHists( const Reco::CalorimeterParticle* clus,
				     const CsTrack *trk, double xtrk, double ytrk,
				     double momentum );
    void FillTrackHists( const CsHelix &helix0, double momentum,
			 bool track_associated );


    //
    // member variables
    //
 
    CsCalorimeter *cscalo;

    /// histogramming levels
    int levelBPHistos;
    int levelCellTimeHistos;

    /// clusters
    std::vector<BPCluster> clusters;

    /// constant: max dist. for association (in multiples of spatial error)
    double Nsigma;

    /// constant: to make some PID(electrons) cut on energy - track momentum
    double NsigmaEcut;
    
    CsHist1D* caloNData;        //!<  Number of cells with data for calorimeter
    CsHist1D* caloN;            //!<  Number of reconstructed clusters for calorimeter
    CsHist1D* clusSize;         //!<  Cluster size
    CsHist1D* caloE;            //!<  Cluster energy for all calorimeter clusters
    CsHist1D* caloEerr;         //!<  Error of cluster energy for all calorimeter clusters
    CsHist1D* caloCellsE;       //!<  Energy for all calorimeter cells
    CsHist1D* caloTime;         //!<  Time Calorimeter for all calorimeter clusters
    CsHist1D* caloTimeS;        //!<  Time Calorimeter/Time Error for all calorimeter clusters
    CsHist1D* caloZ;            //!<  Cluster Z position (assuming charged, interacting particle),
                                //!<  one entry per track
    CsHist1D* deltaX;           //!<  XCalorimeter - XTrack distribution for all clusters
    CsHist1D* deltaXEcut;       //!<  XCalorimeter - XTrack distribution for clusters (dE cut)
    CsHist1D* deltaXS;          //!<  (XCalorimeter - XTrack)/SigmaXCalo distribution for all clusters
    CsHist1D* deltaY;           //!<  YCalorimeter - YTrack distribution for all clusters
    CsHist1D* deltaYEcut;       //!<  YCalorimeter - YTrack distribution for clusters (dE cut)
    CsHist1D* deltaYS;          //!<  (YCalorimeter - YTrack)/SigmaYCalo distribution for all clusters
    CsHist1D* deltaT;           //!<  TCalorimeter - TTrack distribution for all clusters
    CsHist1D* deltaTEcut;       //!<  TCalorimeter - TTrack distribution for all clusters (dE cut)
    CsHist1D* deltaTS;          //!<  (TCalorimeter - TTrack)/(dTcalo+dTtrack) distribution for all clusters
    CsHist1D* deltaTSEcut;      //!<  (TCalorimeter - TTrack)/(dTcalo+dTtrack) distribution for all clusters (dE cut)

    CsHist1D* assocAttempt;     //!<  0: interacting particle, 1: MIP, 2: assoc. failed

    CsHist2D* caloXtrackX;      //!<  XCalorimeter vs XTrack distribution for all clusters
    CsHist2D* caloYtrackY;      //!<  YCalorimeter vs YTrack distribution for all clusters
    CsHist2D* caloXtrackXEcut;  //!<  XCalorimeter vs XTrack distribution for all clusters (dE cut)
    CsHist2D* caloYtrackYEcut;  //!<  YCalorimeter vs YTrack distribution for all clusters (dE cut)
    CsHist2D* caloXY;           //!<  YCalorimeter vs XCalorimeter distribution for all clusters
    CsHist2D* cellsXY;          //!<  YCalorimeter vs XCalorimeter distribution for all clusters cells baised
    CsHist2D* trackXY;          //!<  YTrack vs XTrack distribution for all tracks which hit calorimeter
    CsHist2D* deltaXtrackX;     //!<  XTrack vs (XCalorimeter - XTrack) distribution for all clusters
    CsHist2D* deltaYtrackY;     //!<  YTrack vs (YCalorimeter - YTrack) distribution for all clusters
    CsHist2D* deltaXcaloX;      //!<  XCalo  vs (XCalorimeter - XTrack) distribution for all clusters
    CsHist2D* deltaYcaloY;      //!<  YCalo  vs (YCalorimeter - YTrack) distribution for all clusters
    CsHist2D* deltaXEcuttrackX; //!<  XTrack vs (XCalorimeter - XTrack) distribution for clusters (dE cut)
    CsHist2D* deltaYEcuttrackY; //!<  YTrack vs (YCalorimeter - YTrack) distribution for clusters (dE cut)
    CsHist2D* deltaXEcutcaloX;  //!<  XCalo  vs (XCalorimeter - XTrack) distribution for clusters (dE cut)
    CsHist2D* deltaYEcutcaloY;  //!<  YCalo  vs (YCalorimeter - YTrack) distribution for clusters (dE cut)
    CsHist2D* trackXYas;        //!<  YTrack vs XTrack distribution for all tracks associated with calorimeter cluster
    CsHist2D* trackXYnotas;     //!<  YTrack vs XTrack distribution for all tracks which hit but not associated with calorimeter cluster
    CsHist2D* trackXYas_E;      //!<  YTrack vs XTrack distribution for all tracks with momentum associated with calorimeter cluster
    CsHist2D* trackXYnotas_E;   //!<  YTrack vs XTrack distribution for all tracks with momentum not associated with calorimeter cluster
    CsHist2D* trackXYas_noE;    //!<  YTrack vs XTrack distribution for all tracks without momentum associated with calorimeter cluster
    CsHist2D* trackXYnotas_noE; //!<  YTrack vs XTrack distribution for all tracks without momentum not associated with calorimeter cluster
    CsHist2D* caloEtrackE;      //!<  ETrack vs ECalorimeter for all tracks associated with calorimeter cluster
    CsHist1D* deltaEs;          //!<  (ECalorimeter - ETrack)/ETrack distribution for all clusters
    CsHist2D* deltaEtrackE;     //!<  ETrack vs (Ecalorimeter - Etrack)
    CsHist2D* deltaEtrackX;     //!<  XTrack vs (Ecalorimeter - Etrack)
    CsHist2D* deltaEtrackY;     //!<  YTrack vs (Ecalorimeter - Etrack)
    CsHist2D* caloTtrackT;      //!<  Time Track vs Time Calorimeter for all tracks associated with calorimeter cluster
    CsHist2D* caloTtcsT;        //!<  Time TCS vs Time Calorimeter for all tracks associated with calorimeter cluster
    CsHist2D* caloTcaloE;       //!<  Time Track vs ECalorimeter 
    CsHist2D* dataCHIcaloE;     //!<  Chisq main cell vs ECalorimeter 
    CsHist2D* neutralXY;        //!<  Y vs X distribution of neutral particles produced from not associated calorimeter clusters

    CsHist2D* deltaXYmip;       //!<  XY distance between tracks at "interacting particle" Z and "MIP" Z

    TProfile*   effTrackAss;    //!<  Probability of association of charged track with calorimeter cluster
    TProfile*   effCaloAss;     //!<  Probability of association of calorimeter cluster with charged track   
    TProfile2D* effTrackAss2D;  //!<  Probability of association of charged track with calorimeter cluster
    TProfile2D* effCaloAss2D;   //!<  Probability of association of calorimeter cluster with charged track

    // anonymous enum to set NPdep constant
    enum { NPdep = 8 };

    double      Pdep               [NPdep];  //!< upper limit for momentum bin
    TProfile2D* effTrackAssPdep    [NPdep];
    CsHist1D*   deltaEsPdep        [NPdep];
    CsHist1D*   deltaXPdep         [NPdep];
    CsHist1D*   deltaXEcutPdep     [NPdep];
    CsHist1D*   deltaYPdep         [NPdep];
    CsHist1D*   deltaYEcutPdep     [NPdep];
    CsHist2D*   deltaXcaloXPdep    [NPdep];
    CsHist2D*   deltaYcaloYPdep    [NPdep];
    CsHist2D*   deltaXEcutcaloXPdep[NPdep];
    CsHist2D*   deltaYEcutcaloYPdep[NPdep];
    
    std::vector<CsHist1D*> caloCellsTime;
  };  // end class BPCalo
};  // end class CsBuildParticles

#endif
