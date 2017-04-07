// $Id: CsBuildParticles.cc,v 1.30 2011/02/01 22:05:53 tnagel Exp $

/*!
   \file    BuildParticles.h
   \brief   Track and Calorimeter Cluster Association (Singleton).
   \version $Revision: 1.30 $
   \date    $Date: 2011/02/01 22:05:53 $
*/


#include <string>

#include "CsBuildParticles.h"

#include "CsBeam.h"
#include "CsCalorimeter.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsOpt.h"
#include "CsParticle.h"

#include "Reco/Calorimeter.h"
#include "Reco/CalorimeterParticle.h"
#include "Reco/CellDataRaw.h"
#include "Reco/Shower.h"

using namespace std;


CsBuildParticles *CsBuildParticles::_instance = NULL;


CsBuildParticles *CsBuildParticles::Instance() {
  if ( ! _instance )
    _instance = new CsBuildParticles();
  return _instance;
}


//! Constructor.  Mostly histogram booking.
CsBuildParticles::BPCalo::BPCalo( CsCalorimeter *c, int levelBPHistos_, int levelCellTimeHistos_ ) :
  cscalo(c),
  levelBPHistos(levelBPHistos_),
  levelCellTimeHistos(levelCellTimeHistos_),
  Nsigma( c->GetOptions().GetMiscDoubleOption(1) ),
  NsigmaEcut(3.),
  caloNData(NULL),
  caloN(NULL),
  clusSize(NULL),
  caloE(NULL),
  caloEerr(NULL),
  caloCellsE(NULL),
  caloTime(NULL),
  caloTimeS(NULL),
  caloZ(NULL),
  deltaX(NULL),
  deltaXEcut(NULL),
  deltaXS(NULL),
  deltaY(NULL),
  deltaYEcut(NULL),
  deltaYS(NULL),
  deltaT(NULL),
  deltaTEcut(NULL),
  deltaTS(NULL),
  deltaTSEcut(NULL),
  assocAttempt(NULL),
  caloXtrackX(NULL),
  caloYtrackY(NULL),
  caloXtrackXEcut(NULL),
  caloYtrackYEcut(NULL),
  caloXY(NULL),
  cellsXY(NULL),
  trackXY(NULL),
  deltaXtrackX(NULL),
  deltaYtrackY(NULL),
  deltaXcaloX(NULL),
  deltaYcaloY(NULL),
  deltaXEcuttrackX(NULL),
  deltaYEcuttrackY(NULL),
  deltaXEcutcaloX(NULL),
  deltaYEcutcaloY(NULL),
  trackXYas(NULL),
  trackXYnotas(NULL),
  trackXYas_E(NULL),
  trackXYnotas_E(NULL),
  trackXYas_noE(NULL),
  trackXYnotas_noE(NULL),
  caloEtrackE(NULL),
  deltaEs(NULL),
  deltaEtrackE(NULL),
  deltaEtrackX(NULL),
  deltaEtrackY(NULL),
  caloTtrackT(NULL),
  caloTtcsT(NULL),
  caloTcaloE(NULL),
  dataCHIcaloE(NULL),
  neutralXY(NULL),
  deltaXYmip(NULL),
  effTrackAss(NULL),
  effCaloAss(NULL),
  effTrackAss2D(NULL),
  effCaloAss2D(NULL)
{
  for (int i=0; i<NPdep; i++) {
    effTrackAssPdep    [i] = NULL;
    deltaEsPdep        [i] = NULL;
    deltaXPdep         [i] = NULL;
    deltaXEcutPdep     [i] = NULL;
    deltaYPdep         [i] = NULL;
    deltaYEcutPdep     [i] = NULL;
    deltaXcaloXPdep    [i] = NULL;
    deltaYcaloYPdep    [i] = NULL;
    deltaXEcutcaloXPdep[i] = NULL;
    deltaYEcutcaloYPdep[i] = NULL;
  }
  Pdep[0] =  0.01;
  Pdep[1] =  1.;
  Pdep[2] =  5.;
  Pdep[3] = 10.;
  Pdep[4] = 20.;
  Pdep[5] = 45.;
  Pdep[6] = 90.;
  Pdep[7] = 180.;

  cout << c->GetName() << ": track association Nsigma=" << Nsigma << endl;

  if (!levelBPHistos)
    return;

  const char *cname = c->GetName().c_str();
  string dir_name = string("/CsBPMonitor/Calorimeter_") + cname;
  CsHistograms::SetCurrentPath(dir_name);

  const double xMin = c->GetXmin() + c->GetPositionX();
  const double xMax = c->GetXmax() + c->GetPositionX();
  const double yMin = c->GetYmin() + c->GetPositionY();
  const double yMax = c->GetYmax() + c->GetPositionY();

  caloNData  = new CsHist1D("NData",        string("Number of cells with data (")+cname+")",
			    2500, 0., 2500.);
  caloN      = new CsHist1D("NClusters",    string("Number of clusters (")+cname+")",
			    500, 0., 500.);
  clusSize   = new CsHist1D("Cluster Size", string("Cluster Size (")+cname+")",
			    50, -.5, 49.5);
  caloE      = new CsHist1D("Energy",       string("Cluster Energy (")+cname+")",
			    2500, 0., 250.);
  caloEerr   = new CsHist1D("Energy Error", string("Cluster Energy Error (")+cname+")",
			    2500, 0., 250.);
  caloCellsE = new CsHist1D("CellsEnergy",  string("Cells Energy (")+cname+")",
			    2500, 0., 250.);
  deltaEs    = new CsHist1D("deltaEs",      string("ECalorimeter/ETrack -1 associated (")+cname+")",
			    200, -0.5, 0.5);
  caloTime   = new CsHist1D("Time",         string("Cluster Time (")+cname+")",
			    2000, -1000., 1000.);
  caloTimeS  = new CsHist1D("Time/dTime",   string("Cluster Time (")+cname+")",
			    600, -30., 30.);
  caloZ      = new CsHist1D("Z position",   string("Z Position (")+cname+")",
			    4000, -5000., 35000.);
  deltaX     = new CsHist1D("deltaX",       string("XCalorimeter - XTrack (")+cname+")",
			    1000, -500., 500.);
  deltaXEcut = new CsHist1D("deltaXEcut",   string("XCalorimeter - XTrack with dE cut (")+cname+")",
			    1000, -500., 500.);
  deltaXS    = new CsHist1D("deltaXS",      string("(XCalorimeter - XTrack)/SigmaXCalo (")+cname+")",
			    200, -10., 10.);
  deltaY     = new CsHist1D("deltaY",       string("YCalorimeter - YTrack (")+cname+")",
			    1000, -500., 500.);
  deltaYEcut = new CsHist1D("deltaYEcut",   string("YCalorimeter - YTrack with dE cut (")+cname+")",
			    1000, -500., 500.);
  deltaYS    = new CsHist1D("deltaYS",      string("(YCalorimeter - YTrack)/SigmaYCalo (")+cname+")",
			    200, -10., 10.);
  deltaT     = new CsHist1D("deltaT",       string("CalorimeterT - TrackT (")+cname+")",
			    600, -150., 150.);
  deltaTEcut = new CsHist1D("deltaTEcut",   string("CalorimeterT - TrackT with dE cut (")+cname+")",
			    600, -150., 150.);
  deltaTS    = new CsHist1D("deltaTS",      string("(CalorimeterT - TrackT)/(dcaloT+dtrackT) (")+cname+")",
			    600, -30., 30.);
  deltaTSEcut= new CsHist1D("deltaTSEcut",  string("(CalorimeterT - TrackT)/(dcaloT+dtrackT) with dE cut (")+cname+")",
			    600, -30., 30.);
  assocAttempt=new CsHist1D("assocAttempt", string("0: interacting, 1: MIP, 2: assoc. failed (")+cname+")",
			    3, -.5, 2.5);
  
  if(levelBPHistos > 1) {
    const int fact2Dim = levelBPHistos - 1;
    caloXtrackX = new CsHist2D("caloXtrackX", string("XCalorimeter vs XTrack (")+cname+")",
			       100*fact2Dim, xMin, xMax, 100*fact2Dim, xMin, xMax);
    caloYtrackY = new CsHist2D("caloYtrackY", string("YCalorimeter vs YTrack (")+cname+")",
			       100*fact2Dim, yMin, yMax, 100*fact2Dim, yMin, yMax);
    trackXY     = new CsHist2D("trackXY", string("Y vs X Track (")+cname+")",
			       100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    caloXY      = new CsHist2D("caloXY", string("Y vs X CaloCluster (")+cname+")",
			       100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    if( c->XYRegularGrid() ) {
      const int nx = c->GetNColumns(); 
      const int ny = c->GetNRows(); 
      cellsXY = new CsHist2D("cellsXY", string("Ycell vs Xcell CaloCluster (")+cname+")",
			     nx, -0.5, double(nx)-0.5, ny, -0.5, double(ny)-0.5);
    }
    deltaXtrackX     = new CsHist2D("deltaXtrackX", string("XTrack vs (XCalorimeter - XTrack) (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, xMin, xMax);
    deltaYtrackY     = new CsHist2D("deltaYtrackY", string("YTrack vs (YCalorimeter - YTrack) (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, yMin, yMax);
    deltaXcaloX      = new CsHist2D("deltaXcaloX", string("XCalo vs (XCalorimeter - XTrack) (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, xMin, xMax);
    deltaYcaloY      = new CsHist2D("deltaYcaloY", string("YCalo vs (YCalorimeter - YTrack) (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, yMin, yMax);
    deltaXEcuttrackX = new CsHist2D("deltaXEcuttrackX", string("XTrack vs (XCalorimeter - XTrack) with dE cut (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, xMin, xMax);
    deltaYEcuttrackY = new CsHist2D("deltaYEcuttrackY", string("YTrack vs (YCalorimeter - YTrack) with dE cut (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, yMin, yMax);
    deltaXEcutcaloX  = new CsHist2D("deltaXEcutcaloX", string("XCalo vs (XCalorimeter - XTrack) with dE cut (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, xMin, xMax);
    deltaYEcutcaloY  = new CsHist2D("deltaYEcutcaloY", string("YCalo vs (YCalorimeter - YTrack) with dE cut (")+cname+")",
				    100*fact2Dim, -500., 500., 200*fact2Dim, yMin, yMax);
    trackXYas        = new CsHist2D("trackXYas", string("YTrack vs XTrack associated (")+cname+")",
				    100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    trackXYnotas     = new CsHist2D("trackXYnotas", string("YTrack vs XTrack hit but not associated (")+cname+")",
				    100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    trackXYas_E      = new CsHist2D("trackXYas_E", string("YTrack vs XTrack( with momentum) associated (")+cname+")",
				    100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    trackXYnotas_E   = new CsHist2D("trackXYnotas_E", string("YTrack vs XTrack( with momentum) not associated (")+cname+")",
				    100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    trackXYas_noE    = new CsHist2D("trackXYas_noE", string("YTrack vs XTrack( without momentum) associated (")+cname+")",
				    100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    trackXYnotas_noE = new CsHist2D("trackXYnotas_noE", string("YTrack vs XTrack( without momentum) not associated (")+cname+")",
				    100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);
    deltaEtrackE     = new CsHist2D("deltaEtrackE", string("ETrack vs deltaE (")+cname+")",
				    100*fact2Dim, -50., 50, 100*fact2Dim, 0, 200);
    deltaEtrackX     = new CsHist2D("deltaEtrackX", string("XTrack vs deltaE (")+cname+")",
				    100*fact2Dim, -50., 50, 100*fact2Dim, xMin, xMax);
    deltaEtrackY     = new CsHist2D("deltaEtrackY", string("YTrack vs deltaE (")+cname+")",
				    100*fact2Dim, -50., 50, 100*fact2Dim, yMin, yMax);
    caloEtrackE      = new CsHist2D("caloEtrackE", string("ETrack vs ECalorimeter associated (")+cname+")",
				    100*fact2Dim, 0, 200, 100*fact2Dim, 0, 200);
    caloTtrackT      = new CsHist2D("caloTtrackT", string("TimeTrack vs TimeCalorimeter associated (")+cname+")",
				    100*fact2Dim, -200, 200, 100*fact2Dim, -200, 200);
    caloTtcsT        = new CsHist2D("caloTtcsT", string("TimeTCS vs TimeCalorimeter associated (")+cname+")",
				    100*fact2Dim, -50, 50, 100*fact2Dim, -50, 50);
    caloTcaloE       = new CsHist2D("caloTcaloE", string("ECalorimeter vs Time (")+cname+")",
				    400*fact2Dim, -100, 100, 100*fact2Dim, 0, 200);
    dataCHIcaloE     = new CsHist2D("dataCHIcaloE", string("ECalorimeter vs Chisq in main (")+cname+")",
				    100*fact2Dim, -5000., 5000., 100*fact2Dim, 0, 20);
    neutralXY        = new CsHist2D("neutralXY", string("Y vs X of neutral patricles (")+cname+")",
				    100*fact2Dim, xMin, xMax, 100*fact2Dim, yMin, yMax);

    deltaXYmip       = new CsHist2D("deltaXYmip", string("track XY @ MIP Z - track XY @ interacting particle Z (")+cname+")",
				    50*fact2Dim, -50., 50., 50*fact2Dim, -50., 50.);
    
    effTrackAss2D    = new TProfile2D( "effTrackAss2D",     (string("Probability of track -> cluster assoc (")+cname+")").c_str(),
				       25*fact2Dim, xMin, xMax, 25*fact2Dim, yMin, yMax);
    effCaloAss2D     = new TProfile2D( "effCaloAss2D",      (string("Probability of cluster -> track assoc (")+cname+")").c_str(),
				       25*fact2Dim, xMin, xMax, 25*fact2Dim, yMin, yMax);

    // plots with momentum slices
    for (int i=0; i<NPdep; i++) {
      stringstream s; s << i;
      stringstream srange;
      if ( i == 0 )
	srange << "p < " << Pdep[0];
      else if ( i == NPdep-1 )
	srange << "p > " << Pdep[NPdep-1];
      else
	srange << Pdep[i-1] << " < p < " << Pdep[i];

      effTrackAssPdep[i] = new TProfile2D( (string("effTrackAss3D")+s.str()).c_str(),
					   (string("Probability of track -> cluster assoc, "+srange.str()+" (")+cname+")").c_str(),
					   25*fact2Dim, xMin, xMax, 25*fact2Dim, yMin, yMax);

      deltaEsPdep[i]     = new CsHist1D( string("deltaEs")+s.str(),
					 string("ECalorimeter/ETrack -1 associated, "+srange.str()+" (")+cname+")",
					 200, -0.5, 0.5);
      deltaXPdep[i]      = new CsHist1D( string("deltaX")+s.str(),
					 string("XCalorimeter - XTrack, "+srange.str()+" (")+cname+")",
					 1000, -500., 500.);
      deltaXEcutPdep[i]  = new CsHist1D( string("deltaXEcut")+s.str(),
					 string("XCalorimeter - XTrack with dE cut, "+srange.str()+" (")+cname+")",
					 1000, -500., 500.);
      deltaYPdep[i]      = new CsHist1D( string("deltaY")+s.str(),
					 string("YCalorimeter - YTrack, "+srange.str()+" (")+cname+")",
					 1000, -500., 500.);
      deltaYEcutPdep[i]  = new CsHist1D( string("deltaYEcut")+s.str(),
					 string("YCalorimeter - YTrack with dE cut, "+srange.str()+" (")+cname+")",
					 1000, -500., 500.);

      deltaXcaloXPdep[i]     = new CsHist2D( string("deltaXcaloX")+s.str(),
					     string("XCalo vs (XCalorimeter - XTrack), "+srange.str()+" (")+cname+")",
					     100*fact2Dim, -500., 500., 200*fact2Dim, xMin, xMax);
      deltaYcaloYPdep[i]     = new CsHist2D( string("deltaYcaloY")+s.str(),
					     string("YCalo vs (YCalorimeter - YTrack), "+srange.str()+" (")+cname+")",
					     100*fact2Dim, -500., 500., 200*fact2Dim, yMin, yMax);
      deltaXEcutcaloXPdep[i] = new CsHist2D( string("deltaXEcutcaloX")+s.str(),
					     string("XCalo vs (XCalorimeter - XTrack) with dE cut, "+srange.str()+" (")+cname+")",
					     100*fact2Dim, -500., 500., 200*fact2Dim, xMin, xMax);
      deltaYEcutcaloYPdep[i] = new CsHist2D( string("deltaYEcutcaloY")+s.str(),
					     string("YCalo vs (YCalorimeter - YTrack) with dE cut, "+srange.str()+" (")+cname+")",
					     100*fact2Dim, -500., 500., 200*fact2Dim, yMin, yMax);
    }

    if ( levelCellTimeHistos > 0 ) {
      CsHistograms::SetCurrentPath(string("/CsBPMonitor/Calorimeter_")+cname+"/CellsTime");
      
      for( unsigned int cell = 0; cell < c->NCells(); cell++ ) {
	stringstream sstream; sstream << cell;
	CsHist1D *h1 = new CsHist1D(string("TimeCell_")+cname+"_"+sstream.str(),
				    string("Cluster Time in Cell (")+sstream.str()+") ("+cname+")",
				    200, -100., 100.);
	caloCellsTime.push_back(h1);
      }
      CsHistograms::SetCurrentPath(dir_name);
    }
  
  }
  
  effTrackAss   = new TProfile(   "effTrackAss",   (string("Probability of track -> cluster assoc (")+cname+")").c_str(),
				  100, 0., 200., 0., 1., " ");
  effCaloAss    = new TProfile(   "effCaloAss",    (string("Probability of cluster -> track assoc (")+cname+")").c_str(),
				  100, 0., 200., 0., 1., " ");
} 


/// Fill calo-cluster histograms (loop over clusters).
void CsBuildParticles::BPCalo::FillClustersHists() {
  if ( levelBPHistos <= 0 )
    return;

  caloN->Fill( (double)(clusters.size()) );
  const vector<Reco::CellDataRaw>& signals_ = cscalo->GetSignals();
  caloNData->Fill( (double)(signals_.size()) );
  for( unsigned int jcell = 0; jcell < signals_.size(); jcell++ )
    caloCellsE->Fill( signals_[jcell].GetEnergy() );
    
  for( vector<BPCluster>::const_iterator it = clusters.begin();
       it != clusters.end(); it++) {
    
    const Reco::CalorimeterParticle* clus = it->clus;
    const unsigned imcell = clus->GetMainCells()[0];
    
    clusSize->Fill( clus->GetClusterSize() );
    caloE   ->Fill( clus->GetE()           );
    caloEerr->Fill( clus->GetEerr()        );
    if ( clus->HasTime() ) {
      double err = clus->GetTimeErr();
      if ( err < 1e-10 ) err = 1e100;
      caloTime ->Fill( clus->GetTime() );
      caloTimeS->Fill( clus->GetTime() / err );
      if ( levelCellTimeHistos > 0 && caloCellsTime.size() == cscalo->NCells() )
	caloCellsTime[imcell]->Fill(clus->GetTime());
    }  
    caloXY->Fill(clus->GetX(), clus->GetY());
    if ( cscalo->XYRegularGrid() ) {
      int ix = cscalo->GetColumnOfCell(imcell); 
      int iy = cscalo->GetRowOfCell(imcell); 
      cellsXY->Fill( (double)ix, (double)iy);
    }
    
    if ( !it->associated )
      neutralXY->Fill( clus->GetX(), clus->GetY() );

    effCaloAss  ->Fill( clus->GetE(), it->associated );
    effCaloAss2D->Fill( clus->GetX(), clus->GetY(), it->associated );
    
    if ( clus->HasTime() ) {
      caloTcaloE->Fill(clus->GetTime(), clus->GetE());
      double TCS_T0 = 40.; // maximum of TCS phase distr.
      double tcs_cor = CsEvent::Instance()->getTCSPhaseTime()-TCS_T0;
      caloTtcsT->Fill(clus->GetTime(), tcs_cor );
    }
    pair< bool,double> chisq = clus->GetMiscInfo(Reco::CalorimeterParticle::CHISQ_SADC_SIGNAL );
    if ( chisq.first ) {
      pair< bool,double> ndf = clus->GetMiscInfo(Reco::CalorimeterParticle::NDF_SADC_SIGNAL );
      pair< bool,double> chisq_noise = clus->GetMiscInfo(Reco::CalorimeterParticle::CHISQ_SADC_NOISE );
      pair< bool,double> ndf_noise = clus->GetMiscInfo(Reco::CalorimeterParticle::NDF_SADC_SIGNAL );
      if( ndf.second > 0 && ndf_noise.second > 0 ) {
	if( (chisq_noise.second/ndf_noise.second > chisq.second/ndf.second) ) 
	  dataCHIcaloE->Fill(chisq.second/ndf.second, clus->GetE());
	else   
	  dataCHIcaloE->Fill( -chisq_noise.second/ndf_noise.second-100., clus->GetE());
      }
    }
  }
}


//! Fill cluster and track related histograms
void CsBuildParticles::BPCalo::FillClusterTrackHists( const Reco::CalorimeterParticle* clus,
						      double exclus, double eyclus,
						      double xtrk, double ytrk,
						      double momentum ) {
  if ( levelBPHistos <= 0 )
    return;
  
  caloZ->Fill( clus->GetZ() );

  int i;
  for (i = 0; momentum > Pdep[i] && i < NPdep-1; i++) ;

  if ( fabs( clus->GetY() - ytrk ) < Nsigma*eyclus ) {
    deltaX       ->Fill( clus->GetX() - xtrk );
    deltaXPdep[i]->Fill( clus->GetX() - xtrk );
    if ( fabs( clus->GetE() - momentum ) < NsigmaEcut*clus->GetEerr() ) {
      deltaXEcut       ->Fill( clus->GetX() - xtrk );
      deltaXEcutPdep[i]->Fill( clus->GetX() - xtrk );
    }  
    deltaXS           ->Fill( (clus->GetX()-xtrk)/exclus      );
    caloXtrackX       ->Fill( clus->GetX(),      xtrk         );
    deltaXtrackX      ->Fill( clus->GetX()-xtrk, xtrk         );
    deltaXcaloX       ->Fill( clus->GetX()-xtrk, clus->GetX() );
    deltaXcaloXPdep[i]->Fill( clus->GetX()-xtrk, clus->GetX() );
    if ( fabs( clus->GetE() - momentum ) < NsigmaEcut*clus->GetEerr() ) {
      deltaXEcuttrackX      ->Fill( clus->GetX()-xtrk, xtrk         );
      deltaXEcutcaloX       ->Fill( clus->GetX()-xtrk, clus->GetX() );
      deltaXEcutcaloXPdep[i]->Fill( clus->GetX()-xtrk, clus->GetX() );
    }  
  }  

  if ( fabs( clus->GetX() - xtrk ) < Nsigma*exclus ) {
    deltaY       ->Fill( clus->GetY() - ytrk );
    deltaYPdep[i]->Fill( clus->GetY() - ytrk );
    if ( fabs( clus->GetE() - momentum ) < NsigmaEcut*clus->GetEerr() ) {
      deltaYEcut       ->Fill( clus->GetY() - ytrk );
      deltaYEcutPdep[i]->Fill( clus->GetY() - ytrk );
    }   
    deltaYS           ->Fill( (clus->GetY()-ytrk)/eyclus      );
    caloYtrackY       ->Fill( clus->GetY(),      ytrk         );
    deltaYtrackY      ->Fill( clus->GetY()-ytrk, ytrk         );
    deltaYcaloY       ->Fill( clus->GetY()-ytrk, clus->GetY() );
    deltaYcaloYPdep[i]->Fill( clus->GetY()-ytrk, clus->GetY() );
    if( fabs( clus->GetE() - momentum ) < NsigmaEcut*clus->GetEerr() ) {
      deltaYEcuttrackY      ->Fill( clus->GetY()-ytrk, ytrk         );
      deltaYEcutcaloY       ->Fill( clus->GetY()-ytrk, clus->GetY() );
      deltaYEcutcaloYPdep[i]->Fill( clus->GetY()-ytrk, clus->GetY() );
    }  
  }    
}


/// Fill histograms correlating track and associated cluster, once the
/// association has been made.
void CsBuildParticles::BPCalo::FillClusterTrackAssocHists( const Reco::CalorimeterParticle* clus,
							   const CsTrack *trk, double xtrk, double ytrk,
							   double momentum ) {
  if ( levelBPHistos <= 0 )
    return;
	   
  if ( momentum>0. ) {
    int i;
    for (i = 0; momentum > Pdep[i] && i < NPdep-1; i++) ;
    deltaEs       ->Fill(clus->GetE() / momentum - 1.);
    deltaEsPdep[i]->Fill(clus->GetE() / momentum - 1.);
  }

  caloEtrackE ->Fill( clus->GetE(), momentum );
  deltaEtrackE->Fill( clus->GetE()-momentum, momentum );
  deltaEtrackX->Fill( clus->GetE()-momentum, xtrk );
  deltaEtrackY->Fill( clus->GetE()-momentum, ytrk );

  if ( clus->HasTime() && trk->hasMeanTime() ) {
    caloTtrackT->Fill(clus->GetTime(), trk->getMeanTime() );
    double deltaTime = clus->GetTime() - trk->getMeanTime();
    double dT        = sqrt( pow(clus->GetTimeErr(), 2) + pow(trk->getMeanTimeError(), 2) );
    if ( dT < 1e-10 ) dT = 1e100;
    deltaT->Fill( deltaTime );
    deltaTS->Fill( deltaTime/dT );
    if ( fabs( clus->GetE() - momentum ) < NsigmaEcut*clus->GetEerr() ) {
      deltaTEcut->Fill( deltaTime );
      deltaTSEcut->Fill( deltaTime/dT );
    }  
  } else if ( clus->HasTime() ) {
    deltaT->Fill( -20. );
  } else {
    deltaT->Fill( +20. );
  }
}


/// Fill track-related histograms for all tracks that have passed the general
/// sanity checks.
/// \param helix0 track helix extrapolated to calo z position
void CsBuildParticles::BPCalo::FillTrackHists( const CsHelix &helix0, double momentum,
					       bool track_associated ) {
  if ( levelBPHistos <= 0 )
    return;

  pair<bool,double> in_active_area =
    cscalo->InActiveAreaExternal( helix0.getX(), helix0.getY(), cscalo->GetPositionZ(),
				  helix0.getDXDZ(), helix0.getDYDZ() );
  if( in_active_area.first )
    effTrackAss->Fill(momentum, track_associated);
  effTrackAss2D->Fill(helix0.getX(), helix0.getY(), track_associated);

  {
    int i;
    for (i = 0; momentum > Pdep[i] && i < NPdep-1; i++) ;
    effTrackAssPdep[i]->Fill(helix0.getX(), helix0.getY(), track_associated);
  }

  trackXY->Fill(helix0.getX(), helix0.getY());
  if ( track_associated ) {
    trackXYas->Fill(helix0.getX(), helix0.getY());
    if ( momentum > 0.01 )
      trackXYas_E->Fill(helix0.getX(), helix0.getY());
    else
      trackXYas_noE->Fill(helix0.getX(), helix0.getY());
  } else {
    trackXYnotas->Fill(helix0.getX(), helix0.getY());
    if ( momentum > 0.01 )
      trackXYnotas_E->Fill(helix0.getX(), helix0.getY());
    else 
      trackXYnotas_noE->Fill(helix0.getX(), helix0.getY());
  }
}


CsBuildParticles::CsBuildParticles() {
  vector<CsCalorimeter*> &calorimeters = CsGeom::Instance()->getCalorimeters();

  // stop right here if no calorimeters exist
  if ( ! calorimeters.size() )
    return;
  
  // decide if histograms have to be done or not
  int levelBPHistos = 0;
  CsOpt::Instance()->getOpt( "CsBuildPart", "Hist", levelBPHistos );
  if(levelBPHistos > 0)
    CsErrLog::Instance()->mes( elInfo, "CsBuildParticles Histo ON");
  else
    CsErrLog::Instance()->mes( elInfo, "CsBuildParticles Histo OFF");
  
  int levelCellTimeHistos = 0;
  CsOpt::Instance()->getOpt( "CsBuildPart", "HistCellsTime", levelCellTimeHistos );
  
  CsHistograms::SetCurrentPath("/CsBPMonitor");
  
  for (vector<CsCalorimeter*>::const_iterator c = calorimeters.begin();
       c != calorimeters.end(); c++) {
    _calos.push_back( BPCalo(*c, levelBPHistos, levelCellTimeHistos) );
  }
}


//! Check wether there lies a magnet between zstart and zfinish.
static bool crossMagField( double zstart, double zfinish ) {
  double zs = zstart;
  double zf = zfinish; 
  if( zfinish < zstart ) {
    zs = zfinish;
    zf = zstart; 
  } 
  CsField *field = CsGeom::Instance()->getCsField();
  CsMagInfo* magp = field->getMagInfo();
  for(int i=0; i<field->getNumOfMags(); i++) {
    const double z   = magp[i].zcm;  // z position [mm]
    const double fsc = magp[i].fsc;  // field scale factor
    if( zstart<z && z<zfinish && fsc!=0 ) return true;
  }
  return false;
}

template<class T>
inline T sqr(const T a) {return a*a;}

/*!
  \param recoEvent source of tracks and clusters
  \param parts     output of particles containting of tracks and/or clusters
*/
void CsBuildParticles::Build( const CsRecoEvent& recoEvent, vector<CsParticle*>& parts ) {

  const list<CsBeam*>  beams  = recoEvent.getBeamTracksList();
  const list<CsTrack*> tracks = recoEvent.getTracks();

  /// \todo Could be put in options to enable for special runs like electron calibration
  bool associate_with_beam = false;
  unsigned int parts_min   = 0;
  if ( !associate_with_beam )
    parts_min = beams.size();
    
  // get CalorimeterParticles from recoEvent and sort them into our BPCalos
  {
    const vector<Reco::CalorimeterParticle*> calobjs = recoEvent.getCalObjs();
    
    parts.reserve( beams.size() + tracks.size() + calobjs.size() );
    
    for ( vector<BPCalo>::iterator calo = _calos.begin(); calo != _calos.end(); calo++ ) {
      calo->clusters.clear();

      for ( vector<Reco::CalorimeterParticle*>::const_iterator cobj = calobjs.begin();
	   cobj != calobjs.end(); cobj++ ) {
	if ( (*cobj)->GetCalorimeterName() == calo->cscalo->GetName() )
	  calo->clusters.push_back( BPCalo::BPCluster(*cobj) );
      }
    }
  }

  for( list<CsBeam*>::const_iterator i=beams.begin(); i!=beams.end(); i++ ) {
    CsParticle* part = new CsParticle( *i );
    parts.push_back( part );
  }
  
  for( list<CsTrack*>::const_iterator i=tracks.begin(); i!=tracks.end(); i++ ) {
    CsParticle* part = new CsParticle( *i );
    parts.push_back( part );
  }

  // loop calorimeters
  for ( vector<BPCalo>::iterator calo = _calos.begin(); calo != _calos.end(); calo++ ) {
    const double zcal = calo->cscalo->GetPositionZ();

    // loop particles (tracks)
    for ( unsigned int j=parts_min; j<parts.size(); j++ ) {
      const CsTrack* trk = parts[j]->getTrack();
    
      if ( !trk  )
	continue;

      // get nearest helix upstream of zcal
      const CsHelix *hlx = trk->getHelixUpstreamOf( zcal );

      // track without helices?
      if( hlx == NULL )
	continue;
	
      // skip track if it starts after this detector
      if( hlx->getZ() > zcal )
	continue;

      // take momentum from nearest helix
      double momentum = 0.;
      const double abs_cop = fabs( hlx->getCop() );
      if ( abs_cop > 0.0001 ) momentum = 1. / abs_cop;
      
      // skip track without momentum crossing magnetic field
      if( abs_cop <= 0.0001 && crossMagField( trk->getHelices()[0].getZ(), zcal ) )
	continue;

      // skip track if extrapolation fails
      CsHelix hlx0;
      if( ! hlx->Extrapolate( zcal, hlx0, false ) ) {
	// print error message only once per event
	if ( calo == _calos.begin() )
	  cerr << "Problems with extrapolation???? " << j << endl;
	continue;
      }
      
      pair<bool,double> in_active_area =
        calo->cscalo->InActiveAreaExternal( hlx0.getX(), hlx0.getY(), zcal,
					    hlx0.getDXDZ(), hlx0.getDYDZ() );
      
      // skip tracks not hitting the calorimeter
      if ( !in_active_area.first )
	continue;
      
      // Now let's try to associate: loop clusters
      bool track_associated=false;
      for( vector<BPCalo::BPCluster>::iterator it = calo->clusters.begin();
	   it != calo->clusters.end(); it++) {

	Reco::CalorimeterParticle* clus = it->clus;

	const int mcell = clus->GetMainCells()[0]; 

        /// \todo why not use spatial error ?
	const double exclus = calo->cscalo->GetCells()[mcell].GetCellType().GetSizeX()/4.;
	const double eyclus = calo->cscalo->GetCells()[mcell].GetCellType().GetSizeY()/4.;
	//                 double exclus = clus->GetXerr(); 
	//                 double eyclus = clus->GetYerr();
	CsHelix hlx1;

	// set z position of cluster, assuming that the cluster was caused by
	// a charged, interacting particle (this is the default assumption)
	clus->SetZ( clus->CalcZforTrack( hlx0.getDXDZ(), hlx0.getDYDZ(), false ) );
	if ( ! hlx0.Extrapolate(clus->GetZ(), hlx1, false) ) {
	  cerr << __func__ << ": Cannot extrapolate!" << endl;
	  continue;
	}

	calo->FillClusterTrackHists( clus, exclus, eyclus,
				     hlx1.getX(), hlx1.getY(), momentum );

        // match against error ellipse:  x**2/a**2 + y**2/b**2 == 1
        {
          const double dx = clus->GetX() - hlx1.getX();
          const double dy = clus->GetY() - hlx1.getY();
          if ( sqr( dx / (calo->Nsigma*exclus) ) + sqr( dy / (calo->Nsigma*eyclus) ) <= 1. ) {
            parts[j]->addCalobj( clus );
            it->associated   = true;
            track_associated = true;
            calo->FillClusterTrackAssocHists( clus, trk, hlx1.getX(), hlx1.getY(), momentum );
	    calo->assocAttempt->Fill( 0. );
            continue;
          }
        }
		    
	// use z position of MIP and try to match the track again
	const double int_x = hlx1.getX();  // interacting particle track x
	const double int_y = hlx1.getY();  // interacting particle track y
	const double mip_z = clus->CalcZforTrack( hlx0.getDXDZ(), hlx0.getDYDZ(), true );
	if ( ! hlx0.Extrapolate( mip_z, hlx1, false ) ) {
	  cerr << __func__ <<": Cannot extrapolate!" << endl;
	  continue;
	}

	calo->deltaXYmip->Fill( hlx1.getX() - int_x, hlx1.getY() - int_y );

        // match against error ellipse:  x**2/a**2 + y**2/b**2 == 1
        {
          const double dx = clus->GetX() - hlx1.getX();
          const double dy = clus->GetY() - hlx1.getY();
          if ( sqr( dx / (calo->Nsigma*exclus) ) + sqr( dy / (calo->Nsigma*eyclus) ) <= 1. ) {
            parts[j]->addCalobj( clus );
            it->associated   = true;
            track_associated = true;
            calo->FillClusterTrackAssocHists( clus, trk, hlx1.getX(), hlx1.getY(), momentum );
	    calo->assocAttempt->Fill( 1. );
            continue;
          }
        }
	
	calo->assocAttempt->Fill( 2. );
      }  // end loop clusters

      calo->FillTrackHists( hlx0, momentum, track_associated );
    }  // end loop particles (tracks)
    

    calo->FillClustersHists();
    
    for( vector<BPCalo::BPCluster>::iterator it = calo->clusters.begin();
	 it != calo->clusters.end(); it++) {
            
      if ( it->associated )
        continue;
      
      CsParticle* part = new CsParticle( it->clus );
      parts.push_back( part );
    }
  }   // end loop calorimeters
}


void CsBuildParticles::SetClusterZ( vector<CsParticle*>& parts ) const {
  for ( vector<CsParticle*>::iterator p = parts.begin(); p != parts.end(); p++ ) {
    const CsTrack                     *trk = (*p)->getTrack();
    vector<Reco::CalorimeterParticle*> cps = (*p)->getCalObjects();
    
    // loop calorimeter clusters
    for ( vector<Reco::CalorimeterParticle*>::iterator cp = cps.begin(); cp != cps.end(); cp++ ) {

      // calculate shower midpoint
      double mid;
      const Reco::Calorimeter *calo = (*cp)->GetCalorimeter();
      const unsigned  idx_main_cell = (*cp)->GetMainCells().at(0);
      if ( calo->GetType() == Reco::Calorimeter::ElectroMagnetic ) {
	const double crit_energy = calo->GetCells().at(idx_main_cell).GetCellType().GetCriticalEnergy();
	if ( trk ) {
	  mid = Reco::ZmidShowerElectron( (*cp)->GetE(), crit_energy );
	} else {
	  mid = Reco::ZmidShowerPhoton( (*cp)->GetE(), crit_energy );
	}
	mid *= calo->GetCells().at(idx_main_cell).GetCellType().GetRadiationLength();
      } else {
	mid = Reco::ZmidShowerHadronic( (*cp)->GetE() );
	mid *= calo->GetCells().at(idx_main_cell).GetCellType().GetNuclearLength();
      }

      // apply track angle to z calculation
      double dxdz, dydz;
      const double zcal = calo->GetCells().at(idx_main_cell).GetFront();
      if ( trk ) {

	// find nearest helix upstream of detector
	const CsHelix *nearest_hlx = trk->getHelixUpstreamOf( zcal );
	dxdz = nearest_hlx->getDXDZ();
	dydz = nearest_hlx->getDYDZ();

      } else {
	
	/// \todo: put correct target z position
	dxdz = (*cp)->GetX() / (zcal - 0.);
	dydz = (*cp)->GetY() / (zcal - 0.);
      }

      const double z = zcal + mid / sqrt( dxdz*dxdz + dydz*dydz );
      (*cp)->SetZ(z);
    }
  }
}
