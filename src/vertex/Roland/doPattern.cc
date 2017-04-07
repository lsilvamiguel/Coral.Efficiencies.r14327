/*!
   \file    doPattern.cc
   \brief   Vertex Pattern recognition.
   \author  Alexandre Korzenev
   \version $Revision: 1.5 $
   \date    $Date: 2009/01/09 16:03:46 $ 

*/
#include "CsRolandPattern.h"
#include "CsEvent.h"
#include "CsGeant3.h"

using namespace std;

bool CsRolandPattern::doPattern( vector<CsParticle*> &parts,
				 double *reTrackT0 = 0 )
{
  
  // HISTOGRAM CREATION
  static bool first = true;
  static CsHist1F *hX,*hY,*hZ;
  static CsHist1F *hDRX,*hDRY,*hDRZ;
  if( hist_ && first ) {
    const double Zint = 1500;
    CsHistograms::SetCurrentPath("/CsRolandPattern/");
    hX = new CsHist1F("RX","X coordinate of primary vertex",100,-40,40);
    hY = new CsHist1F("RY","Y coordinate of primary vertex",100,-40,40);
    hZ = new CsHist1F("RZ","Z coordinate of primary vertex",100,-350-Zint,-350+Zint);
    
    hDRX = new CsHist1F("DRX","DeltaX (MC - rec)",100,-1.5,1.5);
    hDRY = new CsHist1F("DRY","DeltaY (MC - rec)",100,-1.5,1.5);
    hDRZ = new CsHist1F("DRZ","DeltaZ (MC - rec)",100,-100,100);

    CsHistograms::SetCurrentPath("/");
    first = false;
  }

  vrts_.clear();
  beams_.clear();
  tracks_.clear();
  
  // LOOP FOR TRACKS SELECTION
  int NSPECIALS = 0;
  vector<CsParticle*>::const_iterator ip = parts.begin();
  for( ip; ip != parts.end(); ip++ ) {

    CsTrack* trk  = const_cast<CsTrack*>( (*ip)->getTrack() );
    if( trk == 0 ) continue;
    
    if( (*ip)->getType() == CsParticle::SPECIAL ) {
      NSPECIALS ++;
      vector<CsHelix> tpar=trk->getHelices();
      if( tpar[0].getZ() < 0 && tpar[1].getZ() < 0 ) {
	beams_.push_back( trk );
	continue;
      }
    }
    
    tracks_.push_back( trk );
  }
  
  if( beams_.size() == 0 || tracks_.size() == 0 ) return false;
  if( NSpec_ && NSPECIALS < 2 ) return false;
  
  // CALL FOR ROLAND ALGORITHM
  float vtx2[3];
  int* nflag=new int [8],bad=0,less=0,muon=0,muon2=0,vertflag=0;
  CsVertex *vertPrim;

  vertPrim=getPrimaryVertex(vtx2,nflag,vertflag,bad,less,muon,muon2);
  double z = vertPrim->getZ();

  if( z < -1500 || z > 1000 ) {
    delete vertPrim;
    return false;
  } else 
    vrts_.push_back(vertPrim);

  // HISTOGRAMING
  if( hist_ ) {
    double x = vertPrim->getX();
    double y = vertPrim->getY();
    double z = vertPrim->getZ();
    hX->Fill( x );
    hY->Fill( y );
    hZ->Fill( z );
    
    if( CsEvent::Instance()->isAMonteCarloEvent() ) {
      list<CsMCVertex*>mcvertices = CsGeant3::Instance()->getMCVertices();
      list<CsMCTrack*> mctracks   = CsGeant3::Instance()->getMCTracks();
      double XXMC=mcvertices.front()->getX();
      double XYMC=mcvertices.front()->getY();
      double XZMC=mcvertices.front()->getZ();
      
      hDRX->Fill( XXMC - x );
      hDRY->Fill( XYMC - y );
      hDRZ->Fill( XZMC - z );      
    }
    
  }

  return true;
}
