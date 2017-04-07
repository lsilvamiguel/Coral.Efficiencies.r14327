/*!
   \file    CsVertexFinder.cc
   \brief   Vertex Reconstruction Procedure
   \author  A.Korzenev, C.Ulvegren 
   \version $Revision: 1.6 $
   \date    $Date: 2004/09/21 15:08:01 $
*/

#include "coral_config.h"

#include "CsKalmanFitting.h"
#include "CsVTrack.h"
#include "CsGeant3.h"
#include "CsEvent.h"
#include "CsMCTrack.h"

using namespace std;

extern TMtx inv3x3(TMtx* M, int& ierr);
void CsHelix2TMtx( CsHelix &hel, TMtx &Vk, TMtx &Pk );
void CsMCHit2TMtx( CsMCHit *hit, TMtx &Pnk );

extern "C" float prob_(float*,int*);


bool CsKalmanFitting::PullsMC( list<CsVTrack>& vVertTrk )
{
  int nhits,i;
  double factor;
  TMtx pulls(5);
  TMtx Pk(5),Pnk(5),Vk(5,5);

  if( ! (CsEvent::Instance()->isAMonteCarloEvent()) ) {
    // MC data only...
    CsErrLog::Instance()->mes( elError, "This method works only on MC data." );
    return false;
  }

  CsTrack  *trk;
  CsMCTrack  *trkMC;

  list<CsVTrack>::iterator iTrack;
  for(iTrack=vVertTrk.begin(); iTrack!=vVertTrk.end(); iTrack++){
    CsVTrack &it = (*iTrack);
    
    if(it.primary() && it.getAssociatedTrk()!=0){
      
      trk=it.getAssociatedTrk();
      
      vector<CsHelix> vH = trk->getHelices();
      CsHelix2TMtx(vH[0],Vk,Pk);

      const CsMCTrack* trkMC = trk->getAssociatedMCTrack();
      if (trkMC){
	list<CsMCHit*>hits = trkMC->getMCHits(); list<CsMCHit*>::iterator ih;
	if(Print_[5]) cout<<setprecision(3)<<"CVF::PullsMC==> Reference plane: "
			  <<"z = "<<vH[0].getZ()/10<<" cm"<<endl;
	for(ih=hits.begin();ih!=hits.end();ih++){
	  if( fabs((*ih)->getZ() - vH[0].getZ()) < 1 && (*ih)->getOrigin()==0 ){

	    CsMCHit2TMtx( *ih, Pnk);

	    if(Print_[5]) cout<<"CVF::PullsMC==> Track momentum:"
	                      <<"  px="<<(*ih)->getP().x()
			      <<", py="<<(*ih)->getP().y()
			      <<", pz="<<(*ih)->getP().z()<<endl;
	    break;
	  }
	}
      }else{
	//pulls(1)=-1; pulls(2)=-1; pulls(3)=-1; pulls(4)=-1; pulls(5)=-1;
	//it.setRes(pulls);
	continue;
      }

      for(i=1;i<6;i++){
	pulls(i) = (Pk(i)-Pnk(i)) / sqrt( Vk(i,i) );
	hPullsMC[i-1]->Fill( pulls(i) );
      }
      
      if(Print_[5]){
	cout<<setprecision(2)<<"CVF::PullsMC==> MC: delta=";
        for(i=1;i<6;i++)
	  cout<<setw(9)<<Pk(i)-Pnk(i);
	cout<<endl;
	cout<<setprecision(2)<<"CVF::PullsMC==> MC: pulls=";
        for(i=1;i<6;i++)
	  cout<<setw(9)<<pulls(i);
	cout<<endl<<endl;
      }
      
    }
    
  }

  return true;
}



bool CsKalmanFitting::Pulls( list<CsVTrack>& vVertTrk, TMtx& Cn, TMtx& Xn )
{
  // HISTOGRAM CREATION 
  static bool first = true;
  static CsHist2F *hPullsSm_Z[5],*hPullsSm_P[5];
  if( hist_ && first ) {
    const double Zint = 1500;
    CsHistograms::SetCurrentPath("/CsKalmanFitting/PullsSm");
    hPullsSm[0] = new CsHist1D("hPullX",   "Smoothed Pulls: x"    ,100,-10,10);
    hPullsSm[1] = new CsHist1D("hPullY",   "Smoothed Pulls: y"    ,100,-10,10);
    hPullsSm[2] = new CsHist1D("hPulldXdZ","Smoothed Pulls: dX/dZ",100,-10,10);
    hPullsSm[3] = new CsHist1D("hPulldYdZ","Smoothed Pulls: dX/dZ",100,-10,10);
    hPullsSm[4] = new CsHist1D("hPullCop", "Smoothed Pulls: q/p"  ,100,-10,10);
    
    hPullsSm_Z[0] = new CsHist2F("hPX_Z",   "Smoothed Pulls: x"    , 50,-350-Zint,-350+Zint ,100,-10,10);
    hPullsSm_Z[1] = new CsHist2F("hPY_Z",   "Smoothed Pulls: y"    , 50,-350-Zint,-350+Zint ,100,-10,10);
    hPullsSm_Z[2] = new CsHist2F("hPdXdZ_Z","Smoothed Pulls: dX/dZ", 50,-350-Zint,-350+Zint ,100,-10,10);
    hPullsSm_Z[3] = new CsHist2F("hPdYdZ_Z","Smoothed Pulls: dX/dZ", 50,-350-Zint,-350+Zint ,100,-10,10);
    hPullsSm_Z[4] = new CsHist2F("hPCop_Z", "Smoothed Pulls: q/p"  , 50,-350-Zint,-350+Zint ,100,-10,10);
    
    hPullsSm_P[0] = new CsHist2F("hPX_P",   "Smoothed Pulls: x"    , 50,0,1 ,100,-10,10);
    hPullsSm_P[1] = new CsHist2F("hPY_P",   "Smoothed Pulls: y"    , 50,0,1 ,100,-10,10);
    hPullsSm_P[2] = new CsHist2F("hPdXdZ_P","Smoothed Pulls: dX/dZ", 50,0,1 ,100,-10,10);
    hPullsSm_P[3] = new CsHist2F("hPdYdZ_P","Smoothed Pulls: dX/dZ", 50,0,1 ,100,-10,10);
    hPullsSm_P[4] = new CsHist2F("hPCop_P", "Smoothed Pulls: q/p"  , 50,0,1 ,100,-10,10);
    
    CsHistograms::SetCurrentPath("/");
    first = false;
  }
  
  int i,ierr5=0;
  TMtx Rnk(5), Vk(5,5), aux(5,5), pulls (5);

  list<CsVTrack>::iterator iTrack = vVertTrk.begin();
  //iTrack++; // beam is excluded
  for(; iTrack != vVertTrk.end(); iTrack++){
    CsVTrack &it = (*iTrack);
    
    if( it.primary() ){
      Vk = it.Gk.i5(ierr5);
      aux = Vk - ( it.Ak * (Cn * it.AkT) +
		   it.Ak * (it.Enk * it.BkT) +
		  (it.Ak * (it.Enk * it.BkT)).t() +
		   it.Bk * (it.Dnk * it.BkT) );
      Rnk = it.Pk - it.Pnk;
      for ( i=1; i < 6; i++ )
	pulls(i) = Rnk(i) / sqrt(fabs(aux(i,i)));
	
      if(Print_[5]){
	cout<<setprecision(2)<<"CVF::Pulls==> Smoothed: delta=";
        for( i=1; i<6; i++ )
	  cout<<setw(9)<<Rnk(i);
	cout<<endl;
	cout<<setprecision(2)<<"CVF::Pulls==> Smoothed: pulls=";
        for( i=1; i<6; i++ )
	  cout<<setw(9)<<pulls(i);
	cout<<endl<<endl;
      }

      for( i=1; i<6; i++ ){
	hPullsSm[i-1]  ->Fill( pulls(i) );
	hPullsSm_Z[i-1]->Fill( 10 * Xn(1), pulls(i) );
	int NDF = 2;
	float Chi2 = it.getChi2();
	hPullsSm_P[i-1]->Fill( prob_( &Chi2, &NDF ), pulls(i) );
      }

    }else 
      for ( i=1; i<6; i++ ) 
         pulls(i)=-1;
	 
    it.setRes(pulls);
  }

  
  if( ierr5!=0 ) return false;
  return true;
}


/////////////////////////////////////////////////////


void CsHelix2TMtx( CsHelix &hel, TMtx &Vk, TMtx &Pk ){
  double factor;
  double* bcov = hel.getCov();

  for(int k=1;k<6;k++){
    for(int l=1;l<=k;l++){
      factor=1;
      if( k<3 && l<3 )                   factor=100;
      if( (k>2 && l<3) || (l>2 && k<3) ) factor=10; 
      Vk(l,k)=Vk(k,l)=bcov[((k-1)*k)/2 +(l-1)]/factor;
    }
  }
  
  Pk(1)=hel.getX()/10.;
  Pk(2)=hel.getY()/10.;
  Pk(3)=hel.getDXDZ();
  Pk(4)=hel.getDYDZ();
  Pk(5)=hel.getCop();
      
  return;
}



void CsMCHit2TMtx( CsMCHit *hit, TMtx &Pnk ){
  
  Pnk(1)=hit->getX()/10.0;
  Pnk(2)=hit->getY()/10.0;
  Pnk(3)=hit->getP().x()/hit->getP().z();
  Pnk(4)=hit->getP().y()/hit->getP().z();
  Pnk(5)=hit->getMCTrack()->getParticle()->getCharge()/
         sqrt((hit->getP().x()) * (hit->getP().x())+ 
	      (hit->getP().y()) * (hit->getP().y())+ 
	      (hit->getP().z()) * (hit->getP().z()));
  
  return;
}
