/*!
   \file    FindSecondaries.cc
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.22 $
   \date    $Date: 2006/07/05 12:29:23 $ 

*/

#include "CsAverPattern.h"
#include "CsGeant3.h"
#include "THlx.h"
#include "CsEvent.h"
#include "CsMCUtils.h"
#include "CsOpt.h"

using namespace std;
using namespace CLHEP;

extern void lindist(THlx *HTr1, THlx *HTr2, double xmin[], double *dist);
//extern void helixdist(THlx *H1, THlx *H2, double a,double b,double x0,
//	       double ver[], double *dist, THlx *H1v, THlx *H2v);
double Mass( THlx&, THlx&, int );
THlx add( THlx &H1, THlx &H2 );




bool CsAverPattern::FindSecondaries( THlx** HTrMom, int* vIact, CsTrack** TrkRef, int tn )
{
  int i,j;

  double v12[3];
  double d12=0;
  THlx HP, HM;
  
  // HISTOGRAM CREATION 
  static bool first = true;
  
  static CsHist1F *hZ;
  static CsHist2F *hXY1,*hXY2;
  
  if( hist_ && first ) {
    
    CsHistograms::SetCurrentPath("/CsAverPattern/Secondaries");
    hZ = new CsHist1F( "hZ", "Z distr. of sec. vertex", 100,-1500,3000 );
    hXY1 = new CsHist2F( "hXY1", "XY distr. of sec. vertex", 100,-100,100, 100,-100,100 );
    hXY2 = new CsHist2F( "hXY2", "XY distr. of sec. vertex", 100,-100,100, 100,-100,100 );
    
    CsHistograms::SetCurrentPath("/");

    first = false;
  }
  
  // if you look for primary particle the 0 track is beam.
  // It must be ignored.
  if( findPrim_ ) i=1;
  else            i=0;
  for( ; i < tn; i++ ){
    
    if( vIact[i] == 5 ) continue;          // skip "special" particles
    if( (*HTrMom[i])(5) < 0 ) continue;        //--- POSITIVE IS ACCEPTED
    if( TrkRef[i]->getXX0() > SecXX0_ ) continue;   // reject muons
    
    // Rejection of mu
    vector<CsHelix> vPl = TrkRef[i]->getHelices();
    if( vPl.size() > 1 && vPl.back().getZ() > 33000 ) continue;
    
    if( Print_[0] ) cout<<"CAP::FindSecondary==> i="<<i<<" mom(+)="<<HTrMom[i]->Mom()<<endl;
    
    if( findPrim_ ) j=1;
    else            j=0;
    for( ; j < tn; j++ ) {

      if( vIact[j] == 5 ) continue;          // skip "special" particles
      if( (*HTrMom[j])(5) > 0 ) continue;        //--- NEGATIVE IS ACCEPTED
      if( TrkRef[j]->getXX0() > SecXX0_ ) continue; // reject muons

      if( Print_[0] ) cout<<"CAP::FindSecondary==>  j="<<j<<" mom(-)="<<HTrMom[j]->Mom()<<endl;
      
      // RELATIVE TIME CHECK
      if( !CsEvent::Instance()->isAMonteCarloEvent() )
	if( TrkRef[i]->hasMeanTime() && TrkRef[j]->hasMeanTime() ) 
	  if( fabs( TrkRef[i]->getMeanTime() - TrkRef[j]->getMeanTime() ) > TimeSecCut_ ) continue;
      
      
      HP = *HTrMom[i];
      HM = *HTrMom[j];
      double ZMAX = (*HTrMom[i])(0) < (*HTrMom[j])(0) ? (*HTrMom[i])(0) : (*HTrMom[j])(0) ;

      //            ********** \param X0 of FindCDA **********
      // - It is the 0-th approximation of CDA X position
      // - Looks like it cannot be outside the search range (3rd and 4th arg's).
      //  This is what I empirically determined when debugging "FindSecondaries"
      //  for the hadron setup, where "THlx::FindCDA", initialised w/ X0=0, as
      //  was originally the case, failed to converge for secondary vertices w/ 
      //  1st helices, and therefore "ZMAX+20", upstream of 0. (Y.B.)
      double X0 = ZMAX-50; //  => Set X0 so that it is always w/in search range.
      if( !HP.FindCDA( HM, X0, -1000, ZMAX + 20 ) ) {
        if( Print_[1] ) cout<<"CAP::FindSecondary==> No intersection"<<endl;
        d12 = 1e10;
        continue;
      }
      
      d12 = HP.Dist( HM );
      v12[0] = ( HP(0) + HM(0) ) / 2;
      v12[1] = ( HP(1) + HM(1) ) / 2;
      v12[2] = ( HP(2) + HM(2) ) / 2;
      
      if( d12 < SecDist_ ) {
	
	// SET VERTEX
	CsVertex *vrt;
	HepMatrix &Cn = *new HepMatrix(3,3);

	vrt = new CsVertex( 10 * v12[1], 10 * v12[2], 10 * v12[0] );

	Cn[0][0] = 4;  Cn[0][1] = 0;  Cn[0][2] = 0;    // sigma=5cm
	Cn[1][0] = 0;  Cn[1][1] = 4;  Cn[1][2] = 0;    // sigma=0.2cm
	Cn[2][0] = 0;  Cn[2][1] = 0;  Cn[2][2] = 2500; // sigma=0.2cm
	
	vrt -> addCov( &Cn );
	vrt -> setType( false );
	
	//TMtx Cov(3,3);
	//Cov *= 0;
	//Cov(1,1) = Cov(2,2)=4; Cov(3,3)=2500;
	//vrt -> setCov( Cov );

	vrt->addTrack( TrkRef[i] );
	vrt->addTrack( TrkRef[j] );
	vrts_.push_back( vrt );
	
	statistics_[7]++; 
	
	if( hist_ ) {
	  
	  double x  = 10 * v12[1];
	  double y  = 10 * v12[2];
	  double z  = 10 * v12[0];
	  hZ->Fill( z );
	  if( z>210 && z<280 ) hXY1->Fill( x, y );
	  if( z>481 && z<533 ) hXY2->Fill( x, y );
	  
        }

	// TEMPORARY
	//delete vrt;

      }
      
    }
  }

  return true;
}


    
THlx add( THlx &H1, THlx &H2 )
{
  THlx H;
    
  if( fabs( H1(0)-H2(0) ) > 1e-4 ) cout<<endl<<"add=> Two THlx have different z coordinates: "<<H1(0)<<",  "<<H2(0)<<"."<<endl;

  H(0) = ( H1(0) + H2(0) ) / 2;
  H(1) = ( H1(1) + H2(1) ) / 2;
  H(2) = ( H1(2) + H2(2) ) / 2;

  double px = H1.Mom() * H1.DirCos(1) + H2.Mom() * H2.DirCos(1);
  double py = H1.Mom() * H1.DirCos(2) + H2.Mom() * H2.DirCos(2);
  double pz = H1.Mom() * H1.DirCos(3) + H2.Mom() * H2.DirCos(3);

  H(3) = py/px;
  H(4) = pz/px;
  H(5) = 1 / sqrt( px*px + py*py + pz*pz );
  //H(5) = 0;

  return H;
}


const double M_Pi = 0.139567;
const double M_P = 0.9382723;

double Mass( THlx& H_PL, THlx& H_MI, int ID )
{
  double M_PL(0),M_MI(0),M;
  double E_PL,E_MI,E;
  double P_PL,P_MI,P;

  if( ID == 0 ) { M_PL=M_P; M_MI=M_Pi; }        //  L 
  else if( ID == 1 ) { M_PL=M_Pi; M_MI=M_P; }   //  AL 
  else if( ID == 2 ) { M_PL=M_Pi; M_MI=M_Pi; }  //  K0 

  P_PL = H_PL.Mom();
  P_MI = H_MI.Mom();

  E_PL = sqrt( M_PL*M_PL + P_PL*P_PL );
  E_MI = sqrt( M_MI*M_MI + P_MI*P_MI );
  E = E_PL + E_MI;

  double PX = P_PL * H_PL.DirCos(1) + P_MI * H_MI.DirCos(1);
  double PY = P_PL * H_PL.DirCos(2) + P_MI * H_MI.DirCos(2);
  double PZ = P_PL * H_PL.DirCos(3) + P_MI * H_MI.DirCos(3);
  P = sqrt( PX*PX + PY*PY + PZ*PZ );

  if( E >= P ) M = sqrt( E*E - P*P );
  else M = -1;

  return M;
}
