/*!
   \file    CsRCnTup.cc
   \brief   Output RICH1 info to an HBOOK ntuple.
   \author  Dmitri Pechekhonov, Benigno Gobbo, Alexandre Korzevev
   \date    $Date: 2010/06/18 10:44:21 $
   \version $Revision: 1.28 $
   \giulia adds root ntuple
*/

//--- recontructed by Paolo 01/10/19.
//    revision of 05/01/10
// giulia 07/07/02
 

#include "CsErrLog.h"
#include "CsMCUtils.h"
#include "CsEvent.h"
#include "CsOpt.h"
#include "CsHbookProto.h"
#include "CsHistograms.h"
#include "CsGeom.h"

#include "CsRCnTup.h"

#include "CsRichOne.h"
#include "CsRCRecConst.h"
#include "CsRCDetectors.h"
#include "CsRCEventParticles.h"
#include "CsRCParticle.h"
#include "CsRCPad.h"
#include "CsRCCluster.h"
#include "CsRCEventPartPhotons.h"
#include "CsRCPartPhotons.h"
#include "CsRCEventRings.h"
#include "CsRCRing.h"
#include "CsRCPhoton.h"
#include "CsRCEventAnalysis.h"

using namespace std;
using namespace CLHEP;

CsRCnTup* CsRCnTup::instance_ = 0;

//===========================================================================
CsRCnTup* CsRCnTup::Instance() {
//------------------------------
  if( instance_ == 0 ) instance_ = new CsRCnTup();
  return( instance_ );
}


//===========================================================================
CsRCnTup::CsRCnTup() {
//--------------------

  cout << endl;
  cout << "RICHONE: CsRCnTup::CsRCnTup() : Ntuple required" << endl;
  cout << "-----------------------------------------------" << endl;


  int i,tmp;
  string msg;
  list<string> block;
  list<string>::iterator Is;
  
  for(i=0;i<4;i++) block_[i]=false;
  
  block.clear();
  CsOpt::Instance()->getOpt("RICHONE","BLOCK",block);
  
  i=0;
  for( Is=block.begin(); Is!=block.end(); Is++, i++ ) {
    block_[i] = false;
    if     ( i == 0  &&  (*Is) == "Ring" ) block_[i] = true;
    else if( i == 1  &&  (*Is) == "Photon" ) block_[i] = true;
    else if( i == 2  &&  (*Is) == "Like" ) block_[i] = true;
    else if( i == 3  &&  (*Is) == "Phys" ) block_[i] = true;
    //cout << block_[i] << endl;
  }

  msg = "Block information: ";
  for( Is=block.begin(); Is!=block.end(); Is++, i++ ) {
    msg += (*Is);
    msg += " ";
  }
  CsErrLog::Instance()->mes( elInfo, msg);

  CsRCnTupNm();
//------------

}


//===========================================================================
void CsRCnTup::CsRCnTupNm() {
//---------------------------

//Paolo  -  January 2005


  CsRCRecConst *cons = CsRCRecConst::Instance();

  static bool firstCall = true;
  if( firstCall ) {
    firstCall = false;

    ntupNum = 110;
  }
  ntupNum++;
  std::cout << "Ntuple  " << ntupNum << "  open" << std::endl;

  switch( CsHistograms::GetImplementation() )
  {
    #if USE_HBOOK
    case CsHistograms::HBOOK:
    {
      if( block_[0] || block_[1] || block_[2] || block_[3] ) {
	CsHistograms::SetCurrentPath( "/RICH/NTUP" );

	pathname = CsHistograms::GetCurrentPath().c_str();

	//ntupNum = 111;
	hbnt  ( ntupNum, pathname, " ");
	hbname( ntupNum, "Header  ", (void*)&(Header.Run), 
		"Run:R*4, Ev:R*4, Index:R*4" );

	std::string str_0 = "nTracks[0,30]:I, ";
	str_0.append( "xIn(nTracks):R*4, yIn(nTracks):R*4, " );
	str_0.append( "zIn(nTracks):R*4, " );
	str_0.append( "tgX(nTracks):R*4, tgY(nTracks):R*4, " );
	str_0.append( "cmom(nTracks):R*4, " );
	str_0.append( "zHx0(nTracks):R*4, zHx1(nTracks):R*4" );
	if( block_[0] ) hbname( ntupNum, "Track", (void*)&(Track.nTracks),
				str_0.c_str() );

	std::string str_00 = "nRings[0,30]:I, ";
	str_00.append( "theta(nRings):R*4, cPaDe(nRings):I, " );
	str_00.append( "xPaDe(nRings):R*4, yPaDe(nRings):R*4, " );
	str_00.append( "nPhoRing(nRings):I, nPhoBkx(nRings):I" );
	if( block_[0] ) hbname( ntupNum, "Ring", (void*)&(Ring.nRings),
				str_00.c_str() );

	std::string str_1 = "nPhotons[0,300]:I, ";
	str_1.append( "the(nPhotons):R*4, phi(nPhotons):R*4, ");
	str_1.append( "PH(nPhotons):R*4, cath(nPhotons):I, " );
	str_1.append( "xc(nPhotons):R*4, yc(nPhotons):R*4, " );
	str_1.append( "nPads(nPhotons):I, paPH(nPhotons):R*4");
	if( block_[1] ) hbname( ntupNum, "Photon", (void*)&(Photon.nPhotons),
				str_1.c_str() );

	std::string str_2 = "nLike[0,30]:I, ";
	str_2.append( "likeBkg(nLike):R*4, likePi(nLike):R*4, " );
	str_2.append( "likeKa(nLike):R*4, likePr(nLike):R*4, " );
	str_2.append( "dLiDnPi(nLike):R*4, dLiDnKa(nLike):R*4, " );
	str_2.append( "dLiDnPr(nLike):R*4, theLike(nLike):R*4, " );
	str_2.append( "theRec(nLike):R*4, nPhoLk(nLike):R*4, " );
	str_2.append( "theFit(nLike):R*4, chiRing(nLike):R*4, " );
	str_2.append( "chiPi(nLike):R*4, chiKa(nLike):R*4, " );
	//str_2.append( "chiPr(nLike):R*4");
	if( cons->outBufferSize() > 15 ) {
	  str_2.append( "chiPr(nLike):R*4, ");
	  str_2.append( "likeEl(nLike):R*4, likeMu(nLike):R*4, " );
	  str_2.append( "dLiDnEl(nLike):R*4, dLiDnMu(nLike):R*4, " );
	  str_2.append( "chiEl(nLike):R*4, chiMu(nLike):R*4" );
	} else {
	   str_2.append( "chiPr(nLike):R*4");
	}
	if( block_[2] ) hbname( ntupNum, "Like", (void*)&(Like.nLike),
				str_2.c_str() );

        if( block_[3] ) hbname( ntupNum, "Phys", (void*)&(Phys.kEvGood),
          "kEvGood:I, kMuPart:I, zVertex:R*4, phiMass:R*4" );
	
	CsHistograms::SetCurrentPath( "/" ); 
      }
      break;
    }
    #endif // USE_HBOOK
    
    case CsHistograms::ROOT:
      {

	if( block_[0] || block_[1] || block_[2] || block_[3] ) {
	  CsHistograms::SetCurrentPath( "/RICH/NTUP" );
	 std::string myTitle="Ntuple of coral results";
	 std::stringstream ss;
	 ss<<ntupNum;
	 ss >> myTitle;

	 MyTree_=new TTree(myTitle.c_str(),"Coral results stored into a tree");
	 MyTree_->Branch("Header",&Header.Run,"Run/F:Ev/F:Index/F");
	 //track
	 if(block_[0]){
	   MyTree_->Branch("nTracks",&Track.nTracks,"nTracks/I");
	   MyTree_->Branch("xIn",Track.xIn,"xIn[nTracks]/F");
	   MyTree_->Branch("yIn",Track.yIn,"yIn[nTracks]/F");
	   MyTree_->Branch("zIn",Track.zIn,"zIn[nTracks]/F");
	   MyTree_->Branch("tgX",Track.tgX,"tgX[nTracks]/F");
	   MyTree_->Branch("tgY",Track.tgY,"tgY[nTracks]/F");
	   MyTree_->Branch("cmom",Track.cmom,"cmom[nTracks]/F");
	   MyTree_->Branch("zHx0",Track.zHx0,"zHx0[nTracks]/F");
	   MyTree_->Branch("zHx1",Track.zHx1,"zHx1[nTracks]/F");
	 }//if(block_[0])
	 
	 //ring
	 if(block_[0]){
	   MyTree_->Branch("nRings",&Ring.nRings,"nRings/I");
	   MyTree_->Branch("theta",Ring.theta,"theta[nRings]/F");
	   MyTree_->Branch("cPaDe",Ring.cPaDe,"cPaDe[nRings]/I");
	   MyTree_->Branch("xPaDe",Ring.xPaDe,"xPaDe[nRings]/F");
	   MyTree_->Branch("yPaDe",Ring.yPaDe,"yPaDe[nRings]/F");
	   MyTree_->Branch("nPhoRing",Ring.nPhoRing,"nPhoRing[nRings]/I");
	   MyTree_->Branch("nPhoBkx",Ring.nPhoBkx,"nPhoBkx[nRings]/I");
	 }//if(block_[0])
	 
	  //photon
	 if(block_[1]){
	   MyTree_->Branch("nPhotons",&Photon.nPhotons,"nPhotons/I");
	   MyTree_->Branch("the",Photon.the,"the[nPhotons]/F");
	   MyTree_->Branch("phi",Photon.phi,"phi[nPhotons]/F");
	   MyTree_->Branch("PH",Photon.PH,"PH[nPhotons]/F");
	   MyTree_->Branch("cath",Photon.cath,"cath[nPhotons]/I");
	   MyTree_->Branch("xc",Photon.xc,"xc[nPhotons]/F");
	   MyTree_->Branch("yc",Photon.yc,"yc[nPhotons]/F");
	   MyTree_->Branch("nPads",Photon.nPads,"nPads[nPhotons]/I");
	   MyTree_->Branch("paPH",Photon.paPH,"paPH[nPhotons]/F");
	 }//if(block_[1])

	 //Like
	 if(block_[2]){
	   MyTree_->Branch("nLike",&Like.nLike,"nLike/I");
	   MyTree_->Branch("likeBkg",Like.likeBkg,"likeBkg[nLike]/F");
	   MyTree_->Branch("likePi",Like.likePi,"likePi[nLike]/F");
	   MyTree_->Branch("likeKa",Like.likeKa,"likeKa[nLike]/F");
	   MyTree_->Branch("likePr",Like.likePr,"likePr[nLike]/F");
	   MyTree_->Branch("dLiDnPi",Like.dLiDnPi,"dLiDnPi[nLike]/F");
	   MyTree_->Branch("dLiDnKa",Like.dLiDnKa,"dLiDnKa[nLike]/F");
	   MyTree_->Branch("dLiDnPr",Like.dLiDnPr,"dLiDnPr[nLike]/F");
	   MyTree_->Branch("theLike",Like.theLike,"theLike[nLike]/F");
	   MyTree_->Branch("theRec",Like.theRec,"theRec[nLike]/F");
	   MyTree_->Branch("nPhoLk",Like.nPhoLk,"nPhoLk[nLike]/F");
	   MyTree_->Branch("theFit",Like.theFit,"theFit[nLike]/F");
	   MyTree_->Branch("chiRing",Like.chiRing,"chiRing[nLike]/F");
	   MyTree_->Branch("chiPi",Like.chiPi,"chiPi[nLike]/F");
	   MyTree_->Branch("chiKa",Like.chiKa,"chiKa[nLike]/F");
	  


	   if( cons->outBufferSize() > 15 ) {
	     MyTree_->Branch("chiPr",Like.chiPr,"chiPr[nLike]/F");
	     MyTree_->Branch("likeEl",Like.likeEl,"likeEl[nLike]/F");
	     MyTree_->Branch("likeMu",Like.likeMu,"likeMu[nLike]/F");
	     MyTree_->Branch("dLiDnEl",Like.dLiDnEl,"dLiDnEl[nLike]/F");
	     MyTree_->Branch("dLiDnMu",Like.dLiDnMu,"dLiDnMu[nLike]/F");
	     MyTree_->Branch("chiEl",Like.chiEl,"chiEl[nLike]/F");
	     MyTree_->Branch("chiMu",Like.chiMu,"chiMu[nLike]/F");
	   } else {
	     MyTree_->Branch("chiPr",Like.chiPr,"chiPr[nLike]/F");
	   }



	 }// if(block_[2])


	 
	 //Phys
	 if( block_[3] ) MyTree_->Branch("Phys", &Phys.kEvGood,"kEvGood/I:kMuPart/I:zVertex/F:phiMass/F" );
	 
	 CsHistograms::SetCurrentPath( "/" ); 

	}//if( block_[0] || block_[1] || block_[2] || block_[3] )
	break;
      }//case ROOT
    
    default:
      CsErrLog::Instance()->mes( elFatal, 
             "Histogram package not set to HBOOK or ROOT" );
  }
}


//===========================================================================
void CsRCnTup::fillOutput() {
//---------------------------


  static int recNum;

  static bool firstCall = true;
  if( firstCall ) {
    firstCall = false;

    recNum = 0;
  }
  recNum++;
  if( recNum%1000000 == 0 ) CsRCnTupNm();
//                          ------------

  switch( CsHistograms::GetImplementation() ) {
#   if USE_HBOOK
    case CsHistograms::HBOOK: {
      fillOutNtuple();
      break;
    }  
#   endif
    case CsHistograms::ROOT: {
 fillOutNtuple();
      fillTree();
      break;
    }
    default: {
      break;
    }
  }
}


//===========================================================================
void CsRCnTup::fillOutNtuple() {
//------------------------------


  CsRCRecConst *cons = CsRCRecConst::Instance();

  if( !block_[0] && !block_[1] && !block_[2] && !block_[3] ) return;

  Header.Run = CsEvent::Instance()->getRunNumber();
  //Header.Ev  = CsEvent::Instance()->getEventNumberInRun();

  Header.Ev  = CsRCRecConst::Instance()->CFRefInd();            //   040512
  //Header.Index  = CsRCRecConst::Instance()->CFRefInd();         //   040519
  Header.Index  = CsRCRecConst::Instance()->CFRefIndVS();         //   070131

//------------------- EVENT PART-PHOTONS or RINGS -----------------
  if( block_[0] ) {

//- loop over rings :
//  -----------------
    //int kRing = 0;
    //int kPhot = 0;
    //list<CsRCRing*> lRings = CsRCEventRings::Instance()->lRings();
    //list<CsRCRing*>::iterator ir;
    //for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
    //  if( !(*ir)->flag() ) continue;
    //  if( kRing >= 30 ) continue;

//- LOOP over PART-PHOTONS or RINGS :
//  ---------------------------------
    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;

    int kPPhot = 0;
    int kRing = 0;
    int kPhot = 0;
    int kPk = 0;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;
      if( kPPhot >= 30 ) continue;
      CsRCRing* ring = papho->pRing();
      bool RING = false;
      if( ring != NULL  &&  ring->flag() ) RING = true;
      //ring->lPhotons().size() >= nPhotMinRing ) RING = true;
//@@  -------------------------------------------------------
      if( cons->likeType() == "RING"  &&  !RING ) {
	kPPhot++;
	continue;
      }

      if( cons->likeType() == "ALL" ) kPk = kPPhot;
      if( cons->likeType() == "RING" ) kPk = kRing;
      //CsRCParticle* part = (*ir)->pPartPhot()->pPart();
      CsRCParticle* part = (*ipf)->pPart();
      float mom = part->mom();
      int charge = part->charge();
      Hep3Vector vPosIn = part->vPosIn();
      Hep3Vector vDirIn = part->vDirIn();
      double zHelix0 = 1.e+06;
      double zHelix1 = 1.e+06;
      double chiSq = 0.;
      int nClus = 0;
      zHelix0 = part->zHelix0();
      zHelix1 = part->zHelix1();
      nClus = part->nClus();
      chiSq = part->chiSq();

      Track.xIn[kPk] = vPosIn.x();
      Track.yIn[kPk] = vPosIn.y();
      Track.zIn[kPk] = vPosIn.z();
      Track.tgX[kPk] = 0.;
      if( vDirIn.z() != 0. ) Track.tgX[kPk] = vDirIn.x()/vDirIn.z();
      Track.tgY[kPk] = 0.;
      if( vDirIn.z() != 0. ) Track.tgY[kPk] = vDirIn.y()/vDirIn.z();
      Track.cmom[kPk] = charge * mom;
      //Track.zHx0[kPk] = zHelix0;
      //Track.zHx1[kPk] = zHelix1;
      Track.zHx0[kPk] = int( chiSq )*100000 + int( zHelix0 );
      Track.zHx1[kPk] = nClus*100000 + int( zHelix1 );
      //if( track ) {
	//cout <<  Track.zHx0[kPk] << "  " << Track.zHx1[kPk] << endl;
	//float chiSqw = int( Track.zHx0[kPk]/100000. );
	//float zHel0w = Track.zHx0[kPk] - chiSqw*100000.;
	//int nClusw = int( Track.zHx1[kPk]/100000. );
	//Float zHel1w = Track.zHx1[kPk] - nClusw*100000;
	//cout << chiSqw << "  " << zHel0w << "  " << nClusw 
	//     << "  " << zHel1w<< endl;
      //}

      //Ring.theta[kPk] = (*ir)->the();
      //int kDetPart = (*ir)->pPartPhot()->kDetPart();
      //Ring.cPaDe[kPk] = (*ir)->pPartPhot()->iCaPa();
      //Ring.xPaDe[kPk] = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].x();
      //Ring.yPaDe[kPk] = (*ir)->pPartPhot()->vPoPaDetW()[kDetPart].y();
      int kDetPart = (*ipf)->kDetPart();
      Ring.cPaDe[kPk] = (*ipf)->iCaPa();
      Ring.xPaDe[kPk] = (*ipf)->vPoPaDetW()[kDetPart].x();
      Ring.yPaDe[kPk] = (*ipf)->vPoPaDetW()[kDetPart].y();
      //cout << kDetPart << "  " 
      //<< (*ipf)->vPoPaDetW()[0].y() << "  "
      //<< (*ipf)->vPoPaDetW()[1].y() << endl;
      //<< (*ipf)->vPoPaDetW()[kDetPart].x() << "  "
      //<< (*ipf)->vPoPaDetW()[kDetPart].y() << "  "
      //<< (*ipf)->pPart()->vPoPaMir0()[kDetPart].x() << "  "
      //<< (*ipf)->pPart()->vPoPaMir0()[kDetPart].y() << "  "
      //<< Ring.xPaDe[kPk] << "  " << Ring.yPaDe[kPk] 
      //<< endl;

      Ring.theta[kPk] = 0.;
      if( cons->likeType() == "ALL"  &&  !RING ) {
	if( (*ipf)->thetaLikeSet() ) Ring.theta[kPk] = (*ipf)->thetaLike();
      }
      else  Ring.theta[kPk] = ring->the();
      //Ring.theta[kPk] = (*ir)->the();
      //list<CsRCPhoton*> lPhotons = (*ir)->lPhotons();
      list<CsRCPhoton*> lPhotons;
      lPhotons.clear();
      list<CsRCPhoton*>::iterator ih;
      if( cons->likeType() == "ALL"  &&  !RING ) lPhotons = papho->lPhotons();
      if( RING ) lPhotons = ring->lPhotons();
      Ring.nPhoRing[kPk] = lPhotons.size();
      CsRCEventAnalysis* ana = CsRCEventAnalysis::Instance();
      Ring.nPhoBkx[kPk] = -1;

//------------------- RING PHOTONS --------------------------
      if( block_[1] ) {

        for( ih=lPhotons.begin(); ih!=lPhotons.end(); ih++ ) {
          if( !(*ih)->flag() ) continue;
	  if( kPhot >= 300 ) break;

	  if( (*ih)->isPMT() ) {
	    Photon.the[kPhot] = (*ih)->the0();
	  } else {
	    Photon.the[kPhot] = - fabs( (*ih)->the0() );
	  }
          Photon.phi[kPhot] = (*ih)->phi();
          Photon.PH[kPhot] = (*ih)->PH();
          Photon.cath[kPhot] = (*ih)->ptToClu()->ic();
	  Photon.xc[kPhot] = (*ih)->ptToClu()->xc();
	  Photon.yc[kPhot] = (*ih)->ptToClu()->yc();

	  list<CsRCPad*> lPads = (*ih)->ptToClu()->lPads();
          Photon.nPads[kPhot] = lPads.size();
	  list<CsRCPad*>::iterator ia;
	  double pPHmx = -1.;
	  for( ia=lPads.begin(); ia!=lPads.end(); ia++ ) {
	    if( (*ia)->PH() > pPHmx ) pPHmx = (*ia)->PH();
	  }
	  Photon.paPH[kPhot] = pPHmx;

	  kPhot++;
	}

      }

//--- 030404
//------------------- RING LIKE'S --------------------------
      if( block_[2] ) {

	CsRCRecConst *cons = CsRCRecConst::Instance();

	static const int nProb = cons->outBufferSize();
	double partProbs0[nProb];
	for( int k=0; k<nProb; k++ ) partProbs0[k] = 0.;
	const double* partProbs;
	partProbs = partProbs0;
	//CsTrack* tk = (*ir)->pPartPhot()->pPart()->pTrack();
	CsTrack* tk = (*ipf)->pPart()->pTrack();
	if( tk ) {
	//--- from raw data
	  if( tk->hasRich1Probs() ) partProbs = tk->getRich1Probs();
	} else {
	//--- from myfile
	  if( cons->likeType() == "RING" ) {
	    if( ring->partProbsSet() ) partProbs = ring->partProbs();
	  }
	  if( cons->likeType() == "ALL" ) {
	    //CsRCPartPhotons* papho = (*ir)->pPartPhot();
	    //if( papho ) {
	    //  if( papho->partProbsSet() ) partProbs = papho->partProbs();
	    //}
	    if( (*ipf)->partProbsSet() ) partProbs = (*ipf)->partProbs();
	  }
	}

	Like.likeBkg[kPk] = partProbs[0];
	Like.likePi[kPk] = partProbs[1];
	Like.likeKa[kPk] = partProbs[2];
	Like.likePr[kPk] = partProbs[3];
        Like.dLiDnPi[kPk] = partProbs[4];
        Like.dLiDnKa[kPk] = partProbs[5];
        Like.dLiDnPr[kPk] = partProbs[6];
	Like.theLike[kPk] = partProbs[7];
	Like.theRec[kPk] = partProbs[8];
	Like.nPhoLk[kPk] = partProbs[9];
	Like.theFit[kPk] = partProbs[10];
	Like.chiRing[kPk] = partProbs[11];
	Like.chiPi[kPk] = partProbs[12];
	Like.chiKa[kPk] = partProbs[13];
	Like.chiPr[kPk] = partProbs[14];
	if( cons->outBufferSize() > 15 ) {
	  Like.likeEl[kPk] = partProbs[15];
	  Like.likeMu[kPk] = partProbs[16];
	  Like.dLiDnEl[kPk] = partProbs[17];
	  Like.dLiDnMu[kPk] = partProbs[18];
	  Like.chiEl[kPk] = partProbs[19];
	  Like.chiMu[kPk] = partProbs[20];
	}
      }
//--- 030404

      kPPhot++;
      if( RING ) kRing++;
    }
    int nPk = 0;
    if( cons->likeType() == "ALL" ) nPk = kPPhot;
    else if( cons->likeType() == "RING" ) nPk = kRing;
    Track.nTracks = nPk;
    Ring.nRings = nPk;
    Like.nLike = nPk;
    Photon.nPhotons = kPhot;
    //Like.nLike = kRing;

  }

// 040519
// WARNING : PHYSICS Ntups writing works ONLY from g-file (type 3) processing!
//------------------- PHYSICS ------------------------------------------------
  if( block_[3] ) {

    CsRichOne* richone = CsRichOne::Instance();
    //std::cout << richone->evGood() << " 1 " << richone->kMuonPart()
    //    << "  " << richone->zVertex()
    //  << "  " << richone->phiMass() << std::endl;
    //int kMuPart = richone->kMuonPart();
    int kPartCode = richone->kMuonPart();
    int kMuPart = 0;
    int kK1Part = 0;
    int kK2Part = 0;
    kK2Part = kPartCode / 10000;
    kK1Part = (kPartCode - kK2Part*10000) / 100;
    kMuPart = kPartCode - kK2Part*10000 - kK1Part*100;
    //std::cout << "------  " << kPartCode << "  ";
    //std::cout << kMuPart << "  " << kK1Part << "  " 
    //	  << kK2Part << std::endl;

    CsRCEventParticles* evparts = CsRCEventParticles::Instance();
    list<CsRCParticle*> lParticles = evparts->lParticles();
    list<CsRCParticle*>::iterator ipk;
    CsRCParticle* pMuPart = NULL;
    CsRCParticle* pK1Part = NULL;
    CsRCParticle* pK2Part = NULL;
    int kPart = 1;
    for( ipk=lParticles.begin(); ipk!=lParticles.end(); ipk++ ) {
      if( kPart == kMuPart ) pMuPart = (*ipk);
      if( kPart == kK1Part ) pK1Part = (*ipk);
      if( kPart == kK2Part ) pK2Part = (*ipk);
      kPart++;
    }
    //std::cout << pMuPart << "  " << pK1Part << "  "
    //  	<< pK2Part << std::endl;

    //list<CsRCRing*> lRings = CsRCEventRings::Instance()->lRings();
    //list<CsRCRing*>::iterator ir;
    //CsRCParticle* pRing = NULL;
    //int kMuRing = 0;
    //int kK1Ring = 0;
    //int kK2Ring = 0;
    //int kRing = 1;
    //for( ir=lRings.begin(); ir!=lRings.end(); ir++ ) {
    //  if( !(*ir)->flag() ) continue;
    //  if( kRing >= 30 ) continue;
    //  pRing = (*ir)->pPartPhot()->pPart();
    //  if( pRing == pMuPart ) kMuRing = kRing;
    //  if( pRing == pK1Part ) kK1Ring = kRing;
    //  if( pRing == pK2Part ) kK2Ring = kRing;
    //  kRing++;
    //}
    int kPPhot = 1;
    int kRing = 1;
    int kPk = 0;
    CsRCEventPartPhotons* evepapho = CsRCEventPartPhotons::Instance();
    list<CsRCPartPhotons*> lPaPhotons = evepapho->lPartPhotons();
    list<CsRCPartPhotons*>::iterator ipf;
    CsRCParticle* pPart = NULL;
    int kMuRing = 0;
    int kK1Ring = 0;
    int kK2Ring = 0;
    for( ipf=lPaPhotons.begin(); ipf!=lPaPhotons.end(); ipf++ ) {
      if( !(*ipf) ) continue;
      CsRCPartPhotons* papho = (*ipf);
      if( !papho->flag() ) continue;
      if( kPPhot >= 30 ) continue;
      CsRCRing* ring = papho->pRing();
      bool RING = false;
      if( ring != NULL  &&  ring->flag() ) RING = true;
      if( cons->likeType() == "RING"  &&  !RING ) {
	kPPhot++;
	continue;
      }
      if( cons->likeType() == "ALL" ) kPk = kPPhot;
      if( cons->likeType() == "RING" ) kPk = kRing;
      pPart = papho->pPart();
      if( pPart == pMuPart ) kMuRing = kPk;
      if( pPart == pK1Part ) kK1Ring = kPk;
      if( pPart == pK2Part ) kK2Ring = kPk;
      kPPhot++;
      if( RING ) kRing++;
    }
    int kRingCode = kMuRing + kK1Ring*100 + kK2Ring*10000;
    //std::cout << kRingCode << "  " << kRing-1 << std::endl;
    //std::cout << richone->evGood() << " 2 " << richone->kMuonPart()
    //    << "  " << richone->zVertex()
    //  << "  " << richone->phiMass() << std::endl;

    Phys.kEvGood = richone->evGood();
    Phys.kMuPart = kRingCode;
    Phys.zVertex = richone->zVertex();
    Phys.phiMass = richone->phiMass();

  }

  //print();                                       // print for check :
//-------
#if USE_HBOOK
  CsHistograms::SetCurrentPath( "/RICH/NTUP" );
  hfnt( ntupNum );                                 // fill  ntuple
  CsHistograms::SetCurrentPath( "/" );      
#endif // USE_HBOOK
}


//===========================================================================
void CsRCnTup::fillTree() {

  CsHistograms::SetCurrentPath( "/RICH/NTUP" );
  MyTree_->Fill();
  CsHistograms::SetCurrentPath( "/" ); 
 }
//---------------------------


//===========================================================================
void CsRCnTup::print() {
//----------------------

  cout << "   " << endl;
  cout << "CsRCnTup::print() : Entuple variables" << endl;
  cout << "-------------------------------------" << endl;

  int kh0 = 0;
  for( int kr=0; kr<Ring.nRings; kr++ ) {

    cout << " ring " << kr << " :" << endl;
    cout << " track pos x " << Track.xIn[kr] << ",  y " << Track.yIn[kr]
	 << ",  z " << Track.zIn[kr] << endl;
    cout << " track dir x " << Track.tgX[kr] << ",  y " << Track.tgY[kr]
	 << endl;
    cout << " mom " << Track.cmom[kr]
	 << ",  track beg " << Track.zHx0[kr]
	 << ",  track end " << Track.zHx1[kr] << endl;

    cout << " thetaC " << Ring.theta[kr] << ",  caC " <<  Ring.cPaDe[kr]
	 << ",  xC " << Ring.xPaDe[kr] << ",  yC " << Ring.yPaDe[kr]
	 << endl;
    cout << " photons " << Ring.nPhoRing[kr]
	 << ",  bk pho " << Ring.nPhoBkx[kr] << endl;
    int kh1 = kh0 + Ring.nPhoRing[kr];
    for( int kh=kh0; kh<kh1; kh++ ) {
      cout << "   the " << Photon.the[kh] << ",  phi " << Photon.phi[kh]
	   << ",  PH " << Photon.PH[kh] << ",  cat " << Photon.cath[kh]
	   << ",  x " << Photon.xc[kh] << ",  y " << Photon.yc[kh]
	   << ",  pads " << Photon.nPads[kh] 
	   << ",  paPH " << Photon.paPH[kh] << endl;
    }
    kh0 = kh1;
  }

  //cout << Phys.n1 << "  -----  " << Phys.f1 << endl;

}
