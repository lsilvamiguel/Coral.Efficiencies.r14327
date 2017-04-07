// $Id: CsEventUtils.cc,v 1.31 2010/06/02 10:56:12 tnagel Exp $

/*!
   \file    CsEventUtils.cc
   \brief   Compass Event Utilities Class.
   \author  Benigno Gobbo
   \version $Revision: 1.31 $
   \date    $Date: 2010/06/02 10:56:12 $
*/

#include "CsEventUtils.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsOpt.h"
#include "CsMCDigit.h"
#include "CsMCUtils.h"
#include "CsHistograms.h"
#include "CsDriftChamberDetector.h"
#include "CsStrawTubesDetector.h"
#include "CsDriftTubeDetector.h"
#include "CsDWDetector.h"
#include "CsRichWallDetector.h"

#ifdef COMPASS_USE_OSPACE_STD
#  include <ospace/std/algorithm>
#else
# include <algorithm>
#endif

using namespace std;

//_____________________________________________________________________________
void CsEventUtils::associateClusters() {

  list<CsDetector*> dets =              // Detectors
    CsGeom::Instance()->getDetectors();
  list<CsDetector*>::iterator idet;
  list<CsCluster*> c0List, c1List;      // List of clusters and iterators
  list<CsCluster*>::iterator ic0, ic1;
  double zT = CsGeom::Instance()->getTargetCenter();
  CsDetector *detPrv; for (idet = dets.begin(), detPrv = 0; idet!=dets.end();
			  idet++) {
    CsDetector *det0 = *idet;    //  ********** LOOP ON DETECTORS **********
    if (det0==detPrv) continue; // Associate detector already porcessed => skip

    int LRMode = 0;     // 0: target pointing; 1: cell overlap; other: nothing
    double dCut = -1;   // Association cut

    CsDriftChamberDetector *dcd = dynamic_cast<CsDriftChamberDetector*>(det0);
    CsStrawTubesDetector   *std = dynamic_cast<CsStrawTubesDetector*>(det0);
    CsDriftTubeDetector    *dtd = dynamic_cast<CsDriftTubeDetector*>(det0);
    CsDWDetector           *dwd = dynamic_cast<CsDWDetector*>(det0);
    CsRichWallDetector     *rwd = dynamic_cast<CsRichWallDetector*>(det0);

    CsDetector* det1 = NULL; //   ***** ASSOCIATE DETECTOR? *****
    // ***** GET Left/Right AMBIGUITY RAISING MODE AND ASSOCIATION CUT *****
    if      (dcd) {
      if (!dcd->hasAssociateDet()) continue; det1 = dcd->getAssociateDet();
      LRMode = dcd->getLRMode(); dCut = dcd->getAssociationCut();
    }
    else if (std) {
      if (!std->hasAssociateDet()) continue; det1 = std->getAssociateDet();
      LRMode = std->getLRMode(); dCut = std->getAssociationCut();
    }
    else if (dtd) {
      if (!dtd->hasAssociateDet()) continue; det1 = dtd->getAssociateDet();
      LRMode = dtd->getLRMode(); dCut = dtd->getAssociationCut();
    }
    else if (dwd) {
      if (!dwd->hasAssociateDet()) continue; det1 = dwd->getAssociateDet();
      LRMode = dwd->getLRMode(); dCut = dwd->getAssociationCut();
    }
    else if (rwd) {
      if (!rwd->hasAssociateDet()) continue; det1 = rwd->getAssociateDet();
      LRMode = rwd->getLRMode(); dCut = rwd->getAssociationCut();
    }
    else continue;
    detPrv = det1;
    double z10T = LRMode==1 ? 1 : (det1->getZcm()-zT)/(det0->getZcm()-zT);
    int wireMn = det0->getNWir()/4, wireMx = wireMn*3;
    double resol0 = det0->getVel()*det0->getSpSli();
    double resol1 = det1->getVel()*det1->getSpSli();
    double dResol = // Resolution w/ which the association is checked. Since
      // we are going to consider a cluster and its exact mirror. Let's
      // divide by sqrt(2). And we take 2 sigmas.
      sqrt((resol0*resol0+resol1*resol1)*2);

    //     ***** GIVE -UP IF aCut < 0 OR UNKNOWN MODE *****
    if (dCut<0 || (LRMode!=0 && LRMode!=1) ) continue;

    //      ***** CHECK CLUSTERIZATION ON "det0" AND "det1"
    if (!(det0->clusterized() && det1->clusterized())) continue;
    //       ***** GET CORRESPONDING CLUSTERS LISTS *****
    c0List = det0->getMyClusters(); c1List = det1->getMyClusters();

    //      ********** LOOP ON det0 CLUSTERS **********
    CsCluster *m0Prv; bool ambigPrv;
    for (ic0 = c0List.begin(), m0Prv = 0, ambigPrv = false; ic0!=c0List.end();
	 ic0++) {
      CsCluster *c0 = *ic0, *m0 = c0->getMirrorCluster();
      if (c0==m0Prv)// ***** c0==PREVIOUS MIRROR => already processed => GIVE UP
	// Note: this condition does not suffice: there can be more than one
	// clusters per cell, and cluster and mirror do not necessarily follow
	// each other.
	continue;
      bool ambiguous = ambigPrv; ambigPrv = false;
      if (m0) m0Prv = m0;
      if (c0->hasAssociateClusters())          // ***** c0 ALREADY ASSOCIATED...
	continue;                              // ***** ...=> GIVE UP
      static double whichSide0; if (std)
	whichSide0 = c0->getDigitsList().front()->getData()[1];
      double z10; if (LRMode==2) {
	// Mode==2: infinity pointing in the center, target pointing outside
	int wire = c0->getDigitsList().front()->getAddress();
	z10 = wireMn<wire && wire<wireMx ? 1 : z10T;
      }
      else z10 = z10T;

      //       ***** LOOK FOR BEST CANDIDATE TO CLUSTER ASSOCIATION *****
      double uc0 = c0->getU()*z10, uMn = uc0-dCut, uMx = uc0+dCut;
      static bool mIsBest; static double um0; if (m0) {
	um0 = m0->getU()*z10; double mMn = um0-dCut, mMx = um0+dCut;
	// "um0" need not be larger than "uc0". Therefore check both limits.
	if (mMn<uMn) uMn = mMn; if (mMx>uMx) uMx = mMx;
      }
      ic1 = c1List.begin(); CsCluster *c1 = 0; double dmin = dCut+1;
      while (ic1!=c1List.end()) {
	if ((*ic1)->hasAssociateClusters()) { ic1++; continue; }
	CsCluster *m1 = (*ic1)->getMirrorCluster();
	if (m1 && m1->hasAssociateClusters()) { ic1++; continue; }
	if (std && whichSide0 &&
	    whichSide0*(*ic1)->getDigitsList().front()->getData()[1]<0) {
	  ic1++; continue;
	}
	double uc1 = (*ic1)->getU();
	if (uc1<uMn) { ic1++; continue; } if (uc1>uMx) break;
	double d = fabs(uc1-uc0); bool mbest = false; if (m0) {
	  double dm = fabs(uc1-um0); if (dm<d) { d = dm; mbest = true; } 
	}
	if (d<dmin) { dmin = d; c1 = *ic1; mIsBest = mbest; }
	ic1++;
      } 
      if (c1) {                                         // ***** CANDIDATE FOUND
	if (!m0 || um0>uc0) {
	  //                   ***** CANDIDATE BETTER ASSOCIATED IN NEXT CELL?
	  // If mirror m0, require it to be larger than c0. Otherwise means
	  // we are reprocessing an already processed case.
	  // Association will be all the better w/ c1's mirror, if exists.
	  double uc1 = c1->getU();
	  CsCluster *m1 = c1->getMirrorCluster(); if (m1) uc1 = m1->getU();
	  list<CsCluster*>::iterator jc0 = ic0; jc0++;
	  if (m0 && /* to be on the safe side...*/ jc0!=c0List.end()) jc0++;
	  bool ok = jc0!=c0List.end(); if (ok && std) {
	    double wSide0 = (*jc0)->getDigitsList().front()->getData()[1];
	    double wSide1 = c1->getDigitsList().front()->getData()[1];
	    if (wSide0*wSide1<0) ok = false;
	  }
	  if (ok) {
	    double dNext = fabs((*jc0)->getU()*z10-uc1);
	    if (dNext<dmin) {
	      if (dNext>dmin-dResol) ambigPrv = true;
	      // Better associated is indeed found. Cancel present one. The case
	      // will be re-assessed in next iteration of loop on c0.
	      continue;
	    }
	    else if (dNext<dmin+dResol) ambiguous = true;
	  }
	}
      }
      else continue;

      //                ***** DC: FILL CORRELATION HISTOGRAMS *****
      // (DC is selected here for historical reasons.)
      if (dcd && det0->GetHistLevel()>=CsDetector::Normal) {
	CsDriftChamberDetector *dcd1 =
	  dynamic_cast<CsDriftChamberDetector*>(det1);  
	if (dcd1) {
	  bool err;
	  dcd->getLRHist()->Fill
	    (dcd1->getDistToWire(c1->getDigitsList().front()->getDatum(),err) + 
	     dcd->getDistToWire(c0->getDigitsList().front()->getDatum(),err));
	  dcd->getLRduMinHist()->Fill(dmin);
	}
      }

      CsCluster *m1 = c1->getMirrorCluster(); // ***** NO MIRRORS => EXIT *****
      if (!m0 && !m1) continue;

      //    ***** LOOK FOR BEST Left/Right COMBINATION *****
      int flag = 0;          // Which combination is the best one
      double d[4];           // Left/Right estimator for each combination
      double uc1 = c1->getU();
      double dLow = ambiguous ? 2*dResol : dResol;

      if (!mIsBest) {                  // (c0,c1) best
	d[0] = dmin; d[2] = m0 ? fabs(uc1-um0) : d[0]; flag = 0;
      }
      else {                           // (m0,c1) best
	d[0] = fabs(uc1-uc0); d[2] = dmin; flag = 2;
      }
      if (d[0]<dLow) d[0] = dLow;
      if (d[2]<dLow) d[2] = dLow;

      if (m1) {                         // (c0,m1)
	d[1] = fabs(m1->getU()-uc0);
	if (d[1]<dmin) { dmin = d[1]; flag=1; }
      }
      else d[1] = d[0];
      if (d[1]<dLow) d[1] = dLow;

      if (m0 && m1) {                  //  (m0,m1)
	d[3] = fabs(m1->getU()-um0);
	if (d[3]<dmin) { dmin = d[3]; flag=3; }
      }
      else if (!m0) d[3] = d[1];
      else /* hence if (!m1) */ d[3] = d[2];
      if (d[3]<dLow) d[3] = dLow;

      // Calculates left/right probability. 
      // Associates clusters properly
   
      double a0 = 0.5, a1 = 0.5; // LR probability for c0/c1 to be true

      // Make proper associations depending on flag and set left/right probabilities
      switch (flag) {

      case 0:        // c0 and c1 are best candidates
	// m0 or m1 might be NULL
	if (d[0]+d[2]>0) a0 = d[0]/(d[0]+d[2]); else a0 = 0.5;
	if (d[0]+d[1]>0) a1 = d[0]/(d[0]+d[1]); else a1 = 0.5;
	c0->addAssociateCluster(*c1); if (m1) c0->addAssociateCluster(*m1);
	if (m0) {
	  m0->addAssociateCluster(*c1); if (m1) m0->addAssociateCluster(*m1);
	}
	c1->addAssociateCluster(*c0); if (m0) c1->addAssociateCluster(*m0);
	if (m1) {
	  m1->addAssociateCluster(*c0); if (m0) m1->addAssociateCluster(*m0);
	}
	break;
    
      case 1:        // c0 and m1 are best candidates
	// WARNING m0 might be NULL
	if (d[1]+d[3]>0) a0 = d[1]/(d[1]+d[3]); else a0 = 0.5;
	if (d[0]+d[1]>0) a1 = d[0]/(d[0]+d[1]); else a1 = 0.5;
	c0->addAssociateCluster(*m1); c0->addAssociateCluster(*c1);
	if (m0) {
	  m0->addAssociateCluster(*m1); m0->addAssociateCluster(*c1);
	}
	c1->addAssociateCluster(*c0); if (m0) c1->addAssociateCluster(*m0);
	m1->addAssociateCluster(*c0); if (m0) m1->addAssociateCluster(*m0);
	break;

      case 2:       // m0 and c1 are best candidates
	// WARNING m1 might be NULL
	if (d[0]+d[2]>0) a0 = d[0]/(d[0]+d[2]); else a0 = 0.5;
	if (d[2]+d[3]>0) a1 = d[2]/(d[2]+d[3]); else a1 = 0.5;
	c0->addAssociateCluster(*c1); if (m1) c0->addAssociateCluster(*m1);
	m0->addAssociateCluster(*c1); if (m1) m0->addAssociateCluster(*m1);
	c1->addAssociateCluster(*m0); c1->addAssociateCluster(*c0);
	if (m1) {
	  m1->addAssociateCluster(*m0);  m1->addAssociateCluster(*c0);
	}
	break;

      case 3:       // m0 and m1 are best candidates
	if (d[1]+d[3]>0) a0 = d[1]/(d[1]+d[3]); else a0 = 0.5;
	if (d[2]+d[3]>0) a1 = d[2]/(d[2]+d[3]); else a1 = 0.5;
	c0->addAssociateCluster(*m1); c0->addAssociateCluster(*c1);
	m0->addAssociateCluster(*m1); m0->addAssociateCluster(*c1);
	c1->addAssociateCluster(*m0); c1->addAssociateCluster(*c0);
	m1->addAssociateCluster(*m0); m1->addAssociateCluster(*c0);
	break;

      default: 
	break;
      }
  
      //          ***** FILL LEFT/RIGHT PROBABILITIES *****  
      if (!m0) c0->setLRProb(1);
      else {
	c0->setLRProb(1-a0); m0->setLRProb(a0);
      }
      if (!m1) c1->setLRProb(1);
      else {
	c1->setLRProb(1-a1); m1->setLRProb(a1);
      }
    }  // End loop on clusters
  }  // End loop on detectors
}

//_____________________________________________________________________________
double CsEventUtils::TPQ( CsCluster* c0, CsCluster* c1, const int LRMode ) {
  
  switch (LRMode) {
  case 0:   // ***** TARGET POINTING *****
//    return fabs( 2.0*( c1->getU()*c0->getW() - c0->getU()*c1->getW() )/( c0->getW() + c1->getW() ) ); 
    
    return fabs( 2.0*
      ( c1->getU()*c0->getW() - c0->getU()*c1->getW() )/
      ( c0->getW() + c1->getW() -2.0*CsGeom::Instance()->getTargetCenter() )
    ); // target center position correction
    break;

  case 1:  // ***** INFINITY POINTING *****
    return fabs(c0->getU()-c1->getU());
    break;

  default:
    CsErrLog::Instance()->mes(elFatal,"Unrecognized LRMode"); return 0;
  }
}

//_____________________________________________________________________________
void CsEventUtils::useLRFilter( list<CsCluster*>& clusters ) {

  list<CsCluster*>::iterator iclus; // iterator on cluster list
  list<CsCluster*>::iterator irem; // iterator on cluster list (locked on clusters to be erased)
  CsDetector* cldet = NULL; // detector associated to cluster
  double cut = -1; // cut on LR probability
  double prob = -1;
    
  for ( iclus=clusters.begin(); iclus!=clusters.end(); iclus++ ) {
    if ( !isALRGoodCluster( *(*iclus) ) ) { 
      irem = iclus;
      iclus--;
      clusters.erase( irem );
    }   
  }
} 

//_____________________________________________________________________________
bool CsEventUtils::isALRGoodCluster( const CsCluster clus ) {

  CsDetector *det = NULL; // Detector associated to cluster
  if (!(clus.hasMirrorCluster() &&    // Return true for non drift-like clusters
	clus.hasAssociateClusters())) // or if association failed
    return true;
  
  // Get detector (we take 1st one...) associated to cluster and its LRProbCut
  if (!( det = clus.getDetsList().front())) return true;

  CsDriftChamberDetector *dcd = dynamic_cast<CsDriftChamberDetector*>(det);
  CsStrawTubesDetector   *std = dynamic_cast<CsStrawTubesDetector*>(det);
  CsDriftTubeDetector    *dtd = dynamic_cast<CsDriftTubeDetector*>(det);
  CsDWDetector           *dwd = dynamic_cast<CsDWDetector*>(det);
  CsRichWallDetector     *rwd = dynamic_cast<CsRichWallDetector*>(det);

  double cut = -1;      // Cut on LR probability
  if( ( dcd != NULL && ( cut = dcd->getLRProbCut() ) < 0.0 ) ||
      ( std != NULL && ( cut = std->getLRProbCut() ) < 0.0 ) ||
      ( dtd != NULL && ( cut = dtd->getLRProbCut() ) < 0.0 ) ||
      ( dwd != NULL && ( cut = dwd->getLRProbCut() ) < 0.0 ) ||
      ( rwd != NULL && ( cut = rwd->getLRProbCut() ) < 0.0 )) return true;

  // Compare clus LRProb and LRProbCut. 
  // return false if 0<=LRProb<LRProbCut
  double prob = clus.getLRProb();
  return ( prob >=0.0 && prob < cut ) ? false : true;
  
}


//_____________________________________________________________________________
#define nLRbin 50
void CsEventUtils::LRmonitor(list<CsCluster*> clus) {

  static bool first_call = true;  // to book histograms at first call
  static bool doLRHistos = false;

  // decide if histograms have to be done or not
  if (first_call) {
    string tag = "LR";
    string key = "make histograms";
    doLRHistos = (CsOpt::Instance()->getOpt( tag, key ) && CsEvent::Instance()->isAMonteCarloEvent()) ? true : false;
    if (!doLRHistos) CsErrLog::Instance()->mes(elInfo, "No Left/Right histogram monitoring." );
  }
  if (!doLRHistos) return;
  
  /* for each kind of histograms, three histos are booked: 
     0 is for drift chambers; 
     1 is for 6mm straw tubes;
     2 is for 10mm straw tubes;
     4 is for drift tubes;
  */
  
  // Book histograms
  static CsHist1D* nLR[6];            //Left/right probability distribution for all clusters
  static CsHist1D* nOKLR[6];          //Left/right probability distribution for true clusters (MC)    
  static CsHist1D* nLRCut[6];         //Left/right probability integrated distribution for true clusters
  static CsHist1D* nOKLRCut[6];       //Left/right probability integrated distribution for all clusters
  static CsHist1D* nDWire[6];         //Distance to wire for true (MC) clusters
  static CsHist1D* nOKDWire[6];       //Distance to wire for true (MC) clusters passing LRCut
    
  static CsHist1D* dAssoc[6];         //association between clusters
  static CsHist2D* tPoint[6];         //target pointing
  
  if (first_call) {
    first_call = false;
    
    char name[20], title[100];
    const char *detName[] = { "DC", "ST 6mm", "ST 10mm", "MB", "WD", "DR", 0 };
    
    double d_max[] = {
      5.0,  // max distance to wire for drift chambers (mm)
      3.0,  // max distance to wire for 6mm straws     (mm)
      5.0,  // max distance to wire for 10mm straws    (mm)
      16.75 // max distance to wire for drift tubes    (mm)
    };
        
    double tet_max[] = {
      0.50, // max angle for drift chambers (mm)
      1.00, // max angle to wire for 6mm straws     (mm)
      0.50, // max angle to wire for 10mm straws    (mm)
      0.15  // max angle to wire for drift tubes    (mm)
    };
		
		// histogram declaration
    CsHistograms::SetCurrentPath("/CsLRMonitor");
    for( int ih = 0 ; ih<6 ; ih++ ) {
      
      // number of clusters OK at a given LRProb
      sprintf(name,"nOKLR%i",ih);
      sprintf(title,"#clusters OK at given LRProb (%s)",detName[ih]); 
      nOKLR[ih] = new CsHist1D(name,title,nLRbin,-0.5/nLRbin,1-0.5/nLRbin);
      
      // number of clusters at a given LRProb
      sprintf(name,"nLR%i",ih);
      sprintf(title,"#clusters at given LRProb (%s)",detName[ih]); 
      nLR[ih] = new CsHist1D(name,title,nLRbin,-0.5/nLRbin,1-0.5/nLRbin);
      
      // number of clusters OK at a LRProb > LRProbCut
      sprintf(name,"nOKLRCut%i",ih);
      sprintf(title,"#clusters OK at LRProb>LRProbCut (%s)",detName[ih]); 
      nOKLRCut[ih] = new CsHist1D(name,title,nLRbin,-0.5/nLRbin,1-0.5/nLRbin);
      
      // number of clusters at a LRProb > LRProbCut
      sprintf(name,"nLRCut%i",ih);
      sprintf(title,"#clusters at LRProb>LRProbCut (%s)",detName[ih]); 
      nLRCut[ih] = new CsHist1D(name,title,nLRbin,-0.5/nLRbin,1-0.5/nLRbin);
        
      // number of clusters OK at a given distance to wire
      sprintf(name,"nOKDWire%i",ih);
      sprintf(title,"#clusters (cut > 0.5) vs distance to wire (%s)",detName[ih]); 
      nOKDWire[ih] = new CsHist1D(name,title,50,0,d_max[ih]);
      
      // number of clusters at a given distance to wire
      sprintf(name,"nDWire%i",ih);
      sprintf(title,"#clusters OK vs distance to wire (%s)",detName[ih]); 
      nDWire[ih] = new CsHist1D(name,title,50,0,d_max[ih]);
       
      // distance between associated clusters
      sprintf(name,"DAssoc%i",ih);
      sprintf(title,"#d between assoc cls (%s)",detName[ih]); 
      dAssoc[ih] = new CsHist1D(name,title,100,-5*d_max[ih],5*d_max[ih]);
      
      // target pointing
      sprintf(name,"TPoint%i",ih);
      sprintf(title,"Target pointing (%s)",detName[ih]); 
      tPoint[ih] = new CsHist2D(name,title,100,-tet_max[ih],tet_max[ih],100,-tet_max[ih],tet_max[ih]);
     }    
        
  } // End of booking block
  
  // Fill histograms

  list<CsCluster*>::iterator Ic;
  for (Ic = clus.begin(); Ic!=clus.end(); Ic++) {
    double prob = -1;

    // Do nothing if cluster isn't drift-like or association failed
    if (!(*Ic)->hasMirrorCluster() || (prob = (*Ic)->getLRProb())<0) continue;

    int ib = (int)(nLRbin*prob), jb = 0;

    CsDetector *det0 = (*Ic)->getDetsList().front();
    CsDigit *digit = (*Ic)->getDigitsList().front();
    if (!det0 || !digit) continue;
    CsMCDigit *mcdigit = dynamic_cast<CsMCDigit*>(digit);
    if (!mcdigit) continue;

    // Reject clusters not corresponding to correlated hits
    CsMCHit *hit = mcdigit->getHits().front();
    if (!hit || !CsMCUtils::isACorrelatedHit(hit)) continue;

    // Choose which row of histogram to fill
    int ih = -1;
    if      (strncmp(det0->getName(),"DC",2)==0) ih = 0;
    else if( strncmp(det0->getName(),"ST",2)==0) {
      if (fabs(det0->getWirP()-6)<1) ih = 1;                // 6mm straw tubes
      else                           ih = 2;                // 10mm straw tubes
    }
    else if (strncmp(det0->getName(),"MB",2)==0) ih = 3;
    else if (strncmp(det0->getName(),"DW",2)==0) ih = 4;
    else if (strncmp(det0->getName(),"DR",2)==0) ih = 5;
    if ( ih<0 ) continue;

    double d = digit->getDatum()*det0->getVel()/det0->getTiSli();

    nLR[ih]->Fill(prob);
    for( jb=0 ; jb<ib ; jb++ ) nLRCut[ih]->Fill((double)jb/nLRbin);
    
    if( CsMCUtils::isAGoodCluster( *Ic ) ) {
      nOKLR[ih]->Fill(prob);
      nDWire[ih]->Fill(d);
      if ( isALRGoodCluster( *(*Ic) ) ) nOKDWire[ih]->Fill(d);
      for( jb=0 ; jb<ib ; jb++ ) nOKLRCut[ih]->Fill((double)jb/nLRbin);
    }

    // get associated detector
    CsDriftChamberDetector *dcd = dynamic_cast<CsDriftChamberDetector*>(det0);
    CsStrawTubesDetector   *std = dynamic_cast<CsStrawTubesDetector*>(det0);
    CsDriftTubeDetector    *dtd = dynamic_cast<CsDriftTubeDetector*>(det0);
    CsDWDetector           *dwd = dynamic_cast<CsDWDetector*>(det0);
    CsRichWallDetector     *rwd = dynamic_cast<CsRichWallDetector*>(det0);

    CsDetector *det1 = NULL; //   ***** ASSOCIATE DETECTOR? *****
    int LRMode;  //    ***** GET Left/Right AMBIGUITY RAISING MODE *****
    if      (dcd) {
      det1 = dcd->getAssociateDet(); LRMode = dcd->getLRMode();
    }
    else if (std) {
      det1 = std->getAssociateDet(); LRMode = std->getLRMode();
    }
    else if (dtd) {
      det1 = dtd->getAssociateDet(); LRMode = dtd->getLRMode();
    }
    else if (dwd) {
      det1 = dwd->getAssociateDet(); LRMode = dwd->getLRMode();
    }
    else if (rwd) {
      det1 = rwd->getAssociateDet(); LRMode = rwd->getLRMode();
    }
    else continue;
    if (!det1) continue;

    // get det1 clusters
    list< CsCluster* > ac = det1->getMyClusters();
    list< CsCluster* >::iterator Ia;
    double dMin = -1;
    bool first = true;
    double uRef = (*Ic)->getU();
 	
    for( Ia = ac.begin(); Ia != ac.end(); Ia++ ) 
    if( first || TPQ( *Ic, *Ia, LRMode ) < fabs(dMin) ) { 
      first = false; 
      dMin = TPQ( *Ic, *Ia, LRMode );
    }
    dAssoc[ih]->Fill( dMin );
    
    
    //get associated cluster
    if( !(*Ic)->hasAssociateClusters() ) continue;
    CsCluster* a;
    if( (*Ic)->getLRProb() >= 0.5 ) a = (*Ic);
    else a = (*Ic)->getMirrorCluster();
    
    CsCluster* b = (*Ic)->getAssociateClusters().front();
    double t0 = ( a->getU() + b->getU() )/( a->getW() + b->getW() );
    double t1 = ( a->getU() - b->getU() )/( a->getW() - b->getW() );
    
    tPoint[ih]->Fill( t0, t1 );        
  }
    
  return;
}

//____________________________________________________________________________
void CsEventUtils::setLRToMatchAngle(CsCluster *c0, CsCluster *c1,
				     const double a) { // Incidence angle

  //   *************** RE-EVALUATE LR PROBAS GIVEN ARG INCIDENCE ***************

  // Rationale: Probability is OO = Optimum Optimorum (match to incidence "a")
  //               compared to o2 = 2nd optimum
  //              P = o2/(OO+o2)            
  // Technically: Let's have Dij = (u0-u1)/dZ-incidence, i,j = c,m
  //                         OO = Optimum Optimorum = max (Dij) = Dkl 
  //                         o2(0/1) = max over all combi's excluding k/l
  // Refinement:  To address the case where c0c1 is almost // to m0m1 and // to
  //            incidence "a" (which can yield an arbitrarily small "O"):
  //             optimum is bound by some lower limit, corresponding to
  //             the maximum angular resolution of the setup. 
  //              In order to keep things simple this lower bound is a const
  //             (=> no need to retrieve the detector's resolution => fast)
  //             = to the DC case = ~0.3 mm resolution / dZ = 7 mm

  CsCluster *m0 = c0->getMirrorCluster(), *m1 = c1->getMirrorCluster();
  if (m0==NULL && m1==NULL) return;  // Probas must have already been evaluated

  double dMin;          // Left/Right estimator for best combination
  double d[4];          // Left/Right estimator for each combination
  static double dLow = 0.3/7;
  double dZ = c0->getW() - c1->getW();

  //                                   1st combination: c0 and c1
  int flag = 0;   // Tells which combination is the best one
  d[0]=dMin = fabs((c0->getU()-c1->getU())/dZ-a);
  if (d[0]<dLow) d[0] = dLow;
  if (m1) {                         // 2nd combination: c0 and m1
    if ((d[1] = fabs((c0->getU()-m1->getU())/dZ-a))<dMin) {
      dMin = d[1]; flag = 1;
    }
    if (d[1]<dLow) d[1] = dLow; 
  }
  else d[1] = d[0];
  if (m0) {                         // 3rd combination: m0 and c1
    if ((d[2] = fabs((m0->getU()-c1->getU())/dZ-a))<dMin) {
      dMin = d[2]; flag = 2;
    }
    if (d[2]<dLow) d[2] = dLow;
  }
  else d[2] = d[0];
  if (m0 && m1) {                  // 4th combination: m0 and m1
    if ((d[3] = fabs((m0->getU()-m1->getU())/dZ-a))<dMin) {
      dMin = d[3]; flag = 3;
    }
    if (d[3]<dLow) d[3] = dLow;
  }
  else if (m1) d[3] = d[1];
  else if (m0) d[3] = d[2];
  else d[3] = d[0];

  // Removes associate clusters for c0, c1, m0, m1
  c0->removeAssociateClusters(); if (m0) m0->removeAssociateClusters();
  c1->removeAssociateClusters(); if (m1) m1->removeAssociateClusters();

  double a0, a1;  // LR probabilities for (c0,c1)
  // Associates clusters properly, so that first association corresponds to
  // the probability attributed to the cluster (i.e. not necessarily the best
  // one, if the cluster does not belong to the optimum optimorum)
  switch (flag) {
  case 0:                   // OO = c0,c1
    c0->addAssociateCluster(*c1); if (m1) c0->addAssociateCluster(*m1);
    c1->addAssociateCluster(*c0); if (m0) c1->addAssociateCluster(*m0);
    if (d[3]<d[2]) {          // o2(0) = m0,m1
      a0 = d[3]/(d[0]+d[3]);
      // Do we need to reshuffle associateship also for the mirrors of
      // argument cluster? So far no use is made of this reassociation.
      // Nevertheless let's do it (temporary decision, which we may want to
      // turn round, 'cause consuming CPU uselessly if not extensively).
      if (m0) {
	if (m1) m0->addAssociateCluster(*m1); m0->addAssociateCluster(*c1);
      }
    }
    else {                    // o2(0) = m0,c1
      a0 = d[2]/(d[0]+d[2]);
      if (m0) {
	m0->addAssociateCluster(*c1); if (m1) m0->addAssociateCluster(*m1);
      }
    }
    if (d[3]<d[1]) {          // o2(1) = m0,m1
      a1 = d[3]/(d[0]+d[3]);
      if (m1) {
	if (m0) m1->addAssociateCluster(*m0); m1->addAssociateCluster(*c0);
      }
    }
    else{                     // o2(1) = c0,m1
      a1 = d[1]/(d[0]+d[1]);
      if (m1) {
	m1->addAssociateCluster(*c0); if (m0) m1->addAssociateCluster(*m0);
      }
    }
    break;

  case 1:                   // OO = c0,m1 (m1 does exist)
    c0->addAssociateCluster(*m1); c0->addAssociateCluster(*c1);
    m1->addAssociateCluster(*c0); if (m0) m1->addAssociateCluster(*m0);
    if (d[2]<d[3]) {          // o2(0) = m0,c1
      a0 = d[2]/(d[1]+d[2]);
      if (m0) {
	m0->addAssociateCluster(*c1); m0->addAssociateCluster(*m1);
      }
    }
    else {                    // o2(0) = m0,m1
      a0 = d[3]/(d[1]+d[3]);
      if (m0) {
	m0->addAssociateCluster(*m1); m0->addAssociateCluster(*c1);
      }
    }
    if (d[2]<d[0]) {          // o2(1) = m0,c1
      a1 = d[1]/(d[1]+d[2]);
      if (m0) c1->addAssociateCluster(*m0); c1->addAssociateCluster(*c0);
    }
    else{                     // o2(1) = c0,c1
      a1 = d[1]/(d[1]+d[0]);
      c1->addAssociateCluster(*c0); if (m0) c1->addAssociateCluster(*m0);
    }
    break;

  case 2:                   // OO = m0,c1 (m0 does exist)
    m0->addAssociateCluster(*c1); if (m1) m0->addAssociateCluster(*m1);
    c1->addAssociateCluster(*m0); c1->addAssociateCluster(*c0);
    if (d[1]<d[0]) {          // o2(0) = c0,m1
      a0 = d[2]/(d[2]+d[1]);
      if (m1) c0->addAssociateCluster(*m1); c0->addAssociateCluster(*c1);
    }
    else {                    // o2(0) = c0,c1
      a0 = d[2]/(d[2]+d[0]);
      c0->addAssociateCluster(*c1); if (m1) c0->addAssociateCluster(*m1);
    }
    if (d[1]<d[3]) {          // o2(1) = c0,m1
      a1 = d[1]/(d[2]+d[1]);
      if (m1) {
	m1->addAssociateCluster(*c0); m1->addAssociateCluster(*m0);
      }
    }
    else{                     // o2(1) = m0,m1
      a1 = d[3]/(d[2]+d[3]);
      if (m1) {
	m1->addAssociateCluster(*m0); m1->addAssociateCluster(*c0);
      }
    }
    break;

  default: // Switch's condition "flag" can only take value 0,1,2 or 3
  case 3:                   // OO = m0,m1 (m0,m1 do exist)
    m0->addAssociateCluster(*m1); m0->addAssociateCluster(*c1);
    m1->addAssociateCluster(*m0); m1->addAssociateCluster(*c0);
    if (d[0]<d[1]) {          // o2(0) = c0,c1
      a0 = d[3]/(d[3]+d[0]);
      c0->addAssociateCluster(*c1); c0->addAssociateCluster(*m1);
    }
    else {                    // o2(0) = c0,m1
      a0 = d[3]/(d[3]+d[1]);
      c0->addAssociateCluster(*m1); c0->addAssociateCluster(*c1);
    }
    if (d[0]<d[2]) {          // o2(1) = c0,c1
      a1 = d[3]/(d[3]+d[0]);
      c1->addAssociateCluster(*c0); c1->addAssociateCluster(*m0);
    }
    else{                     // o2(1) = m0,c1
      a1 = d[3]/(d[3]+d[2]);
      c1->addAssociateCluster(*m0); c1->addAssociateCluster(*c0);
    }
  }
  
  // Fill left/right probabilities
  
  if (m0) {
    c0->setLRProb(a0); m0->setLRProb(1-a0);
  }
  else c0->setLRProb(1);
  if (m1) {
    c1->setLRProb(a1); m1->setLRProb(1-a1);
  }
  else c1->setLRProb(1); 
    
  return;
}  

