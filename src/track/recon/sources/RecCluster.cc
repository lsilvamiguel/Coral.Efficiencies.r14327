#include <iostream>
#include "Coral.h"
#include "CsGeom.h"
#include "CsEventUtils.h"
#include "RecCall.h"
#include "CsCluster.h"
#include "CsMCUtils.h"
#include "RecOpt.h"
#include "CsTrack.h"
#include "CsHelix.h"
#include "CsMCTrkHit.h"
#include "CsMCDigit.h"
#include "CsStopwatch.h"
#include "Recon.h"

using namespace std;

 const int MAX3=300;  //parametr ntrack z PARAM.h
extern "C"
 void recon_(int*, int*, int*, int[], int[],int[], double[], double[],double[],
                   double[], double[],double[],
            int*, int[], int[], int[], int[][MAX3], int[][MAX3],int[][MAX3],
                       float[], float[], float[][5], float[][5] ,
                       float[], float[], float[][4], float[][4], int*);
extern "C"
void  mcini_(int*, int*,int*, int[], int[],int[], double[], double[], double[], double[],double[], double[]
                     , int[], int[],int[], int* , double[], double[][6], int[] );



void RecCall::RecClus(list<CsTrack*>& lCsTrk )
{

  CsGeom*  geom = CsGeom::Instance();  
  CsEvent* ev= CsEvent::Instance();
  CsStopwatch stopwatches;

  list<CsCluster*> listReClusters; listReClusters.clear();   // store ptr  to CsCluster
  list<CsCluster*> RecClus;   RecClus.clear();  // main clusters list
  list<CsCluster> RecObj;          // imported list of clusters  after Traffic (or TRAFDIC) 
  list<CsCluster*> ClusUnused;     // list of uused clusters after reconstruction 

  // import from coral
  const int MAX1=9000;    const int MAX2=5000;
  int beforeSM1IdTrack[MAX1], afterSM1IdTrack[MAX1], afterRICHIdTrack[MAX1];
  int beforeSM1idPlane[MAX2], afterSM1idPlane[MAX2], afterRICHidPlane[MAX2]; 
  double beforeSM1posHit[MAX2], afterSM1posHit[MAX2], afterRICHposHit[MAX2];
  double beforeSM1err[MAX2], afterSM1err[MAX2], afterRICHerr[MAX2];
   CsCluster* beforePointClusSM1[MAX2]; 
   CsCluster* afterPointClusSM1[MAX2]; 
   CsCluster*  afterRICHPointClus[MAX2] ;
   int ic, ic_mirr, ic_good;  // number of all input  clusters,   mirror clusters,  good clusters after  LR 
       ic=ic_mirr=ic_good=0; 
   int ic1,ic2,ic3; // number of clusters in different "recon zones"
       ic1=ic2=ic3=0;

// MC interf
 int ind1,ind2,ind3;    ind1=ind2=ind3=1; 
 int ltruereal,ltruereal1,ltruereal2,lmirrorreal1,lmirrorreal2,lmirror1,lmirror2;
     ltruereal=ltruereal1=ltruereal2=lmirrorreal1=lmirrorreal2=lmirror1=lmirror2=0;
 
  //export to Coral
  int NbRec; NbRec=0;    // nb. of reconstructed tracks
  int NbHitsBeforeSM1[MAX3], NbHitsBehindSM1[MAX3], NbHitsBehindRICH[MAX3];
  int HitNb1[100][MAX3], HitNb2[100][MAX3], HitNb3[100][MAX3];  
  float ChiSq[MAX3], ZposHelix[MAX3], ParamHelix[MAX3][5], DiagCovMatr[MAX3][5]; //first helix
  float ChiSq1[MAX3], ZposHelix1[MAX3], ParamHelix1[MAX3][4], DiagCovMatr1[MAX3][4]; //second helix
  double pTrack[500] ;
  double pHit[500][6];   // momentum in first plane 
 
//  SCHEME 1-RECON, 2-RECON_AFTER_TRAFFIC 
   int chrono3 = stopwatches.start(); // start a chronometer 
   
  const list<CsCluster*> &AllClus=ev->getClusters();   
 
if (Recon::ref().isScheme1==true ) {   
    for(list<CsCluster*>::const_iterator Itc=AllClus.begin(); Itc!=AllClus.end(); Itc++){
     if ( (*Itc)->getW() >  Recon::ref().magType[1] || (*Itc)->getW() < 0.0) continue;
     RecClus.push_back(*Itc);   }
//    printf ("Scheme 1:  AllClus %d RecClus %d \n", AllClus.size(),RecClus.size());
 }
 else  { RecCall::GetClusters(RecObj); 
  for(list<CsCluster*>::const_iterator Itc=AllClus.begin(); Itc!=AllClus.end(); Itc++){
     if ( (*Itc)->getW() >  Recon::ref().magType[1] || (*Itc)->getW() < 0.0) continue;
         for(list<CsCluster>::const_iterator  m=RecObj.begin(); m!=RecObj.end(); m++) {
	  if(  (*Itc)->getW() == m->getW() && (*Itc)->getU()  == m->getU() ) { 
            RecClus.push_back(*Itc);  break;}  
  // better but longer:  if( *(*Itc)==*m) { RecClus.push_back(*Itc); break;}
 	} // loop clus. (RecObj)  
     } // loop clus All
 // printf ("Scheme 2:  RecClus %d  %d \n",RecObj.size(),RecClus.size());  
 } // end  isScheme
 
 float time3 = stopwatches.stop(chrono3); // stop a chronometer
 Recon::ref().countTime3(time3);

// trigger mask
    unsigned trygger_mask=ev->getDaqEvent().GetHeader().GetTrigger();
      int trg=(int) trygger_mask;
      //       printf(" trigger mask %d \n",trg);

//------------------------MAIN LOOP-------------------------------------------
//-----------------------------------------------------------------------------
 for (list<CsCluster*>::iterator Ic=RecClus.begin(); Ic!=RecClus.end(); Ic++) { 
  listReClusters.push_back(*Ic);  // remember pointer to cluster 
   ic++; // counter for all clusters
  	
   list<CsDetector*> lDet = (*Ic)->getDetsList(); // get detector for cluster
      if ( lDet.empty() ) { printf("RecCluster ==> CsCluster has no reference to CsDetector \n"); continue;}
      if ( lDet.size()> 1 ) { printf("RecCluster ==> CsCluster belongs to more than one detector \n"); continue;}
   CsDetector* pDet = *(lDet.begin());
   
   // check LR  for each cluster  
   CsEventUtils LRutil; bool LR=LRutil.isALRGoodCluster(*(*Ic));
   if (! LR) {  // printf(" Mirror Cluster %d  %s \n",ic,pDet->GetTBName().c_str());  
     ic_mirr++; continue; }
   ic_good++;
   HepMatrix sC;  sC=(*Ic)->getCov() ;

   if ( (*Ic)->getW() < Recon::ref().magType[0] &&  (*Ic)->getW() > 0.0){          //zone1   
     ic1++;
     beforePointClusSM1[ic1-1]=(*Ic);                 //array inside C code
     beforeSM1idPlane[ic1-1]=pDet->GetID();           // detector ID
     beforeSM1posHit[ic1-1]=(*Ic)->getU();            // cluster coordinate
     beforeSM1err[ic1-1]=sqrt(sC(1,1));	              // cluster coordinate error 
    //   printf("zone 1:  %d  %s  %5.3f  %4d,  %5.3f,  %5.3f  \n",ic1,pDet->GetTBName().c_str(),(*Ic)->getW(),
      //    beforeSM1idPlane[ic1-1],beforeSM1posHit[ic1-1],beforeSM1err[ic1-1]  );
   } 
   else if ( (*Ic)->getW() > Recon::ref().magType[0] &&  (*Ic)->getW() < 7000.) {  // zone 2
     ic2++;
     afterPointClusSM1[ic2-1]=(*Ic);
     afterSM1idPlane[ic2-1]=pDet->GetID();
     afterSM1posHit[ic2-1]=(*Ic)->getU();
     afterSM1err[ic2-1]=sqrt(sC(1,1));
     // printf("zone 2:  %d  %s  %5.3f %4d,  %5.3f,  %5.3f  \n",ic2,pDet->GetTBName().c_str(),(*Ic)->getW(),
       //    afterSM1idPlane[ic2-1],afterSM1posHit[ic2-1],afterSM1err[ic2-1]  );
     }
   else if ( (*Ic)->getW() > 7000. &&   (*Ic)->getW() < Recon::ref().magType[1] ) {  // zone 3 
     ic3++;
      afterRICHPointClus[ic3-1]=(*Ic);
     afterRICHidPlane[ic3-1]=pDet->GetID();
     afterRICHposHit[ic3-1]=(*Ic)->getU();
     afterRICHerr[ic3-1]=sqrt(sC(1,1));
    //  printf("zone 3:  %d %s  %5.3f  %4d,  %5.3f,  %5.3f  \n",ic3,pDet->GetTBName().c_str(),(*Ic)->getW(),
      //    afterRICHidPlane[ic3-1],afterRICHposHit[ic3-1],afterRICHerr[ic3-1]  );
   } // end of "end if" statemant for different zones


        if (! Recon::ref().isMonteCarlo) continue;
	printf(" MC interface under construction !! \n"); exit(1);
//------------------------------------------------------------------------------------------
//----------------------------- MC INTERFEJS------------------------------------------------

   list<CsDigit*> ldig=(*Ic)->getDigitsList();
   list<CsMCHit*> lhit;

    int ihit=0; int ihit1=0;  int ihit2=0;  int ihit3=0;
  
  for(  list<CsDigit*>::iterator Id=ldig.begin(); Id!=ldig.end(); Id++ ) { //loop on digits
      CsMCDigit* theMCDigit = dynamic_cast<CsMCDigit*>(*Id); 
      if( theMCDigit != NULL )   lhit = theMCDigit->getHits();
      
      for(list<CsMCHit*>::iterator  Ih=lhit.begin(); Ih!=lhit.end(); Ih++ ) { // loop on hits  
	ihit++;
	CsMCTrack* trak=(*Ih)->getMCTrack();
   if ( (*Ic)->getW() < Recon::ref().magType[0] &&  (*Ic)->getW() > 0.0){          //zone1   
    ihit1++;
            if ( (*Ih)->getOrigin() == 0 ){
	             if (CsMCUtils::isAGoodCluster(*Ic)){
                           beforeSM1IdTrack[ind1+ihit1]=trak->getGnum();
                           ltruereal++; ltruereal1++;}
                     else
                        {lmirror1++;
                          beforeSM1IdTrack[ind1+ihit1]=200;}
                     }//endif
           else{  
                   beforeSM1IdTrack[ind1+ihit1]=200;
                   if (CsMCUtils::isAGoodCluster(*Ic)){
                       beforeSM1IdTrack[ind1+ihit1]=trak->getGnum();
                      lmirrorreal1++;}
             }

} // end of zone 1 
   else if ( (*Ic)->getW() > Recon::ref().magType[0] &&  (*Ic)->getW() < 7000.) {  // zone 2
     ihit2++;
       if ( (*Ih)->getOrigin() == 0 ){
             if (CsMCUtils::isAGoodCluster(*Ic)){
               afterSM1IdTrack[ind2+ihit2] =trak->getGnum();
               ltruereal++; ltruereal2++;}
            else
              { lmirror2++;
                afterSM1IdTrack[ind2+ihit2]=200;}
              }//endif
        else
            {
             if (CsMCUtils::isAGoodCluster(*Ic)){
             afterSM1IdTrack[ind2+ihit2]=trak->getGnum();
             lmirrorreal2++;}
             afterSM1IdTrack[ind2+ihit2]=200;
         }
    

  } //end of zone 2
   else if ( (*Ic)->getW() > 7000. &&   (*Ic)->getW() < Recon::ref().magType[1] ) {  // zone 3 
     ihit3++;

     if ( (*Ih)->getOrigin() == 0 )
             afterRICHIdTrack[ind3+ihit3]=trak->getGnum();
           else
             afterRICHIdTrack[ind3+ihit3]=200;
  } //end of zone3
      } //  end of loop on  hits 
    } // end of loop on  digits

  if ( (*Ic)->getW() < Recon::ref().magType[0] &&  (*Ic)->getW() > 0.0){          //zone1   
         beforeSM1IdTrack[ind1]=ihit1;
          beforeSM1IdTrack[ind1-1]=ic1;
          ind1=ind1+2+ihit1; }
  else if ( (*Ic)->getW() > Recon::ref().magType[0] &&  (*Ic)->getW() < 7000.) {  // zone 2
          afterSM1IdTrack[ind2]=ihit2;
          afterSM1IdTrack[ind2-1]=ic2;
          ind2=ind2+2+ihit2;}
  else if ( (*Ic)->getW() > 7000. &&   (*Ic)->getW() < Recon::ref().magType[1] ) {  // zone 3 
          afterRICHIdTrack[ind3]=ihit3;
          afterRICHIdTrack[ind3-1]=ic3;
          ind3=ind3+2+ihit3;   } 
//------------------ END OF MC INTERFEJS-----------------------------------------------------
//-------------------------------------------------------------------------------------------  
 
 }  //end cluster loop
//  printf( " Zone 1 %d, zone 2 %d, zone 3 %d suma %d  \n",ic1,ic2,ic3,ic);
//    printf( "ReconCluster size %d  \n",listReClusters.size());
//---------------------------END OF MAIN LOOP--------------------------------------------------
//---------------------------------------------------------------------------------------------


//---------------------------------NEXT MC PART------------------------------------------------
  double px, py, pz, pp;    px=py=pz=pp=0 ;
   int k=0;
   int iflag[1000];

if (Recon::ref().isMonteCarlo ) {
     CsMCUtils  pileup;
   list<CsMCTrack*> tra = Coral::Instance()->getMCTracks();  list<CsMCTrack*>::iterator  It;
   
   for( It=tra.begin(); It!=tra.end(); It++ ) {
     k++;
     px=(*It)->getPX();
     py=(*It)->getPY();
     pz=(*It)->getPZ();
     pTrack[k-1]=sqrt(px*px+py*py+pz*pz);
         bool status=pileup.isACorrelatedTrack( *It );
     if (status==false)  { pTrack[k-1]=-pTrack[k-1];}
     
     list<CsMCHit*> hits=(*It)->getMCHits();   list<CsMCHit*>::iterator  ih;
     int nhitBefore,nhitAfter;  nhitBefore=nhitAfter=0;
     
     for( ih=hits.begin(); ih!=hits.end(); ih++ ) {
       if ((*ih)->getOrigin() == 0 ) {
	 if ( (*ih)->getZ() <   Recon::ref().magType[0] &&  (*ih)->getZ()  > 0 ) 
	   nhitBefore++;
	 else
	   if ( (*ih)->getZ() < 7000. && (*ih)->getZ()  > 0 )
	     nhitAfter++;
       }
     } //end of hits list
     
     if(hits.empty()){
	pHit[k-1][0]=-1000.;
	pHit[k-1][1]=-1000.;             
	pHit[k-1][2]=-1000.;
	pHit[k-1][3]=pHit[k-1][4]=pHit[k-1][5]=0.;       
      } else {
	pHit[k-1][0]=(hits.front()->getP()).x();
	pHit[k-1][1]=(hits.front()->getP()).y();             
	pHit[k-1][2]=(hits.front()->getP()).z();
	pHit[k-1][3]=hits.front()->getZ();             
	pHit[k-1][4]=hits.front()->getY();
	pHit[k-1][5]=hits.front()->getX();  }
      
      if ( nhitBefore >= RecOpt::McPar[1] ) iflag[k-1]=333;    
      else  iflag[k-1]=999;
   } //end of tracks list
             
  }//end if statement
//-----------------------------END  MC PART-----------------------------------------------------------

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
 int chrono1 = stopwatches.start(); // start a chronometer
   if (Recon::ref().isMonteCarlo) {
       mcini_(&ic1, &ic2, &ic3, beforeSM1idPlane, afterSM1idPlane,afterRICHidPlane,
  	       beforeSM1posHit, afterSM1posHit, afterRICHposHit,
	       beforeSM1err, afterSM1err,afterRICHerr,
               beforeSM1IdTrack, afterSM1IdTrack,afterRICHIdTrack, &k, pTrack, pHit, iflag);
   }
    
    int chrono2 = stopwatches.start(); // start a chronometer
     
      recon_(&ic1, &ic2, &ic3, beforeSM1idPlane, afterSM1idPlane,afterRICHidPlane,
  	 beforeSM1posHit, afterSM1posHit,afterRICHposHit,
         beforeSM1err, afterSM1err,afterRICHerr, 
	 &NbRec, NbHitsBeforeSM1, NbHitsBehindSM1,NbHitsBehindRICH,  HitNb1, HitNb2, HitNb3,
  	 ChiSq, ZposHelix, ParamHelix, DiagCovMatr, ChiSq1, ZposHelix1, ParamHelix1, DiagCovMatr1,&trg );

   float time2 = stopwatches.stop(chrono2); // stop a chronometer
   Recon::ref().countTime2(time2);

   float time1 = stopwatches.stop(chrono1); // stop a chronometer
   Recon::ref().countTime1(time1);
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------

//-------------------------------EXPORT TO CORAL--------------------------------------------------
  CsTrack* pCsTra;
  list<CsCluster*> ReconClusters;
  list<CsZone*> ReconZones;
  vector<CsHelix> ReconHelices;
  list<CsCluster*> listEvUsedClus;   listEvUsedClus.clear();

  for(int trzr=0; trzr<NbRec; trzr++)
    {
      double cov[15], cov1[15];
      for(int ic=0; ic<14; ic++){
	cov[ic]=0;   
        cov1[ic]=0; 
      }
      ReconZones.clear();
      ReconHelices.clear();
      //prepare helix
      double x,y,z,dxdz,dydz,p;
      double x1,y1,z1,dxdz1,dydz1,p1;
      
      x=(double)((ParamHelix[trzr][0])*1000);
      y=(double)((ParamHelix[trzr][1])*1000);
      z=(double)((ZposHelix[trzr])*1000);
      dxdz=(double)ParamHelix[trzr][2];
      dydz=(double)ParamHelix[trzr][3];
      p=(double)ParamHelix[trzr][4];
      cov[0]=(double)((DiagCovMatr[trzr][0])*1000);
      cov[2]=(double)((DiagCovMatr[trzr][1])*1000);
      cov[5]=(double)DiagCovMatr[trzr][2];
      cov[9]=(double)DiagCovMatr[trzr][3];
      cov[14]=(double)DiagCovMatr[trzr][4];
      
      cov[1]=0;
      cov[3]=0; cov[4]=0;
      cov[6]=0; cov[7]=0; cov[8]=0;
      cov[10]=0; cov[11]=0; cov[12]=0; cov[13]=0;
 
      //second helix 
      x1=(double)((ParamHelix1[trzr][0])*1000);
      y1=(double)((ParamHelix1[trzr][1])*1000);
      z1=(double)((ZposHelix1[trzr])*1000);
      dxdz1=(double)ParamHelix1[trzr][2];
      dydz1=(double)ParamHelix1[trzr][3];
      p1=(double)ParamHelix[trzr][4];
      cov1[0]=(double)((DiagCovMatr1[trzr][0])*1000);
      cov1[2]=(double)((DiagCovMatr1[trzr][1])*1000);
      cov1[5]=(double)DiagCovMatr1[trzr][2];
      cov1[9]=(double)DiagCovMatr1[trzr][3];
      cov1[14]=(double)DiagCovMatr[trzr][4];

      cov1[1]=0;
      cov1[3]=0; cov1[4]=0;
      cov1[6]=0; cov1[7]=0; cov1[8]=0;
      cov1[10]=0; cov1[11]=0; cov1[12]=0; cov1[13]=0;     
      
      CsHelix helix(x,y,z,dxdz,dydz,p,cov);
      CsHelix helix1(x1,y1,z1,dxdz1,dydz1,p1,cov1);
      ReconHelices.push_back( helix );
      ReconHelices.push_back( helix1 );
      //prepare chi square
      double  chi2=(double)ChiSq[trzr];
   //prepare zones
      list<CsZone*> zones = geom->getZones();
      list<CsZone*>::iterator Iz;
      for(Iz=zones.begin(); Iz != zones.end(); Iz++)
	{ 
	  CsZone* theZone=*Iz;
	  	
          if ((*Iz)->getZMax() < Recon::ref().magType[1] ){
	    ReconZones.push_back(theZone);
	  }
	}// end of loop over zones
      
      //prepare clusters

      ReconClusters.clear();
      
      for(int nrpl1=0; nrpl1<NbHitsBeforeSM1[trzr]; nrpl1++){
	if(beforePointClusSM1[(HitNb1[nrpl1][trzr])-1] == NULL) {

	  cout<<"RecCall::RecCluster-before SM1 ==> Cluster # "<< beforePointClusSM1[(HitNb1[nrpl1][trzr])-1]
    	      <<" has NULL pointer to CsCluster"<<endl;
	  exit(1);
	}
		
	ReconClusters.push_back(beforePointClusSM1[(HitNb1[nrpl1][trzr])-1]);
	listEvUsedClus.push_back(beforePointClusSM1[(HitNb1[nrpl1][trzr])-1]);
	 if (beforePointClusSM1[(HitNb1[nrpl1][trzr])-1]->hasMirrorCluster()== true){
	    CsCluster* ptrMirrorBef=beforePointClusSM1[(HitNb1[nrpl1][trzr])-1]->getMirrorCluster();
	    listEvUsedClus.push_back(ptrMirrorBef);
	  }
         
	      } // end of hits loop before SM1
      
      for(int nrpl2=0; nrpl2<NbHitsBehindSM1[trzr]; nrpl2++){
        if(afterPointClusSM1[(HitNb2[nrpl2][trzr])-1] == NULL) {
       cout<<"RecCall::RecCluster-behind SM1 ==> Cluster # "<<afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]
	   <<" has NULL pointer to CsCluster"<<endl;
       exit(1);
	}
	
	//for cluster get detector list
	list<CsDetector*> ldet= afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]->getDetsList();
	
	if (ldet.empty()){
	  cout << "RecCall::RecClus ===>>CsCluster has no reference to CsDetector" << endl;
	  break;
	}
	if (ldet.size() > 1 ) {
	  cout << " RecCall::RecClus ===>>CsCluster belong to more then one CsDetector " << endl;
	  continue;
	}

	CsDetector* ptrDet=*(ldet.begin());

	//only for straw  -> 2001 option: 11 changed for 11000 - kk

	if (ptrDet->getType()==11000) {
	  // cout<< ldet.size() << " type: " << ptrDet->getType() << "  " << ptrDet->getName() << endl;
        
	  if (afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]->hasAssociateClusters()==true){
	    
	    list<CsCluster*> assocClus=afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]->getAssociateClusters();
	    list<CsCluster*>::iterator itc;
	    CsCluster* ptrTrue=*(assocClus.begin());
	    ReconClusters.push_back(ptrTrue); 
	    //     cout <<  " assoc. true" << ptrTrue->getU() << "  " << ptrTrue->getLRProb()  << "      " << endl;
	    for ( itc=assocClus.begin(); itc!=assocClus.end(); itc++) {
	      CsCluster* ptrAssocClusters=*itc;
	      if (ptrAssocClusters== NULL) {
		cout << "RecCluster.cc:  NULL pointer to CsCluster"<<endl;
		break;
	      }
	      listEvUsedClus.push_back(ptrAssocClusters);
	    } //end of for assocClus loop
	    
	  } //end if cluster has association
	  
	  else if ( afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]->hasAssociateClusters()==false)
	    {  
	      //  cout << trzr+1 << "    clus " << HitNb2[nrpl2][trzr]  << "  has NO ASSOCATION!! " << endl;
	     
	    } //end of if-else-if
	  //mirror for recon (orginal) cluster
	  if (afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]->hasMirrorCluster()== true){
	    CsCluster* ptrMirror=afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]->getMirrorCluster();
	    listEvUsedClus.push_back(ptrMirror);
	  }
         
	} //end of if straw det.
	
	ReconClusters.push_back(afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]);                          
	listEvUsedClus.push_back(afterPointClusSM1[(HitNb2[nrpl2][trzr])-1]);

      } // end of hits loop before SM2
   
      for(int nrpl3=0; nrpl3< NbHitsBehindRICH[trzr]; nrpl3++){
        if(afterRICHPointClus[(HitNb3[nrpl3][trzr])-1] == NULL) {
       cout<<"RecCall::RecCluster ==> Cluster # "<<afterRICHPointClus[(HitNb3[nrpl3][trzr])-1]
	   <<" has NULL pointer to CsCluster"<<endl;
       exit(1);
	}
	
	ReconClusters.push_back(afterRICHPointClus[(HitNb3[nrpl3][trzr])-1]);
	listEvUsedClus.push_back(afterRICHPointClus[(HitNb3[nrpl3][trzr])-1]);
	

      } // end of hits loop after RICH behind SM2
      
      CsTrack* pCsTra=new CsTrack(ReconHelices,ReconClusters,ReconZones,chi2);
      lCsTrk.push_back(pCsTra);
    } //end of reconstructed tracks loop
   delete  pCsTra;

   //below commented by kk June 25 2002
   /*
  Recon::ref().listUnusedClus.clear();
  for(list<CsCluster*>::iterator RecCl=listReClusters.begin(); RecCl!= listReClusters.end(); RecCl++) {	
    if( find( listEvUsedClus.begin(),listEvUsedClus.end(),*RecCl ) != listEvUsedClus.end()){ //nothing 
    } 
     else { Recon::ref().listUnusedClus.push_back(*RecCl); }
}

*/
  //end of kk comment June 25 2002
/*       
    cout << "RecCluster: size of  input clusters list: " << listReClusters.size() << endl;
    cout << "RecCluster: size of  used clusters list: " <<listEvUsedClus.size() <<  endl; 
    cout << "RecCluster: size of unused clusters list: " << Recon::ref().listUnusedClus.size() << endl;
     cout << "nb. of tracks give back  to coral; " << lCsTrk.size() << endl ;
*/       

} //end all


