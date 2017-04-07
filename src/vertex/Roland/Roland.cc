// $Id: Roland.cc,v 1.4 2003/04/24 08:24:20 benigno Exp $

/*!
   \file    Roland.cc
   \brief   Vertex Reconstruction Procedure
   \author  C.Ulvegren
   \version $Revision: 1.4 $
   \date    $Date: 2003/04/24 08:24:20 $
*/

#include "coral_config.h"

#include "CsRolandPattern.h"
#include "CsMCUtils.h"
#include "CsGeom.h"
#include "CsGeant3.h"
#include "CsEvent.h"

using namespace std;

//For solution of linear system of equations
extern "C" {void reqn_(int &,float *a,int &,int* r,int &,int &,float* bt);}

//-------------------------------------------------------------------
//    Roland Windmolders vertex finder
//____________________________________________________________________

CsVertex* CsRolandPattern::getPrimaryVertex(float* vtx2,int *nflag,int& vertflag,
				   int & bad,int & less,int & muon,int & muon2)
{  
  CsEvent* event = CsEvent::Instance();
  list<CsTrack*> Beam = beams_; //CsBeamRecons::Instance()->getBeam();
  
  list<CsTrack*> tracks= tracks_;  //event->getTracks();
  int mlc=tracks.size();
  
  if( ! event->isAMonteCarloEvent() ) 
    CsErrLog::Instance()->mes( elFatal, "Roland's algorithm work only with MC data.");

  list<CsMCTrack*> Mtracks=event->getMCTracks();
  list<CsMCTrack*>::iterator itM=++Mtracks.begin();
  const CsMCTrack* mctr=(*itM);
  // Muon track is found for 98% of similarity between MCmuon and track
  float minhit=98.0; //Input to getAssociatedTracks()    
  CsMCUtils u; 
  list<CsTrack*> vctracks =u.getAssociatedTracks(mctr,minhit);
  list<CsTrack*>::iterator itv;
  int id=0;
  if(vctracks.size()!=0){
    itv=vctracks.begin();
    id=(*itv)->getId();
    muon++;
  }

  //Temporary beam parameters
  float bvec[6],p=0;

  list<CsTrack*>::iterator itb=Beam.begin();
  vector<CsHelix> beampar=(*itb)->getHelices();
  /*
  cout<<"Beampar :"<<"x= "<<beampar[0].getX()
      <<"y= "<<(beampar[0]).getY()<<flush
      <<"z= "<<(beampar[0]).getZ()<<flush
      <<"dxdz= "<<beampar[0].getDXDZ()<<flush
      <<"dydz= "<<beampar[0].getDYDZ()<<flush
      <<"mom= "<<beampar[0].getCop()<<endl<<flush;
  */
  bvec[5]=(beampar[0]).getZ();
  bvec[0]=(beampar[0]).getX();
  bvec[1]=(beampar[0]).getY();
  bvec[2]=(beampar[0]).getDXDZ();
  bvec[3]=(beampar[0]).getDYDZ();
  bvec[4]=(beampar[0]).getCop();
       
        
  //--------------------------------------------------------------------------
  // Local variables for choosing tracks to be used in the vertexdetermination
  //--------------------------------------------------------------------------
  //Tolerances 
  float tol[3]={40,10,5}; 
  
  tol[0] = c_[0];
  tol[1] = c_[1];
  tol[2] = c_[2];

  //Limit for how close tracks can be
  const float c=1.0;
  
  //Differences between zcalc and mcvertex zcoord
  float difz,difz2;//difsum=0.0;
  
  //For flagsoutput
  int flag;
  int kp=0,jflag=0,kpmu=0,nmu=0,mlassc=0;//flags #tracks and #muons
  
  //Array for trackparameters for test that tracks not the same
  float *testpar=new float[6];
  
  float **par;
  par=new float*[mlc];
  for(int ii=0;ii<mlc;ii++)
    par[ii]=new float[6];
  
  //For vertexvalues
  float *vtxbest=new float[3];
  
  //loop on all tracks and determination of primary vertex for this event
  list<CsTrack*>::iterator it;
  for(it=tracks.begin();it!=tracks.end();it++) {
    //cout<<"I loop on tracks "<<endl<<flush;
    int jtest=0,inat=0;
    float dif1,dif2,dif3,dif4;
    inat=(*it)->getId();//Identification of track
    //cout<<"ID :"<<inat<<endl;
    vector<CsHelix> vetrack=(*it)->getHelices();
    vector<CsHelix> :: iterator v=vetrack.begin();
    //Regard tracks too far away from the target or with momentum = 0
    if((vetrack[0]).getZ()>1400.0 ||
       fabs( (vetrack[0]).getCop())<1.e-7 ||
       (1./fabs((vetrack[0]).getCop())<0.1))
      continue;

    testpar[0]=(vetrack[0]).getX();       // xcoord
    testpar[1]=(vetrack[0]).getY();       // ycoord
    testpar[2]=(vetrack[0]).getDXDZ();    // slopes
    testpar[3]=(vetrack[0]).getDYDZ();
    testpar[4]=((vetrack[0]).getCop());   // momentum
    testpar[5]=(vetrack[0]).getZ();   
    
    //Sort tracks and disregard tracks that are too close
    
    for(int kpp=0;kpp<kp;kpp++) {  
      dif1= testpar[0]-par[kpp][0];
      dif2= testpar[1]-par[kpp][1];
      dif3= testpar[2]-par[kpp][2];
      dif4= testpar[3]-par[kpp][3];
      if(fabs(dif1)<c&&fabs(dif2)<c&&fabs(dif3)<c&&fabs(dif4)<c)
	jtest=1;
    }
    
    if(jtest==0) {
      //cout<<" I if sats 2"<<"\n";
      kp++;
      //cout<<"KP : "<<kp<<endl;
      if(inat==id){
	//cout<<"I muon ifsats"<<endl<<flush;
	kpmu=kp;
	nmu++;
	//cout<<"NMU :"<<nmu<<"kpmu :"<<kpmu<<endl;
      }              
      for(int n=0;n<6;n++) {
	par[kp-1][n]= testpar[n];
	if(n==4){ 
	  par[kp-1][n]=1.0/testpar[n];
	  //cout<<par[kp-1][n]<<' ';
	}
      }
      //cout<<endl;
      //cout<<"Slut pa if sats 2"<<"\n"<<flush;
    }
  }
  mlassc=kp;
  //cout<<endl<<"===============evt,ntr: "<<

  if(mlassc<2){
    less++;
    //cout<<"Less than two tracks "<<endl;
    CsVertex* vertex=new CsVertex(-99999.0,-99999.0,-99999.0);
    vertflag=0;
    nflag[jflag+4]++;
    delete [] par;
    return(vertex);
  }
  
  if(nmu!=1) kpmu=0;
  if(kpmu!=0)muon2++;
  vtxfit(mlassc,par,kpmu,bvec,tol,vtxbest,vtx2,jflag); 
  //cout<<"Jflag"<<' '<<jflag<<endl;
  flag=jflag;
  vertflag=jflag;
  CsVertex* vertex=new CsVertex(vtxbest[0],vtxbest[1],vtxbest[2]);
  
  nflag[(jflag+4)]++;
  delete [] par;
  return (vertex);
  
}


    	 
int CsRolandPattern::vtxfit(int & mlassc,float **par,int & jmu,float* prbeam,
			   float* tol,float* vtxbest,float *vtx2,int & jret) {
  
  //-----------------------------------------------------------------------------
  // Local variables
  //----------------------------------------------------------------------------- 

  const int mxll=50;
  const int mxcc=400;
  const float pmin=0.5;
  int mleft,mml;                // mleft #tracks, mml=#tracks after filtervtx
  float vec[2][6];
  float vecbeam[6];
  float **vtk;
  vtk=new float *[mxll];
  for(int ii=0;ii<mxll;ii++)
    vtk[ii]=new float[6];
  float dist,distmax;
  float zapr,xapr,yapr,disapr;  // For results of quickvtx, first call
  float zapp,xapp,yapp,disapp=0.0;  // For results of quickvtx, second call
  float drad;                   // distance from beam to tracks
  float distfin;                // Average distance for several tracks
                                // to the vertex
  int list[mxcc];
  int iveto[mxcc];
  int nveto=0;
  int jfilerr;                  // Flag to filterfunction, filtervtx
  int iamax,iasec;
  
  //-----------------------------------------------------------------------------
  //              Take away tracks with too low momentum
  //-----------------------------------------------------------------------------
  //cout<<"\n"<<" I vtxfit "<<flush;
  mleft=0;
  for(int i=0; i<mlassc;i++){ 
    if( fabs(par[i][4])>pmin){ 
      for(int k=0;k<6;k++){ 
	par[mleft][k]=par[i][k];
      }
      mleft++;
    }
  }
  mlassc=mleft;
  if(mleft==0) {
    //cout<<"No tracks with momentum >0.5 GeV "<<endl;
    return 0;
  }
  //-----------------------------------------------------------------------------
  //    Select pair of tracks with maximal separation at first detectorplane  
  //-----------------------------------------------------------------------------                                                            
  jret=0;
  
  int w=1;
  while(w){
    distmax=-99999.0;
    vtxbest[0]=vtxbest[1]=vtxbest[2]=-99999.0;
    iamax=1;
    iasec=2;
    //--------------------------------------------------------------------
    //        Case with flagged muon
    //--------------------------------------------------------------------    
    if(jmu>0){ 
      //cout<<"\n"<<"I if sats med muon  "<<flush; 
      iamax=jmu-1;  //the muon track is always used if flagged, choose iasec 
      //with greatest distance to iamax
      for(int j=0;j<mlassc;j++){
	int skip=0;
	if(j==(jmu-1))
	  continue;
	if(nveto>0) {
	  for(int n=0;n<nveto;n++) {
	    if((100*iamax+j)==(iveto[n])) { 
	      skip=1;
	      break;
	    }
	  }
	}
	if(skip)
	  continue;
	dist=pow((par[iamax][0]-par[j][0]),2)+pow((par[iamax][1]-par[j][1]),2);
	if(dist>distmax){
	  iasec=j;
	  distmax=dist;
	}
      }
    }
    // if no muon choose two tracks iamax and iasec furthest apart
    else {
      //cout<<"\n"<<"I else satsen ingen muon "<<flush;
      for(int jj=0;jj<(mlassc-1);jj++){ 
	int skip;
	for (int k=jj+1;k<mlassc;k++){ 
	  skip=0;
	  if(nveto>0){
	    for(int l=0;l<nveto;l++){
	      if((100*jj+k)==iveto[l]){
		skip=1;
		break;
	      }
	    }
	  }
	  if(skip) continue;  
	  
	  dist= pow((par[jj][0]-par[k][0]),2) + pow((par[jj][1]-par[k][1]),2);
	  //cout<<"\n"<<"Dist :"<<dist<<flush;        
	  
	  if(dist>distmax){ 
	    iamax=jj;
	    iasec=k;
	    distmax=dist;
	    
	  }
	  //cout<<"\n"<<"Distmax "<<distmax<<flush;    
	}
      }    
    }
    //-----------Quit if no pair of canditates found-----------------------
    //cout<<"\n"<<"Efter else sats "<<" Distmax "<<distmax<<flush;
    if(distmax<=0){
      jret=-1;
      w=0;
      return 0;
    }
    distmax= sqrt(distmax);
    //---------------------------------------------------------------------------
    //   Selection done.....Now find closest approach by linear extrapolation
    //-------------------------------------------------------------------------- 
    
    quickvtx(par[iamax],par[iasec],zapr,xapr,yapr,disapr);
    //cout<<"\n"<<"Efter quickvtx functionen "<<flush;
    //-----------------------------------------------------------------------
    //        Reject if intersection downstream of zplan1 
    //        Reject if tracks don't intersect
    //-----------------------------------------------------------------------
    
    //cout<<"\n"<<" 2tr-quick :"<<iamax<<' '<<iasec<<' '<<flush; 
    //cout<<"Zapr :"<<zapr<<' '<<" Disapr :"<<disapr<<"\n"<<flush;
    //cout<<"zplan1 :"<<par[iamax][5]<<' '<<"tol1 "<<tol[0]<<"\n"<<flush;  
    if((zapr>par[iamax][5])||(zapr>par[iasec][5])||(disapr>tol[0])){  
      //cout<<"\n"<<"I if zapr >zplan1 eller disapr>tol 0"<<flush;
      nveto++;
      iveto[nveto-1]=100*iamax + iasec;
      continue; //start from beginning of the while loop
    }
    float vect1[6],vect2[6]; 
    for(int lop=0;lop<6;lop++) {
      vect1[lop]=par[iamax][lop];
      vect2[lop]=par[iasec][lop];}
    //-----------------------------------------------------------------------
    //    Propagate tracks in magneticfield form zplan1 to zapr
    //----------------------------------------------------------------------
    
    trackinfield(vect1,zapr,vec[0]);
    trackinfield(vect2,zapr,vec[1]);
    //-----------------------------------------------------------------------
    //    Redo closest approach with parameters at zapr
    //-----------------------------------------------------------------------
    
    quickvtx(vec[0],vec[1],zapp,xapp,yapp,disapp);
    
    //-----------------------------------------------------------------------
    //    Reject if disapp> tolerance 2
    //-----------------------------------------------------------------------
    
    if(disapp>tol[1]){  
      nveto++;
      iveto[nveto-1]=100*iamax +iasec;
      //cout<<"\n"<<"Borja om efter att disapp>tol 2"<<"\n"<<flush; 
      continue;
    }
    
    //-----------------------------------------------------------------------
    //    Check if X,Y position is compatible with beam track
    //-----------------------------------------------------------------------
    vecbeam[0]=prbeam[0]+ prbeam[2]*(zapp-prbeam[5]);
    vecbeam[1]=prbeam[1]+ prbeam[3]*(zapp-prbeam[5]);
    vecbeam[2]=prbeam[2];
    vecbeam[3]=prbeam[3];
    vecbeam[4]=prbeam[4];
    vecbeam[5]=zapp;
    drad=sqrt((vecbeam[0]-xapp)*(vecbeam[0]-xapp) +
	      (vecbeam[1]-yapp)*(vecbeam[1]-yapp));
    //cout<<"\n"<<"Position comp with beam"<<" DRAD; "<<drad<<"\n"<<flush;
    if(drad>tol[2]){
      nveto++;
      iveto[nveto-1] = 100*iamax + iasec;
      //cout<<"\n"<<"Borja om efter drad> tol 3 "<<"\n"<<flush;
      continue; // start from beginning of while loop
    }
    jret=1;
    w=0;
    break;
  }  
  w=0;
  //---------------------------------------------------------------------------
  //    Take care of other tracks (if any)
  //---------------------------------------------------------------------------
  
  //--------first beam in list of tracks-------------------------------
  list[0]=-1;
  for(int m=0;m<mlassc;m++){
    list[m+1]=m;
  }
  mleft=mlassc+1;
  //cout<<"\n"<<"Tracks left :"<<mleft<<flush;
  //-----------Default values-----------------------------------------
  
  vtxbest[0]=xapp;
  vtxbest[1]=yapp;
  vtxbest[2]=zapp;
  for(int kk=0;kk<3;kk++)
    vtx2[kk]=vtxbest[kk];
  
  //-----------------------------------------------------------------------
  //     This is the end if there are no other tracks
  //     If jmu>0 and moun NOT included  error flag is set to -2
  //-----------------------------------------------------------------------
  
  if(mleft<=2){
    if((jmu==1)&&(list[0]!=jmu)&&(list[1]!=jmu))
      jret= -2;
    return 0; 
  }
  
  //------------------------------------------------------------------------
  //   Get all track parameters at Z=zapp and store in vtk
  //------------------------------------------------------------------------
  
  //-----------first copy beam--------------------------------------
  
  for(int b=0;b<6;b++)
    vtk[0][b]=vecbeam[b];
  
  //--------now tracks----------------------------------------------
  for(int t=0;t<mleft-1;t++)
    trackinfield(par[t],zapp,vtk[t+1]);

  //for(int bb=0;bb<mleft;bb++){
  //  for(int tt=0;tt<6;tt++)
  //    cout<<vtk[bb][tt]<<' '<<flush;
  //  cout<<"\n"<<flush;
  //}
  //-----------------------------------------------------------------------
  //   Keep mml tracks (out of mleft) for further work
  //   Selection based on dist. of appr. of present track vs average
  //   dist. of approach
  //-----------------------------------------------------------------------

  //float toldis=0.75*disapp;
  float toldis=c_[3]*disapp;

  //disapp is distance between tracks and toldis is distance to central
  //point i.e s/r filtervtx will accept tracks within a distance wich is
  //less than 1.5*distance found for initial pair
  
  jfilerr=0;
  
  filtervtx(mleft,mml,list,toldis,vtk,jfilerr);
  
  //------------check if beam is still in the list-------------------------
  int binc=0;
  for(int in=0;in<mml;in++){
    if(list[in]==-1) binc=1;
  }
  if(1!=binc){
    jret=-3;
    return 0;
  }
  
  //----if jmu>0, check that muon has not been dropped from list-------
  //cout<<"\n"<<"JMU "<<jmu<<flush;
  if(jmu>0){
    int muinc=0;
    
    for(int jm=0;jm<mml;jm++)
      if(list[jm]==(jmu-1)) muinc=1;
    
    if(1!=muinc){
      jret= -4;
      return 0;       
    }
  }

  
  //-------------------------------------------------------------------------
  //   Quit if jfilerr/=0 ie no additional tracks has been accepted
  //
  //-------------------------------------------------------------------------
  //cout<<"\n"<<"Jfilerr :"<<jfilerr<<" (!=0 means no add tr)"<<flush;
  if(jfilerr!=0)
    return 0;
  
  //-------------------------------------------------------------------------
  //    Now find best vertex for the mml remaining tracks using track
  //    parameters vtk given at Z=zapp, same procedure as in quickvtx 
  //    but written in matrix form for mml tracks.
  //    distfin is the average track distance to the vertex
  //-------------------------------------------------------------------------
  ntrvtx(mml,vtk,zapp,vtxbest,distfin);
  jret=2;
  return 0;
  delete [] vtk;
};



int CsRolandPattern::trackinfield(float* parin,float zend,float* parout)
{
  CsField* M = CsGeom::Instance()->getCsField();
  
  //------------------------------------------------------------------------
  //        Local variables
  //------------------------------------------------------------------------
  //cout<<"I trackinfieldfunc"<<endl<<flush;
  float step=5.0;
  float xx[3],h[3]; 
  float xp,yp,facr,ff1,ff2,x2,y2,dz=0.0;
  int jfin;
  xx[0] = parin[0];
  xx[1] = parin[1];
  xx[2] = parin[5];
  
  xp = parin[2];
  yp = parin[3];
  for(int s=0;s<5;s++)
    jfin = 0;
  
  while(jfin==0){
    int flag=0;
    h[0]=0.0;
    h[1]=0.0;
    h[2]=0.0;
    M->getField(xx[0],xx[1],xx[2],h[0],h[1],h[2]);
    //cout<<"Position : "<<xx[0]<<' '<<xx[1]<<' '<<xx[2]<<' ';
    //cout<<"faltet : "<<h[0]<<' '<<h[1]<<' '<<h[2]<<"\n";
    facr = sqrt(1.0 + xp*xp + yp*yp);
    ff1 = h[2]*yp + h[0]*xp*yp - h[1]*(1.0 + xp*xp);
    
    ff2 = -h[2]*xp - h[1]*xp*yp + h[0]*(1.0 + yp*yp);
    
    x2 = 0.3*facr*ff1/parin[4];
    y2 = 0.3*facr*ff2/parin[4];
    //cout<<"X2,Y2 "<<x2<<' '<<y2<<"\n";
    if((xx[2]-step)<zend){
      dz = zend - xx[2];
      xx[2] = zend;
      flag=1;
      
    }else{
      dz = -step;
      xx[2] -= step;
    }
    
    xp +=(dz*x2/1000.0); //Division by 1000 for units for x,y,z in mm
    yp +=(dz*y2/1000.0);
    
    xx[0] += (xp*dz);
    xx[1] += (yp*dz);
    //cout<<"X,Y;"<<xx[0]<<' '<<xx[1]<<"\n";
    jfin=flag;
  }
  
  parout[0] = xx[0];
  parout[1] = xx[1];
  parout[2] = xp;
  parout[3] = yp;
  parout[4] = parin[4];
  parout[5] = xx[2];
  
  return 0;
};

int CsRolandPattern::trackbeam(float* parin,float zend,float* parout)
{
  CsField* M = CsGeom::Instance()->getCsField();
  
  //------------------------------------------------------------------------
  //        Local variables
  //------------------------------------------------------------------------
  float step=5.0;
  float xx[3],h[3]; 
  float xp,yp,facr,ff1,ff2,x2,y2,dz=0.0;
  int jfin;
  //cout<<"\n"<<"Zstart :"<<parin[5];
  //cout<<"\n"<<"ZEND :"<<zend<<"\n";
  xx[0] = parin[0];
  xx[1] = parin[1];
  xx[2] = parin[5];
  
  xp = parin[2];
  yp = parin[3];
  //for(int s=0;s<5;s++)
  //  cout<<parin[s]<<' '<<xx[s]<<' ';
  //cout<<"\n";
  jfin = 0;
  
  while(jfin==0){
    int flag=0;
    //cout<<"\n"<<"Fore magfield "<<xx[0]<<"\n"<<flush;   
    // Call magnetic field in point xx
    h[0]=0.0;
    h[1]=0.0;
    h[2]=0.0;
    
    M->getField(xx[0],xx[1],xx[2],h[0],h[1],h[2]);
    //cout<<"Position : "<<xx[0]<<' '<<xx[1]<<' '<<xx[2]<<' ';
    //cout<<"faltet : "<<h[0]<<' '<<h[1]<<' '<<h[2]<<"\n";
     
    facr = sqrt(1.0 + xp*xp + yp*yp);
    ff1 = h[2]*yp + h[0]*xp*yp - h[1]*(1.0 + xp*xp);
    ff2 = -h[2]*xp - h[1]*xp*yp + h[0]*(1.0 + yp*yp);
    
    x2 = 0.3*facr*ff1/parin[4];
    y2 = 0.3*facr*ff2/parin[4];
    //cout<<"X2,Y2 "<<x2<<' '<<y2<<"\n";
    if((xx[2]+step)>zend){
      dz = zend - xx[2];
      xx[2] = zend;
      flag=1;
    }else{
      dz = step;
      xx[2] += step;
    }
    
    xp +=(dz*x2/1000.0); //Division by 1000 for units for x,y,z in mm
    yp +=(dz*y2/1000.0);
    //cout<<"Slopes : dx/dz = "<<xp<<' '<<"dy/dz = "<<yp<<endl;
    xx[0] += (xp*dz);
    xx[1] += (yp*dz);
    //cout<<"X,Y;"<<xx[0]<<' '<<xx[1]<<"\n";
    jfin=flag;
  }
  
  parout[0] = xx[0];
  parout[1] = xx[1];
  parout[2] = xp;
  parout[3] = yp;
  parout[4] = parin[4];
  parout[5] = xx[2];
  //for(int q=0;q<6;q++)
  //  cout<<parout[q]<<' ';
  //cout<<"\n";
  
  return 0;
};



int CsRolandPattern::trackinfield2(float* parin1,float *parin2,float & zend,float* parout1,float *parout2,float & dist)
{
  CsField* M = CsGeom::Instance()->getCsField();

  //------------------------------------------------------------------------
  //        Local variables
  //------------------------------------------------------------------------
  dist=9999999.0;
  float step=5.0,zstop=-1000.0,dismin;
  float xx1[3],xx2[3],h[3]; 
  float xp1,xp2,yp1,yp2,facr,ff1,ff2,x1,x2,y1,y2,dz=0.0;
  int jfin;
  //cout<<"\n"<<"Zstart :"<<zstart;
  //cout<<"\n"<<"ZEND :"<<zend<<"\n";
  //cout<<"Parin : ";
  xx1[0] = parin1[0];
  xx1[1] = parin1[1];
  xx1[2] = parin1[5];
  
  xp1 = parin1[2];
  yp1 = parin1[3];
  
  xx2[0] = parin2[0];
  xx2[1] = parin2[1];
  xx2[2] = parin2[5];
  
  xp2 = parin2[2];
  yp2 = parin2[3];
  
  for(int s=0;s<5;s++)
    jfin = 0;
  while(jfin==0){
    int flag=0;
    M->getField(xx1[0],xx1[1],xx1[2],h[0],h[1],h[2]);
    //cout<<"Position : "<<xx[0]<<' '<<xx[1]<<' '<<xx[2]<<' ';
    //cout<<"faltet : "<<h[0]<<' '<<h[1]<<' '<<h[2]<<"\n";
    facr = sqrt(1.0 + xp1*xp1 + yp1*yp1);

    ff1 = h[2]*yp1 + h[0]*xp1*yp1 - h[1]*(1.0 + xp1*xp1);
    ff2 = -h[2]*xp1 - h[1]*xp1*yp1 + h[0]*(1.0 + yp1*yp1);
    
    x1 = 0.3*facr*ff1/parin1[4];
    y1 = 0.3*facr*ff2/parin1[4];
    //cout<<"X2,Y2 "<<x2<<' '<<y2<<"\n";
    
    // Call magnetic field in point xx for track 2
    M->getField(xx2[0],xx2[1],xx2[2],h[0],h[1],h[2]);
    facr = sqrt(1.0 + xp2*xp2 + yp2*yp2);

    ff1 = h[2]*yp2 + h[0]*xp2*yp2 - h[1]*(1.0 + xp2*xp2);
    ff2 = -h[2]*xp2 - h[1]*xp2*yp2 + h[0]*(1.0 + yp2*yp2);
    
    x2 = 0.3*facr*ff1/parin2[4];
    y2 = 0.3*facr*ff2/parin2[4];
    //cout<<"X2,Y2 "<<x2<<' '<<y2<<"\n";
    
    if((xx1[2]-step)<zstop){
      dz = zstop - xx1[2];
      xx1[2] = zstop;
      xx2[2] = zstop;
      flag=1;
    }else{
      dz = -step;
      xx1[2] -= step;
      xx2[2] -= step;
      //cout<<"zpos: "<<xx1[2]<<' '<<xx2[2]<<endl;
    }
    
    xp1 +=(dz*x1/1000.0); //Division by 1000 for units for x,y,z in mm
    yp1 +=(dz*y1/1000.0);
    
    xx1[0] += (xp1*dz);
    xx1[1] += (yp1*dz);
    //cout<<"X,Y;"<<xx[0]<<' '<<xx[1]<<"\n";
    
    xp2 +=(dz*x2/1000.0); //Division by 1000 for units for x,y,z in mm
    yp2 +=(dz*y2/1000.0);
    
    xx2[0] += (xp2*dz);
    xx2[1] += (yp2*dz);
    //cout<<"X,Y;"<<xx[0]<<' '<<xx[1]<<"\n";
    dismin= (xx2[0]-xx1[0])*(xx2[0]-xx1[0])+(xx2[1]-xx1[1])*(xx2[1]-xx1[1]);
    //cout<<"trackdist: "<<dismin<<" zpos: "<<xx1[2]<<endl;
    
    if(dismin<=dist){
      dist=dismin;
      zend=xx1[2];
      //cout<<"trackdist: "<<dist<<" zpos: "<<xx1[2]<<endl;
    }else{
      dist=dismin;
      zend=xx1[2];
      flag=1;
    }
    jfin=flag;
  }
  dist=sqrt(dist);
  parout1[0] = xx1[0];
  parout1[1] = xx1[1];
  parout1[2] = xp1;
  parout1[3] = yp1;
  parout1[4] = parin1[4];
  parout1[5] = xx1[2];
  
  parout2[0] = xx2[0];
  parout2[1] = xx2[1];
  parout2[2] = xp2;
  parout2[3] = yp2;
  parout2[4] = parin2[4];
  parout2[5] = xx2[2];
  
  return 0;
};


int CsRolandPattern::ntrvtx(int nrtr,float **vtk,float zin, float* vtxout,float &disfin)
{
  //----------------------------------------------------------------
  //  Local variables
  //----------------------------------------------------------------

  const int mxll=50;
  int ndim=20,nleft;
  int *rr;
  rr=new int[20];
  float **amat;
  amat=new float *[20];
  for(int p=0;p<20;p++){
    amat[p]=new float[20];}
  float *a;
  a=new float[400];
  float *bt;
  bt=new float[20];
  int ifail,k=1;
  float xr,yr,zr,xaver,yaver,zaver,diskr,disaver;
  //---------------------------------------------------------------
  //cout<<"\n"<<"I ntrvtx functionen "<<flush;
  nleft=nrtr;
  
  if(nleft>20)
    nleft=20;
  
  for(int i=0;i<nleft;i++){
    
    amat[i][i]= (vtk[i][2]*vtk[i][2] + vtk[i][3]*vtk[i][3] + 1.0)*(nleft-1);
    //cout<<"\n"<<"AMAT: "<<amat[i][i]<<flush;
    for(int j=i+1;j<nleft; j++){
      amat[j][i]= -vtk[i][2]*vtk[j][2] - vtk[i][3]*vtk[j][3] -1.0;
    }
    bt[i] = (vtk[i][2]*vtk[i][0] + vtk[i][3]*vtk[i][1])*(1-nleft);     
    
    for(int l=0;l<nleft;l++){
      if(l==i) 
	continue;
      bt[i] += vtk[l][0]*vtk[i][2] +vtk[l][1]*vtk[i][3];   
    }
  }
  //---------------------------------------------------------------
  // Symmetrize amat (matrix)
  //---------------------------------------------------------------
  for(int k=1;k<nleft;k++){
    for(int p=0;p<=k-1;p++){
      amat[p][k] = amat[k][p];
    }  
  }
  
  //-----------solve linear system-----------------------------
  for(int ii=0;ii<20;ii++){
    for(int jj=0;jj<20;jj++)
      a[(ii*20+jj)]=amat[ii][jj];}
  reqn_(nleft,a,ndim,rr,ifail,k,bt);
  xaver = 0.0;
  yaver = 0.0;
  zaver = 0.0;
  
  for(int n=0; n<nleft; n++){
    zr = zin + bt[n];
    xr = vtk[n][0] + vtk[n][2]*(zr - zin);
    yr = vtk[n][1] + vtk[n][3]*(zr - zin);
    
    xaver += xr;
    yaver += yr;
    zaver += zr;  
  }
  
  xaver /=nleft;
  yaver /=nleft;
  zaver /=nleft;
  vtxout[2] = zaver;
  vtxout[1] = yaver;
  vtxout[0] = xaver;
  disaver = 0;
  for(int m=0; m<nleft; m++){
    zr = zin + bt[m];
    xr = vtk[m][0] + vtk[m][2]*(zr - zin);
    yr = vtk[m][1] + vtk[m][3]*(zr - zin);
    diskr = (xaver-xr)*(xaver-xr) + (yaver-yr)*(yaver-yr) + (zaver-zr)*(zaver-zr);
    diskr = sqrt(diskr);
    disaver += diskr;
  } 
  
  disaver /=nleft;
  
  delete [] amat;
  return 0;
};



int CsRolandPattern::filtervtx(int nin, int &nout, int* list,
			      float toldis,float **vtk, int &jflag)
{

  //-----------------------------------------------------------------------   
  //   Local variables
  //-----------------------------------------------------------------------

  const int mxll=50;
  float distr[mxll];
  float xc,yc,xcen,ycen,distmax,max=0.0;
  jflag=0;
  int nleft=nin;
  //cout<<"\n"<<"Toldis "<<toldis<<flush;
  
  //--------initial central position-----------------------------------
  
  int w=1;
  while(w){
     
    xcen=0.0;
    ycen=0.0;
    
    for(int i=0; i<nleft;i++){
      xcen+=(vtk[i][0]/nleft);
      ycen+=(vtk[i][1]/nleft);
    }
    //cout <<"\n"<<"XCEN,YCEN :"<<xcen<<ycen<<flush;    

    //------distance to center for each track-----------------------------
    
    for(int j=0;j<nleft;j++){
      distr[j]=(vtk[j][0]-xcen)*(vtk[j][0]-xcen)+(vtk[j][1]-ycen)*(vtk[j][1]-ycen);
      //cout<<distr[j]<<' '<<flush;
    }
    
    //---sort according to distance and compare largest dist. with average
    
    max=distr[0];
    for(int k=1;k<nleft;k++){
      if(max<distr[k])
	max=distr[k];
      //cout<<"\n"<<"max :"<<max<<' '<<flush;
    }
    
    distmax=sqrt(max);
    //--------------------------------------------------------------------
    //  Normal end: All tracks accepted
    //--------------------------------------------------------------------
   
    if(distmax<=toldis) w=0;

    //--------------------------------------------------------------------
    //  Quit if only two tracks are left
    //--------------------------------------------------------------------
  
    if(nleft==3){
      //cout<<"\n"<<"Only two tracksleft "<<flush;
      nout=2;
      jflag=1;
      w=0; 
      return 0;
    }
    
    //--------------------------------------------------------------------
    //  Drop track with largest distance and recalculate
    //--------------------------------------------------------------------
    
    int dt=0;
    for(int d=0;d<nleft;d++){
      if(distr[d]==max)
	continue;
      dt++;
      list[dt-1]=list[d];
      distr[dt-1]=distr[d];
      for(int i=0;i<5;i++)
	vtk[dt-1][i]=vtk[d][i];   
    }
    nleft--;
    //cout<<"\n"<<"Nleft :"<<nleft<<flush;
  }
  nout=nleft;
  return 0; 
};

int CsRolandPattern::quickvtx(float* vec1,float* vec2,float &zapp,float & xapp, float &yapp,float &disapp)
{
  float tk=0.0,tl=0.0,c1=0.0,c2=0.0,rho1=0.0,rho2=0.0,xvt1=0.0,xvt2=0.0;
  float tm=0.0,x1=0.0,x2=0.0,y1=0.0,y2=0.0,dist=0.0,zstart1=0.0,zstart2=0.0;
  
  disapp=0.0;     
  zstart1=vec1[5];
  zstart2=vec2[5];
  tk=vec1[2]*vec1[2]+vec1[3]*vec1[3] +1.0;
  tl=vec1[2]*vec2[2] + vec1[3]*vec2[3] + 1.0;
  tm=vec2[2]*vec2[2] + vec2[3]*vec2[3] + 1.0;
  
  c1= vec1[2]*(vec1[0]-vec2[0]) + vec1[3]*(vec1[1]-vec2[1]);
  c2= vec2[2]*(vec1[0]-vec2[0]) + vec2[3]*(vec1[1]-vec2[1]);
  c1= -c1;
  c2= -c2;
  
  rho1 = (c1*tm - c2*tl)/(tk*tm - tl*tl);
  rho2 = (c1*tl - c2*tk)/(tk*tm - tl*tl);
  
  xvt1 = rho1 + zstart1;
  xvt2 = rho2 + zstart2;
  
  zapp = 0.5*(xvt1 + xvt2);
  x1 = vec1[0] + vec1[2]*(xvt1 - zstart1);
  y1 = vec1[1] + vec1[3]*(xvt1 - zstart1);

  x2 = vec2[0] + vec2[2]*(xvt2 - zstart2);
  y2 = vec2[1] + vec2[3]*(xvt2 - zstart2);
  
  xapp = 0.5*(x1 + x2);
  yapp = 0.5*(y1 + y2);
  dist =(((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2)) +((xvt1-xvt2)*(xvt1-xvt2)));
  
  //cout<<"Dist i quick"<<dist<<endl<<flush;
  disapp = sqrt(dist);
  return 0;
  
};

