// $Id: TDisplayDraw.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include "CsGeom.h"
#include "CsEvent.h"
#include "TDisplay.h"
#include "TAlgo.h"
#include "TSetup.h"
#include "TOpt.h"
#include "THlx.h"
#include "TEv.h"
#include "TConstants.h"
#include "higz.h"

using namespace std;

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TDisplayDrawText":
   i) Add "Write Back" feature.
*/

void Wait(float);

/*!
  Method for
  drawing of different parts of the event
//
// control flag key == 0  - Draw everything
// control flag key == 1  - Draw static objects
// control flag key == 2  - Draw hits
// control flag key == 3  - Draw tracks and write event info
// control flag key == 4  - Draw menu and set key = 0;
//
*/

int TDisplay::Draw(int key)
{ 
  
  if(TOpt::Graph[0] <= 0) return(0); // do nothing if the graphics is OFF

  const TEv& ev       = TEv::Ref();
  const TSetup& setup = TSetup::Ref();

  int im=0, nt=0, istat=0, iret=0;
  float x1,y1,x2,y2,dx,dy;
  int ip=-1, ih=-1, id = -1;


 draw:
  
  if(key == 0 || key == 1){

    // Clear worksation
    Clear();
    
    // Draw target, RICH, muon walls etc.
    DrawMisc();
    //Draw material map (if requesred). Draw planes
    DrawMaterialMap(); 
    //Draw detectors
    DrawDet();
    //Draw magnets (or mag. field)
    DrawMag(); 
  }

  if(key == 0 || key == 2){

    // Draw hits
    
    if(DrOpt[0] > 0) DrawHitsMC(); // draw MC  hits
    if(DrOpt[1] > 0) DrawHits();   // draw clusters 

  }    

  if(key == 0 || key == 3){

    // Draw reconstructed tracks
    if(DrOpt[9]  > 0) DrawMCTracks();
    if(DrOpt[10] > 0) DrawTracks();
  }

 menu:
  
  if(key == 0 || key == 4){

    key = 0;

    // Call to user drawing 
    UserDraw();

    // Write some general event info and draw the frame 
    DrawText();
    DrawFrame();
    
    // Draw menu and get user's choice
    
    if(TOpt::Graph[7] > 0){ // "movie" mode
      Wait(float(TOpt::Graph[7])); 
      im = 15; // next event
    } else { 
      im=DrawMenu();
    }

    // process menu eveny
    switch(im) 
      {
      case 1:   // Zoom
	stXvf0.push(Xvf[0]); stXvf1.push(Xvf[1]); // save current view filed
	stYvf0.push(Yvf[0]); stYvf1.push(Yvf[1]); // into stacks
      zoom_again:
	IGLOC2(1,nt,x1,y1,x2,y2,istat," ");
	if(istat==0) goto menu;
	if((x1 > 0 && x1 < 1.0 && y1 > 0 && y1 < 1.0) ||
	   (x2 > 0 && x2 < 1.0 && y2 > 0 && y2 < 1.0)) {
	  cout<<"Your rubber box must be inside the view field. Zoom again"<<endl;
	  goto zoom_again; // rubber box goes out of view field  
	}
	if(x1 == x2 || y1 == y2) {
	  cout<<"One side of rubber box == 0. Zoom again"<<endl;
	  goto zoom_again;   
	}
	// set new view field
	Xvf[0]=(x1 < x2 ? x1 : x2); 
	Xvf[1]=(x1 < x2 ? x2 : x1); 
	Yvf[0]=(y1 < y2 ? y1 : y2); 
	Yvf[1]=(y1 < y2 ? y2 : y1);
	goto draw;
      case 2:   // undo zoom
	if(stXvf0.empty()) {
	  cout<<"'Unzoom' steak is empty"<<endl;
	  goto menu;
	} else {
	  Xvf[0]=stXvf0.top(); Xvf[1]=stXvf1.top();
	  Yvf[0]=stYvf0.top(); Yvf[1]=stYvf1.top();
	  stXvf0.pop(); stXvf1.pop();
	  stYvf0.pop(); stYvf1.pop();
	  goto draw;
	}
      case 3:  // zoom out
	stXvf0.push(Xvf[0]); stXvf1.push(Xvf[1]); // save current view filed
	stYvf0.push(Yvf[0]); stYvf1.push(Yvf[1]); // into stacks
	double dx,dy;
	dx=0.5*(Xvf[1]-Xvf[0]);
	dy=0.5*(Yvf[1]-Yvf[0]);
	Xvf[0]-=dx; Xvf[1]+=dx;
	Yvf[0]-=dy; Yvf[1]+=dy;
	goto draw;
      case 4:   // Move
	IRQLC(1,20,istat,nt,x1,y1);
	if(istat==0) goto menu;
	IRQLC(1,40,istat,nt,x2,y2);
	if(istat==0) goto menu;     
	Xvf[0]-=x2-x1; Xvf[1]-=x2-x1;
	Yvf[0]-=y2-y1; Yvf[1]-=y2-y1;
	goto draw;
      case 5:  // Locator
	do {
	  IRQLC(1,20,istat,nt,x1,y1);
	  double xx[3], b[3], bb;
	  xx[0] = x1;
	  if      (Proj==0) { xx[1]=y1; xx[2]=0.; }
	  else if (Proj==1) { xx[1]=0.; xx[2]=y1; }
	  else              { xx[1]=y1; xx[2]=0.; }// As if it were horiz. proj.
	  bb = TAlgo::Field(xx,b);
	  cout.fill(' ');
	  cout<<"X = "
  	      <<setprecision(7)<<setw(10)<<xx[0]
	      <<"   Y = "
	      <<setprecision(7)<<setw(10)<<xx[1]
	      <<"   Z = "
	      <<setprecision(7)<<setw(10)<<xx[2];
	  cout<<"   Bx = "
	      <<setprecision(4)<<setw(10)<<b[0]
	      <<"   By = "
	      <<setprecision(4)<<setw(10)<<b[1]
	      <<"   Bz = "
	      <<setprecision(4)<<setw(10)<<b[2]<<endl;
	  if(TOpt::ReMode[20] > 0) { // material map is ON
	    THlx H; H(0)=xx[0]; H(1)=xx[1]; H(2)=xx[2]; 
	    float radl(0), step(0);
	    CsMaterialMap *map = CsGeom::Instance()->getCsMaterialMap();
	    map->getRadLength(H,true,radl,step);
	    H(5)=0.01; // assume 100 GeV to calculate dE/dX
	    cout<<"Rad. Length = "<<radl<<" cm. 10*log(X0air/X0) = "<<10*log(30420/radl)
		<<"  dE/dX = "<<map->getdE(H, 1.)<<endl;
	  }
	  cout<<endl;
	} while(istat!=0);
	goto menu;
      case 6:  // Ruler
	do {
	  IRQLC(1,20,istat,nt,x1,y1);
	  if(istat==0) break;
	  IRQLC(1,40,istat,nt,x2,y2);
	  cout<<"Delta X = ";
	  cout<<setprecision(8)<<setw(11)<<x2-x1;
	  if(Proj==0) cout<<"   Delta Y = ";
	  if(Proj==1) cout<<"   Delya Z = ";
	  cout<<setprecision(8)<<setw(12)<<y2-y1;
	  cout<<"   Dist = "<<sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	  if((x2-x1) != 0.) cout<<"   Slope = "<<(y2-y1)/(x2-x1)<<"  ("<<180.*atan((y2-y1)/(x2-x1))/M_PI<<" deg.)";
	  cout<<endl;
	} while(istat!=0);
	goto menu;
      case 7: // Track info
	cout<<endl;
	while((id = DrawTracks(1)) != -1){
	  TrackInfo(id);
	}
	goto menu;
      case 8: // Plane info
	cout<<endl;
	while((ip = DrawDet(1)) != -1){
	  DetInfo(ip);
	}
	goto menu;
      case 9: // Hit info
	cout<<endl;
	while((ih = DrawHits(1)) != -1){
	  HitInfo(ih);
	}
	goto menu;
      case 10: // MC Hit info
	if(ev.IsMC() && DrOpt[0] > 0){
	  cout<<endl;
	  while((id = DrawHitsMC(1)) != -1){
	    ev.vHitMC(id).Print();
	  }
	}
	goto menu;
      case 11: // MC track info
	cout<<endl;
	while((id = DrawMCTracks(1)) != -1){
	  MCTrackInfo(id);
	}
	goto menu;
      case 12: // Reset view
	// clear zoom stacks
	while(!stXvf0.empty()){
	  stXvf0.pop(); stXvf1.pop();
	  stYvf0.pop(); stYvf1.pop();
	}
	Xvf[0]=TOpt::DefView[0]; Xvf[1]=TOpt::DefView[1];
	Yvf[0]=TOpt::DefView[2]; Yvf[1]=TOpt::DefView[3];
	Proj = 0;
	goto draw;
      case 13:   // Open/Close PS file
	if(DrOpt[2]==-1){
	  char num[5];
	  const char* psdir = TOpt::PSdir.c_str(); // conversion to C-like string
	  sprintf(num,"%u",DrOpt[3]);
	  char filnam[sizeof(psdir)+10];
	  strcpy(filnam,psdir);
	  strcat(filnam,"Pic"); strcat(filnam,num); strcat(filnam,".ps");
	  cout<<"Open  PS file : "<<filnam<<endl;
	  OPPS(filnam);
	  DrOpt[2]*=-1;
	}
	else{
	  CLOPS();
	  DrOpt[2]*=-1;
	  DrOpt[3]++;
	  cout<<"Close PS file"<<endl;
	}
	goto draw;
      case 14:   // Drawing modes
	DrawModeMenu();
	goto draw;
      case 15:   // Next event
	this->ev_count++;
	
	break;
      case 16: 
	goto menu;
      case 17:   // Projection toggle
	if(++Proj == int(setup.vProj().size())) Proj=0;
	goto draw;
      case 18:  // Write back
	static int iStream = -1;
	if (iStream<0) {
	  if (CsEvent::Instance()->isAMonteCarloEvent()) {
	    cout<<"\n\a** No writing back MC **\a\n\n"; goto menu;
	  }
	  string s("TraFDic.WriteBack.raw");
	  if ((iStream = CsEvent::Instance()->openRawEventOutputStream(s))<0)
	    cout<<"\n\a** Cannot open WriteBack stream \""<<s<<"\" **\a\n\n";
          else
	    cout<<"\nOpen WriteBack output file \""<<s<<"\"\n\n";
	}
	if (iStream>=0) {
	  CsEvent::Instance()->outputRawEventToStream(iStream);
	  cout<<"\nEvent "<<CsEvent::Instance()->getEventNumberInRun()<<
	    " written back\n";
	}
	goto menu;
      case 19:  //Exit
	iret = 1;
	break;
      case -1000: // Right mouse button (redraw without menu)
	goto draw;
      default:
	goto menu;
      }
  }
  return(iret);
}
