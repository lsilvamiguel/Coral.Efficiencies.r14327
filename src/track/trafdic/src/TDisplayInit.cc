// $Id: TDisplayInit.cc 13647 2013-02-01 09:16:15Z ybedfer $

#include <iostream>
#include "TDisplay.h"
#include "TOpt.h"
#include "TConstants.h"
#include <stdio.h>
#include "higz.h"

using namespace std;

/*!
  Private method for
  initialization of all datamembers of the class TDisplay. 
  Called by TDisplay constructor.
*/

bool TDisplay::Init(){

  // Default menu position and width
  Xmenu[0]=0.01;   Xmenu[1]=1.-0.01;        // XY of right or left (if X > 0.5) upper corner
  Wmenu[0]=0.006;  Wmenu[1]=0.02;           // Width/character, width/line

  ifstream DrawOpt("drawing.options");
  if(!DrawOpt){ // file do not exists. Use defaults
    Proj=0; // default projection

    // Default view field from option file
    Xvf[0]=TOpt::DefView[0]; Xvf[1]=TOpt::DefView[1];
    Yvf[0]=TOpt::DefView[2]; Yvf[1]=TOpt::DefView[3];

    for(int i=0; i < 20; i++) DrOpt[i]=-1; // reset

    // Default mode flags
    DrOpt[0]=-1;            // MC        hits     toggle
    DrOpt[1]= 1;            // Real hits (digits) toggle
    DrOpt[2]=-1;            // PS file status toggle (-1 = closed)
    DrOpt[3]= 0;            // PS file number
    DrOpt[4]=-1;            // Draw/Not draw planes, not in current projections
    DrOpt[5]= 1;            // Draw tracks like joined hits or trajectory
    DrOpt[6]=-1;            // short/long track info printout
    DrOpt[7]=-1;            // detailed/simple detector drawing
    if(TOpt::Graph[3] == 0) DrOpt[8]=-1; // field OFF
    else                    DrOpt[8]= 1; // fieled ON

    DrOpt[9] =-1;           // MC track drawing toggle
    DrOpt[10]= 1;           // Reconstructed track drawing toggle
    DrOpt[11]= 1;           // Filled graphycs / hollow graphycs toggle
    DrOpt[12]=-1;           // Material map drawing OFF
    DrOpt[13]=-1;           // Extrapolation to MC vtx OFF
    DrOpt[14]=-1;           // Draw text tags
    DrOpt[15]= 0;           // Image distorsion mode. 0 means "no distortion"
      
  } else {
    DrawOpt>>Proj;
    DrawOpt>>Xvf[0]>>Xvf[1]>>Yvf[0]>>Yvf[1];
    for(int i = 0; i < 20; i++) DrawOpt>>DrOpt[i];
    string dummy; DrawOpt>>dummy;
    if(!DrawOpt.eof()){
      cout<<"TDisplay::Init() ==> ./drawing.options file looks corrupted"<<endl
	  <<"Remove it and start again."<<endl;
      assert(false);
    }
    DrawOpt.close();
    cout<<"TDisplay::Init() ==> previously saved drawing options will be used"<<endl;
  }
  // Display modes, depended from reconstruction modes
  if(TOpt::ReMode[5] != 0) DrOpt[5]=-1;

  ev_count=1;      // init event counter

  // Write higz_windows.dat file in working directory
  ofstream HigzWindow;
  HigzWindow.open("higz_windows.dat",ios::out|ios::trunc);
  if(!HigzWindow) {
    cout<<"TDisplay::Init() ==> Cannot create \"higz_windows.dat\" file.\n"
	<<"Check your write permission for current directory"<<endl;
    assert(false);
  }

  HigzWindow<<"   -1    1 1350 1000"<<endl; // 1 big
  HigzWindow<<"   -1    1 1100 0800"<<endl; // 2 medium
  HigzWindow<<"   -1    1 0850 0603"<<endl; // 3 small
  HigzWindow<<"   -1    1 0950 1000"<<endl; // 4 different aspect ratio 
  HigzWindow<<"   -1    1 0800 1000"<<endl; // 5 closing line


  HigzWindow.close();
    
  // Init HIGZ
  int ityp=1, konid=1, kwkid=1;
  IGINIT(0);
  if(TOpt::Graph[7] <= 0 ) { // NOT a "movie" mode
    cout<<endl;
    cout<<"     As HIGZ graphic window can not be resized,"<<endl;
    cout<<"specify workstation type, corresponding to desired window size"<<endl<<endl;
    cout<<" Workstation type           Window size (pixels)\n";
    cout<<"           1                  1350 x1000 (large)\n";
    cout<<"           2                  1100 x 800 (medium)\n";
    cout<<"           3                   850 x 600 (small)\n";
    cout<<"           4                   950 x1000 (widescreen)\n";
    IGWKTY(ityp);  // request workstation type
  }
  if(ityp < 1 || ityp > 4){
    cout<<"Workstation type is out of range. Event display is OFF"<<endl<<endl;
    TOpt::Graph[0] = 0;
    return(true);
  }

  IGSSE(6,ityp); // open workstation
  ISELNT(0);
  IGSET("CHHE",0.03);
  IGSET("TXAL",23);
  IGSET("TXFP",-132);
  char str[50];
  int ic[]={1,0,4,2};
  for(int i = 0; i < 4; i++){
    ISTXCI(ic[i]);
    float a = 0.0015;
    ITX(0.5+a*i, 0.50+a*i, "Track Finding and Fit in COMPASS");
    ITX(0.5+a*i, 0.40+a*i, "( TRAFFIC )");
    sprintf(str,"Version %4.2f",TConstants_TraFFiC_Version);
    ITX(0.5+a*i, 0.30+a*i,str);
    IUWK(0,1);
  }

  // define new colors
  ISCR(1,101, 1.0,  0.94,  0.87); // for background
    
  ISCR(1,102, 0.7,  0.9, 1.0 );   // for RICH
  
  ISCR(1,103, 0.8,  1.0 ,  0.8  ); // for PT
  ISCR(1,104, 1.00, 0.93 , 1.00 ); // for magnets
  ISCR(1,105, 0.5,  0.5 ,  0.5  ); // grey for outlines
  ISCR(1,106, 0.4,  0.4 ,  0.44 ); // for targets
  
  ISCR(1,107, 0.0,   0.0 , 0.4 );  // for the field lines
  
  ISCR(1,108, 0.8,  0.4 ,  0. );   // for the other-proj-planes
  ISCR(1,109, 0.,   0.4 ,  0.  );  // muon walls
  ISCR(1,110, 0.,   0.7 ,  0.7 );  // for MC hits
  
  ISCR(1,111, 0.6,   0. ,  0.0 );  // for tracks not associated with MC track
  ISCR(1,112, 1.0,   0.45, 0.15);  // for "golden" tracks
  
  /*
  IOPKS(6);
  IOPWK(kwkid,konid,ityp);
  IACWK(kwkid);
  */

  return(true);
 
}




