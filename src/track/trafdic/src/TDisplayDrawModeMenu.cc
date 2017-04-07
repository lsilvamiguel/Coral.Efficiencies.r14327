/*!

  Draw menu of visualisation modes
  Change modes according to user's choice

*/

#include <iostream>
#include <stdio.h>
#include "TDisplay.h"
#include "TOpt.h"
#include "higz.h"

void TDisplay::DrawModeMenu()
{
  int ich=0,im=0;
  int nbu=0;
  int n=20;
  const int len = 16;
  float x1,x2,y1,y2;

  char* chu[nbu];

  // Memu items text
  char  chi[n][len];
  char  chd[n][len];
  char  chv[n][len];

  strcpy(chi[0] , " hollow/filled "); // DrOpt[11]
  strcpy(chi[1] , "    hit/helix  "); // DrOpt[5]
  strcpy(chi[2] , "    prt. mode  "); // DrOpt[6]
  strcpy(chi[3] , "   field mode  "); // DrOpt[8]
  strcpy(chi[4] , "   material    "); // DrOpt[12]
  strcpy(chi[5] , "   MC tracks   "); // DrOpt[9]
  strcpy(chi[6] , "   Reco tracks "); // DrOpt[10]
  strcpy(chi[7] , "  xtrap to vtx "); // DrOpt[13]
  strcpy(chi[8] , "     Details   "); // DrOpt[7]
  strcpy(chi[9] , "    Other det. "); // DrOpt[4]
  strcpy(chi[10], "    Clusters   "); // DrOpt[1]
  strcpy(chi[11], "    MC hits    "); // DrOpt[0]
  strcpy(chi[12], "    Tags       "); // DrOpt[14]
  strcpy(chi[13], "               ");
  strcpy(chi[14], "     Save      ");
  strcpy(chi[15], "-     Quit     ");
  n=16;

  if(DrOpt[9] >0) chi[ 5][0]='-';
  if(DrOpt[10]>0) chi[ 6][0]='-';
  if(DrOpt[7] >0) chi[ 8][0]='-';
  if(DrOpt[4] >0) chi[ 9][0]='-';
  if(DrOpt[1] >0) chi[10][0]='-';
  if(DrOpt[0] >0) chi[11][0]='-';

  if(DrOpt[11] < 0){
    strcpy(chi[0], "    filled     ");
  } else {
    strcpy(chi[0], "    hollow     ");
  }
  if(DrOpt[5] == -1){
    strcpy(chi[1], "   Smooth      ");
  }
  if(DrOpt[5] ==  0){
    strcpy(chi[1], "   Extrap.     ");
  }
  if(DrOpt[5] == +1){
    strcpy(chi[1], "   Polyline    ");
  }

  if(DrOpt[6] > 0){
    strcpy(chi[2], "  short print  ");
  } else {
    strcpy(chi[2], "  long print   ");
  }
  if(DrOpt[8] > 0){
    strcpy(chi[3], " field OFF     ");
  } else {
    strcpy(chi[3], " field ON      ");
  }
  if(DrOpt[12] > 0){
    strcpy(chi[4], " m.map OFF     ");
  } else {
    strcpy(chi[4], " m.map ON      ");
  }
  if(DrOpt[13] > 0){
    strcpy(chi[7], " vertex OFF    ");
  } else {
    strcpy(chi[7], " vertex ON     ");
  }
  if(DrOpt[14] < 0){
    strcpy(chi[12],"  Tags ON      ");
  } else {
    strcpy(chi[12],"  Tags OFF     ");
  }


  for(int i=0; i<n; i++) strcpy(chd[i],"               ");

  ISELNT(0);  // Default NT
  if(Xmenu[0] < 0.5){
    x1=Xmenu[0]+20*Wmenu[0];
    x2=x1+len*Wmenu[0];
  } else {
    x2=Xmenu[0]-20*Wmenu[0];
    x1=x2-len*Wmenu[0];
  }
  y1=Xmenu[1]-0.05;
  y2=y1- n *Wmenu[1];
  ich=0;
  while(ich==0){
    char s1[] = "Modes", s2[] = "STCHD";
    IGMENU(im,s1,x1,x2,y1,y2,nbu,chu,n,chi,chd,chv,ich,s2);
  }

  switch(ich-1)
    {
    case 0:
      DrOpt[11]*=-1;
      break;
    case 1:
      DrOpt[5]+=1;
      if(DrOpt[5] == 2) DrOpt[5]=-1;
      break;
    case 2:
      DrOpt[6]*=-1; 
      break;
    case 3:
      if(TOpt::Graph[3] == 0 && DrOpt[8] == -1) {
	TOpt::Graph[3]=3; // field was OFF in options, but let's switch them ON
	std::cout<<"TDisplay::DrawModeMenu() ==> Redefine 'TraF Graph [3] 0' option !"<<std::endl;   
      }
      DrOpt[8]*=-1;
      break;
    case 4:
      DrOpt[12]*=-1;
      break;
    case 5:
      DrOpt[9]*=-1;
      break;
    case 6:
      DrOpt[10]*=-1;
      break;
    case 7:
      DrOpt[13]*=-1;
      break;
    case 8:
      DrOpt[7]*=-1;
      break;
    case 9:
      DrOpt[4]*=-1;
      break;
    case 10:
      DrOpt[1]*=-1;
      break;
    case 11:
      DrOpt[0]*=-1;
      break;
    case 12:
      DrOpt[14]*=-1;
      break;
    case 13:
      break;

    case 14: // save drawing modes
      std::ofstream DrawOpt("drawing.options");
      if(!DrawOpt) {
	std::cout<<"TDisplay::DrawModeMenu() ==> Can not crerate drawing.options file. "<<std::endl
	    <<"Check your write permission for current directory"<<std::endl;
	assert(false);
      }
      DrawOpt<<Proj<<" ";
      DrawOpt<<Xvf[0]<<" "<<Xvf[1]<<" "<<Yvf[0]<<" "<<Yvf[1]<<" ";
      for(int i = 0; i < 20; i++) DrawOpt<<DrOpt[i]<<" ";
      DrawOpt.close();

      std::cout<<std::endl<<"Drawing options are saved"
	  <<std::endl<<"To use defaults, erase 'drawing.options' file in the currect directory"<<std::endl;
      break;
    }


  ISELNT(10); // NT defined by view field size

}












