#include "GroupScifiJ.h"
#include "TThread.h"
#include "GroupScifiJPanel.h"

ClassImp(GroupScifiJ);


void GroupScifiJ::Init() {
  fHistList.push_back(new TH2F((fName+"_y1VSx1").c_str(),"y1 VS x1",96,0,96,96,0,96));
  fHistList.push_back(new TH2F((fName+"_y2VSx2").c_str(),"y2 VS x2",96,0,96,96,0,96));
  fHistList[0]->SetOption("box");
  fHistList[1]->SetOption("box");
}

void GroupScifiJ::EndEvent(const CS::DaqEvent &event) {

  if(fPlanes.size()<4) return;

  std::vector<Variable*>& vx1=((Plane*)fPlanes[0])->GetVariables(); // Plane* cast is needed by gcc 2.95.2
  float nhitsx1=((Plane*)fPlanes[0])->GetNhits();
  float* chx1=vx1[0]->GetValues();
  //float* tx1=vx1[1]->GetValues();
	
  std::vector<Variable*>& vy1=((Plane*)fPlanes[1])->GetVariables();
  float nhitsy1=((Plane*)fPlanes[1])->GetNhits();
  float* chy1=vy1[0]->GetValues();
  //float* ty1=vy1[1]->GetValues();

  std::vector<Variable*>& vx2=((Plane*)fPlanes[2])->GetVariables();
  float nhitsx2=((Plane*)fPlanes[2])->GetNhits();
  float* chx2=vx2[0]->GetValues();
  //float* tx2=vx2[1]->GetValues();
	
  std::vector<Variable*>& vy2=((Plane*)fPlanes[3])->GetVariables();
  float nhitsy2=((Plane*)fPlanes[3])->GetNhits();
  float* chy2=vy2[0]->GetValues();
  //  float* ty2=vy2[1]->GetValues();

  if (thr_flag) TThread::Lock();
  for(int i=0; i<nhitsx1;i++) {
    for(int j=0; j<nhitsy1;j++) {
      fHistList[0]->Fill(chx1[i],chy1[j]);
    }
  }      

  for(int i=0; i<nhitsx2;i++) {
    for(int j=0; j<nhitsy2;j++) {
      fHistList[1]->Fill(chx2[i],chy2[j]);
    }
  }
  if (thr_flag) TThread::UnLock();
}

void GroupScifiJ::ControlPanel(const TGWindow *p, const TGWindow *main) {
  
  if (!fControlPanel) fControlPanel = new GroupScifiJPanel(p, main, 100, 100, this);
}









