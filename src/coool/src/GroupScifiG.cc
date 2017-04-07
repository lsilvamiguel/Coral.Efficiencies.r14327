#include "GroupScifiG.h"
#include "TThread.h"
#include "GroupPanel.h"

ClassImp(GroupScifiG);

void GroupScifiG::Init() {
  fHistList.push_back(new TH2F((fName+"_y1VSx1").c_str(),"y1 VS x1",36,0,144,36,0,144));
  fHistList[0]->SetOption("box");
}

void GroupScifiG::EndEvent(const CS::DaqEvent &event) {

  if(fPlanes.size()<2) return;

  std::vector<Variable*>& vx1=((Plane*)fPlanes[0])->GetVariables(); // Plane* cast is needed by gcc 2.95.2
  int nhitsx1=((Plane*)fPlanes[0])->GetNhits();
  float* chx1=vx1[0]->GetValues();
  //float* tx1=vx1[1]->GetValues();
	
  std::vector<Variable*>& vy1=((Plane*)fPlanes[1])->GetVariables();
  int nhitsy1=((Plane*)fPlanes[1])->GetNhits();
  float* chy1=vy1[0]->GetValues();
  //  float* ty1=vy1[1]->GetValues();

  if (thr_flag) TThread::Lock();
  for(int i=0; i<nhitsx1;i++) {
    for(int j=0; j<nhitsy1;j++) {
      fHistList[0]->Fill(chx1[i],chy1[j]);
    }
  }      
  if (thr_flag) TThread::UnLock();
}

void GroupScifiG::ControlPanel(const TGWindow *p, const TGWindow *main) {
  
  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}









