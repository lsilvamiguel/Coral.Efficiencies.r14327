#include "GroupTrigHodo.h"
#include "TThread.h"
#include "GroupPanel.h"


TriggerMatrix::TriggerMatrix(Plane1V* plane4,
			     Plane1V* plane5,
			     unsigned int trigger_mask,
			     std::vector<TH1*> *hist_list,
			     int index)
{
  p4=plane4;
  p5=plane5;
  if (!p4 || !p5) return;
  trigger=trigger_mask;

  char key[256];
  std::string namebase,name;
  namebase += p4->GetName();
  namebase += " VS ";
  namebase += p5->GetName();
  sprintf(key,"matrix%d_good",index);
  name = namebase + " (This trigger, in time)";
  hist_good = new TH2F(key,name.c_str(),
		  p5->GetNchannels(),0,p5->GetNchannels(),
		  p4->GetNchannels(),0,p4->GetNchannels());
  hist_list->push_back(hist_good);
  hist_good->SetOption("box");

//    sprintf(key,"matrix%d_all",index);
//    name = namebase + " (all combinations)";
//    hist_all = new TH2F(key,name.c_str(),
//  		  p5->GetNchannels(),0,p5->GetNchannels(),
//  		  p4->GetNchannels(),0,p4->GetNchannels());
//    hist_list->push_back(hist_all);
//    hist_all->SetOption("box");

//    sprintf(key,"matrix%d_intime",index);
//    name = namebase + " (all triggers, in time)";
//    hist_intime = new TH2F(key,name.c_str(),
//  		  p5->GetNchannels(),0,p5->GetNchannels(),
//  		  p4->GetNchannels(),0,p4->GetNchannels());
//    hist_list->push_back(hist_intime);
//    hist_intime->SetOption("box");

  sprintf(key,"matrix%d_mult1",index);
  name = namebase + " (this trigger, in time, mult1)";
  hist_mult1 = new TH2F(key,name.c_str(),
		  p5->GetNchannels(),0,p5->GetNchannels(),
		  p4->GetNchannels(),0,p4->GetNchannels());
  hist_list->push_back(hist_mult1);
  hist_mult1->SetOption("box");
   
}
void
TriggerMatrix::Fill(const CS::DaqEvent &event)
{
  if (!p4 || !p5) return;
//   CS::DaqEvent::Header head=event.GetHeader();
//   unsigned int real_trigger=head.typeAttribute[1];
  register unsigned int real_trigger=event.GetTrigger();

  std::vector<Variable*>& vh4=p4->GetVariables();
  int nh4=p4->GetNhits();;
  float* ch4=vh4[0]->GetValues();
  float* t4=vh4[1]->GetValues();
  Variable *t_on_trig_4=p4->GetVariable("_t_on_trigger");

  std::vector<Variable*>& vh5=p5->GetVariables();
  int nh5=p5->GetNhits();
  float* ch5=vh5[0]->GetValues();
  float* t5=vh5[1]->GetValues();
  Variable *t_on_trig_5=p5->GetVariable("_t_on_trigger");
	
  if ((real_trigger & ~trigger)==0) {
    int mult_4=0;
    int mult_5=0;
    float hit4=0,hit5=0;
    for(int i=0; i<nh4;i++) {
      if (t_on_trig_4->Test(t4[i])) {
	mult_4++;
	hit4=ch4[i];
	for(int j=0; j<nh5;j++) {
	  if (t_on_trig_5->Test(t5[j])) { 
	    hist_good->Fill(ch5[j],ch4[i]);
	    mult_5++;
	    hit5=ch5[j];
	  }
	}
      }
    }
    if (mult_4==1 && mult_5==1) 
      hist_mult1->Fill(hit5,hit4);
  }
}



ClassImp(GroupTrigHodo);


void GroupTrigHodo::Init() {
  TriggerMatrix *TM;
  if (fPlanes.size() < 17) return;
  TM = new TriggerMatrix((Plane1V*)fPlanes[0],(Plane1V*)fPlanes[1], // I u 
			 1, &fHistList, 0);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[2],(Plane1V*)fPlanes[3], // I d
			 1, &fHistList, 1);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[4],(Plane1V*)fPlanes[5], // L
			 4, &fHistList, 2);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[6],(Plane1V*)fPlanes[7], // MX u
			 0x102, &fHistList, 3);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[8],(Plane1V*)fPlanes[9], // MX d
			 0x102, &fHistList, 4);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[10],(Plane1V*)fPlanes[11], // MY u
			 0x1021, &fHistList, 5);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[12],(Plane1V*)fPlanes[13], // MY d
			 0x102, &fHistList, 6);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[11],(Plane1V*)fPlanes[7], // MXY u
			 0x102, &fHistList, 7);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[13],(Plane1V*)fPlanes[9], // MXY d
			 0x102, &fHistList, 8);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[14],(Plane1V*)fPlanes[15], // O1 u
			 8, &fHistList, 9);
  TriggerMatrixlist.push_back(TM);
  TM = new TriggerMatrix((Plane1V*)fPlanes[14],(Plane1V*)fPlanes[16], // O2 u
			 8, &fHistList, 10);
  TriggerMatrixlist.push_back(TM);
}

void GroupTrigHodo::EndEvent(const CS::DaqEvent &event) {
  if (thr_flag) TThread::Lock();
  for (std::vector<TriggerMatrix*>::iterator it = TriggerMatrixlist.begin();
       it!=TriggerMatrixlist.end(); it++) {
    TriggerMatrix *mat;
    mat=*it;
    mat->Fill(event);
  }
  if (thr_flag) TThread::UnLock();
}


void GroupTrigHodo::ControlPanel(const TGWindow *p, const TGWindow *main) {
  
  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}










