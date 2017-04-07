#ifndef __GroupTrigHodo__
#define __GroupTrigHodo__

#include "Group.h"
#include "Plane1V.h"
class TriggerMatrix {
 private:
  Plane1V* p4;
  Plane1V* p5;
  unsigned int trigger;
  TH2F *hist_good;
  TH2F *hist_all;
  TH2F *hist_intime;
  TH2F *hist_mult1;
 public:
  TriggerMatrix(Plane1V* plane4,
		Plane1V* plane5,
		unsigned int trigger_mask,
		std::vector<TH1*> *hist_list,
		int index);
#ifndef __CINT__
  void Fill(const CS::DaqEvent &event);
#endif
};

class GroupTrigHodo : public Group {
 private:
  std::vector<TriggerMatrix*> TriggerMatrixlist;  
 public:
  GroupTrigHodo(const char* name) : Group(name) {}
  void Init();
 
#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif  
  void ControlPanel(const TGWindow *p, const TGWindow *main);
  

  ClassDef(GroupTrigHodo,0)
};

#endif
