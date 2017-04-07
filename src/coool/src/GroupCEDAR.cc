#include "GroupCEDAR.h"
#include "TThread.h"
#include "TStyle.h"
#include "GroupPanel.h"
#include "PlaneCEDAR.h"
#include "TriggerTime.h"
#include "DaqEvent.h"

ClassImp(GroupCEDAR);

void GroupCEDAR::Init() {

}


void GroupCEDAR::EndEvent(const CS::DaqEvent &event) {
  // fRateCounter++; removed, done by Monitor now
  Group::EndEvent();
}


void GroupCEDAR::ControlPanel(const TGWindow *p, const TGWindow *main) {
  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}

