#include "GroupBeamKiller.h"
#include "TThread.h"
#include "TStyle.h"
#include "GroupPanel.h"
#include "PlaneBeamKiller.h"
#include "TriggerTime.h"
#include "DaqEvent.h"

ClassImp(GroupBeamKiller);

void GroupBeamKiller::Init() {
  const unsigned int numPlanes = fPlanes.size();

  if (numPlanes>0) for (unsigned int i=0; i<numPlanes-1; i++)
    for (unsigned int j=i+1; j<numPlanes; j++) {
      std::ostringstream nametmp;
      nametmp<<fName << "_tc_" << i << "-" << j << std::ends;
      std::string name = nametmp.str();
      nametmp<<fName << "_tc_" << fPlanes[i]->GetName() << "-" << fPlanes[j]->GetName() << std::ends;
      std::string title = nametmp.str();
      fHistList.push_back(new TH1F_Ref(name.c_str(), title.c_str(),
		  (dynamic_cast<const Plane1V*>(fPlanes[i]))->GetTimeVariable().GetNbins(),
		  (dynamic_cast<const Plane1V*>(fPlanes[i]))->GetTimeVariable().GetMin(),
		  (dynamic_cast<const Plane1V*>(fPlanes[i]))->GetTimeVariable().GetMax(),
		  fRateCounter));
    }

#if USE_DATABASE == 1
  setDBpt(fDataBase);
#endif

  OpenReference();
  if (fReferenceDirectory) {
    for (unsigned int i=0; i < fHistList.size() ; i++) {
      (dynamic_cast<TH1F_Ref*>(fHistList[i]))->SetReference(fReferenceDirectory);
    }
  }

}


void GroupBeamKiller::EndEvent(const CS::DaqEvent &event) {
  // fRateCounter++; removed, done by Monitor now
  Group::EndEvent();

  if (thr_flag) TThread::Lock();

  register unsigned int histIndex = 0;
  register unsigned int numPlanes = fPlanes.size();
  if (numPlanes>0) for (register unsigned int i=0; i<numPlanes-1; i++) {
    for (register unsigned int j=i+1; j<numPlanes; j++) {

      register int numTime1 = (dynamic_cast<const Plane1V*>(fPlanes[i]))->GetTimeVariable().GetNvalues();
      register int numTime2 = (dynamic_cast<const Plane1V*>(fPlanes[j]))->GetTimeVariable().GetNvalues();

      register float* valuesTime1 =(dynamic_cast<const Plane1V*>(fPlanes[i]))->GetTimeValues();
      register float* valuesTime2 =(dynamic_cast<const Plane1V*>(fPlanes[j]))->GetTimeValues();

      for (register int t1=0; t1 < numTime1 ; t1++)
        for (register int t2=0; t2 < numTime2 ; t2++) {
	  (dynamic_cast<TH1F_Ref*>(fHistList[histIndex]))->Fill(valuesTime2[t2]-valuesTime1[t1]);
	}
      histIndex++;
    }
  }

  if (thr_flag) TThread::UnLock();
}


void GroupBeamKiller::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}

