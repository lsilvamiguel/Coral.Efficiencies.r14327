#include "PlaneBMS.h"
#include "Plane1V.h"
#include "Plane1VPanel.h"
#include "ChipF1.h"
#include "TriggerTime.h"
#include "Reference.h"

#include "TF1.h"

ClassImp(PlaneBMS);

// implement your member functions here


void PlaneBMS::Init(TTree* tree) {

  Plane1V::Init(tree);

  // Book your histograms here :

  {
    for (int i=0; i<fMaxTriggerNumber; i++)
      {
	char istring[12];
	sprintf(istring,"%d",i);
	std::string name = fName + "_t_" + istring;
	std::string title = fName + " trigger " + istring;
	fHtrigger_spec_time[i] = new TH1F_Ref(name.c_str(),title.c_str(),	
					      fVt->GetNbins(),
					      fVt->GetMin(),
					      fVt->GetMax(),
					      fRateCounter);
	((TH1F_Ref*)fHtrigger_spec_time[i])->SetReference(fReferenceDirectory);
	AddHistogram(fHtrigger_spec_time[i]);
      }
  }
  if (strncmp(fName.c_str(),"BM03P1",6)==0) {
    std::string tname = fName + "_t_corr";
    fBMS3_repaired_timing=new TH1F_Ref(tname.c_str(),tname.c_str(),
				       fVt->GetNbins(),
				       fVt->GetMin(),
				       fVt->GetMax(),
				       fRateCounter);
    ((TH1F_Ref*)fBMS3_repaired_timing)->SetReference(fReferenceDirectory);
    AddHistogram(fBMS3_repaired_timing);
  } else {
    fBMS3_repaired_timing=NULL;
  }
}





void PlaneBMS::EndEvent(const CS::DaqEvent &event) {

  // Std EndEvent
  Plane1V::EndEvent(event);

  // Class Customized EndEvent

//   const float reftime = fIsHR ? event.GetTT().GetTimeHigh() :
//                                 event.GetTT().GetTimeNorm();

  if (thr_flag) TThread::Lock();

//   CS::DaqEvent::Header head=event.GetHeader();
  register unsigned int triggerwd = event.GetTrigger();

//   for (int i=0; i<fNhits; i++) {
//     double time = CS::ChipF1::TimeDifference(fTime[i],reftime);
//     int channel=fChannel[i];
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"PlaneBMS::EndEvent: a digit is not a F1 one, strange...\n";
      continue;
    }
    register double time = (double) (iii->GetTimeDecoded() / iii->GetTimeUnit());
    register int channel = iii->GetChannel();


    if( fVch->Test(channel) ) {
      if( fVt->Test(time) ) {
	if (fBMS3_repaired_timing) {
	    double t = time;
	    if (channel >=40) t+=5.165728e+03;
	    fBMS3_repaired_timing->Fill(t);
	}
	unsigned int m=1;
        for (int i=0; i<fMaxTriggerNumber; ++i, m<<=1)
// 	  if ((head.typeAttribute[1] & 0x0fffffff) == m)
// 	  if ((triggerwd & 0x0fffffff) == m)
          if ((triggerwd & 0x0fff) == m)
	    fHtrigger_spec_time[i]->Fill(time);
      }
    }
  }
  if (thr_flag) TThread::UnLock();
}


