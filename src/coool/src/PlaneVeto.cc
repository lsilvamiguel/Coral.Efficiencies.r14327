#include "PlaneVeto.h"
#include "Plane1V.h"
#include "Plane1VPanel.h"
#include "ChipF1.h"
#include "TriggerTime.h"

#include "TF1.h"

ClassImp(PlaneVeto);

// implement your member functions here


void PlaneVeto::Init(TTree* tree) {

  Plane1V::Init(tree);

  // Book your histograms here :

  {
    for (int i=0; i<fMaxTriggerNumber; i++)
      {
	char istring[12];
	sprintf(istring,"%d",i);
	std::string name = fName + "_t_" + istring;
	std::string title = fName + " trigger " + istring;
	fHtrigger_spec_time[i] = new TH1F(name.c_str(),title.c_str(),	
					  fVt->GetNbins(),
					  fVt->GetMin(),
					  fVt->GetMax());
	AddHistogram(fHtrigger_spec_time[i]);
	name = fName + "_tVSch_" + istring;
	title = fName + " time vs channel for trigger " + istring;
	fHtrigger_spec_timeVSch[i] = new TH2F(name.c_str(),title.c_str(),	
					      fVt->GetNbins(),
					      fVt->GetMin(),
					      fVt->GetMax(),
					      fVt->GetNbins(),
					      fVt->GetMin(),
					      fVt->GetMax());
	AddHistogram(fHtrigger_spec_timeVSch[i]);
      }
  }
}





void PlaneVeto::EndEvent(const CS::DaqEvent &event) {

  // Std EndEvent
  Plane1V::EndEvent(event);

  // Class Customized EndEvent

//   const float reftime = fIsHR ? event.GetTT().GetTimeHigh() :
//                                 event.GetTT().GetTimeNorm();


  if (thr_flag) TThread::Lock();
//   CS::DaqEvent::Header head=event.GetHeader();
  register unsigned int triggerwd = event.GetTrigger() ;

//   for (int i=0; i<fNhits; i++) {
//     double time = CS::ChipF1::TimeDifference(fTime[i],reftime);
//     int channel=fChannel[i];
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"PlaneVeto::EndEvent: a digit is not a F1 one, strange...\n";
      continue;
    }
    register double time = (double) (iii->GetTimeDecoded() / iii->GetTimeUnit());
    register int channel = iii->GetChannel();


    if( fVch->Test(channel) ) {
      if( fVt->Test(time) ) {
	unsigned int m=1;
	for (int i=0; i<fMaxTriggerNumber; ++i, m<<=1)
// 	  if ((head.typeAttribute[1] & 0x0fffffff)== m)
	  if ((triggerwd & 0x0fffffff)== m)
	    {
	      fHtrigger_spec_time[i]->Fill(time);
	      fHtrigger_spec_timeVSch[i]->Fill(channel,time);
	    }
      }
    }
  }
  if (thr_flag) TThread::UnLock();
}



