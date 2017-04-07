#include "GroupDAQ.h"
#include "TThread.h"
#include "TStyle.h"
#include "DaqError.h"
#include "Reference.h"
#include "Plane.h"

#include <vector>

using std::vector;

ClassImp(GroupDAQ);


void GroupDAQ::Init() {
  fRateCounter = 0;
  fHistList.push_back(new TH1F("spill_number","spill number",200,0,200));
  fHistList.push_back(new TH1F_Ref("trigger_bits","trigger bits",12,0,12,fRateCounter));
  fHistList.push_back(new TH1F_Ref("eventsize","Event Size (Bytes)",100,0,1E5,fRateCounter));
  fHistList.push_back(new TH1F_Ref("number_of_errors_few","DAQ Errors",100,0,100,fRateCounter));
  fHistList.push_back(new TH1F_Ref("number_of_errors_many","DAQ Error levels",4,0,4,fRateCounter));
  fHistList.push_back(new TH1F_Ref("trigger_rates","trigger rates",12,0,12,fRateCounter));
  fHistList.push_back(new TH1F_Ref("trigger_type","trigger type",12,0,12,fRateCounter));
  fHistList.push_back(new TH1F("evtnum_in_spill","event number in spill",250,0,25000));
  fHistList.push_back(new TH2I("err_sourceID","errors per source ID",100,0,100,1012,0,1012));
  fHistList.push_back(new TH1I("errors_in_spill","errors in spill",5000,0,25000));
  gStyle -> SetPalette(1);
  fHistList[8]->SetOption("COLZ");
#if USE_DATABASE == 1
  setDBpt(fDataBase);
#endif
  OpenReference();

  if (fReferenceDirectory) {
    for (int i=1; i<=6;i++) {
      ((TH1F_Ref*)fHistList[i])->SetReference(fReferenceDirectory);
    }
  }
}

void GroupDAQ::EndEvent(const CS::DaqEvent &event)  {
  // fRateCounter++; removed, done by Monitor now
  fHistList[0]->Fill(event.GetBurstNumber());
  int errors_in_src[1012];
  int sourceID;
  for(int i=0;i<1012;i++)
  {
      errors_in_src[i]=0;
  }

//   CS::DaqEvent::Header head=event.GetHeader();
  register unsigned int triggerwd = event.GetTrigger() ;
  
  //cout<<"***************************** "
  //    <<head.typeAttribute[0]<<" "<<head.typeAttribute[1]<<dec<<endl;

  int m=1;
  for(size_t i=0; i<sizeof(int)*8; i++) {

//     if(int bit = ((head.typeAttribute[1] & m) >> i)) {
//     if ((head.typeAttribute[1] & m) >> i) {
    if ((triggerwd & m) >> i) {
      fHistList[1]->Fill(i);
      //cout<<i<<endl;
    }
    m = m << 1;    
  }

  fHistList[2]->Fill(event.GetLength());
  fHistList[7]->Fill(event.GetEventNumberInBurst());

    Int_t nof_errors=0;
    for( vector<CS::DaqError>::const_iterator it=event.GetDaqErrors().errors.begin(); it!=event.GetDaqErrors().errors.end(); it++ )
    {
        CS::DaqError::SeverityLevel level;
        level=it->GetDaqErrorType().GetSeverityLevel();
        fHistList[4]->Fill(level);
        if( level==CS::DaqError::MINOR_PROBLEM ||
                level==CS::DaqError::WARNING ||
                level==CS::DaqError::SEVERE_PROBLEM )
        {
          nof_errors++;

          //sourceID = it->GetArg(CS::DaqError::SOURCE_ID);
          errors_in_src[0]++;
        }
    }

// Commented out part due to changes in decoding lib (BTW, I don't know what this part did...) DN 5/9/2007
//
//         for( std::vector<CS::Chip*>::const_iterator chip=event.GetChips().begin();
//              chip!=event.GetChips().end(); chip++ )
//         {
// 	  CS::DaqError::SeverityLevel level;
// 	  for( vector<CS::DaqError>::const_iterator e=(*chip)->GetDaqErrors().errors.begin();
// 	       e!=(*chip)->GetDaqErrors().errors.end(); e++ )
//             {
//               level=e->second.GetDaqErrorType().GetSeverityLevel();
// 	      fHistList[4]->Fill(level);
// 	      if (level==CS::DaqError::MINOR_PROBLEM || level==CS::DaqError::SEVERE_PROBLEM)
// 		nof_errors++;
// 	    }
// 	}

	fHistList[3]->Fill(nof_errors);
    for(int i=0;i<1012;i++)
    {
        fHistList[8]->Fill(errors_in_src[i],i);
    }
    fHistList[9]->SetBinContent(event.GetEventNumberInBurst(),nof_errors);


  fHistList[6]->Fill(event.GetType());
 
};








