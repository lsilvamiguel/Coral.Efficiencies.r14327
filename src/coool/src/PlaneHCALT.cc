
#include "PlaneHCALT.h"
#include "TriggerTime.h"

ClassImp(PlaneHCALT);
using namespace std;

void PlaneHCALT::Init(TTree* tree) {

  Plane2V::Init(tree);

  // Book your histograms here :
}


void PlaneHCALT::StoreDigit(CS::Chip::Digit* digit) {

  lDigits.push_back(digit);
  fNhits++;
//cerr<<fNhits<<" "<<GetName()<<" PlaneHCALT stored digits ////////////////////////////////////////////////////////////////////////////////////////////////////////"<<endl;
//     CS::ChipF1::Digit* f1digit = dynamic_cast<CS::ChipF1::Digit*>(digit);
//         cerr<<"calling storedigit"<<std::endl;
//     if(f1digit) {
//              cerr<<f1digit->GetX()<<" "<<f1digit->GetY()<<" "<<f1digit->GetTimeDecoded()/f1digit->GetTimeUnit()<<std::endl;
//     Plane2V::StoreDigit(f1digit->GetX(),f1digit->GetY(),f1digit->GetTime());
//     }
}


void PlaneHCALT::EndEvent(const CS::DaqEvent &event) {

//   const float reftime = event.GetTT().GetTimeNorm();
 
  if (thr_flag) TThread::Lock();
  fNrows = 2; //this is the value needed for HC02 Trigger
  if (fName == "HC01P1T1" || fName == "HC01P1T2" || fName == "HC01P3T1" || fName == "HC01P3T2")  fNrows = 5;
  if (fName == "HC01P2T1" || fName == "HC01P2T2" || fName == "HC01P4T1" || fName == "HC01P4T2") fNrows = 4; 



  //std::cout<<"EndEvent"<<std::endl;
//   for (int i=0; i<fNhits; i++) {
//     int row=fRow[i];
//     int col=fCol[i];
//     int time = (int) CS::ChipF1::TimeDifference((double)fAmp[i], (double)reftime);

//         const CS::DaqEvent::Header head = event.GetHeader();
//         const CS::DaqEvent::EventType eventType = head.GetEventType();

  int TM = 0xFFFE&event.GetTrigger() ;

  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"Plane1V::EndEvent ("<<GetName()<<"): a digit is not a F1 one, strange...\n";
      continue;
    }
    register int time = (int) (iii->GetTimeDecoded() / iii->GetTimeUnit());
    register int col = iii->GetX();
    register int row = iii->GetY();

    if(fVrow->Test(row) &&
       fVcol->Test(col) &&
       fVamp->Test(time)) {
//    std::cerr<<"PlaneHCALT::EndEvent: "<<GetName()<<" raw= "<<row<<" col= "<<col<<" time= "<<time<<std::endl;

      int adr=row + col*fNrows;

      fHa->Fill(time);
      fHavsadr->Fill(adr,time);
      fHrc1->Fill(col,row);
//      fHrca->Fill(col,row,-time);

      fVrow->Store(row);
      fVcol->Store(col);
      fVamp->Store(time);

      fNhitsKept++;
    }
  }

//cout<<"fNhits in PlaneHCALT= "<<fNhits<<"           fNhitsKept= "<<fNhitsKept<<" Plane Name= "<<GetName()<< "------------------------------------"<<endl;


  fHhit->Fill(fNhitsKept);
  if (thr_flag) TThread::UnLock();
}


