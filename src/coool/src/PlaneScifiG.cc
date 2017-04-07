#include "PlaneScifiG.h"
#include "TProfile.h"
#include "ChipF1.h"
#include "TriggerTime.h"

ClassImp(PlaneScifiG);

const int PlaneScifiG::fMAX_MULT = 8;
// const int PlaneScifiG::fUP_BOUND = 30000;
// const int PlaneScifiG::fLOW_BOUND = -30000;
// const int PlaneScifiG::fMAX_CLOCK = 64239;

PlaneScifiG::PlaneScifiG(const char *detname,int nchan, int center, int width) :
  Plane(detname),fNchan(nchan) {

  INeedGeom();


  fVch = AddVariable("_ch",fNchan,0,fNchan,fNchan*fMAX_MULT);
  fVtl = AddVariable("_tl",100,center-width,center+width,fNchan*fMAX_MULT);
  fVth = AddVariable("_th",100,center-width,center+width,fNchan*fMAX_MULT);
  fVstat = AddVariable("_stat", 3, -1, 2, fNchan*fMAX_MULT);
 fVtc = AddVariable("_tcorr",100,-1000,1000,fNchan*fMAX_MULT);
}

void PlaneScifiG::Init(TTree* tree) {

  // 0 multiplicity within the time cut on the whole detector
  std::string hitsname = fName + "_hits";
  fHhits=new TH1F_Ref(hitsname.c_str(),hitsname.c_str(),50,0,50,fRateCounter);
  ((TH1F_Ref*)fHhits)->SetReference(fReferenceDirectory);
  AddHistogram(fHhits);
  std::string hitsleavlist = hitsname + "/I";

  // 1 histogram : hit map
  std::string chname = fName + fVch->GetName();
  std::string chleavlist = chname + "[" + hitsname +"]/F";

  std::string chname_h = chname + "_h";
  fHchh=new TH1F_Ref(chname_h.c_str(),chname_h.c_str(),
		     fVch->GetNbins(),
		     fVch->GetMin(),
		     fVch->GetMax(),fRateCounter);
  ((TH1F_Ref*)fHchh)->SetReference(fReferenceDirectory);
  AddHistogram(fHchh);

  // 2 histogram : hit map
  std::string chname_l = chname  + "_l";
  fHchl=new TH1F_Ref(chname_l.c_str(),chname_l.c_str(),
		 fVch->GetNbins(),
		 fVch->GetMin(),
		 fVch->GetMax(),fRateCounter);
  ((TH1F_Ref*)fHchl)->SetReference(fReferenceDirectory);
  AddHistogram(fHchl);

  // 3 histogram : tl
  std::string d1name = fName + fVtl->GetName();
  fHtl=new TH1F_Ref(d1name.c_str(),d1name.c_str(),
		fVtl->GetNbins(),
		fVtl->GetMin(),
		fVtl->GetMax(),fRateCounter);
  ((TH1F_Ref*)fHtl)->SetReference(fReferenceDirectory);
  AddHistogram(fHtl);
  std::string d1leavlist = d1name + "[" + hitsname +"]/F";

  // 4 histogram : th
  std::string d2name = fName + fVth->GetName();
  fHth=new TH1F(d2name.c_str(),d2name.c_str(),
		fVth->GetNbins(),
		fVth->GetMin(),
		fVth->GetMax());
  AddHistogram(fHth);
  std::string d2leavlist = d2name + "[" + hitsname +"]/F";

  // 5 histogram : Data2-Data1 % channel
  std::string d21name = fName + fVth->GetName() + "-" +
    fName + fVtl->GetName();
  fHdtvsch=new TH2F(d21name.c_str(),d21name.c_str(),
		    fVch->GetNbins(),
		    fVch->GetMin(),
		    fVch->GetMax(),
		    100,-1000,1000);
  fHdtvsch->SetOption("col");
  AddHistogram(fHdtvsch);

  // 6 histogram: low threshold trigger time correlation
  std::string tlname = fName + "_tlVSch";
  fHtlvsch=new TH2F(tlname.c_str(),tlname.c_str(),
		  fVch->GetNbins(),
		  fVch->GetMin(),
		  fVch->GetMax(),
		  fVtl->GetNbins(),
		  -10000,-6000);
  fHtlvsch->SetOption("col");
  AddHistogram(fHtlvsch);

  // 7 histogram: high threshold trigger time correlation
  std::string thname = fName + "_thVSch";
  fHthvsch=new TH2F(thname.c_str(),thname.c_str(),
                    fVch->GetNbins(),
                    fVch->GetMin(),
                    fVch->GetMax(),
                    fVth->GetNbins(),
                    -10000,-6000);
  AddHistogram(fHthvsch);

 // 8 histogram : pos
  std::string d3name = fName + fVstat->GetName();
  fHstat=new TH1F(d3name.c_str(),d3name.c_str(),
		   fVstat->GetNbins(),
		   fVstat->GetMin(),
		   fVstat->GetMax());
  AddHistogram(fHstat);
  std::string d3leavlist = d3name + "[" + hitsname +"]/F";


  std::string tcname = fName+"_tcorr"; 
  fHtc=new TH1F(tcname.c_str(),tcname.c_str(), 400, -80, 80);
  AddHistogram(fHtc);
  fHtc->GetXaxis()->SetTitle("corrected time (ns)");
  


  // Cluster histograms
  float posmax = 50;

  if(fGeom)
    posmax = fGeom->GetHalfSize();

  fVcch = AddVariable("_clusterPosition",100,-posmax,posmax,fNchan*fMAX_MULT);
  fVcs = AddVariable("_clusterSize",30, -0.5, 29.5, fNchan*fMAX_MULT);

  // clusters ----------------------------------------------------------

  //multiplicity within the time cut on the whole detector
  std::string chitsname = fName + "_clusterMultiplicity";
  fHchits=new TH1F(chitsname.c_str(),chitsname.c_str(),70,-0.5,69.5);
  AddHistogram(fHchits);
  std::string chitsleavlist = chitsname + "/I";

  //branch : cluster map
  std::string cchname = fName + fVcch->GetName();
  fHcch=new TH1F(cchname.c_str(),cchname.c_str(),
                 fVcch->GetNbins(), fVcch->GetMin(), fVcch->GetMax());
  AddHistogram(fHcch);
  std::string cchleavlist = cchname + "[" + chitsname +"]/F";

  //Cluster size
  std::string csname = fName+"_clusterSize";
  fHcs=new TH1F(csname.c_str(),csname.c_str(),
                fVcs->GetNbins(), fVcs->GetMin(), fVcs->GetMax());
  AddHistogram(fHcs);
  std::string csleavlist = csname + "[" + chitsname +"]/F";


  fHchhVsTis=new TH2F((fName+"_ch_h_VStis").c_str(), (fName+" channel (high thr.) versus time in spill").c_str(), 
		                  fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),
                      100, 0., 15.);
  fHchhVsTis->SetOption("COLZ");
  AddHistogram(fHchhVsTis);
  fHchlVsTis=new TH2F((fName+"_ch_l_VStis").c_str(), (fName+" channel (low thr.) versus time in spill").c_str(),
		                  fVch->GetNbins(), fVch->GetMin(), fVch->GetMax(),
                      100, 0., 15.);
  fHchlVsTis->SetOption("COLZ");
  AddHistogram(fHchlVsTis);


  if(tree) {
    fIsInTree = true;
    tree->Branch(hitsname.c_str(),&fNhitsKept,
		 hitsleavlist.c_str(),32000);
    tree->Branch(chname.c_str(),fVch->GetValues(),
		 chleavlist.c_str(),32000);
    tree->Branch(d1name.c_str(),fVtl->GetValues(),
		 d1leavlist.c_str(),32000);
    tree->Branch(d2name.c_str(),fVth->GetValues(),
		 d2leavlist.c_str(),32000);
    tree->Branch(d3name.c_str(),fVstat->GetValues(),
		 d3leavlist.c_str(),32000);

  }
}

PlaneScifiG::~PlaneScifiG() {
//   delete fChannel; delete fData;
//   delete fPos;
}

void PlaneScifiG::StoreDigit(int channel, int data, int pos) {
  std::cerr<<"PlaneScifiG::StoreDigit(int, int, int) ("<<GetName()<<"): Error ! this method is obsolete !!!\n";
  // if (fNhits==0) fDigits.clear();
  // std::vector<double> dataX;
  //  dataX.push_back(data);
  // fDigits.insert(make_pair(channel,CDigit1(channel,dataX)));
//   if (fNhits < fNchan*fMAX_MULT) {
//     fChannel[fNhits]=channel;
//     fData[fNhits]=data;
//     fPos[fNhits]=pos;
//     fNhits++;
//   }
}

void PlaneScifiG::StoreDigit(CS::Chip::Digit* digit) {
//   std::vector<float> data=digit->GetNtupleData();
//   if(data.size()>2);
//   this->StoreDigit(data[0],
// 		   data[1],
// 		   data[2]);
  lDigits.push_back(digit);
  fNhits++;
}




void PlaneScifiG::EndEvent(const CS::DaqEvent &event) {
  if (thr_flag) TThread::Lock();

//   const float reftime = fIsHR ? event.GetTT().GetTimeHigh() :
//                                 event.GetTT().GetTimeNorm();
  
//   const float factor = fIsHR ? CS::ChipF1::GetUnitHigh() :  
//     CS::ChipF1::GetUnitNorm();
  
  
  
  fDigits.clear();
  
  //    bool haslow=false;
//    bool lasthit = false;
//    int lasttime=0;
//    int lastchan=0;
//    int last_hightime=0;
//    int last_highchan=0;

  //int fNlowhits = 0;
  //int fNhighhits = 0;

  int fLow = 0, fHigh= 0;
//   int time, chan, pos;

  //std::cout << "\nNEW EVENT (or plane in event)" << std::endl;

//    for (int i=0; i<fNhits; i++) {
//      if(IsLowThreshold(pos)) fNlowhits++;
//      else fNhighhits++;
//    }

//    int low_chan[fNlowhits];
//    int high_chan[fNhighhits];
//    int low_time[fNlowhits];
//    int high_time[fNhighhits];

  int low_chan[fNhits];
  int high_chan[fNhits];
  int low_time[fNhits];
  int high_time[fNhits];

//   for (register int i=0; i<fNhits; i++) {
//     time = (int)CS::ChipF1::TimeDifference(fData[i],reftime);
//     chan=fChannel[i];
  typedef std::list<CS::Chip::Digit*>::iterator lDIter;
  for (lDIter ii = lDigits.begin(); ii != lDigits.end(); ii++) {
    CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
    if (!iii) {
      std::cerr<<"Plane1V::EndEvent ("<<GetName()<<"): a digit is not a F1 one, strange...\n";
      continue;
    }
    register float factor = iii->GetTimeUnit();
    register int time = (int) (iii->GetTimeDecoded() / factor);
    register int chan = iii->GetChannel();
    register int pos = iii->GetChannelPos();

    float timeT0=-8000;
    if (fUseCalib && chan<(int)calib_data.size()) timeT0 = calib_data[chan].t0;
    //if (fUseCalib) std::cout <<"use"<<std::endl;
    double tcorr = (time - timeT0)*factor;
    //   std::cout <<"time"<<time<<"calib"<<timeT0<<"tcorr"<<tcorr<<"reftime"<<reftime<<std::endl;

    std::vector<double> dataX;
    dataX.push_back(tcorr);
    fDigits.insert(make_pair(chan,CDigit1(chan,dataX)));


  
    //  chan=fChannel[i];
//     pos =fPos[i];
    fHtc->Fill(tcorr);
    if(IsLowThreshold(pos)) {
      
      fHchl->Fill(chan);	
      fHtl->Fill(time);	
      fHtlvsch->Fill(chan,time);
      fHchlVsTis->Fill(chan, event.GetTT().GetTimeInSpill());	

      low_chan[fLow] = chan;
      low_time[fLow] = time;
      //std::cout << "low: " << fLow << " " << low_chan[fLow] << " " << low_time[fLow]<< std::endl;
      fLow++;
    }

    else {

      fHchh->Fill(chan);	
      fHth->Fill(time);	
      fHthvsch->Fill(chan,time);
      fHchhVsTis->Fill(chan, event.GetTT().GetTimeInSpill());	

      high_chan[fHigh] = chan;
      high_time[fHigh] = time;
      //std::cout << "high: " << fHigh << " " << high_chan[fHigh] << " " << high_time[fHigh] << std::endl;
      fHigh++;

    }

  }

  //std::cout << "That were the Rohdaten that was, now for the Daten wot is: " << std::endl;;
  for (register int i=0; i<fLow; i++) {
    for (register int j=0; j<fHigh; j++) {

      if (low_chan[i] != high_chan[j]) continue;
      //std::cout << "they're equal!" << std::endl;

// These two histograms have not been defined anywhere -> commented out !!
//       fHtimediff->Fill(high_time[j]-low_time[i]);
//       fHtdvsch->Fill(low_chan[i], high_time[j]-low_time[i]);

      if ( high_time[j] < low_time[i] + 50 && high_time[j] > low_time[i] - 25){
	//std::cout << "we have a winner!" << std::endl;
	if (low_time[i] == 0) continue;
	fVch->Store(low_chan[i]);
	fVtl->Store(low_time[i]);
	fVth->Store(high_time[j]);
	fVstat->Store(0);
	fHstat->Fill(0);
	
	//std::cout << fNhits << ": " << low_chan[i] << " " << high_chan[j] << " " <<low_time[i] << " " << high_time[j]  << std::endl;

	fNhitsKept++;
	low_time[i] = 0;
	high_time[j] = 0;
	
      }
    }
  }

  for (register int i=0; i<fLow; i++) {
    if (low_time[i] != 0){

	fVch->Store(low_chan[i]);
	fVtl->Store(low_time[i]);
	fVth->Store(-30000);
	fVstat->Store(-1);
	fHstat->Fill(-1);
	//std::cout << fNhits << ": " << low_chan[i] << " " <<low_time[i] << " " << "-30000"  << std::endl;
	fNhitsKept++;

    }
  }


  for (register int j=0; j<fHigh; j++) {
    if (high_time[j] != 0){

	fVch->Store(high_chan[j]);
	fVtl->Store(-30000);
	fVth->Store(high_time[j]);
	fVstat->Store(1);
	fHstat->Fill(1);
	//std::cout << fNhits << ": " << high_chan[j] << " " << -30000 << " " << high_time[j]  << std::endl;
	fNhitsKept++;

    }
  }

  fHhits->Fill(fNhitsKept);
  if (thr_flag) TThread::UnLock();

}

//  void PlaneScifiG::EndEvent(const CS::DaqEvent &event) {

//    const float reftime = fIsHR ? event.GetTT().GetTimeHigh() :
//                                  event.GetTT().GetTimeNorm();

//    bool haslow=false;
//    int lasttime=0;
//    int lastchan=0;
//    int last_hightime=0;
//    int last_highchan=0;

//    if (thr_flag) TThread::Lock();

//   for (int i=0; i<fNhits; i++) {
//      int time = (int) CS::ChipF1::TimeDifference(fData[i], reftime);
//  //     int time=CorrectRollOver(fData[i]-reftime);

//      int chan=fChannel[i];
//      int pos =fPos[i];

//      if(IsLowThreshold(pos)) {

//             fHchl->Fill(chan);	
//        fHtl->Fill(time);	
//        fHtlvsch->Fill(chan,time);

//        if (i == fNhits - 1){ // for the last hit
//  	fVch->Store(chan);
//  	fVtl->Store(time);
//  	fVth->Store(-30000);
//  	fVstat->Store(pos);
//  	fHstat->Fill(pos);
//  	fNhitsKept ++;
//        }

//        if (haslow == true && chan != lastchan){
//  	fVch->Store(lastchan);
//  	fVtl->Store(lasttime);
//  	fVth->Store(-30000);
//  	fVstat->Store(-1);
//  	fHstat->Fill(-1);
//  	fNhitsKept ++;
//        }

//        if (haslow == true && chan == lastchan){
//  	fVch->Store(lastchan);
//  	fVtl->Store(lasttime);
//  	fVth->Store(-30000);
//  	fVstat->Store(-1);
//  	fHstat->Fill(-1);
//  	fNhitsKept ++;
//        }

//        if (haslow == false && chan != lastchan && fNhits != 0){
//  	fVch->Store(last_highchan);
//  	fVtl->Store(-30000);
//  	fVth->Store(last_hightime);
//  	fVstat->Store(-1);
//  	fHstat->Fill(-1);
//  	fNhitsKept ++;
//        }

//        if (haslow == false && chan == lastchan && fNhits != 0){
//  	fVch->Store(last_highchan);
//  	fVtl->Store(-30000);
//  	fVth->Store(last_hightime);
//  	fVstat->Store(-1);
//  	fHstat->Fill(-1);
//  	fNhitsKept ++;
//        }

//        lasttime = time;
//        lastchan=chan;
//        haslow = true;
//      }

//      else{

//        fHchh->Fill(chan);	
//        fHth->Fill(time);	

//        if (i == fNhits - 1){ // for the last hit
//  	fVch->Store(chan);
//  	fVtl->Store(-30000);
//  	fVth->Store(time);
//  	fVstat->Store(pos);
//  	fHstat->Fill(pos);
//  	fNhitsKept ++;
//        }

//        if( haslow == true && chan == lastchan){

//  	fHdtvsch->Fill(chan,lasttime-time);	
//  	fHthvsch->Fill(chan,time);	
	
//  	fVch->Store(lastchan);
//  	fVtl->Store(lasttime);
//  	fVth->Store(time);
//  	fVstat->Store(0);
//  	fHstat->Fill(0);
//  	fNhitsKept++;
//  	}


//        if(haslow == false && chan == lastchan){
//  	fVch->Store(last_highchan);
//  	fVtl->Store(-30000);
//  	fVth->Store(last_hightime);
//  	fVstat->Store(-1);
//  	fHstat->Fill(-1);
//  	fNhitsKept++;
//  	}


//        if(haslow == true && chan != lastchan){
//  	fVch->Store(lastchan);
//  	fVtl->Store(lasttime);
//  	fVth->Store(-30000);
//  	fVstat->Store(-1);
//  	fHstat->Fill(-1);
//  	fNhitsKept++;
//  	}

//        if (haslow == false && chan != lastchan){
//  	fVch->Store(chan);
//  	fVtl->Store(-30000);
//  	fVth->Store(time);
//  	fVstat->Store(1);
//  	fHstat->Fill(1);
//  	fNhitsKept++;
//  	}
//        last_hightime = time;
//        last_highchan=chan;
//        haslow = false;
//      }

//    }

//    fHhits->Fill(fNhitsKept);

//    if (thr_flag) TThread::UnLock();
//  }

void PlaneScifiG::ControlPanel(const TGWindow *p, const TGWindow *main) {

  if (!fControlPanel) fControlPanel = new PlanePanel(p, main, 100, 100, this);

}

void PlaneScifiG::Clusterize() {
  if (!fGeom) return;
  std::set<CCluster1>& clusters = fGeom->Clusterize(fDigits);

  if (thr_flag) TThread::Lock();
  int nclusters=0;
  typedef std::set<CCluster1>::iterator IC;
  for (IC i = clusters.begin(); i!=clusters.end(); i++) {

    float cpos = i->pos;
    float size = i->size;

    if (fVcch->Test(cpos) &&
        fVcs->Test(size)) {

      fHcs->Fill(size);
      fHcch->Fill(cpos);

      fVcch->Store(cpos);
      nclusters++;
    }
    else {
      //update set of clusters
//      clusters.erase(i);
    }
  }
  fNclustKept = nclusters;
  fHchits->Fill(nclusters);

  if (thr_flag) TThread::UnLock();
}





#if USE_DATABASE == 1
void PlaneScifiG::ReadCalib(const tm &t)
{
  // read-in corresponding calibration constants
  try{
    ReadFromDataBase(calib_data,t);
    fUseCalib=0;
    //!!!!
    // fUseCalib=true;
    //for(int i=0; i<calib_data.size();i++)
    // std::cout <<calib_data[i].ch<<" "<<calib_data[i].t0<<" "<<calib_data[i].flag<<std::endl;
    
    //    std::cout<<"PlaneScifiG::ReadCalib() ==> "
    //	<<this->GetName()<<" calibrations are found !"<<std::endl;
    //    std::cout <<"fN"<<fNchan<<"max"<<fMAX_MULT<<std::endl;
    if(calib_data.size() < (unsigned) fNchan){
      //   std::cerr<<"Size of Calibration File is not correct ! Should be : "
      //	  <<fNchan<<" Is "<<calib_data.size()<<" "
      //	  <<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
      //	  <<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec<<std::endl;
    }
    else
      fUseCalib = true;

    //is not ok!!!!!!!
    fUseCalib=true;
  }
  
  catch(CS::Exception& e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(const std::exception &e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(...) {
    std::cout<<"PlaneScifiG::ReadCalib() ==> "<<GetName()
	<<" calibrations, valid for ";
    std::cout<<t.tm_mday<<"."<<t.tm_mon+1<<"."<<t.tm_year+1900<<" "
	<<t.tm_hour<<":"<<t.tm_min<<":"<<t.tm_sec
	<<", not found in DB"<<std::endl;
  }
  
}
#endif //USE_DATABASE

















































































