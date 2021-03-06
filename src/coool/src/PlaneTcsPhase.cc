#include "PlaneTcsPhase.h"

#include "ChipF1.h"
#include "TriggerTime.h"

#include "TThread.h"
#include <vector>
#include <stdexcept>

ClassImp(PlaneTcsPhase)

// implement your member functions here
PlaneTcsPhase::PlaneTcsPhase(const char *detname) : Plane1V(detname, 1, 0, 1000) {
}

void PlaneTcsPhase::Init(TTree *tree) {
    fRateCounter = 0;

    // hits --------------------------------------------------------

    //multiplicity within the time cut on the whole detector
    std::string hitsname = fName + "_hits";
    fHhits=new TH1F_Ref(hitsname.c_str(),hitsname.c_str(),10,0,10,fRateCounter);
    AddHistogram(fHhits);
    std::string hitsleavlist = hitsname + "/I";

    //branch : hit map
    std::string chname = fName + fVch->GetName();
    fHch=new TH1F_Ref(chname.c_str(),chname.c_str(),
                      fVch->GetNbins(),
                      fVch->GetMin(),
                      fVch->GetMax(),
                      fRateCounter);
    AddHistogram(fHch);
    std::string chleavlist = chname + "[" + hitsname +"]/F";

    //branch : time
    std::string tname = fName + fVt->GetName();
    fHt=new TH1F_Ref(tname.c_str(),tname.c_str(),
                     fVt->GetNbins(),
                     fVt->GetMin(),
                     fVt->GetMax(),
                     fRateCounter);
    AddHistogram(fHt);
    std::string tleavlist = tname + "[" + hitsname +"]/F";

    // time vs channel histo
    std::string tcname = fName + "_tVSch";
    fHtvsch=new TH2F(tcname.c_str(),tcname.c_str(),
                     fVch->GetNbins(),
                     fVch->GetMin(),
                     fVch->GetMax(),
                     fVt->GetNbins(),
                     fVt->GetMin(),
                     fVt->GetMax());
    fHtvsch->SetOption("col");
    AddHistogram(fHtvsch);

    // On-Trigger Time Distribution
    std::string onttname = fName + "_t_on_trig";
    fHonTT=new TH1F_Ref(onttname.c_str(),onttname.c_str(),   //onTT is "on Trigger Time"
                        fOnTrigTVarID->GetNbins(),
                        fOnTrigTVarID->GetMin(),
                        fOnTrigTVarID->GetMax(),
                        fRateCounter);
    AddHistogram(fHonTT);

    // Off-Trigger Time Distribution
    std::string offttname = fName + "_t_off_trig";
    fHoffTT=new TH1F(offttname.c_str(),offttname.c_str(),   //offTP is "off Trigger Time"
                     fOffTrigTVarID->GetNbins(),
                     fOffTrigTVarID->GetMin(),
                     fOffTrigTVarID->GetMax());
    AddHistogram(fHoffTT);


    // On-Trigger Profile
    std::string offtcname = fName + "_ch_on_trig";
    fHonTP=new TH1F_Ref(offtcname.c_str(),offtcname.c_str(),   //onTP is "on Trigger Profile"
                        fVch->GetNbins(),
                        fVch->GetMin(),
                        fVch->GetMax(),fRateCounter);
    AddHistogram(fHonTP);

    // Off-Trigger Profile
    std::string offtpname = fName + "_ch_off_trig";
    fHoffTP=new TH1F(offtpname.c_str(),offtpname.c_str(),   //offTP is "off Trigger Profile"
                     fVch->GetNbins(),
                     fVch->GetMin(),
                     fVch->GetMax());
    AddHistogram(fHoffTP);

    // Corrected On-Trigger Profile
    std::string ontcorrname = fName + "_ch_on_corr_trig";
    fHoncorrTP=new TH1F(ontcorrname.c_str(),ontcorrname.c_str(),   
                        fVch->GetNbins(),
                        fVch->GetMin(),
                        fVch->GetMax());
    AddHistogram(fHoncorrTP);

    // Rates HISTOGRAMS
    std::string rname = fName + "_rates";
    fHrates=new TH1F_Ref(rname.c_str(),rname.c_str(),
                         fVch->GetNbins(),
                         fVch->GetMin(),
                         fVch->GetMax(),
                         fRateCounter,true);
    fHrates->SetYTitle("rates per channel (kHz)");
    fHrates->SetTitleOffset(1.2,"Y");
    fHrates->SetOption("hist");
    AddHistogram(fHrates);

    //TCSPhase
    std::string phname = fName + "_TCSphase";
    fTCSPhase=new TH1F(phname.c_str(), "TCS phase from TDC digits", 35, -5, 30);
    AddHistogram(fTCSPhase);

    phname = fName + "_TCSphaseDDD";
    fTCSPhaseDDD=new TH1F(phname.c_str(), "TCS phase from DDD", 35, -5, 30);
    AddHistogram(fTCSPhaseDDD);

    if(tree) {
        fIsInTree = true;
        tree->Branch(hitsname.c_str(),&fNhitsKept,
                     hitsleavlist.c_str(),32000);
        tree->Branch(chname.c_str(),fVch->GetValues(),
                     chleavlist.c_str(),32000);
        tree->Branch(tname.c_str(),fVt->GetValues(),
                     tleavlist.c_str(),32000);
        tree->Branch(fName.c_str(), &fTCSPhaseVal,
                     (fName+"/F").c_str(), 32000);
    }

    if (fReferenceDirectory) {
        ((TH1F_Ref*)fHhits)->SetReference(fReferenceDirectory);
        ((TH1F_Ref*)fHch)->SetReference(fReferenceDirectory);
        ((TH1F_Ref*)fHt)->SetReference(fReferenceDirectory);
        ((TH1F_Ref*)fHonTT)->SetReference(fReferenceDirectory);
        ((TH1F_Ref*)fHonTP)->SetReference(fReferenceDirectory);
        ((TH1F_Ref*)fHrates)->SetReference(fReferenceDirectory);
    }
}

#ifndef __CINT__

std::vector<float> PlaneTcsPhase::GetTCSPhase(const CS::DaqEvent &event) const {
    float min=1000;
    std::vector<float> phase;

    typedef std::list<CS::Chip::Digit*>::const_iterator lDIter;
    for (lDIter ii=lDigits.begin(); ii!=lDigits.end(); ii++) {
        const CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*>(*ii);
        if (!iii) {
            std::cerr<<"PlaneTcsPhase::GetTCSPhase ("<<GetName()<<"): a digit is not a F1 one, strange..."<<std::endl;
            continue;
        }

        int channel = iii->GetChannel();
        if (channel!=0) {
            std::cerr<<"PlaneTcsPhase::GetTCSPhase ("<<GetName()<<"): unexpected channel for TCS phase "<<channel<<std::endl;
            continue;
        }

        double time = iii->GetTimeDecoded();
        phase.push_back(time);
        if (time < min) min = time;
    }

    phase.push_back(min);

    return phase;
}

void PlaneTcsPhase::Reset() {
    Plane1V::Reset();
    fTCSPhaseVal = 0;
}

void PlaneTcsPhase::EndEvent(const CS::DaqEvent &event) {
    if (thr_flag) TThread::Lock();

    typedef std::list<CS::Chip::Digit*>::const_iterator lDIter;
    for (lDIter ii=lDigits.begin(); ii!=lDigits.end(); ii++) {
        CS::ChipF1::Digit* iii = dynamic_cast<CS::ChipF1::Digit*> (*ii);
        if (!iii) {
            std::cerr<<"PlaneTcsPhase::EndEvent ("<<GetName()<<"): a digit is not a F1 one, strange..."<<std::endl;
            continue;
        }

        double time = (double) (iii->GetTimeDecoded() / iii->GetTimeUnit());
        int channel = iii->GetChannel();

        if( fVch->Test(channel) ) { 
            if( fVt->Test(time) ) {
                fVch->Store(channel);
                fVt->Store(time);
                fHch->Fill(channel);
                fHt->Fill(time);
                fHtvsch->Fill(channel,time);
                fNhitsKept++;
            }

            if( fOnTrigTVarID->Test(time) )  {
                fOnTrigTVarID->Store(time);
                fHonTT->Fill(time);        // on trigger time distribution
                fHonTP->Fill(channel);  // on trigger profile
                fNhitsKeptOnTrig++;
            } 

            if( fOffTrigTVarID->Test(time) )  {
                fOffTrigTVarID->Store(time);
                fHoffTT->Fill(time);        // off trigger time distribution
                fHoffTP->Fill(channel);  // off trigger profile
                fNhitsKeptOffTrig++;
            } 
        }
    }
    fHhits->Fill(fNhitsKept);

    if((fRateCounter%Plane1V::fRATE_UPDATE)==0) {
        // corrected on-trigger profile
        if( fOffTrigTVarID->GetMax() != fOffTrigTVarID->GetMin() ) {
            float scale = (float) ( fOnTrigTVarID->GetMax() - fOnTrigTVarID->GetMin() )/( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() );
            for(int i=1; i<=fHonTP->GetNbinsX(); i++)
                fHoncorrTP->SetBinContent(i, (float) fHonTP->GetBinContent(i) - scale*fHoffTP->GetBinContent(i) );
        }  

        // rates histogram
        float timewin = ( fOffTrigTVarID->GetMax() - fOffTrigTVarID->GetMin() )*fF1_TICK*fRateCounter;
        if( timewin > 0 )
            for(int i=1; i<=fHonTP->GetNbinsX(); i++) {
                fHrates->SetBinContent(i,(float) fHoffTP->GetBinContent(i)/timewin);
                fHrates->SetBinError(i,(float) fHoffTP->GetBinError(i)/timewin);
            }
    }

    // fill histogram with TCS phase
    std::vector<float> sphase = GetTCSPhase(event);
    for(unsigned int i=0;i<sphase.size()-1;i++){ 
        fTCSPhase->Fill(sphase[i]);
    } 

    try {
        fTCSPhaseDDD->Fill(event.GetTT().GetPhaseTCS());
        fTCSPhaseVal = event.GetTT().GetPhaseTCS();
    } catch (std::logic_error& e) {
        std::cerr << "PlaneTcsPhase::EndEvent: TCS phase not available: " << e.what() << std::endl;
        fTCSPhaseVal = -1e6;
    }

    if (thr_flag) TThread::UnLock();
}

#endif // __CINT__

