#ifndef CsGEMPlane_H
#define CsGEMPlane_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store/analyze events from COMPASS GEM
 detectors

--------------------------------------------------------------------------

class CsGEMPlane :

                 Object representing a plane of hits in one event.
                 The data member fHits is a list of pointers to CsGEMHit
                 objects.
                 The data member fClusters is a list of pointers to
                 CsGEMCluster objects, each containing a list of references
                 to the CsGEMHit objects it is made of.

--------------------------------------------------------------------------

v1.0     19/11/2002    by    Bernhard Ketzer

v2.0     27/05/2008    by    Bernhard Ketzer
--------------------------------------------------------------------------
*/

// GEM headers
#include "CsGEMChan.h"
#include "CsGEMCluster.h"
#include "CsGEMHit.h"
#include "CsGEMPlaneHeader.h"
#include "CsGEMPlanePar.h"
#include "CsGEMTimeCal.h"

// C++ headers
#include <map>
#include <string>

// ROOT headers
#include "TCanvas.h"
#include "TH1F.h"

//-----------------------------------------------------------------------------
// CsGEMPlane class declaration
//-----------------------------------------------------------------------------
class CsGEMPlane {

 private:

    CsGEMPlaneHeader   fHeader;                  //  plane header
    CsGEMPlanePar      fPar;                     //  plane parameters
    CsGEMTimeCals      fTimeCals;                //  time calibrations
    std::map<CsGEMChanId,CsGEMChan*> fChannels;  //! map of all channels including calibrations and active neighbours
    std::list<CsGEMHit*>             fHits;      //  list with pointers to hits
    std::list<CsGEMCluster*>         fClusters;  //  list with pointers to clusters
    bool               fIsHitSorted;             //  list of hits sorted by channel Id
    bool               fIsHitAmpSorted;          //  list of hits sorted by amplitude
    bool               fIsClusterSorted;         //  cluster sorted by channel Id
    bool               fIsClusterized;
    std::string        fName;                    //  Name of plane
    struct             fCompareHits;             //  to compare list of pointers
    struct             fCompareHitAmps0;         //  to compare list of pointers
    struct             fCompareHitAmps1;         //  to compare list of pointers
    struct             fCompareHitAmps2;         //  to compare list of pointers
    struct             fCompareClusters;         //  to compare list of pointers
    TCanvas*           fDisplayCanvas;           //! ROOT canvas to display single event
    TH1F*              fDisplayHist;             //! ROOT historam to display single event
    int                fKey;                     //! keystroke in fDisplayCanvas

 public:

                        CsGEMPlane();
                        CsGEMPlane(std::string name);

    virtual            ~CsGEMPlane();

    void                Clear();
    void                ClearChannels();
    void                ClearHits();
    void                ClearClusters();
    void                SetHeader(unsigned int date, unsigned int run, unsigned int evtNrRun, unsigned int evtNrSpill);
    void                AddChan(int _detchan, int _hem, int _flag, float _ped, float _sigma, int _apvchip=-1, int _apvchan=-1 );
    void                AddChan(CsGEMChan* _chan) { fChannels[*(_chan->GetId())] = _chan; }
    void                SetTimeCals(const CsGEMTimeCals& _timecals) { fTimeCals.Clear(); fTimeCals = _timecals; }
    void                SetTimeCals(const CsGEMTimeCals* _timecals) { fTimeCals.Clear(); fTimeCals = *_timecals; }
    void                AddHit(int _detchan, int _hem, float _amp0, float _amp1, float _amp3);
    void                AddHit(int _detchan, int _hem, const std::vector<float>& _amp);
    void                Clusterize();
    void                SimpleClustering();
    void                FullClustering();
    void                FindXTalk();
    void                SelectClusters();
    void                SortHits();         // sort hits by channel Id
    void                SortHitAmps();      // sort hits by amplitude
    void                SortClusters();
    bool                IsHitSorted() const { return fIsHitSorted; }
    bool                IsClusterSorted() const { return fIsClusterSorted; }
    bool                IsClusterized() const { return fIsClusterized; }
    int                 GetNchannel() const { return fChannels.size(); }
    int                 GetNhit() const { return fHits.size(); }
    int                 GetNcluster() const { return fClusters.size(); }
    const std::string&  GetName() const     { return fName; }
    CsGEMPlaneHeader*   GetHeader()         { return &fHeader; }
    CsGEMPlanePar*      GetPar()            { return &fPar; }
    const CsGEMTimeCals* GetTimeCals() const { return &fTimeCals; }
    const std::map<CsGEMChanId,CsGEMChan*>  &GetChannels() const { return fChannels; }
    const std::list<CsGEMHit*>              &GetHits()     const { return fHits; }
    const std::list<CsGEMCluster*>          &GetClusters() const { return fClusters; }
    const CsGEMChan* GetChan(const CsGEMChanId& _id)              const;
    void          Display(bool wait=true);
    void          GetKey(int _event, int _key, int _keysym, TObject* _selected);
    void          PrintChannels();
    void          PrintClusters();

    ClassDef(CsGEMPlane, 3)
};

#endif

