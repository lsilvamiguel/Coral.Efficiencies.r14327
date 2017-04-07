#ifndef CsPixelGEMPlane_H
#define CsPixelGEMPlane_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store/analyze events from COMPASS PixelGEM/Sil
 detectors

 Author: Markus Kr�mer / Thiemo Nagel     March 2007
         Bernhard Ketzer                  June 2008
--------------------------------------------------------------------------

class CsPixelGEMPlane :

                 Object representing a plane of hit strips in one event.
                 The data member fHits is a list of pointers to CsGEMHit
                 objects.
                 The data member fClusters is a list of pointers to
                 CsPixelGEMCluster objects, each containing a list of references
                 to the CsGEMHit objects it is made of.

--------------------------------------------------------------------------

class CsPixelGEMPlaneHeader

                 Header object.
*/

// GEM headers
#include "CsGEMPlaneHeader.h"
#include "CsGEMTimeCal.h"
#include "CsPixelGEMCluster.h"
#include "CsPixelGEMPlanePar.h"

// C++ headers
#include <map>
#include <string>

// ROOT headers
#include "TCanvas.h"
#include "TH2F.h"

//-----------------------------------------------------------------------------
// CsPixelGEMPlane class declaration
//-----------------------------------------------------------------------------
class CsPixelGEMPlane {

private:

  CsGEMPlaneHeader                  fHeader;            //  plane header
  CsPixelGEMPlanePar                fPar;               //  plane parameters
  CsGEMTimeCals                     fTimeCals;          //  time calibrations
  std::map<CsGEMChanId,CsGEMChan*>  fChannels;          //! map of all channels including calibrations and active neighbours
  std::list<CsGEMHit*>              fHits;              //  list with pointers to hits
  std::list<CsPixelGEMCluster*>     fClusters;          //  list with pointers to clusters
  bool                              fIsHitAmpSorted;
  bool                              fIsClusterSorted;
  bool                              fIsClusterized;
  std::string                       fName;              //  Name of plane
  struct                            fCompareHitAmps0;   //  Needed to sort hits
  struct                            fCompareHitAmps1;   //  Needed to sort hits
  struct                            fCompareHitAmps;    //  Needed to sort hits
  TCanvas*                          fDisplayCanvas;     //! ROOT canvas to display single event
  TH2F*                             fDisplayHist;       //! ROOT histogram to display single event
  int                               fKey;               //! keystroke in fDisplayCanvas

public:

  CsPixelGEMPlane();
  CsPixelGEMPlane(std::string name);

  virtual                             ~CsPixelGEMPlane();

  virtual void                        Clear();  // Clear all entries of plane
  void                                ClearChannels();
  void                                ClearHits();
  void                                ClearClusters(); // Clear all clusters
  void                                SetHeader(unsigned int date, unsigned int run, unsigned int evtNrRun, unsigned int evtNrSpill);
  void                                AddChan(int _channr, int _xpos, int _ypos, int _flag, float _ped, float _sigma, int _apvchip=-1, int _apvchan=-1 );
  void                                AddChan(CsGEMChan* _chan) { fChannels[*(_chan->GetId())] = _chan; }
  void                                SetTimeCals(const CsGEMTimeCals& _timecals) { fTimeCals.Clear(); fTimeCals = _timecals; }
  void                                SetTimeCals(const CsGEMTimeCals* _timecals) { fTimeCals.Clear(); fTimeCals = *_timecals; }
  void                                AddHit(int _detchan, int _px, int _py, const std::vector<float>& _amp);
  void                                AddHit(int _detchan, int _px, int _py, float _amp0, float _amp1, float _amp2);
  void                                Clusterize();
  void                                SimpleClustering();
  void                                FullClustering();
  void                                SelectClusters();
   // Add a cluster to plane, default is no flagged pad
  void                                SortHitAmps();      // sort hits by amplitude
  bool                                IsHitAmpSorted() const { return fIsHitAmpSorted; };
  bool                                IsClusterSorted() const { return fIsClusterSorted; };
  bool                                IsClusterized() const { return fIsClusterized; };
  int                                 GetNchannel() const { return fChannels.size(); }
  int                                 GetNhit() const { return fHits.size(); };
  int                                 GetNcluster() const { return fClusters.size(); };
  const std::string&                  GetName() const { return fName; };
  CsGEMPlaneHeader*                   GetHeader() { return &fHeader; };
  CsPixelGEMPlanePar*                 GetPar()            { return &fPar; }
  const CsGEMTimeCals* GetTimeCals() const { return &fTimeCals; }
  const std::map<CsGEMChanId,CsGEMChan*>  &GetChannels() const { return fChannels; }
  const std::list<CsGEMHit*>         &GetHits() const { return fHits; };
  const std::list<CsPixelGEMCluster*> &GetClusters() const { return fClusters; };
  const CsGEMChan* GetChan(const CsGEMChanId& _id)              const;
  void                                Display(bool wait=true);
  void                                GetKey(int event, int key, int keysym, TObject* _selected);
  void                                PrintChannels();
  void                                PrintClusters();
  virtual void                        FindXTalk();

  ClassDef(CsPixelGEMPlane, 3)
};

#endif

