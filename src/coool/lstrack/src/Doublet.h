#ifndef __Doublet__
#define __Doublet__

#include <string>
#include <vector>

#include "TNode.h"

#include "Structures.h"
#include "Tracker.h"

class Doublet : public TObject {

 private:
  
  std::string fName;
  float fPos;

  std::vector<Tracker*> fTrackers;
  //  TClonesArray *fPoints;
  
  std::vector<Point3*> fPoints;

  //Int_t fNpoints;
  Bool_t InActiveZone(Float_t y, Float_t z);

  static const unsigned int fCnofplanes=2;
  
 public:
  Doublet(const char* name, Tracker* p1, Tracker* p2);
  ~Doublet();
  
  //  void AddTracker(Tracker* plane);
  void Reset();
  const std::vector<Point3*>& Point();

  const std::vector<Point3*>& GetPoints() {return fPoints;}
  std::vector<Tracker*>& GetTrackers() {return fTrackers;}
  float GetPos() const {return fPos;}

  void Draw(TNode *worldnode);
  void Dump();
  
  static unsigned int GetNTrackers() {return Doublet::fCnofplanes;}

  ClassDef(Doublet,0)
};

#endif
