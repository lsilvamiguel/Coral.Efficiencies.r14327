#ifndef __Station__
#define __Station__

#include <string>
#include <vector>

#include "TObject.h"
#include "TH2.h"
#include "TNode.h"

#include "Structures.h"
#include "Doublet.h"

class Station  : public TObject {

 private:

  std::string fName;
  
  static const unsigned int fNdoublets;
  std::vector<Doublet*> fDoublets;
  std::vector<Point3*> fPoints;

  float fPos;
  TH2F *fProfile;
  
 public:
  Station(const char* name,Doublet* d1,Doublet* d2);
  
  const std::vector<Point3*>& Coinc3(int plane);
  
  void Draw(TNode *worldnode);  
  void Reset();
  TH2F* GetProfile() {return fProfile;}
  
  std::vector<Doublet*>& GetDoublets() {return fDoublets;}
  
  static unsigned int GetNDoublets() {return Station::fNdoublets;}

  ClassDef(Station,0)
};



#endif







