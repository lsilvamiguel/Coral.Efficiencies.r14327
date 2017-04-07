#ifndef __Variable__
#define __Variable__

#include <string>
#include "TObject.h"

class Variable : public TObject {
  
 private:
  std::string fName;
  int fNbins;
  float fMin;
  float fMax;
  float* fValues;
  int fNvalues;
  
 public:
  Variable(const char* name, int nbins, float min, float max, int maxsize)
    : TObject(),fName(name), fNbins(nbins), fMin(min), fMax(max), 
    fValues(new float[maxsize]), fNvalues(0) {}
  
  virtual ~Variable() {delete[] fValues;}
  
  void SetRange(float min, float max) {
    fMin=min;
    fMax=max;
  }

  void Reset() {fNvalues=0;}
  void MyDump();
  
  bool Test(float value) const {
    if(fMin<=value && value<=fMax) return true;
    else return false;
  }

  void Store(float value) {fValues[fNvalues]=value; fNvalues++;}

  float* GetValues() {return fValues;}    
  std::string& GetName() {return fName;}
  int   GetNbins() const {return fNbins;}
  float GetMin() const {return fMin;}
  float GetMax() const {return fMax;}
  int GetNvalues() const {return fNvalues;} 

  void  SetNbins(int nbins) { fNbins = nbins; }


  ClassDef(Variable,0)
};


#endif
