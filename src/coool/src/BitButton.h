#ifndef __BitButton__
#define __BitButton__

#include "TGButton.h"

class BitButton : public TGTextButton {

 private:
  int fValue;

 public:
  BitButton(const TGWindow* p,int id=-1) : 
    TGTextButton(p,"0",id), fValue(0) {}

  Bool_t HandleButton(Event_t *event); 
  
  int GetValue() const {return fValue;}
  void SetValue(int value);

  ClassDef(BitButton,0)

};


#endif
