#include <cassert>
#include <iostream>

#include "BitButton.h"

ClassImp(BitButton)

Bool_t BitButton::HandleButton(Event_t *event) {

  TString label;

  if(event->fType == kButtonPress) {
    
    if(fValue) {
      fValue=0;
      label = "0";
    }
    else {
      fValue=1;
      label = "1";    
    }
  
    SetText(label);
    DoRedraw();
  }
  return TGButton::HandleButton(event);
}

void BitButton::SetValue(int value) {
  
  assert(value==1 || !value);

  fValue = value;
  TString label;
  if(value)
    label="1";
  else 
    label="0";

  SetText(label);
  DoRedraw();  
}
