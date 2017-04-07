#ifndef __MWPC_EVENT_WINDOW_H__
#define __MWPC_EVENT_WINDOW_H__

#include <TH1.h>

class TH1;
class EventDisplay;

class EventDisplay : public TH1
{
  float   width;
  float   height;
  struct MwpcEvent* data;
  static float uslope;
  static float vslope;
  bool   isModified;
 public:
  EventDisplay(char* name);
  ~EventDisplay();
  void Paint(Option_t* option);
  void Fill(MwpcEvent* event);
 private:
  void PaintPicture();
};

#endif
