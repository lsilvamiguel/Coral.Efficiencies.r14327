#include "GroupRW.h"
#include "TThread.h"
#include "PlaneRW.h"
#include <TH1F.h>

ClassImp(GroupRichWall);

GroupRichWall::GroupRichWall(const char* name) : Group(name)
{
}

void GroupRichWall::Init() 
{
  Group::Init();
  std::string hMultName = fName + "__mult";
  fHHitMult = new TH1F(hMultName.c_str(),"Hit multiplicity",64,0,64);
  AddHistogram(fHHitMult);
}

#ifndef __CINT__

void GroupRichWall::EndEvent(const CS::DaqEvent &event)
{
  if (thr_flag) TThread::Lock();
  int nHits = 0;
  for(unsigned i=0; i<fPlanes.size(); i++) 
    {
      PlaneRichWall* p = (PlaneRichWall*)(fPlanes[i]);
      int nh = p->GetNhits();
      nHits += nh;
    }
  fHHitMult->Fill(nHits);
  if (thr_flag) TThread::UnLock();
}

#endif
