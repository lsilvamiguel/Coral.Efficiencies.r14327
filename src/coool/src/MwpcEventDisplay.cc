#include <math.h>
#include <stdio.h>
#include <memory.h>
#include <TVirtualPad.h>
#include <TArc.h>
#include <TLine.h>
#include <TROOT.h>
#include <TPad.h>
#include "MwpcEventDisplay.h"
#include "MwpcReconst.h"

float EventDisplay::uslope = UPlane.slope;
float EventDisplay::vslope = VPlane.slope;

const int MAX_N_WIRE = 800;

static float uvshift = tan(UPlane.slope)*XPlane.height/2.; 

void EventDisplay::Fill(MwpcEvent* event)
{
  for(int i = 0; i < 4; i++)
    {
      data->fNhits[i] = event->fNhits[i];
      if(data->fNhits[i] > MAX_N_WIRE)
	data->fNhits[i] = MAX_N_WIRE;
      memcpy(data->fHits[i],event->fHits[i],sizeof(int)*data->fNhits[i]);
    }
  data->fNpoints = event->fNpoints;
  memcpy(data->fX,event->fX,sizeof(float)*data->fNpoints);
  memcpy(data->fY,event->fY,sizeof(float)*data->fNpoints);
}

void EventDisplay::PaintPicture()
{
  TLine line;
  line.SetLineColor(4);
  float cm2unitx = 1./XPlane.width;
  float cm2unity = 1./XPlane.height;
  float sy = cm2unity*YPlane.pitch;
  float sx = cm2unitx*XPlane.pitch;
  float x0,y0,x1,y1;
  // X-plane hits:
  if(data->fNhits[X_PLANE])
    {
      y0 = 0.0;
      y1 = 1.0;
      float x = 0.5+sx*(0.5-XPlane.nWires/2.);
      for(int i = 0; i < data->fNhits[X_PLANE]; i++)
	{
	  x0 = x1 = x+sx*data->fHits[X_PLANE][i];
	  line.PaintLine(x0,y0,x1,y1);
	}
    }
  // Y-plane hits:
  if(data->fNhits[Y_PLANE])
    {
      x0 = 0.0;
      x1 = 1.0;
      float y = 0.5+sy*(0.5-YPlane.nWires/2.);
      for(int i = 0; i < data->fNhits[Y_PLANE]; i++)
	{
	  y0 = y1 = y+sy*data->fHits[Y_PLANE][i];
	  line.PaintLine(x0,y0,x1,y1);
	}
    }
  // U-plane hits:
  if(data->fNhits[U_PLANE])
    {
      y0 = 0.0;
      y1 = 1.0;
      for(int i = 0; i < data->fNhits[U_PLANE]; i++)
	{
	  x0  = x1 = sx*(data->fHits[U_PLANE][i]-UPlane.nWires/2.+.5);
	  x0  += cm2unitx*uvshift;
	  x1  -= cm2unitx*uvshift;
	  line.PaintLine(0.5+x0,y0,0.5+x1,y1);
	}
    }
  // V-plane hits:
  if(data->fNhits[V_PLANE])
    {
      y0 = 0.0;
      y1 = 1.0;
      for(int i = 0; i < data->fNhits[V_PLANE]; i++)
	{
	  x0  = x1 = sx*(data->fHits[V_PLANE][i]-VPlane.nWires/2.+.5);
	  x0  -= cm2unitx*uvshift;
	  x1  += cm2unitx*uvshift;
	  line.PaintLine(0.5+x0,y0,0.5+x1,y1);
	}
    }
  for(int i = 0; i < data->fNpoints; i++)
    {
      float x = data->fX[i]*cm2unitx+0.5;
      float y = data->fY[i]*cm2unity+0.5;
      TArc arc(x,y,0.01);
      arc.SetLineWidth(1);
      arc.SetLineColor(2);
      arc.Paint();
    }
  line.SetLineWidth(3);
  line.SetLineColor(2);
  line.PaintLine(0,0,1,0);
  line.PaintLine(1,0,1,1);
  line.PaintLine(1,1,0,1);
  line.PaintLine(0,1,0,0);
}

EventDisplay::EventDisplay(char* name)
{
  SetName(name);
  data = new MwpcEvent;
  memset(data,0,sizeof(MwpcEvent));
  for(int i = 0; i < 4; i++)
    data->fHits[i] = new int[800];
}

EventDisplay::~EventDisplay() 
{
  for(int i = 0; i < 4; i++)
    delete data->fHits[i];
  delete data;
}

void EventDisplay::Paint(Option_t* option)
{
  ((TPad*)gPad)->Range(0.00,0.00,1.00,1.00);
  PaintPicture();
}
