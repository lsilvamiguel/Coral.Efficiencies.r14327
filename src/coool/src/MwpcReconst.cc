#include "MwpcReconst.h"
#include <math.h>
#include <memory.h>

const float SLOPE  = 0.17700547497;
const float PITCH  = 0.19687507788;
const float XPITCH = 0.2;
const float YPITCH = 0.2;
const float WIDTH  = 152.0;
const float HEIGHT = 120.0;
const int   MAX_DRIFT_TIME = 120; 

PlaneGeometry XPlane = { X_PLANE, 760, XPITCH,   0.0,    WIDTH, HEIGHT };
PlaneGeometry UPlane = { U_PLANE, 760, PITCH,  SLOPE,    WIDTH, HEIGHT };
PlaneGeometry VPlane = { V_PLANE, 760, PITCH, -SLOPE,    WIDTH, HEIGHT };
PlaneGeometry YPlane = { Y_PLANE, 520, YPITCH, M_PI/2.,  WIDTH, HEIGHT };

const float MAX_ERROR = 2.; // Maximum reconstruction error is 5 wires

//*********************** class CAbstractArray ***************************//

void CAbstractArray::AddObject(CObject* object)
{
//   CObject** tmp = new (CObject*)[size+1];
//   for(int i = 0; i < size; i++) tmp[i] = array[i];
//   tmp[size] = object;
//   size++;
//   delete array;
//   array = tmp;
  array.push_back(object);
  size = array.size();
  isEmpty = false;
}

void CAbstractArray::Destroy()
{
  for(int i = 0; i < size; i++) 
    delete array[i];
  ClearArray();
}

void CAbstractArray::ClearArray()
{
//   delete array;
  array.clear();
  size = 0;
//   array = NULL;
  isEmpty = true;
}

void CAbstractArray::RemoveObject(int n)
{
  if(n >= size) return;
  delete array[n];

//   CObject** tmp = new (CObject*)[size-1];
//   int j = 0;
//   for(int i = 0; i < size; i++) 
//     {
//       if(i == n) continue;
//       tmp[j++] = array[i];
//     }
//   size--;
//   delete array;
//   array = tmp;
//   if(size == 0) isEmpty = true;

  array.erase(array.begin() + n);
  size = array.size();
  if (size == 0) isEmpty = true;
}

//*************************** CWire class ****************************//

float CWire::xcenter = XPlane.nWires/2.-0.5;
float CWire::ycenter = YPlane.nWires/2.-0.5;
float CWire::ucenter = UPlane.nWires/2.-0.5;
float CWire::vcenter = VPlane.nWires/2.-0.5;
float CWire::uvclaim = HEIGHT*tan(SLOPE);

float cosalpha = cos(SLOPE);
float sinalpha = sin(SLOPE);
float tanalpha = tan(SLOPE);
float uvshift  = HEIGHT*tanalpha/XPITCH;

void CWire::SetWire(int planeType, int _number)
{
  plane = planeType;
  number = _number;
  ClearStatistic();
}

CIntersection* CWire::GetIntersection(CWire* wire)
{
  if(wire->plane == plane) return NULL;
  CWire* X = NULL;
  CWire* Y = NULL;
  CWire* U = NULL;
  CWire* V = NULL;
  switch(plane)
    {
    case X_PLANE:
      X = this;
      break;
    case Y_PLANE:
      Y = this;
      break;
    case U_PLANE:
      U = this;
      break;
    case V_PLANE:
      V = this;
      break;
    }
  switch(wire->plane)
    {
    case X_PLANE:
      X = wire;
      break;
    case Y_PLANE:
      Y = wire;
      break;
    case U_PLANE:
      U = wire;
      break;
    case V_PLANE:
      V = wire;
      break;
    }
  float x = 0;
  float y = 0;
  if(X && Y)
    {
      x = (X->number-xcenter)*XPITCH;
      y = (Y->number-ycenter)*YPITCH;
    }
  else if(X && V)
    {
      x = (X->number-xcenter)*XPITCH;
      y = (X->number-V->number)*XPITCH/tanalpha;
    }
  else if(X && U)
    {
      x = (X->number-xcenter)*XPITCH;
      y = (U->number-X->number)*XPITCH/tanalpha;
    }
  else if(U && V)
    {
      float diff = (U->number - V->number)*XPITCH;
      x = ((U->number + V->number)/2.-xcenter)*XPITCH;
      y = diff/(2.*tanalpha);
    }
  else if(U && Y)
    {
      y = (Y->number-ycenter)*YPITCH;
      x = (U->number-xcenter)*XPITCH-y*tanalpha;
    }
  else if(V && Y)
    {
      y = (Y->number-ycenter)*YPITCH;
      x = (V->number-xcenter)*XPITCH+y*tanalpha;
    }
//    if(x < -WIDTH/2. || x > WIDTH/2. || y < -HEIGHT/2. || y > HEIGHT/2.)
//      return NULL;
  CIntersection* ret = new CIntersection;
  ret->SetX(x);
  ret->SetY(y);
  ret->wire[0] = this;
  ret->wire[1] = wire;
  return ret;
}

void CWire::ClearStatistic()
{
  nTrig = nXtalk = nNoise = nEff = nIneff = 0.;
}

bool CWire::HaveIntersection(CWire* wire)
{
  CWire* X = NULL;
  CWire* Y = NULL;
  CWire* U = NULL;
  CWire* V = NULL;
  switch(plane)
    {
    case X_PLANE:
      X = this;
      break;
    case Y_PLANE:
      Y = this;
      break;
    case U_PLANE:
      U = this;
      break;
    case V_PLANE:
      V = this;
      break;
    }
  switch(wire->plane)
    {
    case X_PLANE:
      X = wire;
      break;
    case Y_PLANE:
      Y = wire;
      break;
    case U_PLANE:
      U = wire;
      break;
    case V_PLANE:
      V = wire;
      break;
    }
  if((X || U || V) && Y) return true;
  if(U && V)
    {
      if(fabs(U->GetNumber()-V->GetNumber()) < uvshift)
	return true;
    }
  if(U && X)
    {
      if(fabs(U->GetNumber()-X->GetNumber()) < uvshift/2.)
	return true;
    }
  if(V && X)
    {
      if(fabs(V->GetNumber()-X->GetNumber()) < uvshift/2.)
	return true;
    }
  return false;
}

bool CWire::IsMultipleIntersected(CWireArray* array, int nIntersections)
{
  if(!array || array->IsEmpty()) return false;
  if(IsTriggered()) return false;
  int nInt = 0;
  for(int i = 0; i < array->GetSize(); i++)
    {
      if(((*array)[i])->IsTriggered()) continue;
      if(HaveIntersection((*array)[i]))
	{
	  nInt++;
	  if(nInt >= nIntersections)
	    return true;
	}
    }
  return false;
}


//*************************** CTrPoint class ******************************//

CTrPoint::CTrPoint(float _x, float _y) : x(_x), y(_y), error(-1.0)
{
}

float CTrPoint::GetDistance(CWire* w)
{
  float n = w->GetNumber();
  if(w->plane == X_PLANE)
    return fabs((n-CWire::xcenter)*XPITCH-x);
  else if(w->plane == Y_PLANE)
    return fabs((n-CWire::ycenter)*YPITCH-y);
  else if(w->plane == U_PLANE)
    {
      return fabs(y*sinalpha-((n-CWire::ucenter)*XPITCH-x)*cosalpha);
    }
  else if(w->plane == V_PLANE)
    {
      return fabs(y*sinalpha-((n-CWire::vcenter)*XPITCH-x)*cosalpha);
    }
  return 0.;
}

//*************************** CMwpcCluster class ******************************//

CMwpcCluster::CMwpcCluster(int clusterSize, CWire* hits)
{
  size = clusterSize;
//  cluster = new (CWire*)[size];
  cluster = (CWire**) malloc(size * sizeof(CWire*));
  float n = 0;
  for(register int i = 0; i < size; i++)
    {
      n += hits[i].GetNumber();
      cluster[i] = &hits[i];
    }
  SetWire(hits[0].plane,(int)(n/size));
}

void CMwpcCluster::TriggeredEvent()
{
  for(int i = 0; i < size; i++)
    cluster[i]->TriggeredEvent();
  CWire::TriggeredEvent();
}

void CMwpcCluster::EffectiveEvent(float weight)
{
  for(int i = 0; i < size; i++)
    cluster[i]->EffectiveEvent(weight/size);
  CWire::EffectiveEvent(weight);
}

void CMwpcCluster::NoiseEvent(float weight)
{
  for(int i = 0; i < size; i++)
    cluster[i]->NoiseEvent(weight);
  CWire::NoiseEvent(weight);
}

//************************ class CIntersection   ** **********************//

float CIntersection::GetX()
// returns "virtual" number of "hot" wire of X-plane
{
  return (float(XPlane.nWires-1)/2. + x/XPITCH);
}

float CIntersection::GetY()
// returns "virtual" number of "hot" wire of Y-plane
{
  return (float(YPlane.nWires-1)/2. + y/YPITCH);
}

float CIntersection::GetU()
// returns "virtual" number of "hot" wire of U-plane
{
  float u = x*cosalpha+y*sinalpha;
  return (float(UPlane.nWires-1)/2. + u/PITCH);
}

float CIntersection::GetV()
// returns "virtual" number of "hot" wire of V-plane
{
  float v = x*cosalpha-y*sinalpha;
  return (float(VPlane.nWires-1)/2. + v/PITCH);
}

//************************ class CIntersectionArray **********************//

void CIntersectionArray::operator+=(CIntersectionArray* arr)
{
//   if(!arr) return;
//   CObject** tmp = new (CObject*)[size+arr->size];
//   memcpy((void*)tmp,(void*)array,size*sizeof(void*));
//   memcpy((void*)(&tmp[size]),(void*)arr->array,arr->size*sizeof(void*));
//   size += arr->size;
//   delete array;
//   array = tmp;

  array.reserve(size + arr->GetSize());
  for (register int ii = 0; ii < arr->GetSize(); ii++) {
    CObject* tmpcobj = arr->GetObject(ii);
    array.push_back(tmpcobj);
  }
  size = array.size();
}

//************************ class CTrPointArray  ***************************//

void CTrPointArray::Print(FILE* outstream)
{
  fprintf(outstream,"The array of track points:\n");
  for(int i = 0; i < size; i++)
    {
      CTrPoint* p = (CTrPoint*)GetObject(i);
      p->Print(outstream);
    }
  fprintf(outstream,"\n");
}
//******************************* CWireArray class ************************//

CIntersectionArray* CWireArray::GetIntersections(CWireArray* wires)
{
  if(wires->GetSize() == 0 || GetSize() == 0) return NULL;
  CIntersectionArray* array = new CIntersectionArray;
  for(int i = 0; i < GetSize(); i++)
    {
      CWire* wi = (*this)[i];
      for(int j = 0; j < wires->GetSize(); j++)
	{
	  CWire* wj = (*wires)[j];
	  CIntersection* intersection = wi->GetIntersection(wj);
	  if(intersection)
	    array->AddObject(intersection);
	}
    }
  if(!array->IsEmpty())
    return array;
  else
    {
      delete array;
      return NULL;
    }
}

void CWireArray::Print(FILE* outstream)
{
  if(!size)
    {
      fprintf(outstream,"Empty\n");
      return;
    }
  fprintf(outstream,"%d  -  ",(*this)[0]->GetPlane());
  for(int i = 0; i < size; i++)
    {
      fprintf(outstream,"%.2f ",(*this)[i]->GetNumber());
    }
  fprintf(outstream,"\n");
}

//************************* CPlane class *********************************//

CPlane::CPlane(PlaneGeometry& geometry) : planeType(geometry.planeType)
{
  nWires = geometry.nWires;
  wires = new CWire[nWires];
  for(int i = 0; i < nWires; i++)
    {
      wires[i].SetWire(planeType,i);
    }
  hits = new int[nWires];
  hitTime = new int[nWires];
  nReconstrEvents = 0;
  ClearHits();
}

CPlane::~CPlane()
{
  delete[] wires;
  delete hits;
}

void CPlane::PutHit(int nWire, int time)
{ 
  if(nWire < 0 || nWire >= nWires) return;
  hits[nWire]++;
  hitTime[nWire] = time;
  if(minNum > nWire) minNum = nWire;
  if(maxNum < nWire) maxNum = nWire;
  nHits++;
}

void CPlane::ClearHits()
{ 
  memset(hits,0,sizeof(int)*nWires);
  memset(hitTime,0,sizeof(int)*nWires);
  nHits  = 0;
  minNum = 0x7FFF;
  maxNum = -1;
}

void CPlane::ClearStatistic()
{
  for(int i = 0; i < nWires; i++)
    {
      wires[i].ClearStatistic();
    }
  nReconstrEvents = 0;
}

CMwpcClusterArray* CPlane::GetClusters()
{
  CMwpcClusterArray* array = new CMwpcClusterArray;
  if(!nHits) return array;
  int clusterStart = -1;
  int clusterSize  =  0;
  int minTime = 0;
  int maxTime = -100000;
  for(int i = minNum; i <= maxNum; i++)
    {
      if(hits[i])
	{
	  if(clusterSize == 0)
	    clusterStart = i;
	  if(hitTime[i] < minTime)
	    minTime = hitTime[i];
	  if(hitTime[i] > maxTime)
	    maxTime = hitTime[i];
	  if(maxTime-minTime < MAX_DRIFT_TIME)
	    clusterSize++;
	  else if(clusterSize)
	    {
	      array->AddObject(new CMwpcCluster(clusterSize,&wires[clusterStart]));
	      clusterSize = 1;
	      clusterStart = i;
	      minTime = maxTime = hitTime[i];
	    }
	}
      else
	{
	  if(clusterSize)
	    {
	      array->AddObject(new CMwpcCluster(clusterSize,&wires[clusterStart]));
	      clusterSize = 0;
	      minTime = 0;
	      maxTime = -100000;
	    }
	}
    }
  if(clusterSize)
    {
      array->AddObject(new CMwpcCluster(clusterSize,&wires[clusterStart]));
    }
  return array;
}

bool CPlane::UneffectiveCluster(float n)
{
  if(n < 0 || n >= nWires) return false;
  if(n < 0 || n >= nWires) return false;
  int n0 = (int)floor(n);
  int n1 = (int)ceil(n);
  if(n0==n1 || n1 >= nWires)
    {
      wires[n0].TriggeredEvent();
      wires[n0].IneffectiveEvent(1.);
    }
  else
    {
      wires[n0].TriggeredEvent();
      wires[n1].TriggeredEvent();
      wires[n0].IneffectiveEvent(0.5);
      wires[n1].IneffectiveEvent(0.5);
    }
  return true;
}

float CPlane::GetEfficiency()
{
  float eff = 0.;
  float ntr = 0.;
  for(int i = 0; i < nWires; i++)
    {
      float e = wires[i].GetNEff();
      float a = e + wires[i].GetNIneff();
      if(a > 0.01)
	{
	  eff += e/a;
	  ntr += 1.;
	}
    }
  if(ntr)
    return eff/ntr;
  else
    return 0.;
}


float CPlane::GetNoise()
{
  float noise = 0.;
  for(int i = 0; i < nWires; i++)
    {
      noise += wires[i].GetNNoise();
    }
  return noise;
}

int CPlane::GetHits(int* data)
{
  int nHits = 0;
  for(int i = minNum; i <= maxNum; i++)
    {
      if(hits[i])
	{
	  data[nHits] = i;
	  nHits++;
	}
    }
  return nHits;
}

//**************************** CChamber class *****************************//

CChamber::CChamber(PlaneGeometry* config[4])
{
  nReconstrEvents = 0;
  fReconData = new MwpcEvent;
  memset(fReconData,0,sizeof(MwpcEvent));
  for(int i = 0; i < 4; i ++)
    {
      if(config[i])
	{
	  planes[i] = new CPlane(*config[i]);
	  fReconData->fHits[i] = new int[config[i]->nWires];
	}
      else
	planes[i] = NULL;
    }
}

CChamber::~CChamber()
{
  for(int i = 0; i < 4; i++)
    {
      delete planes[i];
      delete fReconData->fHits[i];
    }
  delete fReconData;
}

void CChamber::ClearStatistic()
{
  for(int i = 0; i < 4; i++)
    {
      if(planes[i]) planes[i]->ClearStatistic();
    }
  nReconstrEvents = 0;
}

void CChamber::ClearHits()
{
  for(int i = 0; i < 4; i++)
    {
      if(planes[i])
	planes[i]->ClearHits();
    }
}

CTrPointArray* CChamber::GetTrackPoints()
{
  CWireArray* uhits = planes[U_PLANE]->GetClusters();
  CWireArray* vhits = planes[V_PLANE]->GetClusters();
  CWireArray* xhits = planes[X_PLANE]->GetClusters();
  CWireArray* yhits = NULL;
  if(planes[Y_PLANE])
    yhits = planes[Y_PLANE]->GetClusters();
  CTrPointArray* array = new CTrPointArray;
  // Rejection of ambiguos events:
  for(int i = 0; i <  uhits->GetSize(); i++)
    {
      CMwpcCluster* ucluster = (CMwpcCluster*)(*uhits)[i];
      if(ucluster->IsMultipleIntersected(vhits) || ucluster->IsMultipleIntersected(xhits))
	{
	  goto STOP;
	}
    }
  for(int i = 0; i <  vhits->GetSize(); i++)
    {
      CMwpcCluster* vcluster = (CMwpcCluster*)(*vhits)[i];
      if(vcluster->IsMultipleIntersected(uhits) || vcluster->IsMultipleIntersected(xhits))
	{
	  goto STOP;
	}
    }
  for(int i = 0; i <  xhits->GetSize(); i++)
    {
      CMwpcCluster* xcluster = (CMwpcCluster*)(*xhits)[i];
      if(xcluster->IsMultipleIntersected(uhits) || xcluster->IsMultipleIntersected(vhits))
	{
	  goto STOP;
	}
    }
  // Track points reconstruction block:
  {
    planes[U_PLANE]->nReconstrEvents++;
    planes[V_PLANE]->nReconstrEvents++;
    planes[X_PLANE]->nReconstrEvents++;
    if(planes[Y_PLANE]) planes[Y_PLANE]->nReconstrEvents++;
    // Looking for 3-wires track points:
    for(int i = 0; i < uhits->GetSize(); i++)
      {
	CMwpcCluster* u = (CMwpcCluster*)(*uhits)[i];
	if(u->IsTriggered()) continue;
	for(int j = 0; j < vhits->GetSize(); j++)
	  {
	    CMwpcCluster* v = (CMwpcCluster*)(*vhits)[j];
	    if(v->IsTriggered()) continue;
	    for(int k = 0; k < xhits->GetSize(); k++)
	      {
		CMwpcCluster* x = (CMwpcCluster*)(*xhits)[k];
		if(x->IsTriggered()) continue;
		if(FindTrackPoint(u,v,x,array))
		  {
		    u->TriggeredEvent();
		    v->TriggeredEvent();
		    x->TriggeredEvent();
		    u->EffectiveEvent(1);
		    v->EffectiveEvent(1);
		    x->EffectiveEvent(1);
		  }
		else if(IsTriplet(u,v,x)) // Reject "ambigues" triplet intersections
		  {
		    u->TriggeredEvent();
		    v->TriggeredEvent();
		    x->TriggeredEvent();
		  }
	      }
	  }
      }
    // Looking for 2-wires track points:
    // U and V wires intersections:
    for(int i = 0; i < uhits->GetSize(); i++)
      {
	CMwpcCluster* u = (CMwpcCluster*)(*uhits)[i];
	if(u->IsTriggered()) continue;
	for(int j = 0; j < vhits->GetSize(); j++)
	  {
	    CMwpcCluster* v = (CMwpcCluster*)(*vhits)[j];
	    if(v->IsTriggered()) continue;
	    if(u->HaveIntersection(v))
	      {
		CIntersection* intersec = u->GetIntersection(v);
		u->EffectiveEvent(1);
		v->EffectiveEvent(1);
		u->TriggeredEvent();
		v->TriggeredEvent();
		planes[X_PLANE]->UneffectiveCluster(intersec->GetX());
		array->AddObject(new CTrPoint(intersec->x,intersec->y));
		delete intersec;
	      }
	  }
      }
    // U and X wires intersections:
    for(int i = 0; i < uhits->GetSize(); i++)
      {
	CMwpcCluster* u = (CMwpcCluster*)(*uhits)[i];
	if(u->IsTriggered()) continue;
	for(int j = 0; j < xhits->GetSize(); j++)
	  {
	    CMwpcCluster* x = (CMwpcCluster*)(*xhits)[j];
	    if(x->IsTriggered()) continue;
	    if(u->HaveIntersection(x))
	      {
		CIntersection* intersec = u->GetIntersection(x);
		u->EffectiveEvent(1);
		x->EffectiveEvent(1);
		u->TriggeredEvent();
		x->TriggeredEvent();
		planes[V_PLANE]->UneffectiveCluster(intersec->GetV());
		array->AddObject(new CTrPoint(intersec->x,intersec->y));
		delete intersec;
	      }
	  }
      }
    // V and X wires intersections:
    for(int i = 0; i < vhits->GetSize(); i++)
      {
	CMwpcCluster* v = (CMwpcCluster*)(*vhits)[i];
	if(v->IsTriggered()) continue;
	for(int j = 0; j < xhits->GetSize(); j++)
	  {
	    CMwpcCluster* x = (CMwpcCluster*)(*xhits)[j];
	    if(x->IsTriggered()) continue;
	    if(v->HaveIntersection(x))
	      {
		CIntersection* intersec = v->GetIntersection(x);
		v->EffectiveEvent(1);
		x->EffectiveEvent(1);
		v->TriggeredEvent();
		x->TriggeredEvent();
		planes[U_PLANE]->UneffectiveCluster(intersec->GetU());
		array->AddObject(new CTrPoint(intersec->x,intersec->y));
		delete intersec;
	      }
	  }
      }
    // Looking for noising channels:
    // U-wires:
    for(int i = 0; i < uhits->GetSize(); i++)
      {
	CMwpcCluster* u = (CMwpcCluster*)(*uhits)[i];
	if(u->IsTriggered()) continue;
	u->TriggeredEvent();
	u->NoiseEvent(1);
      }
    // V-wires:
    for(int i = 0; i < vhits->GetSize(); i++)
      {
	CMwpcCluster* v = (CMwpcCluster*)(*vhits)[i];
	if(v->IsTriggered()) continue;
	v->TriggeredEvent();
	v->NoiseEvent(1);
      }
    // X-wires:
    for(int i = 0; i < xhits->GetSize(); i++)
      {
	CMwpcCluster* x = (CMwpcCluster*)(*xhits)[i];
	if(x->IsTriggered()) continue;
	x->TriggeredEvent();
	x->NoiseEvent(1);
      }
  }
 STOP:
  delete uhits;
  delete vhits;
  delete xhits;
  delete yhits;
  return array;
}

bool CChamber::IsTriplet(CMwpcCluster* u, CMwpcCluster* v, CMwpcCluster* x)
{
  if(u->IsTriggered() || v->IsTriggered() || x->IsTriggered())
    return false;
  if((u->HaveIntersection(x) && u->HaveIntersection(v)) ||
     (v->HaveIntersection(x) && v->HaveIntersection(u)) ||
     (x->HaveIntersection(u) && x->HaveIntersection(v)))
    return true;
  return false;
}

bool CChamber::FindTrackPoint(CMwpcCluster* u, CMwpcCluster* v, CMwpcCluster* x, 
		      CTrPointArray* array)
{
  if(!IsTriplet(u,v,x)) return false;
  float minerr = MAX_ERROR+1.;
  int nu = 0;
  int nv = 0;
  int nx = 0;
  for(int i = -1; i < u->GetSize(); i++)
    {
      CWire* wu;
      if(i == -1)
	wu = (CWire*) u;
      else
	wu = (*u)[i];
      for(int j = -1; j < v->GetSize(); j++)
	{
	  CWire* wv;
	  if(j == -1)
	    wv = (CWire*) v;
	  else
	    wv = (*v)[j];
	  for(int k = -1; k < x->GetSize(); k++)
	    {
	      CWire* wx;
	      if(k == -1)
		wx = (CWire*) x;
	      else
		wx = (*x)[k];
	      if(wv->HaveIntersection(wu) && wu->HaveIntersection(wx) &&
		 wv->HaveIntersection(wx))
		{
		  float nx = wx->GetNumber();
		  float nu = wu->GetNumber();
		  float nv = wv->GetNumber();
		  float tt = fabs(nx-(nu+nv)/2.);
		  if(tt < minerr)
		    {
		      minerr = tt;
		      nu = i;
		      nv = j;
		      nx = k;
		    }
		}
	    }
	}
    }
  if(minerr <  MAX_ERROR)
    {
      CWire* xw = NULL;
      CWire* uw = NULL;
      CWire* vw = NULL;
      if(nx == -1)
	xw = x;
      else
	xw = (CWire*)(*x)[nx];
      if(nu == -1)
	uw = u;
      else
	uw = (CWire*)(*u)[nu];
      if(nv == -1)
	vw = v;
      else
	vw = (CWire*)(*v)[nv];
      CIntersection* int_vx = vw->GetIntersection(xw);
      CIntersection* int_ux = uw->GetIntersection(xw);
      CIntersection* int_uv = uw->GetIntersection(vw);
      float x0 = (int_vx->x+int_uv->x+int_ux->x)/3.;
      float y0 = (int_vx->y+int_uv->y+int_ux->y)/3.;
      CTrPoint* p = new CTrPoint(x0,y0);
      p->error = minerr;
      array->AddObject(p);
      delete int_vx;
      delete int_uv;
      delete int_ux;
      return true;
    }
  return false;
}

MwpcEvent* CChamber::GetReconPoints()
{
  CTrPointArray* parr = GetTrackPoints();
  if(parr && !parr->IsEmpty())
    {
      fReconData->fNpoints = 
	((parr->GetSize() < MAX_POINTS_NUM) ? parr->GetSize() : MAX_POINTS_NUM);
      for(int i = 0; i < fReconData->fNpoints; i++)
	{
	  CTrPoint* p = (CTrPoint*)parr->GetObject(i);
	  fReconData->fX[i] = p->GetX();
	  fReconData->fY[i] = p->GetY(); 
	  fReconData->fErr[i] = p->GetError();
	}
    }
  else
    {
      fReconData->fNpoints = 0;
    }
  delete parr;
  return fReconData;
}

MwpcEvent* CChamber::GetHits()
{
  for(int i = 0; i < 4; i++)
    {
      if(planes[i])
	fReconData->fNhits[i] = planes[i]->GetHits(fReconData->fHits[i]);
    }
  return fReconData;
}

bool CChamber::CalculateEfficiency()
{
  float eff[4];
  bool retVal = true;
  memset(efficiency,0,sizeof(efficiency));
  for(int i = 0; i < 4; i++)
    {
      if(!planes[i]) continue;
      if(planes[i]->GetNRecEvents())
	{
	  eff[i] = planes[i]->GetEfficiency();
	  if(eff[i] == 0. && i != Y_PLANE)
	    {
	      retVal = false;
	      break;
	    }
	}
      else if(!planes[i]->GetNRecEvents())
	{
	  retVal = false;
	  break;
	}
    }
  if(retVal)
    {
      float a = eff[X_PLANE]+eff[U_PLANE]+eff[V_PLANE]-2.;
      float b = eff[X_PLANE]+eff[U_PLANE]-1.;
      float c = eff[X_PLANE]+eff[V_PLANE]-1.;
      float d = eff[U_PLANE]+eff[V_PLANE]-1.;
      if(a <= 0. || b <= 0. || c <= 0. || d <= 0.) retVal = false;
      else
	{
	  efficiency[X_PLANE] = a/d;
	  efficiency[U_PLANE] = a/c;
	  efficiency[V_PLANE] = a/b;
	}
    }
  return retVal;
}
