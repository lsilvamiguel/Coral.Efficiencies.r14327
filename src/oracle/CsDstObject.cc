#include <memory.h>
#include "CsDstObject.h"

using namespace std;
using namespace CLHEP;

bool CsDstObject::Save(ofstream& file) const
{
  vector<uint8> buffer;
  GetBuffer(buffer);
  file.write((const char*)&buffer[0],buffer.size());
  if(file.fail()) return false;
  return true;
}

bool CsDstObject::Save(FlatFile& file) const
{
  vector<uint8> buffer;
  GetBuffer(buffer);
  if(!file.Write((void*)&buffer[0],buffer.size()))
    return false;
  return true;
}

CsDst3Vector::CsDst3Vector()
{
  type = VECTOR3_TYPE;
  memset(vec,0,sizeof(vec));
}

void CsDst3Vector::Init(const float64* v)
{
  memcpy(vec,v,sizeof(vec));
}

void CsDst3Vector::Init(const Cs3Vector& v)
{
  Cs3Vector& vv = const_cast<Cs3Vector&>(v);
  for(int i = 0; i < 3; i++)
    vec[i] = vv(i);
}

void CsDst3Vector::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer, &vec[0], sizeof(vec));  
}

bool CsDst3Vector::Load(std::istream& file)
{
  file.read((char*)&vec[0], sizeof(vec));
  if(file.fail()) return false;
  return true;    
}

bool CsDst3Vector::Load(FlatFile& file)
{
  if(!file.Read((char*)&vec[0], sizeof(vec)))
    return false;
  return true;    
}

bool CsDst3Vector::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDst3Vector* v = (const CsDst3Vector*)obj;
  for(int i = 0; i < 3; i++)
    {
      if(vec[i] != v->vec[i]) return false;
    }
  return true;
}

void CsDst3Vector::Print(std::ostream& msgout) const
{
  char MSG_STR[256];
  sprintf(MSG_STR,"\t(%.2f; %.2f; %.2f)\n", vec[0],vec[1],vec[2]);
  msgout << MSG_STR;
}

CsDst3Matrix::CsDst3Matrix()
{
  type = MATRIX3_TYPE;
  memset(cov,0,sizeof(cov));
}

void CsDst3Matrix::Init(const float32* m)
{
  memcpy(cov,m,sizeof(cov));
}

void CsDst3Matrix::Init(const HepMatrix& m)
{
  for(int i = 0; i < 9; i++)
    {
      cov[i] = m(1+i/3, 1+i%3);
    }
}

HepMatrix* CsDst3Matrix::GetHepMatrix() const
{
  HepMatrix* m = new HepMatrix(3,3);
  for(int i = 0; i < 9; i++)
    (*m)(1+i/3, 1+i%3) = cov[i];
  return m;
}

void CsDst3Matrix::GetBuffer(std::vector<uint8>& buffer) const
{
  AddValueToVector(buffer, &cov[0], sizeof(cov));
}

bool CsDst3Matrix::Load(std::istream& file)
{
  file.read((char*)&cov[0], sizeof(cov));
  if(file.fail()) return false;
  return true;    
}

bool CsDst3Matrix::Load(FlatFile& file)
{
  if(!file.Read((char*)&cov[0], sizeof(cov)))
    return false;
  return true;    
}

bool CsDst3Matrix::IsEquil(const CsDstObject* obj) const
{
  if(type != obj->type) return false;
  const CsDst3Matrix* m = (const CsDst3Matrix*)obj;
  for(int i = 0; i < 9; i++)
    if(cov[i] != m->cov[i]) return false;
  return true;
}

void CsDst3Matrix::Print(std::ostream& msgout) const
{
  char MSG_STR[256];
  for(int i = 0; i < 3; i++)
    {
      sprintf(MSG_STR,"\t%15f %15f %15f\n", cov[3*i], cov[3*i+1], cov[3*i+2]);
      msgout << MSG_STR;
    }
}

CsDstEvent::CsDstEvent()
{
  eventNumber = eventInBurst = burstNumber = runNumber =
    timeInSec = timeInUSec = errorCode = triggerMask = 0;
  isExist = false;
}

bool CsDstRun::Load(std::istream&) { return false; }
bool CsDstRun::Load(FlatFile&) { return false; }
bool CsDstRun::IsEquil(const CsDstObject*) const { return false; }
void CsDstRun::GetBuffer(std::vector<uint8>&) const {;}
