// $Id: CsDet.cc,v 1.37 2010/06/08 19:53:19 suhl Exp $

/*!
   \file    CsDet.cc
   \brief   General properties of any COMPASS detector.
   \version $Revision: 1.37 $
   \author  COMPASS Software Group
   \date    $Date: 2010/06/08 19:53:19 $
*/

#include <cassert>
#include <ctime>

#include "CsInit.h"
#include "CsEvent.h"
#include "CsDigit.h"
#include "CsDet.h"
#include "CsHistograms.h"
#include "CDB.h"
#include "FileDB.h"
#if USE_MySQL
#include "MySQLDB.h"
#endif
#include "CsOpt.h"

#include "DaqDataDecoding/DaqEvent.h"
#include "DaqDataDecoding/ObjectXML.h"
#include "DaqDataDecoding/Chip.h"
#include "DaqDataDecoding/ChipF1.h"

using namespace std;
using CS::DetID;

////////////////////////////////////////////////////////////////////////////////

map<string,CsDet*> CsDet::all_detectors;

////////////////////////////////////////////////////////////////////////////////

CsDet::~CsDet(void)
{
  string name = GetTBName();

  if( name.empty() )
    throw CS::Exception("CsDet::CsDet(): TBname is empty");
  if( all_detectors[name] == this)
  {
    all_detectors.erase(name);
  }
}

////////////////////////////////////////////////////////////////////////////////

CsDet::CsDet(const DetID &id_,const string &_TBname) :
  id(id_),
  TBname(_TBname),
  cdb_(NULL),
  mccdb_(NULL),
  clusteringDone_(false)
{
  string name = GetTBName();

  if( GetTBName().empty() )
    throw CS::Exception("CsDet::CsDet(): TBname is empty");
  if( all_detectors.count(GetTBName())!=0 )
    for( int i=2; i<10000; i++ )
    {
      char s[TBname.length()+11];
      sprintf(s,"%s____%d",TBname.c_str(),i);
      if( all_detectors.count(s)==0 )
      {
        CS::Exception e("CsDet::CsDet():WW: detector number %d with the TBname=\"%s\" name=\"%s\" id=%d",
                         i,GetTBName().c_str(),id.GetName().c_str(),id.GetNumber());
        e.Print();
        name = s;
        break;
      }
    }
  all_detectors[name]=this;

  std::string filedbloc;

  if (CsOpt::Instance()->getOpt(name,"FileDB",filedbloc) ||
      CsOpt::Instance()->getOpt(name.substr(0,2),"FileDB",filedbloc)) {
    cerr<<"Warning on "<<name<<": using private FileDB for calib, location: "<<filedbloc<<endl;
    cdb_ = new FileDB(filedbloc);
  } else {
    cdb_ = CsInit::Instance()->getDB();
  }

  std::string mcfiledbloc;

  if (CsOpt::Instance()->getOpt(name,"MCFileDB",mcfiledbloc) ||
      CsOpt::Instance()->getOpt(name.substr(0,2),"MCFileDB",mcfiledbloc)) {
    cerr<<"Warning on "<<name<<": using private MCFileDB for MC calib, location: "<<mcfiledbloc<<endl;
    mccdb_ = new FileDB(mcfiledbloc);
  }

#if USE_MySQL
  if (CsInit::Instance()->useMySQLDB()) {
    std::string entrytime;
    if (CsOpt::Instance()->getOpt(name,"CDBentrytime",entrytime) ||
        CsOpt::Instance()->getOpt(name.substr(0,2),"CDBentrytime",entrytime)) {

      MySQLDB* mysqldb_ = dynamic_cast<MySQLDB*>(cdb_);
      if (mysqldb_) {
        std::string entrymysqltime = MySQLInterface::toMySQLtime(entrytime.c_str());
        mysqldb_->setSpecificEntryTime(name.c_str(), entrymysqltime.c_str());
      } else {
        std::cerr<<"Warning in CsDet::CsDet: CDBentrytime entry for "<<name<<" works only with MySQLDB"<<endl;
      }
    }
  }
#endif


  string tbname = this->GetTBName();
  folderSet_ = "/COMPASS/";
  for(int i=0;i<2;i++) folderSet_+=tbname[i];folderSet_+="/";
  for(int i=0;i<8;i++) folderSet_+=tbname[i];
}

////////////////////////////////////////////////////////////////////////////////

CsDet *CsDet::FindDetector(const string &name)
{
  for( map<string,CsDet*>::iterator it=all_detectors.begin(); it!=all_detectors.end(); it++ )
    if( it->second->GetTBName() == name )
      return it->second;
  return NULL;
}

////////////////////////////////////////////////////////////////////////////////

void CsDet::Clear(void)
{
  clusteringDone_ = false;
  myClusters_.clear();
  myDigits_.clear();
}

////////////////////////////////////////////////////////////////////////////////

unsigned CsDet::AddMCHitAll(int detector_number,const void *data)
{
  unsigned n=0; // counter of detectors that accepted the hit

  for( map<string,CsDet*>::iterator it=all_detectors.begin(); it!=all_detectors.end(); it++ )
    if( it->second->AddMCHit(detector_number,data) )
      n++;

  return n;
}

////////////////////////////////////////////////////////////////////////////////

#warning first parameter is wrong
bool CsDet::AddMCHit(int detector_number,const void *data)
{
  //cerr << "CsDet::AddMCHit() was not re-implemented for detector " << getName() << endl;
  return false;
}

////////////////////////////////////////////////////////////////////////////////

void CsDet::DecodeChipDigits( const CS::Chip::Digits &digits )
{
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  // get all digits for the detector
  std::pair<m_it,m_it> m_range = digits.equal_range(GetTBName());

  // loop on all found digits
  for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
    DecodeChipDigit(*d_it->second);
}

////////////////////////////////////////////////////////////////////////////////

void CsDet::getAssociatedClusters(list<CsCluster*> &myclus) const
{
  myclus.clear();
  const list<CsCluster*> &allclus = CsEvent::Instance()->getClusters();
  for( list<CsCluster*>::const_iterator Ic=allclus.begin(); Ic!=allclus.end(); Ic++ )
  {
    const list<CsDetector*> &mydets = (*Ic)->getDetsList();
    for( list<CsDetector*>::const_iterator Id=mydets.begin(); Id!=mydets.end(); Id++ )
      if( (*Id)->GetID() == GetID() )
        myclus.push_back( (*Ic) );
  }
}

////////////////////////////////////////////////////////////////////////////////

struct CsDet::sortClusters_ : public binary_function<CsCluster*, CsCluster*, bool>
{
  bool operator() ( CsCluster* c1, CsCluster* c2 )
  {
    if( c1->getU() < c2->getU() )
      return true;
    else
      return false;
  }
};

////////////////////////////////////////////////////////////////////////////////

void CsDet::sortClusters( void )
{
  myClusters_.sort( sortClusters_() );
}

////////////////////////////////////////////////////////////////////////////////

void CsDet::readCalibration(time_t timePoint){

  CDB::Time tp(timePoint,0);
  string data;
  string tmp;
  float calib;
  cdb_->read(GetTBName(),data,tp);
  if(data.size()!=0){
    ;  //   cout << data << endl;
  }else{
    tm *t = localtime(&tp.first);
    cout << GetTBName() << ", no calibration for local time "
	 << t <<" in CDB"<< endl;
  }
}

/////////////////////////////////////////////////////////////////////////////////

void CsDet::readMCCalibration(time_t timePoint)
{
//   CDB::Time tp(timePoint,0);
//   string data;
//   string tmp;
//   float calib;
//   mccdb_->read(GetTBName(),data,tp);
//   if(data.size()!=0){
//     cout << data << endl;
//   }else{
//     tm *t = localtime(&tp.first);
//     cout << GetTBName() << ", no calibration for local time "
// 	 << t <<" in MC CDB"<< endl;
//   }
}

//------------------------------------------------------------------------------


////////////////////////////////////////////////////////////////////////////////
//LS Eff
void CsDet::readMCEffMaps(time_t timePoint){

}
////////////////////////////////////////////////////////////////////////////////

ostream& operator<<(ostream &o, const tm* t) {

  if( t != NULL ) {
    // save current setting of fill character
    char       prevC(o.fill('0'));

    o << setw(2)<<t->tm_mday<<"."
      << setw(2)<<right<<t->tm_mon+1<<"."
      << setw(4)<<t->tm_year+1900<<" ";
    o << setw(2)<<t->tm_hour<<":"
      << setw(2)<<t->tm_min<<":"
      << setw(2)<<t->tm_sec;

    // restore old setting of fill character
    o.fill(prevC);
  }
  return o;
}

