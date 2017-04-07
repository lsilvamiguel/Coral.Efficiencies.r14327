#include "CsPPI.h"
#include <sstream>
#include <cstdio>

using CS::DetID;

/***********************
 ** CsPPImsg section  **
 ***********************/


CsPPImsg::CsPPImsg(const error_level level, const int Code,
                   const std::string& what, const std::string& file,
                   const unsigned int line) :
  fLevel(level),fErrCode(Code),fWhat(what),fFile(file),fLine(line)
{
  
  time(&fWhen);

};

CsPPImsg::~CsPPImsg() {};

error_level
CsPPImsg::Level() const { return fLevel; };

int
CsPPImsg::errcode() const { return fErrCode; };

std::string
CsPPImsg::what() const {

  const char *level;
  switch( fLevel ) {
    case DEBUG3 :
        level = "[DEBUG3]";
        break;
    case DEBUG2 :
        level = "[DEBUG2]";
        break;
    case DEBUG1 :
        level = "[DEBUG1]";
        break;
    case WARNING : 
         level = "[WARNING]";
         break;
    case ERROR :
         level = "[ERROR]";
         break;
    case FATAL :
         level = "[FATAL]";
         break;
    case INFO :
         level = "[INF0]";
         break;
    default :
         level = "[UKNOWN]";
         break;
  }
  std::stringstream ret;
  char _time[33];
  strftime(_time, 33,"[%c (GMT)]",gmtime(&fWhen));
  ret << "CsPPI::MSG " << level << " on " << _time << " : in '" << fFile << "'@line'" << fLine << "' : " << fWhat << "\n";


  return ret.str();

};

/*******************
 ** CsPPI section **
 *******************/

CsPPI::CsPPI(const DetID &id, const std::string &TBname, const error_level thr) :
  fId(id),fTBName(TBname), fMsgThr(thr)
{};


CsPPI::~CsPPI() {
  
  // clear MSG without output
  ClearMSG(1);

};

const DetID&
CsPPI::GetID() const { return fId; };

const char*
CsPPI::getName() const { return fId.GetName().c_str(); };

const std::string&
CsPPI::GetTBName() const { return fTBName; };

const std::vector<CsPPImsg*>&
CsPPI::GetMSG() const { return fMSG; };

unsigned int 
CsPPI::GetNBfatal() const { return fNBfatal; };

unsigned int 
CsPPI::GetNBerror() const { return fNBerror; };

unsigned int 
CsPPI::GetNBwarning() const { return fNBwarning; };

const CsPPImsg*
CsPPI::GetLastError(int code) const {

  // Are there errors at all?
  if ( !fNBerror && !fNBfatal ) 
    return NULL;

  // Backward loop over fMSG
  for ( int i = fMSG.size() - 1; i >= 0; --i) {

    // Check if MSG is error or fatal
    if ( fMSG[i]->Level() >= ERROR ) {

      // Check error code
      if ( code == 0 || code == fMSG[i]->errcode() ) {
        return fMSG[i];
      }

    }
    
  }

  // empty list
  return NULL;

};

const CsPPImsg*
CsPPI::GetLastFatal() const {

  // Are there fatals at all?
  if ( !fNBfatal ) 
    return NULL;

  // Backward loop over fMSG
  for ( int i = fMSG.size() - 1; i >= 0; --i) {

    // Check if MSG is error or fatal
    if ( fMSG[i]->Level() == FATAL ) {

        return fMSG[i];

    }
    
  }

  // empty list
  return NULL;

}; 

const CsPPImsg*
CsPPI::GetLastMSG(int thr) const {
  
  // Backward loop over fMSG
  for ( int i = fMSG.size() - 1; i >= 0; --i) {
    
    // Check if MSG is error or fatal
    if ( fMSG[i]->Level() >= thr ) {

        return fMSG[i];

    }
    
  }

  return NULL;

};

void
CsPPI::ClearMSG(int thr) {

  if( !fMSG.size() )
    return;

  std::cout << "****************************************************************************************************\n"
            << "CsPPI Message Summary:\n";
  while ( fMSG.size() ) {

    std::vector<CsPPImsg*>::iterator it = fMSG.begin();
    if ( *it ) {
      if ( (*it)->Level() >= thr )
        std::cout << (*it)->what();
      delete (*it);
      *it = NULL;
    }
    fMSG.erase(it);
    
  }

  std::cout << "****************************************************************************************************\n";

  fNBfatal=0;
  fNBerror=0;
  fNBwarning=0;


};

void
CsPPI::MSG(CsPPImsg *MSG) {
  
  if( MSG->Level() < fMsgThr ) {
      delete MSG;
      return;
  }

  std::cout << MSG->what();

  fMSG.push_back(MSG);
  
  if ( MSG->errcode() == ERROR )
    fNBerror++;
  if ( MSG->errcode() == WARNING )
    fNBwarning++;
  if ( MSG->errcode() == FATAL ) {
    fNBfatal++;
  }

};


void
CsPPI::MSG(const error_level level, const int Code,
           const std::string& what, const std::string& file,
           const unsigned int line) {

  CsPPImsg *msg = new CsPPImsg(level,Code,what,file,line);
  MSG(msg);  

};

int
CsPPI::process() {

  //This prototype function should never be called
  MSG(FATAL, -1, "Called CsPPI::process(), which should be overloaded", __FILE__, __LINE__);

  return -1;

};

int
CsPPI::end() {


  //This prototype function should never be called
  MSG(FATAL, -1, "Called CsPPI::end(), which should be overloaded", __FILE__, __LINE__);

  return -1;
}


