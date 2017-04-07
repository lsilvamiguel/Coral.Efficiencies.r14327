// -5xx for configuration errors
#define  NO_CSPP_CONFIGFILE -500
#include "CsPP.h"
#include "CsOpt.h"
#include "CsEvent.h"
#include "DaqDataDecoding/Chip.h"
#include "CsErrLog.h"

using namespace std;
using CS::DetID;

// static data members
CsPP* CsPP::_instance = NULL;

/*****************
 ** Constructor **
 *****************/

CsPP::CsPP() {
  // CsPP::CsPP(const string* ppiTableName) {

  if ( _instance ) exit(1);

  _instance = this;
  string ppiTableName("not used yet");
  /*
     if ( !CsOpt::Instance()->getOpt( " " , "CsPP_table",ppiTableName)) { 
     char* message;
     asprintf(&message,"FATAL: No CsPP table specified \n");
     CsErrLog::Instance()->mes(elFatal, message);
     free(message);
     }

     ifstream ppi_file;
     ppi_file.open( ppiTableName.c_str(), ios::in );
     if( !ppi_file.good() )
     {
     char* message;
     asprintf(&message,"FATAL: cannot open CsPP_table:%s \n",ppiTableName.c_str());

     CsErrLog::Instance()->mes(elFatal, message);
     free(message);
     exit(NO_CSPP_CONFIGFILE);
     }

     const int lineSize = 16384; 
     char   line[lineSize];
     string opt;


     do {
     ppi_file.getline( line, lineSize, '\n' );
     if( ppi_file.eof() ) continue;
     if( line[0] == '\0' || line[0] != ' ' ) continue;
#   if __GNUC__ >= 3
istringstream s(line);
#   else
istrstream s(line);
#   endif
s >> opt;
if( ! s ) continue;
   */
  //     if( opt == "ppi" ){
  int id=400;
  string TBname("EC02P1__");
  string name("EC02P1__");
  CsEvent* event = CsEvent::Instance();
  DetID *detid=NULL;
  //   const CS::Chip::Digits& Digs=  event->getChipDigits();
  //   CS::ChipSADC::Digits::const_iterator it;    
  //   for(it = Digs.begin(); it != Digs.end(); ++it) { 
  //     if(it->first.GetName()==name) {
  //       detid =new DetID(it->first);
  //       break;
  //     }
  //   }
  for( map<string,CsDet*>::iterator it=CsDet::GetAllDetectors().begin();
      it!=CsDet::GetAllDetectors().end(); it++ ){
    if(it->second->GetID().GetName()==name) {
      detid =new DetID(it->second->GetID());
      cout<<       it->second->GetID()<<endl<<endl<<"---------------------\n\n\n";

    }
  }




  if(detid==NULL) {
    char* message;
    asprintf(&message,"FATAL: Class CsPP was called, but no detid for ECAL02 was found.\n");
    CsErrLog::Instance()->mes(elFatal, message);
    free(message);
  }




  if(  strncmp( TBname.c_str(), "EC02P1__", 8 ) == 0)
  {
    CsPPI_EC02time* ppi = new CsPPI_EC02time(*detid, TBname);
    _CsPPI.push_back(ppi);
  }
  else {
    char* message;
    asprintf(&message,"FATAL: Class CsPP was called, but no module was set in the file %s\n",ppiTableName.c_str());
    CsErrLog::Instance()->mes(elFatal, message);
    free(message);
  }
  //     }
  //     else {
  //       char* message;
  //       asprintf(&message,"FATAL: Class CsPP was called, but no correct entry in the file %s\n",ppiTableName.c_str());
  //       CsErrLog::Instance()->mes(elFatal, message);
  //     }

  //} while( ppi_file.good() );


};

/****************************************
 ** Instance function to get singelton **
 ****************************************/

CsPP*
CsPP::Instance() {
  // CsPP::Instance(const string* ppiTableName) {

  if ( !_instance ) 
    //     _instance =new  CsPP(ppiTableName);
    _instance =new  CsPP();

  return _instance;

};

/**********************************************
 ** Do the proceessing for the current event **
 **********************************************/

int
CsPP::ProcessEvent() {


  for ( list<CsPPI*>::iterator it = _CsPPI.begin();
      it != _CsPPI.end(); it++) {

    (*it)->process();

  }

  return 0;

};
/**********************************************
 ** End the proceessing of the current data **
 **********************************************/
int
CsPP::end() {


  for ( list<CsPPI*>::iterator it = _CsPPI.begin();
      it != _CsPPI.end(); it++) {

    (*it)->end();

  }

  return 0;

};





/**********************************
 ** Access to the list of CsPPIs **
 **********************************/

list<CsPPI*>&
CsPP::GetPPI() {
  return _CsPPI;
};

