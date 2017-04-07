#ifndef MySQLDB_h
#define MySQLDB_h

#include <sstream>
#include <string>

#include "CDB.h"
#include "MySQLDBInterface.h"

// #define MySQLDBSERVER "pccoeb03"
// #define MySQLDBSERVER "wwwcompass2.cern.ch"
// #define MySQLDBREADER "anonymous"
#define MySQLDBNAME   "runlb"

class MySQLDB : public CDB, public MySQLDBInterface{

 public:

  /// Constructor
  MySQLDB(const char* server, const char* username,
          const char* passwd, const char* dbname = MySQLDBNAME);

  /// Destructor
  ~MySQLDB()    { DisconnectDB(); };

  void read(const std::string &folder, std::string &data, Time tPoint, const char* keyword=0);

  bool   ConnectDB(void);
  void   DisconnectDB(void)     { disconnect(); }
  bool   isConnected(void)      { return MySQLInterface::isConnected(); }

  void   setEntryTime(const char* entrytime)  { fEntrytime = entrytime; }
  void   setEntryTime(std::string& entrytime)      { fEntrytime = entrytime; }
  const std::string&   getEntryTime(void) const    { return fEntrytime; }

  void setSpecificEntryTime(const char* tbname, const char* entrytime)
    { fSpecificEntrytime.insert(std::map<std::string, std::string>::value_type(tbname, entrytime));
      fSpecETfg = true; }
  const std::map<std::string, std::string>& getSpecificEntryTime(void)
    { return fSpecificEntrytime; }

  /// gives the reference run file corresponding to a detector type (used in coool)
  const char* getRootRefFile(const char* dettype, const char* wkenv = "AFS")
                { return MySQLDBInterface::getRootRefFile(dettype, wkenv); }
  /// obsolete
  const char* RootRefFile(const char* dettype, const char* wkenv = "AFS")
                { return MySQLDBInterface::getRootRefFile(dettype, wkenv); }
  /// give averaged values of the target polarization offline calculated
  /// the string values of the map are: runnb, uppol, downpol, solenoid
  std::map<std::string,double> getTgtOfflinePolar (int runnb)
                { return MySQLDBInterface::getTgtOfflinePolar(runnb); }

  /// give the averaged value over the run of the SM2 NMR measurements
  /// calculated by E. Weise. Return -1 in case of errors or 0 if no measurement found
  double getSM2NMR(int runnb)
                { return MySQLDBInterface::getSM2NMR(runnb); }

  /// give the path (on shift dir in tbed014d) of the 1st run event file
  char* get1stEvtPath(int runnb)
                { return MySQLDBInterface::get1stEvtPath(runnb); }


 protected:

  std::string decodeTBNameFromFolder(const std::string& folder);

  std::string fEntrytime;    //!< entries must have entry time before fEntrytime (MySQL format)

  std::map<std::string, std::string> fSpecificEntrytime;  // entry times for specific TBNames
  bool fSpecETfg;  // true if fSpecificEntrytime is not empty

 private:

  /// Don't use that contructor
  MySQLDB(const MySQLDB &);

};





//------------------------------------------------------------------------------


#endif
