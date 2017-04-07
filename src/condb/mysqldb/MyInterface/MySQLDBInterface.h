
#ifndef __MySQLDBInterface__
#define __MySQLDBInterface__

#include <map>
#include "MySQLInterface.h"

// #define DBNAME "runlb"
// #define DBTABLE "tb_calibDB"


class MySQLDBInterface: public MySQLInterface {


 public:

  /// Destructor
  virtual  ~MySQLDBInterface (void);

  /// Constructor
  MySQLDBInterface (const char* server, const char* username,
                    const char* passwd, const char* dbname);


  /// return entire file path from detector tbname, date and
  /// optional calibration type and DB entry time
  std::string giveFilepath (const char* detname, const char* datetime,
                       const char* typecalib=0, const char* entrytime=0);


  /// add an entry in the calib database (this method does not copy the file)
  int addEntry (const char* detname, unsigned int authorID,
                 unsigned int dirID, const char* filename,
                 const char* starttime, const char* endtime,
                 const char* typecalib=0, const char* creattime=0,
                 const char* dettype=0);

  /// return ID of entries with these particular detname, start and end time, type and entry time
  int checkEntry (const char* detname,
                  const char* starttime, const char* endtime,
                  const char* typecalib=0, const char* creattime=0);

  void updateEntry (int entryID, int dirID, const char* filename);
  std::string getEntryTime(int entryID);
  std::string getCreatTime(int entryID);
  int getUserID(const char* name);
  std::string getDateFromRunNb(int runnb);
  std::string convertNameToTBName(const char* detname) {
           std::string tbname = detname; if(tbname.size()<8) tbname.resize(8, '_'); return tbname; }

  /// return an entry identified by its name from the cfg_directories table
  std::pair<int,std::string> getCfgDirectory(const char* entryname);

  /// gives the reference run file corresponding to a detector type and a environnement(used in coool)
  const char* getRootRefFile(const char* dettype, const char* wkenv = "AFS");
  /// obsolete
  const char* RootRefFile(const char* dettype, const char* wkenv = "AFS") { return getRootRefFile(dettype, wkenv); }

  /// give averaged values of the target polarization offline calculated
  /// the string values of the map are: runnb, uppol, downpol, solenoid, status
  /// Status meanings are: 0 no result, 1 offline result, 2 online result, 3 old online result
  std::map<std::string,double> getTgtOfflinePolar (int runnb);

  /// give averaged values of the target magnets currents for a given run
  /// the first member of the pair is the solenoid one, the second the dipole
  std::pair<double,double> getTgtCurrents (int runnb);

  /// give values of the target magnets currents for a given time (interpolated between 2 measurements)
  /// the first member of the pair is the solenoid one, the second the dipole
  std::pair<double,double> getTgtCurrentsTime (time_t mestime);

  /// give the averaged value over the run of the SM2 NMR measurements
  /// calculated by E. Weise. Return -1 in case of errors or 0 if no measurement found
  double getSM2NMR(int runnb);

  /// give SM1 and 2 current measurements taken at 1st spill
  /// the first member of the pair is SM1, the second is SM2
  /// return 1e9 if no value found
  std::pair<double,double> getSMcurrents(int runnb);

  /// give the path (on shift dir in tbed014d) of the 1st run event file
  char* get1stEvtPath(int runnb);

  /// set and give the keyword to use with tb_directories if AFS is not usable
  void  setSpecialPlace(std::string& spplace)      { fSpecialPlace = spplace; }
  const std::string&   getSpecialPlace(void) const    { return fSpecialPlace; }

protected:
  std::string fSpecialPlace;    // calibration directory are not in AFS but elsewhere, looking at tb_directories workstation fSpecialPlace
  time_t tgtcurtime1, tgtcurtime2;
  double tgtsolen1, tgtsolen2;
  double tgtdipol1, tgtdipol2;


};


#endif
