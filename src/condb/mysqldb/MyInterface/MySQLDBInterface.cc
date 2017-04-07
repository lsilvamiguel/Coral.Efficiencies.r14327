#include <unistd.h>
#include <stdio.h>
#include <pwd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

#include <cstdlib>

#include "MySQLDBInterface.h"

#define MAXONLPOLARTIMEDIFF 3600

using namespace std;

namespace {
  char qrystr[1024];
}

//----------------------------------------------------------------------------

MySQLDBInterface::MySQLDBInterface (const char* server, const char* username,
                                    const char* password, const char* dbname)
  : MySQLInterface (server, username, password, dbname),
    fSpecialPlace("AFS"), tgtcurtime1(0), tgtcurtime2(0),
    tgtsolen1(0), tgtsolen2(0),
    tgtdipol1(0), tgtdipol2(0)
  { }


//----------------------------------------------------------------------------

MySQLDBInterface::~MySQLDBInterface() { }


//----------------------------------------------------------------------------

string MySQLDBInterface::giveFilepath
                    (const char* detname, const char* datetime,
                     const char* typecalib, const char* entrytime) {
  char strtmp1[500], strtmp2[400], strtmp3[400];

  string tbname = convertNameToTBName(detname);

  if (typecalib && (string(typecalib) != "")) {
      sprintf(strtmp2, " and (typecalib = '%s') ", typecalib); }
    else { sprintf(strtmp2, " and (typecalib = 'default') "); }
  if (entrytime && (string(entrytime) != ""))
      { sprintf(strtmp3, " and (entrytime <= '%s') ", entrytime); }
    else strtmp3[0] = 0;
  if (datetime) { sprintf(strtmp1, " and (starttime <= '%s') and (endtime >= '%s') ",
                          datetime, datetime); }
    else strtmp1[0] = 0;

  sprintf(qrystr, "select filename,Directory from tb_calibDB left join "
                  "tb_directories on tb_calibDB.dirID=tb_directories.ID "
                  "where (detname = '%s' ) %s %s %s "
                  "and (tb_directories.workstation = '%s') "
                  "order by creattime desc limit 1",
                   tbname.c_str(), strtmp1, strtmp2, strtmp3,
                   getSpecialPlace().c_str());

  if ( query(qrystr) && getNRows() ) {
    string res(getCol(1));
    res += "/";
    res += getCol(0);
    endQuery();
    return res;
  }

  cerr << "MySQLDBInterface::giveFilepath() failed, no answer found for "
       << detname << ", type " << ( typecalib!=NULL ? typecalib : "default" )
       << " at " << datetime << endl;

  endQuery();
  return "";
}


//----------------------------------------------------------------------------


int MySQLDBInterface::addEntry
                    (const char* detname, unsigned int authorID,
                     unsigned int dirID, const char* filename,
                     const char* starttime, const char* endtime,
                     const char* typecalib, const char* creattime,
                     const char* dettype) {

  char dettypestr[3];
  const char* deftypecalib = "default";
  string creattime_2;
  unsigned int authentID;

  if ((typecalib == 0) || (typecalib[0] == 0)) { typecalib = deftypecalib; }

  if (creattime == 0) {
    creattime_2 = "NOW()";
  } else {
    creattime_2 = string("'") + toMySQLtime(creattime) + "'";
  }

  if (dettype == 0) {
    dettypestr[0] = detname[0];
    dettypestr[1] = detname[1];
    dettypestr[2] = 0;
    dettype = dettypestr;
  }

  authentID = getUserID(getenv("LOGNAME"));
  if (authentID == 0) { authentID = 8; }  // 8 correspond to unknown in tb_authors table

  sprintf(qrystr, "insert into tb_calibDB (creattime, dettype, detname, "
                  "typecalib, authorID, dirID, filename, starttime, endtime, "
                  "authentID, entrytime) "
                  "values (%s, '%s', '%s', '%s', %d, %d, '%s', '%s', '%s', "
                  "%d, NOW() ) ",
                   creattime_2.c_str(), dettype, detname, typecalib, authorID,
                   dirID, filename, toMySQLtime(starttime).c_str(),
                   toMySQLtime(endtime).c_str(), authentID);

// fprintf(stderr, "addEntry: starttime %s endtime %s\n", starttime, endtime);
// fprintf(stderr, "addEntry: mysqltime starttime %s endtime %s\n",
//         toMySQLtime(starttime).c_str(), toMySQLtime(endtime).c_str());

  if (query(qrystr)) {
    endQuery();
    return getLastInsertID();
  }

  endQuery();
  return -1;
}


//----------------------------------------------------------------------------


int MySQLDBInterface::checkEntry
                    (const char* detname,
                     const char* starttime, const char* endtime,
                     const char* typecalib, const char* creattime) {

  const char* deftypecalib = "default";
  string creattime_2;

  if ((typecalib == 0) || (typecalib[0] == 0)) { typecalib = deftypecalib; }

  if ((creattime == 0) || (creattime[0] == 0)) {
    creattime_2 = "";
  } else {
    creattime_2 = string("and (creattime = '") + toMySQLtime(creattime) + "')";
  }


  sprintf(qrystr, "select calibID from tb_calibDB "
                  "where (detname = '%s') and (starttime = '%s') "
                  "  and (endtime = '%s') and (typecalib = '%s') %s ",
                   detname, toMySQLtime(starttime).c_str(),
                   toMySQLtime(endtime).c_str(), typecalib, creattime_2.c_str());


  if (query(qrystr)) {
    if (getNRows()) {
      endQuery();
      return atoi(getCol(0));
    }
  }

  endQuery();
  return 0;
}


//----------------------------------------------------------------------------

void MySQLDBInterface::updateEntry (int entryID, int dirID,
                                    const char* filename) {

  sprintf(qrystr, "update tb_calibDB set dirID=%d,filename='%s' "
                  "where calibID=%d limit 1" ,
                   dirID, filename, entryID);

  query(qrystr);
  endQuery();
}


//----------------------------------------------------------------------------


string MySQLDBInterface::getEntryTime(int entryID) {

  sprintf(qrystr, "select entrytime from tb_calibDB where calibID=%d limit 1",
          entryID);

  if (query(qrystr)) {
    if (getNRows()) {
      endQuery();
      return toCDBtime(getCol(0));
    }
  }

  cerr <<"MySQLDBInterface::getEntryTime: can't find entry ID "<<entryID<<std::endl;
  return "";
}




//----------------------------------------------------------------------------


string MySQLDBInterface::getCreatTime(int entryID) {

  sprintf(qrystr, "select creattime from tb_calibDB where calibID=%d limit 1",
          entryID);

  if (query(qrystr)) {
    if (getNRows()) {
      endQuery();
      return toCDBtime(getCol(0));
    }
  }

  cerr <<"MySQLDBInterface::getCreatTime: can't find entry ID "<<entryID<<std::endl;
  return "";
}






//----------------------------------------------------------------------------

int MySQLDBInterface::getUserID(const char* name) {

  if (name == 0) {
    cerr <<"MySQLDBInterface::getUserID: no user name given\n";
    return -1;
  }
  sprintf(qrystr, "select ID from tb_authors where Login='%s' limit 1", name);
  if (query(qrystr)) {
    if (getNRows()) {
      endQuery();
      return atoi(getCol(0));
    }
  }
  endQuery();

  struct passwd *pwname = 0;
  pwname = getpwnam(name);
  if (pwname) {
    char *firstn, *lastn, *ctmp;
//     firstn = strdup(pwname->pw_gecos);
    lastn = firstn = (char*) malloc(strlen(pwname->pw_gecos) + 1);
    strcpy(firstn, pwname->pw_gecos);
    if ((ctmp = strchr(firstn, ' '))) {
      ctmp[0] = 0;
      lastn = ctmp + 1;
    }
    if ((ctmp = strrchr(lastn, ' '))) lastn = ctmp + 1;
    sprintf(qrystr, "insert into tb_authors (login,FirstName,LastName) "
                    "values ('%s','%s','%s')", name, firstn, lastn);
    free(firstn);
  } else
  {
    sprintf(qrystr, "insert into tb_authors (login) values ('%s')", name);
  }

  query(qrystr);
  endQuery();
  return getLastInsertID();
}

//----------------------------------------------------------------------------

std::pair<int,std::string> MySQLDBInterface::getCfgDirectory(const char* entryname) {

  std::pair<int,std::string> result;
  result.first = -1;
  result.second = "";

  sprintf(qrystr, "select dirid,path from cfg_directories where name='%s' limit 1", entryname);
  if (query(qrystr)) {
    if (getNRows()) {
      endQuery();
      result.first = atoi(getCol(0));
      result.second = getCol(1);
      return result;
    }
  }
  endQuery();

  cerr <<"Warning: MySQLDBInterface::getCfgDirectory: entry name "<<entryname
       <<" not found in cfg_directories table\n";
  return result;
}


//----------------------------------------------------------------------------

string MySQLDBInterface::getDateFromRunNb(int runnb) {

  sprintf(qrystr, "select starttime from tb_run where runnb='%d' limit 1", runnb);
  if (query(qrystr)) {
    if (getNRows()) {
      endQuery();
      return toCDBtime(getCol(0));
    }
  }
  endQuery();

  cerr <<"MySQLDBInterface::getDateFromRunNb: run "<<runnb<<" not found in logbook database\n";
  return "";
}


//----------------------------------------------------------------------------

const char* MySQLDBInterface::getRootRefFile(const char* dettype, const char* wkenv) {

  if (!dettype) {
    cerr << "MySQLDBInterface::getRootRefFile: warning, dettype char* is NULL\n";
    return 0;
  }
  string sqlreq;
  sqlreq = "select runnb,filename from tb_cooolref where tb_cooolref.dettype like '";
  sqlreq += dettype;
  sqlreq += "'";
  if(!query(sqlreq.c_str())) {
    cout << " No SQL database answer about reference root file for detector type "<<dettype<<std::endl;
    return 0;
  }
  if (!getCol(0)) {
    cout << " No reference root file registered for detector type "<<dettype<<" in SQL database\n";
    cout << "   taking default reference file\n";
    return 0;
  }

  string runnb(getCol(0));
  string filename = "";
  if (getCol(1)) filename = getCol(1);
  endQuery();

  sqlreq = "select dirID,Directory from tb_directories where tb_directories.workstation='";
  sqlreq += wkenv;
  sqlreq += "' and tb_directories.key='reffiles' limit 1";
  if(!query(sqlreq.c_str())) {
    cout << " No SQL database answer about reference root file directory for work enviromment "<<wkenv<<std::endl;
    endQuery();
    return 0;
  }
  string dirname = "/afs/cern.ch/compass/detector/monitor/References/";
  string dirID;
  if (!getCol(0)) {
    cout << " No reference root file directory registered for work environment "<<wkenv<<" in SQL database\n";
    cout << "   taking default directory\n";
  } else {
    dirID = getCol(0);
    if (getCol(1)) dirname = getCol(1);
  }
  endQuery();

  string name = "";
  name += dirname;
  name += "/";
  if ( filename != "" ) {
    name += filename;
  } else {
    name += "coool_";
    name += runnb;
    name += ".root";
    struct stat statbuf;
    if ( ! stat(name.c_str(), &statbuf) ) {
      // stat succedeed
      if ( S_ISREG(statbuf.st_mode) ) {
        // regular file, we use it instead the file on web server (faster !)
        return name.c_str();
      }
    }

//     name = "http://pccoeb03.cern.ch/rootfile.php/runnb=";
    sqlreq = "select dirID,Directory from tb_directories where tb_directories.workstation='";
    sqlreq += wkenv;
    sqlreq += "' and tb_directories.key='refrunweb' limit 1";
    if(!query(sqlreq.c_str())) {
      cout << " No SQL database answer about reference run web page for work enviromment "<<wkenv<<std::endl;
      endQuery();
      return 0;
    }
    name = "http://wwwcompass.cern.ch/rootfile.php/runnb=";
    if (!getCol(0)) {
      cout << " No reference run web page registered for work environment "<<wkenv<<" in SQL database\n";
      cout << "   taking default web page on wwwcompass\n";
    } else {
      dirID = getCol(0);
      if (getCol(1)) name = getCol(1);
    }
    endQuery();
    name += runnb;
  }

  return name.c_str();
}



//----------------------------------------------------------------------------


std::map<std::string,double> MySQLDBInterface::getTgtOfflinePolar(int runnb) {
  std::map<std::string,double> result;

  result["status"] = 0;
  sprintf(qrystr, "select starttime,stoptime from tb_run where runnb=%d", runnb);
  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getTgtOfflinePolar: can't query SQL command "<<qrystr<<std::endl;
    return result;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getTgtOfflinePolar: can't access to run info for run "<<runnb<<std::endl;
    return result;
  }

  if (getCol(1)) {
    sprintf(qrystr, "select runnb,avg(tb_offlinepolar.Ups_D_polar), "
                    "  avg(tb_offlinepolar.Dwn_D_polar), "
                    "  avg(tb_offlinepolar.Centr_D_polar), "
		    "  avg(tb_offlinepolar.Solenoid_current) "
                    "from tb_run left join tb_offlinepolar "
                    "  on (tb_offlinepolar.Timing between tb_run.starttime and tb_run.stoptime) "
                    "where runnb=%d group by runnb", runnb);
  } else {
    cerr<<"\nWarning in MySQLDBInterface::getTgtOfflinePolar: No end time for run "<<runnb<<std::endl;
    cerr<<"     30 minutes length taken by default \n"<<std::endl;
    sprintf(qrystr, "select runnb,avg(tb_offlinepolar.Ups_D_polar), "
                    "  avg(tb_offlinepolar.Dwn_D_polar), "
                    "  avg(tb_offlinepolar.Centr_D_polar), "
                    "  avg(tb_offlinepolar.Solenoid_current) "
                    "from tb_run left join tb_offlinepolar "
                    "  on (tb_offlinepolar.Timing between tb_run.starttime and DATE_ADD(tb_run.starttime,INTERVAL 30 MINUTE)) "
                    "where runnb=%d group by runnb", runnb);
  }
  endQuery();

  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getTgtOfflinePolar: can't query SQL command "<<qrystr<<std::endl;
    return result;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getTgtOfflinePolar: can't access to offline polar info for run "<<runnb<<std::endl;
    return result;
  }

  result["runnb"] = atof(getCol(0));
  if (getCol(1)) {
    result["uppol"] = atof(getCol(1));
    result["downpol"] = atof(getCol(2));
    result["centrpol"] = atof(getCol(3));
    result["solenoid"] = atof(getCol(4));
    result["status"] = 1;
    endQuery();
    return result;
  }

  // no offline polar data, searching for online one
  cerr<<"\nWarning in MySQLDBInterface::getTgtOfflinePolar: No offline polar values found for run "<<runnb<<std::endl;
  cerr<<"     giving online ones \n"<<std::endl;
  sprintf(qrystr, "select runnb,tb_target.uppol,tb_target.downpol, "
                  "  tb_target.centralpol,tb_target.solencur, "
                  "  UNIX_TIMESTAMP(starttime)-UNIX_TIMESTAMP(magtime) "
                  "from tb_run left join tb_target "
                  "  on (tb_run.targetid=tb_target.targetid) "
                  "where runnb=%d", runnb);

  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getTgtOfflinePolar: can't query SQL command "<<qrystr<<std::endl;
    return result;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getTgtOfflinePolar: can't access to online polar info for run "<<runnb<<std::endl;
    return result;
  }
  result["runnb"] = atof(getCol(0));
  if (getCol(1)) {
    float timediff = atof(getCol(5));
    result["status"] = 2;
    if (timediff > MAXONLPOLARTIMEDIFF) {
      cerr<<"\nWARNING in MySQLDBInterface::getTgtOfflinePolar: the online polar found is very old for run "<<runnb<<std::endl;
      cerr<<" difference between polar measurement and start of run is "<<timediff<<" seconds !!"<<std::endl;
      cerr<<"     returning that measurement anyway, but be careful...\n"<<std::endl;
      result["status"] = 3;
    }

    // looks good, returning online values
    cerr<<"\nMySQLDBInterface::getTgtOfflinePolar: online polar found, returning it for run "<<runnb<<std::endl;
    result["uppol"] = atof(getCol(1));
    result["downpol"] = atof(getCol(2));
    result["centrpol"] = atof(getCol(3));
    result["solenoid"] = atof(getCol(4));
    endQuery();
    return result;
  }

  // nothing found, returning null values
  endQuery();
  cerr<<"\nERROR in MySQLDBInterface::getTgtOfflinePolar: No online polar values found for run "<<runnb<<std::endl;
  cerr<<"     returning null values !!!!!!! \n"<<std::endl;
  result["uppol"] = 0;
  result["downpol"] = 0;
  result["centrpol"] = 0;
  result["solenoid"] = 0;
  result["status"] = 0;
  return result;
}



//----------------------------------------------------------------------------


std::pair<double,double> MySQLDBInterface::getTgtCurrents(int runnb) {
  std::pair<double,double> result;

  result.first = 0;  // return 0 if no value found
  result.second = 0;

  if ( runnb < 32833 ) {  // runs before 1/1/2004
    cerr<<"Warning: MySQLDBInterface::getTgtCurrents not suited for runs before 2004"<<endl;
    cerr<<"     ... using getTgtOfflinePolar instead (no dipole info then)"<<endl;
    std::map<std::string,double> tgtoff = getTgtOfflinePolar(runnb);
    if ( tgtoff.size() > 0 ) {
      result.first = tgtoff["solenoid"];
      result.second = 0; // no dipole information
    }
    return result;
  }


  sprintf(qrystr, "select starttime,stoptime from tb_run where runnb=%d", runnb);
  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getTgtCurrents: can't query SQL command "<<qrystr<<std::endl;
    cerr << "... returning 0 "<<std::endl;
    return result;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getTgtCurrents: can't access to run info for run "<<runnb<<std::endl;
    cerr << "... returning 0 "<<std::endl;
    return result;
  }

  if (getCol(1)) {
    sprintf(qrystr, "select runnb,avg(tb_tgtcurrent.solencur), "
                    "  avg(tb_tgtcurrent.dipolcur), avg(tb_tgtcurrent.`SOLN-DCCT_ai`) "
                    "from tb_run left join tb_tgtcurrent "
                    "  on (tb_tgtcurrent.mestime between tb_run.starttime and tb_run.stoptime) "
                    "where runnb=%d group by runnb", runnb);
  } else {
    cerr<<"\nWarning in MySQLDBInterface::getTgtCurrents: No end time for run "<<runnb<<std::endl;
    cerr<<"     30 minutes length taken by default \n"<<std::endl;
    sprintf(qrystr, "select runnb,avg(tb_tgtcurrent.solencur), "
                    "  avg(tb_tgtcurrent.dipolcur), avg(tb_tgtcurrent.`SOLN-DCCT_ai`) "
                    "from tb_run left join tb_tgtcurrent "
                    "  on (tb_tgtcurrent.mestime between tb_run.starttime and DATE_ADD(tb_run.starttime,INTERVAL 30 MINUTE)) "
                    "where runnb=%d group by runnb", runnb);
  }
  endQuery();


  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getTgtCurrents: can't query SQL command "<<qrystr<<std::endl;
    cerr << "... returning 0 "<<std::endl;
    return result;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getTgtCurrents: can't access to offline polar info for run "<<runnb<<std::endl;
    cerr << "... returning 0 "<<std::endl;
    return result;
  }
  if (getCol(1) == 0) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getTgtCurrents: no target current measurement taken during the run "<<runnb<<std::endl;
    cerr << "... returning 0 "<<std::endl;
    return result;
  }

  result.first = atof(getCol(1));
  result.second = atof(getCol(2));
  // from 2006 the solencur column is signless, the good sign is to be taken from the SOLN-DCCT_ai column
  if (getCol(3) != 0) { result.first = fabs(result.first) * (atof(getCol(3))<0 ? -1 : 1); }
  endQuery();
  return result;
}



//----------------------------------------------------------------------------


std::pair<double,double> MySQLDBInterface::getTgtCurrentsTime(time_t mestime) {
  std::pair<double,double> result;

  result.first = 0;  // return 0 if no value found
  result.second = 0;

  if ( mestime < 1072911600 ) {  // before 1/1/2004
    cerr<<"Warning: MySQLDBInterface::getTgtCurrentsTime not suited for time before 2004"<<endl;
    cerr<<"     ... be careful, the result could be meanless"<<endl;
  }


  if ((mestime < tgtcurtime1) || (mestime > tgtcurtime2)) {

    bool connectedfg = isConnected();
    if (!connectedfg) {
      register bool tmpflag;
      register int retcode = 0;
      tmpflag = connect();
      if (!tmpflag) {
        cerr << "Error in MySQLDBInterface::getTgtCurrentsTime: can't connect to DB"<<endl;
        return result;
      }
      if ((retcode = selectDB())) {
        cerr << "Error in MySQLDBInterface::getTgtCurrentsTime: can't access DB, error code "<<retcode<<endl;
         return result;
     }
    }

    sprintf(qrystr, "select UNIX_TIMESTAMP(mestime), solencur, dipolcur, `SOLN-DCCT_ai` "
                    "from tb_tgtcurrent "
                    "where UNIX_TIMESTAMP(mestime)<=%ld and solencur is not null "
                    "order by mestime desc limit 1", mestime);

    if (!query(qrystr)) {
      cerr << "Error in MySQLDBInterface::getTgtCurrentsTime: can't query SQL command "<<qrystr<<std::endl;
      cerr << "... returning 0 "<<std::endl;
      return result;
    }
    if (!getNRows()) {
      endQuery();
      cerr << "Error in MySQLDBInterface::getTgtCurrentsTime: can't access to any target current measurement for time before "<<ctime(&mestime)<<std::endl;
      cerr << "... returning 0 "<<std::endl;
      return result;
    }
    tgtcurtime1 = atoi(getCol(0));
    tgtsolen1 = atof(getCol(1));
    // dipole info can be null, return 0 then
    if (getCol(2) != 0) { tgtdipol1 = atof(getCol(2)); } else { tgtdipol1 = 0; }
    // from 2006 the solencur column is signless, the good sign is to be taken from the SOLN-DCCT_ai column
    if (getCol(3) != 0) { tgtsolen1 = fabs(tgtsolen1) * (atof(getCol(3))<0 ? -1 : 1); }
    endQuery();

    sprintf(qrystr, "select UNIX_TIMESTAMP(mestime), solencur, dipolcur, `SOLN-DCCT_ai` "
                    "from tb_tgtcurrent "
                    "where UNIX_TIMESTAMP(mestime)>%ld and solencur is not null "
                    "order by mestime asc limit 1", mestime);

    if (!query(qrystr)) {
      cerr << "Error in MySQLDBInterface::getTgtCurrentsTime: can't query SQL command "<<qrystr<<std::endl;
      cerr << "... returning 0 "<<std::endl;
      return result;
    }
    if (!getNRows()) {
      endQuery();
      cerr << "Error in MySQLDBInterface::getTgtCurrentsTime: can't access to any target current measurement for time after "<<ctime(&mestime)<<std::endl;
      cerr << "... returning 0 "<<std::endl;
      return result;
    }
    tgtcurtime2 = atoi(getCol(0));
    tgtsolen2 = atof(getCol(1));
    if (getCol(2) != 0) { tgtdipol2 = atof(getCol(2)); } else { tgtdipol2 = 0; }
    // from 2006 the solencur column is signless, the good sign is to be taken from the SOLN-DCCT_ai column
    if (getCol(3) != 0) { tgtsolen2 = fabs(tgtsolen2) * (atof(getCol(3))<0 ? -1 : 1); }
    endQuery();

    if (!connectedfg) disconnect();
  }

  register double t2t1 = tgtcurtime2 - tgtcurtime1;
  if ( t2t1 == 0 ) {
    result.first = tgtsolen1;
    result.second = tgtdipol1;
    return result;
  }

  register double t2t = tgtcurtime2 - mestime;
  register double fa1 = t2t / t2t1;
  register double tt1 = mestime - tgtcurtime1;
  register double fa2 = tt1 / t2t1;

  if ( t2t > 600 || tt1 > 600 ) {
    cerr<<"Warning in MySQLDBInterface::getTgtCurrentsTime: measurements taken ";
    cerr<<"to give target magnets currents are at more than 10 minutes from "<<ctime(&mestime);
    cerr<<"  measurement time 1: "<<ctime(&tgtcurtime1);
    cerr<<"  measurement time 2: "<<ctime(&tgtcurtime2);
  }


  result.first = fa1 * tgtsolen1 + fa2 * tgtsolen2;
  result.second = fa1 * tgtdipol1 + fa2 * tgtdipol2;
  return result;
}



//----------------------------------------------------------------------------


double MySQLDBInterface::getSM2NMR(int runnb) {

  sprintf(qrystr, "select starttime,stoptime from tb_run where runnb=%d", runnb);
  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getSM2NMR: can't query SQL command "<<qrystr<<std::endl;
    return -1;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getSM2NMR: can't access to run info for run "<<runnb<<std::endl;
    return -1;
  }

  if (getCol(1)) {
    sprintf(qrystr, "select runnb,avg(tb_NMRSM2.nmr) "
                    "from tb_run left join tb_NMRSM2 "
                    "  on (tb_NMRSM2.date between tb_run.starttime and tb_run.stoptime) "
                    "where runnb=%d group by runnb", runnb);
  } else {
    cerr<<"\nWarning in MySQLDBInterface::getSM2NMR: No end time for run "<<runnb<<std::endl;
    cerr<<"     30 minutes length taken by default \n"<<std::endl;
    sprintf(qrystr, "select runnb,avg(tb_NMRSM2.nmr) "
                    "from tb_run left join tb_NMRSM2 "
                    "  on (tb_NMRSM2.date between tb_run.starttime and DATE_ADD(tb_run.starttime,INTERVAL 30 MINUTE)) "
                    "where runnb=%d group by runnb", runnb);
  }
  endQuery();

  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getSM2NMR: can't query SQL command "<<qrystr<<std::endl;
    return -1;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getSM2NMR: can't access to SM2 NMR info for run "<<runnb<<std::endl;
    return -1;
  }

  double result;
  if (getCol(1)) {
    result = strtod(getCol(1),0);
  } else {
    cerr<<"\nWarning in MySQLDBInterface::getSM2NMR: No SM2 NMR values found for run "<<runnb<<std::endl;
    cerr<<"     returning null value \n"<<std::endl;
    result = 0;
  }
  endQuery();
  return result;
}



//----------------------------------------------------------------------------


pair<double,double> MySQLDBInterface::getSMcurrents(int runnb) {
  std::pair<double,double> result;

  result.first = 1e9;  // return 1e9 if no value found
  result.second = 1e9;


  sprintf(qrystr, "select tb_run.runnb,tb_beamvalues.sm1,tb_beamvalues.sm2 "
                  "from tb_run left join tb_beamvalues "
                  "  using (bvalueid) "
                  "where runnb=%d", runnb);

  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::getSMcurrents: can't query SQL command "<<qrystr<<std::endl;
    cerr << "... returning 0 "<<std::endl;
    return result;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::getSMcurrents: can't access to SM1 and SM2 info for run "<<runnb<<std::endl;
    cerr << "... returning 0 "<<std::endl;
    return result;
  }

  result.first = atof(getCol(1));
  result.second = atof(getCol(2));
  endQuery();
  return result;
}



//----------------------------------------------------------------------------


char* MySQLDBInterface::get1stEvtPath(int runnb) {
  static char strtmp[300];

  if (fServer.find("lxfs1657") == std::string::npos) {
    cerr<<"Warning in MySQLDBInterface::get1stEvtPath: 1st event files are only available on lxfs1657 server\n";
    cerr<<"  will return /shift/lxfs1657/... path name"<<std::endl;
  }

  sprintf(qrystr, "select tb_1stevtfiles.runnb,tb_1stevtfiles.filename,tb_directories.Directory "
                  "from tb_1stevtfiles left join tb_directories "
                  "  on (tb_1stevtfiles.dirid = tb_directories.ID) "
                  "where tb_1stevtfiles.runnb=%d and tb_directories.workstation='SHIFT' limit 1", runnb);
  if (!query(qrystr)) {
    cerr << "Error in MySQLDBInterface::get1stEvtPath: can't query SQL command "<<qrystr<<std::endl;
    return 0;
  }
  if (!getNRows()) {
    endQuery();
    cerr << "Error in MySQLDBInterface::get1stEvtPath: can't access to run 1st event info for run "<<runnb<<std::endl;
    return 0;
  }

  if (getCol(1) && getCol(2)) {
    sprintf(strtmp, "%s/%s", getCol(2), getCol(1));
    cerr <<"MySQLDBInterface::get1stEvtPath: found 1st event file for run ";
    cerr <<getCol(0)<<": "<<strtmp<<std::endl;
  } else {
    cerr <<"Error in MySQLDBInterface::get1stEvtPath: no 1st event file found for run ";
    cerr <<runnb<<std::endl;
    endQuery();
    return 0;
  }

  endQuery();
  return strtmp;
}


//----------------------------------------------------------------------------




