
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <pwd.h>

#include "configMySQLDB.h"

MySQLDBInterface* db;


char *startrun=0, *endrun=0;
const char *starttime=0, *endtime=0;
char *creattime=0, *detname=0;
char *dettype=0, *passwd=0, *typecalib=0, *filename=0;


void usage(char* pgmname) {
  fprintf(stderr, "usage: %s [OPTIONS] <calibration file>\n\n", pgmname);
  fprintf(stderr, "  OPTIONS are:\n");
  fprintf(stderr, "    -begin <date>      : begin of validity period by date (one begin or beginrun mandatory)\n");
  fprintf(stderr, "    -beginrun <runnb>  : begin of validity period by run nb (see above)\n");
  fprintf(stderr, "    -end <date>        : end of validity period by date (one end or endrun mandatory)\n");
  fprintf(stderr, "    -endrun <runnb>    : end of validity period by run nb (see above)\n");
  fprintf(stderr, "    -detname <string>  : detector name (mandatory)\n");
  fprintf(stderr, "    -dettype <string>  : detector type (2 chars, deduced from detname if not present)\n");
  fprintf(stderr, "    -typecalib <string>: type of calib file (like \"timing\", default if not present)\n");
  fprintf(stderr, "    -creattime <date>  : gives a specific creation time\n");
  fprintf(stderr, "    -passwd <string>   : MySQL write account password\n");
  fprintf(stderr, "  <date> format is: YYYY-MM-DD-hh:mm:ss\n");
  exit(1);
}


void decode_opt(int argc, char** argv) {

  if (argc < 2) { usage(argv[0]); }
  int parg = 1;
  while (parg < argc) {
    std::string sargv(argv[parg]);

    if ((sargv == "-begin") && (parg + 1 < argc)) {
      starttime = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-beginrun") && (parg + 1 < argc)) {
      startrun = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-end") && (parg + 1 < argc)) {
      endtime = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-endrun") && (parg + 1 < argc)) {
      endrun = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-creattime") && (parg + 1 < argc)) {
      creattime = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-passwd") && (parg + 1 < argc)) {
      passwd = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-detname") && (parg + 1 < argc)) {
      detname = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-dettype") && (parg + 1 < argc)) {
      dettype = argv[parg+1];
      parg += 2;
      continue;
    }

    if ((sargv == "-typecalib") && (parg + 1 < argc)) {
      typecalib = argv[parg+1];
      parg += 2;
      continue;
    }

    if (parg+1 == argc) {
      filename = argv[parg];
      parg++;
      continue;
    } else {
      fprintf(stderr, "Error: unknow option %s or more than one calib file given\n\n", argv[parg]);
      usage(argv[0]);
    }
  }
}


int main(int argc, char** argv) {

  decode_opt(argc, argv);

  if (!startrun && !starttime) {
    std::cerr <<"Error: no -begin or -beginrun present\n\n"; usage(argv[0]);
  }
  if (startrun && starttime) {
    std::cerr <<"Error: -begin and -beginrun present together, use only one\n\n"; usage(argv[0]);
  }
  if (!endrun && !endtime) {
    std::cerr <<"Error: no -end or -endrun present\n\n"; usage(argv[0]);
  }
  if (endrun && endtime) {
    std::cerr <<"Error: -end and -endrun present together, use only one\n\n"; usage(argv[0]);
  }

  if (!filename) { std::cerr <<"Error: no calib file name given\n\n"; usage(argv[0]);}
  if (!detname) { std::cerr <<"Error: no detector name given\n\n"; usage(argv[0]);}

  if (strlen(detname) != 8) {
    fprintf(stderr, "Error: the detector name (%s) must have 8 chars\n\n", detname);
    usage(argv[0]);
  }

  if (!passwd) {
    passwd = getpass("Enter DB write password: ");
  }

  if (!dettype) {
    dettype = (char*) malloc(3);
    dettype[0] = detname[0];
    dettype[1] = detname[1];
    dettype[2] = 0;
  } else {
    if (strlen(dettype) > 2) { dettype[2] = 0; }
  }

  if (!creattime) {
    creattime = (char*) malloc(1024);
    time_t tnow = time(0);
    strftime(creattime, 1024, "%Y-%m-%d-%H:%M:%S", localtime(&tnow));
  }


  // connection to DB
  db = new MySQLDBInterface(DBSERVERWRITE, DBWRITER, passwd, DBNAME);
  if (!db->connect()) {
    std::cerr << "can't connect to database, exiting...\n";
    exit(1);
  }
  db->selectDB();

  // if valid period given with runnb, transform it in date
  std::string starttimestr("");
  std::string endtimestr("");
  if (startrun) {
    register int runnb = atoi(startrun);
    starttimestr = db->getDateFromRunNb(runnb);
    if (starttimestr == "") { std::cerr<<"Error: run nb "<<runnb<<" not found !"<<std::endl; exit(1); }
    std::cout << " start of valid period: run nb "<<runnb<<" date "<<starttimestr<<std::endl;
    starttime = starttimestr.c_str();
  }
  if (endrun) {
    register int runnb = atoi(endrun);
    endtimestr = db->getDateFromRunNb(runnb);
    if (endtimestr == "") { std::cerr<<"Error: run nb "<<runnb<<" not found !"<<std::endl; exit(1); }
    std::cout << " end of valid period: run nb "<<runnb<<" date "<<endtimestr<<std::endl;
    endtime = endtimestr.c_str();
  }


  int authorID = db->getUserID(getenv("LOGNAME"));
  if (authorID == 0) { authorID = DBNOAUTHORID; }

  // check if the entry is not yet in the DB
  int entryID = db->checkEntry(detname, starttime, endtime, typecalib,
                               creattime);

  if (entryID) {        // entry already there
    std::cout <<"ERROR: Entry for "<<detname;
    if (typecalib) { std::cout <<" typecalib "<<typecalib; }
    std::cout << " starting "<<starttime;
    std::cout << " ending "<<endtime<<" already registered !!\n";
    std::cout <<"  entry ID "<<entryID<<" creation time "<<db->getCreatTime(entryID);
    std::cout <<"  entered in MySQLDB "<<db->getEntryTime(entryID);
    std::cout <<std::endl;
    db->disconnect();
    exit(1);
  }

// cerr<<"arret pour test"<<std::endl;
// exit(10);

  // add the entry to the DB
  entryID = db->addEntry(detname, authorID, 0, "", starttime, endtime,
                             typecalib, creattime, dettype);

  if (entryID < 0) {  // entry failed
    std::cerr <<"can't add database entry for "<<filename<<" detector "<<detname<<std::endl;
    db->disconnect();
    exit(1);
  }

  char* typecalib_2 = "";
  if (typecalib) typecalib_2 = typecalib;
  std::string creattime_2 = db->getCreatTime(entryID);


  // determine directory where to put the calib file
  std::string dbdirectory_str;
  int dbdir_id;
  std::pair<int, std::string> dbdir_entry = db->getCfgDirectory("calibdbfiles");
  dbdir_id = dbdir_entry.first;
  dbdirectory_str = dbdir_entry.second;

  // copy the file to the MySQLDB directory with a new file name
  char newfilename[1024], command[2048];
  sprintf(newfilename,"entry_%d~%s~%s~%s~%s~%s",
          entryID, creattime_2.c_str(), detname, typecalib_2, starttime, endtime);
  sprintf(command,"cp %s %s/%s", filename, dbdirectory_str.c_str(), newfilename);
  register int retcp = system(command);
  if (retcp) {
    std::cerr <<"Error: can't copy "<<filename<<" to "<<dbdirectory_str<<"/"<<newfilename<<std::endl;
  } else {
    db->updateEntry(entryID, dbdir_id, newfilename);
  }

  std::cout <<"Inserted entry "<<entryID<<" for "<<filename<<" in database\n";
  std::cout <<"  with detname: "<<detname<<", authorID: "<<authorID;
  std::cout <<", starttime: "<<starttime<<", endtime: "<<endtime;
  std::cout <<", creattime: "<<creattime_2;
  if (typecalib) { std::cout <<", typecalib: "<<typecalib; }
  std::cout <<std::endl;
  std::cout <<"Insertion in MySQLDB time: "<<db->getEntryTime(entryID)<<std::endl;;

  db->disconnect();
}
















