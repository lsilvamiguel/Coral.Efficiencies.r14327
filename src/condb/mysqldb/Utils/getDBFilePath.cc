// $Id: getDBFilePath.cc 14069 2015-09-17 20:44:46Z lsilva $

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "configMySQLDB.h"

MySQLDBInterface* db;

void usage(const char* pgmname, const char* errmsg);
void getDBFilePath(const char* detname, int runnb, const char* datestr,
                   const char* typecalib, const char* entrydate);

char* servername = DBSERVER;


int main(int argc, char** argv) {
  int runnb = -1;
  char *datestr = 0, *typecalib = 0, *detname = 0, *entrydate = 0, *workenv = 0;
  int port = 0;

  if (argc < 2) { usage(argv[0], "not enough arguments"); }

  register int ii = 1;
  while (ii < argc) {
    register char *argvii, *argvii1;
    argvii = argv[ii];
    if ((ii+1) < argc) argvii1 = argv[ii+1];
      else argvii1 = 0;

    if (argvii[0] == '-') {     // it is an option
      if (strcmp(argvii, "-r") == 0) {
        if (argvii1 && (argvii1[0] != '-')) runnb = atoi(argvii1);
          else usage(argv[0], "badly formed -r option");
        ii += 2;
        continue;
      }
      if (strcmp(argvii, "-d") == 0) {
        if (argvii1 && (argvii1[0] != '-')) datestr = argvii1;
          else usage(argv[0], "badly formed -d option");
        ii += 2;
        continue;
      }
      if (strcmp(argvii, "-t") == 0) {
        if (argvii1 && (argvii1[0] != '-')) typecalib = argvii1;
          else usage(argv[0], "badly formed -t option");
        ii += 2;
        continue;
      }
      if (strcmp(argvii, "-e") == 0) {
        if (argvii1 && (argvii1[0] != '-')) entrydate = argvii1;
          else usage(argv[0], "badly formed -e option");
        ii += 2;
        continue;
      }
      if (strcmp(argvii, "-s") == 0) {
        if (argvii1 && (argvii1[0] != '-')) servername = argvii1;
          else usage(argv[0], "badly formed -s option");
        ii += 2;
        continue;
      }
      if (strcmp(argvii, "-w") == 0) {
        if (argvii1 && (argvii1[0] != '-')) workenv = argvii1;
          else usage(argv[0], "badly formed -w option");
        ii += 2;
        continue;
      }
      if (strcmp(argvii, "-P") == 0) {
        if (argvii1 && (argvii1[0] != '-')) {
	  char *end, **endptr = &end; port = strtol(argvii1,endptr,10);
	  if (*endptr!=argvii1+strlen(argvii1))
	    usage(argv[0],"bad arg. to -P option: not an integer");
	}
	else usage(argv[0], "badly formed -P option");
        ii += 2;
        continue;
      }
      if (strcmp(argvii, "--help") == 0) {
        usage(argv[0], 0);
        ii++;
        continue;
      }
      // unknown option found
      char stmp[200];
      sprintf(stmp, "unknown option %s", argvii);
      usage(argv[0], stmp);
      ii++;
      continue;
    }
    else {
      if ((ii+1) == argc) detname = argvii;
       else usage(argv[0], "tbname is not the last argument or more than one tbname");
      ii++;
      continue;
    }
  }

  if (!detname) usage(argv[0], "no tbname given");
  if (datestr && (runnb > -1))
      usage(argv[0], "-d and -r option can't be given at the same time");
  if (!datestr && (runnb == -1)) {
    std::cerr <<"Warning: no -d or -r option given !\n";
  }

//   char* strpasswd = getpass("Enter DB write password: ");

  db = new MySQLDBInterface(servername, DBREADER, "", DBNAME);
  if (port) db->setNumPort(port);
  if (!db->connect()) {
    std::cerr << "can't connect to database, exiting...\n";
    exit(1);
  }
  db->selectDB();
  if (workenv) {
    std::string workstr(workenv);
    db->setSpecialPlace(workstr);
  }

  getDBFilePath(detname, runnb, datestr, typecalib, entrydate);

  db->disconnect();
}



void getDBFilePath(const char* detname, int runnb, const char* datestr,
                   const char* typecalib, const char* entrydate) {

  std::string datestring("");
  if (runnb != -1) {
    datestring = db->getDateFromRunNb(runnb);
  }
  if (!datestr && (datestring != "")) datestr = datestring.c_str();
  if (!datestr) {
    std::cerr << "getDBFilePath: run "<<runnb<<" not found in database, exiting..."<<std::endl;
    exit(2);
  }

fprintf(stderr, "detname: %s, datestr: %s\n", detname, datestr);


  std::string filepath = db->giveFilepath(detname, datestr, typecalib, entrydate);
  std::cout << filepath<<std::endl;

}





void usage(const char* pgmname, const char* errmsg) {
  if (errmsg) fprintf(stderr, "Error: %s\n\n", errmsg);
  fprintf(stderr, "usage: %s {-r <run number> | -d <date>} [-t <calibration type>] [-e <max entry date>] [-s <server_name>] [-w <work env>] [-P <port#>] <tbname>\n", pgmname);
  exit(1);
}

