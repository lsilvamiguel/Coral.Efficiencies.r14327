
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "configMySQLDB.h"

MySQLDBInterface* db;

void usage(const char* pgmname, const char* errmsg);
bool get1stEvtPath(int runnb);




int main(int argc, char** argv) {
  int runnb = -1;

  if (argc < 2) { usage(argv[0], "not enough arguments"); }
  runnb = atoi(argv[1]);


//   char* strpasswd = getpass("Enter DB write password: ");

  db = new MySQLDBInterface(DBSERVER, DBREADER, "", DBNAME);
  if (!db->connect()) {
    std::cerr << "can't connect to database, exiting...\n";
    exit(1);
  }
  db->selectDB();

  if (get1stEvtPath(runnb)) {
    db->disconnect();
    return 0;
  } else {
    db->disconnect();
    exit(1);
  }
}



bool get1stEvtPath(int runnb) {

  char* filepath = db->get1stEvtPath(runnb);
  if (filepath) { std::cout << filepath<<std::endl; return true; }
  else { std::cerr<<"No first event file found for run "<<runnb<<std::endl; return false; }
}





void usage(const char* pgmname, const char* errmsg) {
  if (errmsg) fprintf(stderr, "Error: %s\n\n", errmsg);
  fprintf(stderr, "usage: %s run_number\n", pgmname);
  exit(1);
}

