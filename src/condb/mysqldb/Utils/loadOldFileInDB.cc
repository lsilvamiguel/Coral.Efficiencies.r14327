
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <regex.h>
#include <pwd.h>

#include "configMySQLDB.h"

MySQLDBInterface* db;

void usage(char* pgmname);
int filetype(const char* filename);
void treat_dir_or_file(const char* filename);
int treat_file(const char* filename);

bool actualtimefg = false;


int main(int argc, char** argv) {

  if (argc < 2) { usage(argv[0]); }

  register int deb_arg = 1;
  while ( (deb_arg < argc) && (argv[deb_arg][0] == '-') ) {
    std::string sargv(argv[deb_arg]);
    if ( sargv == "-actualtime" ) {
      actualtimefg = true;
      std::cout << "Using actual time instead of file last modification time as entry time\n";
    }
    deb_arg++;
  }

  if (argc <= deb_arg) { usage(argv[0]); }

  char* strpasswd = getpass("Enter DB write password: ");

  db = new MySQLDBInterface(DBSERVERWRITE, DBWRITER, strpasswd, DBNAME);
  if (!db->connect()) {
    std::cerr << "can't connect to database, exiting...\n";
    exit(1);
  }
  db->selectDB();

  for (register int ii = deb_arg; ii < argc; ii++) {
    treat_dir_or_file(argv[ii]);
  }

  db->disconnect();
}



int treat_file(const char* filename) {
  register const char* fnameidx;
  char *starttime, *endtime, *detname, *typecalib;
  char newfilename[1024], command[2048], strtmp[1024];
  struct stat info;

  if (filename == 0) return 1;    // no file name
  fnameidx = strrchr(filename, '/');
  if (fnameidx == 0) { fnameidx = filename; }   // no '/' found
    else { fnameidx++; }        // placing just after the '/'
  if (fnameidx[0] == 0) return 2; // if '/' is the last char

  regex_t rx_calibname;
  // building regexp command
  register int rx_comperror = regcomp(&rx_calibname,
                            "([^~]+)(~~start-)([^~]+)(~~finish-)([^~]+)", REG_EXTENDED);
  if (rx_comperror) {
    std::cerr <<"treat_file: Error in regexp expression\n";
    regerror(rx_comperror, &rx_calibname, strtmp, 1024);
    std::cerr <<"  regexp error string: "<<strtmp<<std::endl;
    exit(1);
  }

  regmatch_t pmatch[8];
  // apply on filename
  register int rx_matcherror = regexec(&rx_calibname, fnameidx, 8, pmatch, 0);
  if (rx_matcherror) {
    std::cerr <<"file "<<fnameidx<<" has not a good name format\n";
    regerror(rx_matcherror, &rx_calibname, strtmp, 1024);
    std::cerr <<"  regexp error string: "<<strtmp<<std::endl;
    return 4;
  }

  // first part (tbname) should have at least 8 chars
  if ((pmatch[1].rm_eo - pmatch[1].rm_so) < 8) {
    std::cerr <<"file "<<fnameidx<<" has not a good TB name\n";
    return 3;
  }

  detname = strdup(fnameidx + pmatch[1].rm_so);
  detname[8] = 0;        // to get only the TBname

  if ((pmatch[1].rm_eo - pmatch[1].rm_so) > 8) {  // there is a type of calib
    typecalib = strdup(fnameidx + pmatch[1].rm_so + 8);
    typecalib[pmatch[1].rm_eo - pmatch[1].rm_so - 8] = 0;
  } else {
    typecalib = 0;
  }

  starttime = strdup(fnameidx + pmatch[3].rm_so);
  starttime[pmatch[3].rm_eo - pmatch[3].rm_so] = 0;
  endtime = strdup(fnameidx + pmatch[5].rm_so);
  endtime[pmatch[5].rm_eo - pmatch[5].rm_so] = 0;

  if (stat(filename, &info)) {
    std::cerr <<"Can't stat file "<<filename<<"\n  it probably does not exist\n";
    return -1;
  }

  // looks for user name in the DB
  struct passwd *pwstruct = 0;
  int authorID = 0;
  pwstruct = getpwuid(info.st_uid);
  if (pwstruct) authorID = db->getUserID(pwstruct->pw_name);
//   authorID = db->getUserID(getenv("LOGNAME"));
  if (authorID == 0) { authorID = DBNOAUTHORID; }

  // take entry time from the last modification date of the file
  char* creattimefile = "";
  if (strftime(strtmp, 1024, "%Y-%m-%d-%H:%M:%S", localtime(&info.st_mtime))) {
    creattimefile = strtmp;
  }
  if (actualtimefg) {  // using real time instead of file time if flag true
    const time_t actualtime = time(0);
    if (strftime(strtmp, 1024, "%Y-%m-%d-%H:%M:%S", localtime(&actualtime))) {
      creattimefile = strtmp;
    }
  }

  // check if the entry is not yet in the DB
  int entryID = db->checkEntry(detname, starttime, endtime, typecalib,
                               creattimefile);

  if (entryID) {        // entry already there
    std::cout <<"Entry for "<<detname;
    if (typecalib) { std::cout <<" typecalib "<<typecalib; }
    std::cout << " starting "<<starttime;
    std::cout << " ending "<<endtime<<" already registered\n";
    std::cout <<"  entry ID "<<entryID<<" creation time "<<db->getCreatTime(entryID);
    std::cout <<std::endl;
    std::cout <<"  inserted in MySQLDB the "<<db->getEntryTime(entryID);
    std::cout <<std::endl;

    free(starttime);
    free(endtime);
    free(detname);
    if (typecalib) free(typecalib);
    return 0;
  }

  // add the entry to the DB
  entryID = db->addEntry(detname, authorID, 0, "", starttime, endtime,
                             typecalib, creattimefile);

  if (entryID < 0) {  // entry failed
    std::cerr <<"can't add database entry for "<<fnameidx<<std::endl;
    free(starttime);
    free(endtime);
    free(detname);
    if (typecalib) free(typecalib);
    return 5;
  }

  char* typecalib_2 = "";
  if (typecalib) typecalib_2 = typecalib;
  std::string creattime = db->getCreatTime(entryID);


  // determine directory where to put the calib file
  std::string dbdirectory_str;
  int dbdir_id;
  std::pair<int, std::string> dbdir_entry = db->getCfgDirectory("calibdbfiles");
  dbdir_id = dbdir_entry.first;
  dbdirectory_str = dbdir_entry.second;

  // copy the file to the MySQLDB directory with a new file name
  sprintf(newfilename,"entry_%d~%s~%s~%s~%s~%s",
          entryID, creattime.c_str(), detname, typecalib_2, starttime, endtime);
  sprintf(command,"cp %s %s/%s", filename, dbdirectory_str.c_str(), newfilename);
  register int retcp = system(command);
  if (retcp) {
    std::cerr <<"Error: can't copy "<<filename<<" to "<<dbdirectory_str<<"/"<<newfilename<<std::endl;
  } else {
    db->updateEntry(entryID, dbdir_id, newfilename);
  }

  std::cout <<"inserted entry "<<entryID<<" for "<<fnameidx<<" in database\n";
  std::cout <<"  with detname: "<<detname<<", authorID: "<<authorID;
  std::cout <<", starttime: "<<starttime<<", endtime:"<<endtime;
  std::cout <<", creattime: "<<creattime;
  if (typecalib) { std::cout <<", typecalib: "<<typecalib; }
  std::cout <<std::endl;
  std::cout <<"Insertion in MySQLDB time: "<<db->getEntryTime(entryID)<<std::endl;;

  free(starttime);
  free(endtime);
  free(detname);
  if (typecalib) free(typecalib);
  return 0;
}


void treat_dir_or_file(const char* filename) {
  register int ftype;

  ftype = filetype(filename);
  if (ftype < 1) return;

  if (ftype == 2) {
    struct dirent *filedir;
    DIR* dirfd;
    std::string filedirname;

    dirfd = opendir(filename);
    while ((filedir = readdir(dirfd))) {
      std::string subfilename(filedir->d_name);
      if ((subfilename == ".") || (subfilename == "..")) continue;
      filedirname = std::string(filename) + "/" + subfilename;
      treat_dir_or_file(filedirname.c_str());
    }
    closedir(dirfd);
    return;
  }

  if ((ftype = treat_file(filename)) != 0) {
    std::cerr <<"failed to add "<<filename<<" to MySQLDB calib database (code "<<ftype<<")\n";
  }
}


int filetype(const char* filename) {
  struct stat info;

  if (stat(filename, &info)) {
    std::cerr <<"Can't stat file "<<filename<<"\n  it probably does not exist\n";
    return -1;
  }

  if (S_ISREG(info.st_mode)) return 1;
  if (S_ISDIR(info.st_mode)) return 2;
  std::cerr <<filename<<" is not a regular file or a directory\n";
  return 0;
}


void usage(char* pgmname) {
  fprintf(stderr, "usage: %s <file or directory> [<file or dir> ...]\n", pgmname);
  exit(1);
}


