#ifndef DTPIMP_H

#include "Dtp.h"
#include "../MyInterface/MyInterface.h"
#include <iostream>
#include <string>

class DtpImp : public Dtp {
  Q_OBJECT

 private:
  // mysql instance
  MyInterface *db_;

  // time format for mysql
  static const char* timeFormat_;

  // time format used in filenames
  static const char* timeFormatSasha_;

  // file for user personal information
  static const char* persoParams_;

  // pixmaps for connect/disconnect button
  QPixmap* window_destroyXpm_;
  QPixmap* emptyXpm_;

  // user information
  string login_;
  string firstName_;
  string lastName_;
  string email_;

  // user id (from DB)
  string authorID_;

 public:
  DtpImp();
  ~DtpImp();
  
  // get user personal parameters. 
  bool registerUser();
  
  // connect/disconnect to/from DB
  bool connectDB();

  // select files from DB and fills table
  int  selectDB();

  // download file(s) from DB
  void downloadDB();

  // finds all DB folders, and updates locationCombo_ 
  void locationDB() const;

  // warns user that there is no connection
  void NotConnected() const;

  // puts file(s) in the DB
  void uploadDB();

  // removes files from DB (must belong to user)
  void removeDB();

  // updates file information (not implemented yet)
  void updateDB();
  
  // directory writeable ? (replace by the similar function in MyInterface)
  bool dbDirWritable(const char*) const;

  // checks that tbname is in the DB 
  bool tbnameInDB(const char* tbname) const;

  // not used 
  bool string_match(const string &str, const string &pattern);    
};

#endif


