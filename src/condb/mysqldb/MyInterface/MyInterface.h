#ifndef __MyInterface__
#define __MyInterface__

#include "MySQLInterface.h"


class MyInterface: public MySQLInterface {

 public:
  
  /// Destructor
  virtual  ~MyInterface (void);

  MyInterface (const char* server, const char* username,
	       const char* passwd, const char* dbname);
 

  // some tools ---------------------------------------------------
  
  // is dirname writeable
  bool   dirWritable(const std::string& dirname) const; 
  
  // forms a cdb filename
  std::string fileName(const std::string& tbname,const std::string& start, 
		       const std::string& stop);
  
  // convenient functions for cdb ---------------------------------
  
  //puts one file to DB
  std::string   cdbUpload(const std::string& start, const std::string& stop, 
			  const std::string& tbname,
			  const std::string& file, const std::string& dir, 
			  const std::string& authorid);

  // puts one file to DB
  std::string   cdbUpload(const std::string& file,const std::string& dir,
			  const std::string& authorid);

  // returns  directory id 
  std::string cdbDirID(const std::string& dir);

  // returns author id
  std::string authorID(const std::string& login,
		       const std::string& firstname,
		       const std::string& lastname,
		       const std::string& email);
  
  // checks that time has the correct format
  bool checkFormat(const std::string& time);

  // checks that tbname is ok (8 letters)
  bool checkTBName(const std::string& tbname);

  // checks local file
  bool checkLocalFile(const std::string& filename);

  // checks that id is in the DB
  bool cdbCheckAuthorID(const std::string& id);

  // returns "" under failure, "myself" if authorid is the owner,
  // and otherwise the name of the owner 
  std::string fileOwner(const std::string& start, const std::string& stop, 
			const std::string& file, const std::string& authorid);
};

#endif 


