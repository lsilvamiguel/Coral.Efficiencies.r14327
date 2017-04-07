#ifndef __Parameters__
#define __Parameters__

#include "mysql.h"
#include <string>

// #define HOST "pccoeb03"
#define HOST "pccodb00"
#define READER "anonymous"
#define WRITER "toto"
#define ADMINISTRATOR "toto"
#define PASSWORD "toto"
#define DBNAME "runlb"

enum  OpenModes {READ,WRITE,ADMIN};

extern MYSQL *connection;
extern MYSQL mysql;

extern bool DBOpen(MYSQL *conx, MYSQL* sql,int mode);
extern std::string MakeTime(const std::string & in);
extern std::string MakeSQLTime(const std::string & in);
extern bool DBEntry(const std::string& starttime, const std::string& stoptime, 
		    const std::string& dettype, const std::string& detname, 
		    const std::string& name, int authorid);
extern std::string GetLocation();
extern const char* MakeFileName(const std::string& detname, const std::string& starttime, 
				const std::string& stoptime);
extern int DBEntryAuthor(); 

#endif






