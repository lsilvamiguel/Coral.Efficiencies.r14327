#ifndef __Parameters__
#define __Parameters__

#include "mysql.h"
#include <string>

#define HOST "pccoeb03"
#define READER "anonymous"
#define WRITER "onl"
#define ADMINISTRATOR "onl"
#define DBNAME "runlb"

enum  OpenModes {READ,WRITE,ADMIN};

extern MYSQL *connection;
extern MYSQL mysql;

extern bool DBOpen(MYSQL *conx, MYSQL* sql,int mode);
extern string MakeTime(const string & in);
extern string MakeSQLTime(const string & in);
extern bool DBEntry(const string& starttime, const string& stoptime, 
		    const string& dettype, const string& detname, 
		    const string& name, int authorid);
extern string GetLocation();
extern const char* MakeFileName(const string& detname, const string& starttime, 
				const string& stoptime);
extern int DBEntryAuthor(); 

#endif






