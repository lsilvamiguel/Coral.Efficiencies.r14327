#include "monDB.h"

#define DBSERVER "pccoeb03"
#define DBUSER "anonymous"
#define DBNAME "runlb"

int main() {

  CSMon::monDB mondb(DBSERVER,DBUSER,DBNAME);
  mondb.connectDB();

  cout<<mondb.RootRefFile("MM")<<endl;

  return 0;
}









