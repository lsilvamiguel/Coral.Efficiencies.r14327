#include "MyInterface.h"

#define DBSERVER "pccoeb03"
#define DBUSER "toto"
#define DBNAME "runlb"

int main() {

  MyInterface db(DBSERVER,DBUSER,"toto",DBNAME);
  db.connect();
 
  string errmsg = db.cdbUpload("/2002/PS01Y1__~~start-2002-06-04-00:00:00~~finish-2002-06-05-00:00:00","5");
  if(!errmsg.empty()) cout<<errmsg<<endl;
  
  return 0;
}









