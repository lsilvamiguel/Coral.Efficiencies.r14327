#include "MyInterface.h"

#define DBSERVER "pccoeb03"
#define DBUSER "toto"
#define DBNAME "runlb"

int main() {

  MyInterface db(DBSERVER,DBUSER,"toto",DBNAME);
  db.connect();
 

  return 0;
}









