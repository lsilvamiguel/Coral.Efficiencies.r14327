#include "monDB.h"
#include <stdio.h>
#include <unistd.h>

#define DBSERVER "pccoeb03"
#define DBUSER "anonymous"
#define DBNAME "runlb"

#define RUNTOREAD 10000

int main(int argc, char** argv) {
  long sumrunnb = 0;
  int retrunnb, nbtrouv = 0, runnb, runtoread = RUNTOREAD;
  char sqlreq[200];
  char* filename = 0;
  time_t debut, fin;
  FILE* fd;

  srand(getpid());

  if (argc >= 2) { runtoread = atoi(argv[1]); }
  if (argc >= 3) { filename = argv[2]; }

  CSMon::monDB mondb(DBSERVER,DBUSER,DBNAME);
  mondb.connectDB();

  debut = time(0);

  for (register int ii=0; ii < runtoread; ii++) {
    runnb = (int) (rand()/(1.*RAND_MAX) * 8000 + 16000);
    sprintf(sqlreq, "select runnb,starttime from tb_run where runnb=%d", runnb);
    if(mondb.queryDB(sqlreq) && mondb.getColDB(0)) {
      retrunnb = atoi(mondb.getColDB(0));
      sumrunnb += retrunnb;
      nbtrouv++;
    }
    mondb.endQueryDB();
  }
  fin = time(0);

  if (filename) {
    fd = fopen(filename, "a");
  } else {
    fd = stdout;
  }

  fprintf(fd, "%d runs found on %d in %d s, sum of the run numbers: %ld\n",
          nbtrouv, runtoread, (int)(fin - debut), sumrunnb);
  fclose(fd);

  return 0;
}









