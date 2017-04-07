#include "DataBase.h"

using namespace Reco;
using std::cout;
using std::cerr;
using std::endl;

//         int     tm_sec;         /* seconds */
//         int     tm_min;         /* minutes */
//         int     tm_hour;        /* hours */
//         int     tm_mday;        /* day of the month */
//         int     tm_mon;         /* month */
//         int     tm_year;        /* year */
//         int     tm_wday;        /* day of the week */
//         int     tm_yday;        /* day in the year */
//         int     tm_isdst;       /* daylight saving time */

int main(void)
{
  try
  {
    system("rm -rf reco.db"); // remove previous data base

    DataBase db("reco.db",DataBase::WRITE|DataBase::READ|DataBase::CREATE);

    // ---------------------------------------------
    // Write data to data base.

    tm t[2];  // Start,Finish

    // Start time
    t[0].tm_year = 2000 - 1900;
    t[0].tm_mon  = 4    - 1;
    t[0].tm_mday = 25;
    t[0].tm_hour = 0;
    t[0].tm_min  = 0;
    t[0].tm_sec  = 0;
    
    // Finish time
    t[1].tm_year = 2000 - 1900;
    t[1].tm_mon  = 4    - 1;
    t[1].tm_mday = 25;
    t[1].tm_hour = 23;
    t[1].tm_min  = 59;
    t[1].tm_sec  = 59;

    DataBase::Element<int> a("int",1,&t[0],&t[1]);
    db.Write(a);

    DataBase::Element<int> b("int",2);
    db.Write(b);

    t[0].tm_hour = 10;
    b.SetTimeStart(t[0]);
    b.SetData(3);
    db.Write(b);

    t[0].tm_hour = 20;
    b.SetTimeFinish(t[0]);
    b.SetData(4);
    db.Write(b);

    // ---------------------------------------------

    {
      db.Write("container",1.1,t[0],t[1]);
      double v=0;
      t[0].tm_mday = 24;
      db.Read ("container",v,t[0]);
      cout << v << endl;
    }

    // ---------------------------------------------
    
    // Read elements of data base
    
    DataBase::Element<int> c("int");
    
    db.Read(c);
    
    std::vector<DataBase::ElementBase> elements;
    db.FindAllVersions("int",elements);
    for( size_t i=0; i<elements.size(); i++ )
    {
      elements[i].Print();
      cout << endl;
    }
    
  }
  catch( const std::exception &e )
  {
    cerr << e.what() << endl;
  }
  catch( const char *s )
  {
    cerr << s << endl;
  }
  catch( ... )
  {
    cerr << "Unknown exception\n";
  }
}

////////////////////////////////////////////////////////////////////////////////
