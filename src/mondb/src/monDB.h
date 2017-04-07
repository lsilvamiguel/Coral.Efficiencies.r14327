/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/mondb/src/monDB.h,v $
   $Date: 2002/10/24 15:32:20 $
   $Revision: 1.6 $
   -------------------------------------------------------------------------

   This file is part of Compass monitoring program.

   Authors:
     based on original program of Alexander Zvyagin ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     adapted by Colin Bernet & Damien Neyret

   Copyright(C): 1999-2000  A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef monDB___include
#define monDB___include

#define  NO_CLIENT_LONG_LONG

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <typeinfo>
#include <mysql.h>

namespace CSMon {

////////////////////////////////////////////////////////////////////////////////


/*! \brief Data base.

*/
class monDB
{
  //============================================================================
  // Types
  //============================================================================



  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
  virtual  ~monDB (void);


    /*! \brief Create data base
     */
              monDB (const char* server, const char* username, const char* dbname);


  //============================================================================
  // Operators
  //============================================================================


  //============================================================================
  // Methods
  //============================================================================

  public:

    template <class T>
    void                Write           (const string &container,const T& v, const tm &start, const tm &end);

    template <class T>
    void                Read            (const string &container,T& v, const tm &point);

    static void         PrintTime       (const tm &t,ostream &o=cout);

    const char*         RootRefFile(const char* dettype);

    bool                connectDB(void);
    void                disconnectDB(void);
    bool                isConnected(void)   { return fConnected; }
    bool                queryDB(string& request) { return queryDB(request.c_str()); }
    bool                queryDB(const char* request);
    char*               getColDB(int i)  { if (fRow) return fRow[i]; else return 0; }
    void                endQueryDB(void);
    bool                getNextRowDB(void);


  //============================================================================
  // Attributes, data
  //============================================================================

  private:

    const string    fServer;
    const string    fUserName;
    const string    fDBName;
    MYSQL           fMysql;
    MYSQL          *fConnection;
    bool            fConnected;
    MYSQL_RES      *fResult;
    bool            fResInMem;
    MYSQL_ROW       fRow;
    int             fNRows;
};



////////////////////////////////////////////////////////////////////////////////

template <class T>
void monDB::Write(const string &container,const T& v, const tm &start, const tm &end)
{
  //    throw CS::Exception("monDB::Write():  no Write implementation actually... sorry !");


//   if( !(open_mode&WRITE) )
//     throw CS::Exception("monDB::Write():  monDB \"%s\" was opened in read-only mode = %d",
//                      path.c_str(),open_mode);
// 
//   const string file_name = path+'/'+Element<T>(container,&start,&end).CreateFileName();
//   ofstream f(file_name.c_str());
// 
//   if( !f.is_open() )
//     throw CS::Exception("monDB::Write():  container \"%s\"  can not open file \"%s\"",
//                      container.c_str(), file_name.c_str());
// 
//   if( !(f<<v) )
//     throw CS::Exception("monDB::Write():  container \"%s\"  error in writing to file \"%s\"",
//                      container.c_str(), file_name.c_str());
}


////////////////////////////////////////////////////////////////////////////////

template <class T>
void monDB::Read(const string &tbname, T& v, const tm &point)
{
  if(!fConnected)
    if(!connectDB()) {
      cerr<<"WARNING monDB::Read. Not connected to database."<<endl;
      return;
  }

  string sqlreq;
  sqlreq = "select Directory from tb_directories where tb_directories.key like 'calibration'";
  if(!queryDB(sqlreq)) 
    return;
  string directory(getColDB(0));
  endQueryDB();

  sqlreq = "";
  sqlreq += "select FileName,StartTime,EndTime from tb_calibration where StartTime <= ";
  sqlreq += TimeStrSQL(point);
  sqlreq = sqlreq + " and DetName like \"" + tbname + "\" order by StartTime desc";
  queryDB(sqlreq);
  if (!getColDB(0)) {
    endQueryDB();
    disconnectDB();
    cerr<<"monDB::Read(): "<<tbname<<", can't get calib file for "<<TimeStr(point)<<endl;
    return;
    //throw CS::Exception("monDB::Read():  TBname \"%s\",  can't get calib file for %s",
    //tbname.c_str(), TimeStr(point).c_str());
  }
  const string filename = directory + '/' + getColDB(0);

/*    do { */
/*      cerr << " plane "<<tbname <<" file "<<getColDB(0)<<endl; */
/*    } while (getNextRowDB()); */

  endQueryDB();
  ifstream f(filename.c_str());
  if( !f.is_open() ) {
    cerr<<"monDB::Read(): "<<tbname<<", can't open file "<<filename<<endl;
    return;
  }
  f>>v;
  if( !f.good() && !f.eof()) {
    cerr<<"monDB::Read(): "<<tbname<<", can't read from file "<<filename<<endl;
    return;
  }
  //disconnectDB();
}


////////////////////////////////////////////////////////////////////////////////
        string TimeStr(const tm &t);
        string TimeStrSQL(const tm &t);

        bool operator == (const tm &t1, const tm &t2);
inline  bool operator != (const tm &t1, const tm &t2) {return !(t1==t2);}
        bool operator >  (const tm &t1, const tm &t2);
inline  bool operator <  (const tm &t1, const tm &t2) {return t2>t1;}
inline  bool operator >= (const tm &t1, const tm &t2) {return t1==t2 || t1>t2;}
inline  bool operator <= (const tm &t1, const tm &t2) {return t1==t2 || t1<t2;}

////////////////////////////////////////////////////////////////////////////////

template <class T>
ostream& operator<<(ostream &o,const vector<T> &v)
{
  for( typename vector<T>::const_iterator it=v.begin(); it!=v.end(); it++ )
    if( !(o<<(*it)<<"\n") ) {
      cerr<<"monDB::operator<< ()  error in writing to the stream."<<endl;
      return o;
    }
  return o;
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
istream& operator>>(istream &in,vector<T> &v)
{
  v.clear();

  do
    v.push_back(T());
  while( in>>v.back() );

  v.pop_back();  // Remove last element becase it was not read.

  if( in.eof() )
    in.clear();

  return in;
}

////////////////////////////////////////////////////////////////////////////////

} /// using namespace CSMon
#endif // monDB___include
