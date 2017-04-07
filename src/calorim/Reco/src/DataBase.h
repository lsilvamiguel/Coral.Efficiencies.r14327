/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/DataBase.h,v $
   $Date: 2010/09/20 07:12:29 $
   $Revision: 1.14 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Author:
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

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

#ifndef RecoDataBase___include
#define RecoDataBase___include

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <typeinfo>

#include "Exception.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////


/*! \brief Simple filesystem-based database to store structured data together
  with validity ranges.

  One database corresponds to one directory in the file system.

*/

class DataBase
{
  //============================================================================
  // Types
  //============================================================================

  public:

    /// Open mode.
    enum OpenMode {READ=1,WRITE=2,CREATE=4};

    class ElementBase
    {
      public:
                       ~ElementBase             (void) {delete time_start; delete time_finish;}

                       ElementBase              (void) : time_start(NULL), time_finish(NULL) {}

                        ElementBase             (const std::string &name_,const tm *time_start_=NULL, const tm *time_finish_=NULL) :
                                                  name          (name_),
                                                  time_start    (NULL),
                                                  time_finish   (NULL)
                                                {
                                                  if( std::string::npos!=name.find('~') )
                                                    throw Exception("DataBase::Element::Element():  name can not contain character '~'  name=\"%s\"",
                                                                     name.c_str());

                                                  if( NULL!=time_start_ )
                                                    time_start = new tm(*time_start_);

                                                  if( NULL!=time_finish_ )
                                                    time_finish = new tm(*time_finish_);
                                                }

                        ElementBase             (const ElementBase &e) : time_start(NULL), time_finish(NULL) {*this=e;}

        ElementBase    &operator=               (const ElementBase &e)
                                                {
                                                  if( &e!=this )
                                                  {
                                                    name = e.name;
                                                    delete time_start;  time_start =NULL;
                                                    delete time_finish; time_finish=NULL;

                                                    if( NULL!=e.time_start )
                                                      time_start = new tm(*e.time_start);

                                                    if( NULL!=e.time_finish )
                                                      time_finish = new tm(*e.time_finish);
                                                  }
                                                  return *this;
                                                }

        void            Print                   (std::ostream &o=std::cout) const;

        const std::string   &GetName                 (void) const {return name;}

        std::string          CreateFileName          (void) const;

        static
        ElementBase     DecodeFullName          (const std::string &full_name);

        static void     DecodeTime              (const std::string &s,tm &t);

        void            ClearTimeStart          (void) {delete time_start;  time_start =NULL;}

        void            ClearTimeFinish         (void) {delete time_finish; time_finish=NULL;}

        tm             *GetTimeStart            (void)       {return time_start;}
        const tm       *GetTimeStart            (void) const {return time_start;}

        tm             *GetTimeFinish           (void)       {return time_finish;}
        const tm       *GetTimeFinish           (void) const {return time_finish;}

        void            SetTimeStart            (const tm &t);

        void            SetTimeFinish           (const tm &t) { if(time_finish==NULL) time_finish=new tm(t); else *time_finish=t; }

      protected:

        std::string          name;
        tm             *time_start;
        tm             *time_finish;
    };

  ////////////////////////////////////////////////////////////////////////////////

  public:

    template <class T>
    class Element : public ElementBase
    {
      public:
        typedef T ElementType;

        virtual        ~Element                 (void) {delete data;}

                        Element                 (const std::string &name,const tm *time_start=NULL, const tm *time_finish=NULL) :
                                                 ElementBase(name,time_start,time_finish), data(NULL) {}

                        Element                 (const std::string &name,const T &v,tm *time_start=NULL, tm *time_finish=NULL) :
                                                 ElementBase(name,time_start,time_finish), data(new T(v)) {}

                        Element                 (const Element &e) : data(NULL) {*this=e;}

                        Element                 (const T &v) : data(new T(v)) {}

          Element      &operator=               (const Element &e)
                                                {
                                                  if( &e!=this )
                                                  {
                                                    ElementBase::operator=(e);

                                                    if( e.data!=NULL )
                                                      if( data!=NULL )
                                                        *data = *e.data;
                                                      else
                                                        data = new T(*e.data);
                                                    else
                                                    {
                                                      delete data;
                                                      data = NULL;
                                                    }
                                                  }
                                                  return *this;
                                                }

          Element      &operator=               (const T &v) {SetData(v);return *this;}

          T            *GetData                 (void)       {return data;}

          const T      *GetData                 (void) const {return data;}

          void          SetData                 (const T &v) { if(data==NULL) data=new T(v); else *data=v; }

        protected:

          T              *data;
    };

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
    virtual            ~DataBase                (void) {}


    /*! \brief Create data base

       open_mode=DataBase::WRITE|DataBase::READ open existed DB, or to create new DB
       open_mode=DataBase::WRITE|DataBase::READ|DataBase::CREATE to create new DB
           (and remove old one)
    */
                        DataBase                (const std::string &path, int open_mode=READ);

    /// Copy constructor.
                        DataBase                (const DataBase &db) : open_mode(db.open_mode) {*this=db;}

  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator
    DataBase           &operator =              (const DataBase &db);

  //============================================================================
  // Methods
  //============================================================================

  public:

    bool                WasCreated              (void) const {return creation_flag;}

    /// \return Name of this cell type.
    const std::string        GetPath                 (void) const {return path;}

    template <class T>
    void                Write                   (const Element<T> &e);

    void                FindAllVersions         (const std::string &name,std::vector<DataBase::ElementBase> &versions) const;

    template <class T>
    void                Read                    (Element<T> &e) const;

    template <class T>
    void                Write                   (const std::string &container,const T& v, const tm &start, const tm &end);

    template <class T>
    void                Read                    (const std::string &container,T& v, const tm &point);

    std::string              CalibFilePath           (const std::string &container, const tm &point);

    static void         PrintTime               (const tm &t,std::ostream &o=std::cout);

    static std::string       TimeStr                 (const tm &t);

  //============================================================================
  // Attributes, data
  //============================================================================

  private:

    const std::string        path;

    const int           open_mode;

    bool                creation_flag;
};

////////////////////////////////////////////////////////////////////////////////

template <class T>
void DataBase::Read(Element<T> &e) const
{
  if( !(open_mode&READ) )
    throw Exception("DataBase::Read(): data base \"%s\" is not opened for reading",
                     path.c_str());

  // Create file name.
  const std::string file_name = path+'/'+e.CreateFileName();
  std::ifstream f(file_name.c_str());
  if( !f.is_open() )
    throw Exception("DataBase::Read():  element \"%s\",  can not open file \"%s\"",
                     e.GetName().c_str(), file_name.c_str());

//   // Read variable type.
//   std::string type_name;
//   if( !(f>>type_name) )
//     throw Exception("DataBase::Read():  element \"%s\"  error(1) in reading from file \"%s\"",
//                      e.GetName().c_str(), file_name.c_str());
//
//   // Check that types are matched.
//   if( type_name!=typeid(T).name() )
//     throw Exception("DataBase::Read():  element \"%s\"  error in reading from file \"%s\"\n"
//                     "                   types are mismatched:  element type is %s     type in file is %s\n",
//                      e.GetName().c_str(), file_name.c_str(),typeid(T).name(),type_name.c_str() );

  // Create space for data, if needed.
  if( NULL==e.GetData() )
    e.SetData(T());

  // Read data.
  f >> (*e.GetData());
//   if( !(f>>(*e.GetData())) )
//     throw Exception("DataBase::Read():  element \"%s\"  error(2) in reading from file \"%s\"",
//                      e.GetName().c_str(), file_name.c_str());
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void DataBase::Write(const Element<T> &e)
{
  if( !(open_mode&WRITE) )
    throw Exception("DataBase::Write():  database \"%s\" was opened in read-only mode = %d",
                     path.c_str(),open_mode);

  if( NULL==e.GetData() )
    throw Exception("DataBase::Write():  element \"%s\"  has no data",e.GetName().c_str());

  const std::string file_name = path+'/'+e.CreateFileName();
  std::ofstream f(file_name.c_str());

  if( !f.is_open() )
    throw Exception("DataBase::Write():  element \"%s\"  can not open file \"%s\"",
                     e.GetName().c_str(), file_name.c_str());

//   if( !( f << typeid(T).name() << "\n" << (*e.GetData()) ) )
  if( !( f << (*e.GetData()) ) )
    throw Exception("DataBase::Write():  element \"%s\"  error in writing to file \"%s\"",
                     e.GetName().c_str(), file_name.c_str());
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void DataBase::Write(const std::string &container,const T& v, const tm &start, const tm &end)
{
  if( !(open_mode&WRITE) )
    throw Exception("DataBase::Write():  database \"%s\" was opened in read-only mode = %d",
                     path.c_str(),open_mode);

  const std::string file_name = path+'/'+Element<T>(container,&start,&end).CreateFileName();
  std::ofstream f(file_name.c_str());

  if( !f.is_open() )
    throw Exception("DataBase::Write():  container \"%s\"  can not open file \"%s\"",
                     container.c_str(), file_name.c_str());

  if( !(f<<v) )
    throw Exception("DataBase::Write():  container \"%s\"  error in writing to file \"%s\"",
                     container.c_str(), file_name.c_str());
}

////////////////////////////////////////////////////////////////////////////////

        bool operator == (const tm &t1, const tm &t2);
inline  bool operator != (const tm &t1, const tm &t2) {return !(t1==t2);}
        bool operator >  (const tm &t1, const tm &t2);
inline  bool operator <  (const tm &t1, const tm &t2) {return t2>t1;}
inline  bool operator >= (const tm &t1, const tm &t2) {return t1==t2 || t1>t2;}
inline  bool operator <= (const tm &t1, const tm &t2) {return t1==t2 || t1<t2;}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void DataBase::Read(const std::string &container,T& v, const tm &point)
{
  if( !(open_mode&READ) )
    throw Exception("DataBase::Read(): data base \"%s\" is not opened for reading",
                     path.c_str());

  std::vector<ElementBase> vers;
  FindAllVersions(container,vers);
  for( size_t i=0; i<vers.size(); i++ )
  {
    ElementBase &e = vers[i];
    if( (e.GetTimeStart ()==NULL || *e.GetTimeStart ()<=point) &&
        (e.GetTimeFinish()==NULL || *e.GetTimeFinish()>=point) )
    {
      const std::string file_name = path+'/'+Element<T>(container,e.GetTimeStart(),e.GetTimeFinish()).CreateFileName();
      std::ifstream f(file_name.c_str());
      if( !f.is_open() )
        throw Exception("DataBase::Read():  container \"%s\",  can not open file \"%s\"",
                         container.c_str(), file_name.c_str());
      f>>v;
      if( !f.good() )
        throw Exception("DataBase::Read():  conatiner \"%s\":  can not read from file \"%s\"",
                         container.c_str(), file_name.c_str());
      return;
    }
  }

  throw Exception("DataBase::Read():  can not find container \"%s\" in given time point.",
                   container.c_str());
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
std::ostream& operator<<(std::ostream &o,const std::vector<T> &v)
{
  for( typename std::vector<T>::const_iterator it=v.begin(); it!=v.end(); it++ )
    if( !(o<<(*it)<<"\n") )
      throw Exception("operator<<()  error in writing to the stream.");
  return o;
}

////////////////////////////////////////////////////////////////////////////////

#ifndef OPERATOR_ISTR_VECTOR_T
#define OPERATOR_ISTR_VECTOR_T

template <class T>
std::istream& operator>>(std::istream &in,std::vector<T> &v)
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

#endif


////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco
#endif // RecoDataBase___include
