/*
   $Source: RunDataBase.h,v $ 
   $Date: $ 
   $Revision:  $ 
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

     Modified by V.Kolosov ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     from time to run DataBase   Mon Feb 11 2002 

*/

#ifndef MNRunDataBase___include
#define MNRunDataBase___include

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <typeinfo>

#include "Reco/Exception.h"

namespace MN {

////////////////////////////////////////////////////////////////////////////////


/*! \brief Data base.

*/
class RunDataBase
{
  //============================================================================
  // Types
  //============================================================================

  public:

    /// Open mode.
    enum OpenMode {READ=1,WRITE=2,CREATE=4};
    /// DB structure.
    enum StructureDB {PLANE=1,DETSUBDIRS=2,DETCALIBSUBDIRS=3};
    
    const static unsigned RUN_MAX = 9999999;
  
  public:
  
   
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor
    virtual            ~RunDataBase                (void) {}


    /*! \brief Create data base
    
       open_mode=RunDataBase::WRITE|RunDataBase::READ open existed DB, or to create new DB
       open_mode=RunDataBase::WRITE|RunDataBase::READ|RunDataBase::CREATE to create new DB 
           (and remove old one)
    */
                        RunDataBase                (const std::string &path, int open_mode=READ);

    /// Copy constructor.
                        RunDataBase                (const RunDataBase &db) {*this=db;}

  //============================================================================
  // Operators
  //============================================================================

  private:

    /// Assignment operator
    RunDataBase           &operator =              (const RunDataBase &db);

  //============================================================================
  // Methods
  //============================================================================

  public:

    bool                WasCreated              (void) const {return creation_flag;}

    /// \return Name of this cell type.
    const std::string        GetPath                 (void) const {return path;}
     
    void                SetStructureID ( StructureDB structure_id ) { structure_id_ = structure_id; }

    StructureDB         GetStructureID ( void ) const { return structure_id_; }
    
    template <class T>
    void                Write                   (const std::string &container, const std::string &tag, const T& v, const size_t &start, const size_t &end);

    template <class T>
    bool                Read                    (const std::string &container, const std::string &tag, T& v, const size_t &point);

  private:
    
    void                FindAllGoodVersions         (const std::string &container, const std::string &tag, std::vector< std::string > &versions, const size_t &point ) const;

    void                CreateStructure          (const std::string &container, const std::string &tag, const size_t &start, const size_t &end) const;

    std::string          CreateFileName           (const std::string &container, const std::string &tag, const size_t &start, const size_t &end) const;

    std::string          CreateFileNameLast       (const std::string &container, const std::string &tag, const size_t &start, const size_t &end) const;

    std::string          CreateTimeStamp         ( const size_t &start, const size_t &end) const;

    std::string          CreatePathName           (const std::string &container, const std::string &tag) const;

    bool                 FileNameSyntaxOK         (const std::string &file_name, const std::string &container, const std::string &tag, std::pair <size_t, size_t> &time_stamps ) const;
    
    static     bool     DecodeRun                (const std::string &s,size_t &t);
    
    static void         PrintRun                (const size_t &t,std::ostream &o=std::cout);
    
    static std::string       RunStr                  (const size_t &t);

  //============================================================================
  // Attributes, data
  //============================================================================

  private:
    
    const std::string        path;
    
    int                  open_mode;

    bool                 creation_flag;
    
    StructureDB          structure_id_;
};

////////////////////////////////////////////////////////////////////////////////
    
template <class T>
void RunDataBase::Write(const std::string &container,const std::string &tag, const T& v, const size_t &start, const size_t &end)
{
       try
       {
  if( !(open_mode&WRITE) )
    throw Reco::Exception("RunDataBase::Write():  RunDataBase \"%s\" was opened in read-only mode = %d",
                          path.c_str(),open_mode);

  CreateStructure(container, tag, start, end);
  const std::string file_name = CreateFileName(container, tag, start, end);

  std::ofstream f(file_name.c_str());

  if( !f.is_open() )
  {
    throw Reco::Exception("RunDataBase::Write():  container \"%s\"  can not open file \"%s\"",
                          container.c_str(), file_name.c_str());
  }

  if( !(f<<v) )
    throw Reco::Exception("RunDataBase::Write():  container \"%s\"  error in writing to file \"%s\"",
                          container.c_str(), file_name.c_str());
       }
       catch( std::exception &e )
       {
        std::cerr << "Exception:\n" << e.what() << "\n";
       }

// Keep last information in the file

       try
       {
  if( !(open_mode&WRITE) )
    throw Reco::Exception("RunDataBase::Write():  RunDataBase \"%s\" was opened in read-only mode = %d",
                          path.c_str(),open_mode);

  const std::string file_name = CreateFileNameLast(container, tag, start, end);

  std::ofstream f(file_name.c_str());

  if( !f.is_open() )
    throw Reco::Exception("RunDataBase::Write():  container \"%s\"  can not open file \"%s\"",
                          container.c_str(), file_name.c_str());

  if( !(f<<v) )
    throw Reco::Exception("RunDataBase::Write():  container \"%s\"  error in writing to file \"%s\"",
                          container.c_str(), file_name.c_str());
       }
       catch( std::exception &e )
       {
        std::cerr << "Exception:\n" << e.what() << "\n";
       }
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
bool RunDataBase::Read(const std::string &container,const std::string &tag,T& v, const size_t &point)
{
       try
       {
  if( !(open_mode&READ) )
    throw Reco::Exception("RunDataBase::Read(): data base \"%s\" is not opened for reading",
                          path.c_str());

  std::vector< std::string > vers;
//  std::cout << " RunDataBase::Read() Try to find All versions " <<  std::endl;
  FindAllGoodVersions(container, tag, vers, point);
//  std::cout << " RunDataBase::Read() versions " <<  vers.size() << std::endl;
  if( vers.size() > 1 ) 
    throw Reco::Exception("RunDataBase::Read(): Ambigious and so not consistent  \"%s\" in given time point.",
                          container.c_str());
  for( size_t i=0; i<vers.size(); i++ )
  {
    const std::string file_name = vers[i];
//      const std::string file_name = CreateFileName(container, tag, e.GetRunStart(),e.GetRunFinish());
//  std::cout << " Try to open file " <<  file_name << std::endl;
      std::ifstream f(file_name.c_str());
      if( !f.is_open() )
        throw Reco::Exception("RunDataBase::Read():  container \"%s\",  can not open file \"%s\"",
                              container.c_str(), file_name.c_str());
      f>>v;
      if( !f.good() )
        throw Reco::Exception("RunDataBase::Read():  conatiner \"%s\":  can not read from file \"%s\"",
                              container.c_str(), file_name.c_str());
      return true;
  }
  
  throw Reco::Exception("RunDataBase::Read():  can not find container \"%s\" tag \"%s\" in given time point.",
                        container.c_str(),tag.c_str());
       }
       catch( std::exception &e )
       {
         
         std::cerr << "Exception:\n" << e.what() << "\n";
         return false;
       }
}

////////////////////////////////////////////////////////////////////////////////
// 
//         bool operator == (const tm &t1, const tm &t2);
// inline  bool operator != (const tm &t1, const tm &t2) {return !(t1==t2);}
//         bool operator >  (const tm &t1, const tm &t2);
// inline  bool operator <  (const tm &t1, const tm &t2) {return t2>t1;}
// inline  bool operator >= (const tm &t1, const tm &t2) {return t1==t2 || t1>t2;}
// inline  bool operator <= (const tm &t1, const tm &t2) {return t1==t2 || t1<t2;}
// 
////////////////////////////////////////////////////////////////////////////////

template <class T>
std::ostream& operator<<(std::ostream &o,const std::vector<T> &v)
{
  for( typename std::vector<T>::const_iterator it=v.begin(); it!=v.end(); it++ )
    if( !(o<<(*it)<<"\n") )
      throw Reco::Exception("operator<<()  error in writing to the stream.");
  return o;
}

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

} /// using namespace MN
#endif // MNRunDataBase___include
