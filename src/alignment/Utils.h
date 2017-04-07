// $Id: Utils.h,v 1.16 2010/08/27 08:46:14 tnagel Exp $

/*!
  \file    Utils.h
  \brief   some utilities used everywhere
  \author  Hugo Pereira
  \version $Revision: 1.16 $
  \date    $Date: 2010/08/27 08:46:14 $
*/
#ifndef Utils_h
#define Utils_h
#include <stdio.h>
#include <math.h>
#include<string>
#include<vector>
#include<TROOT.h>
#include<TObject.h>
#include<TMatrix.h>

//class ostream;
//class TMatrix;
class TCanvas;
class TH1;
class TH1D;

class TH2;
class TH2D;
/*! 
  \class Utils
  \brief some utilities used everywhere. All string manipulation methods are
  not accessible via root cint
*/

class Utils:public TObject {
public:
   
  /*! \fn static void  DumpMatrix( TMatrix m, std::ostream &out = cout )
    \brief dump root matrix to ostream in way I prefer to root way
    \param m   the matrix to be dumped
    \param out the output stream 
  */
  //  static void  DumpMatrix( TMatrix m, std::ostream &out = cout );
   static void  DumpMatrix( TMatrix m, std::ostream &out);
  /*! \fn static TMatrix  InvertMatrix( TMatrix m )
    \brief invert 2x2 root matrix (because of bug in root version 
    \param m   the matrix to be dumped
  */
  static TMatrix InvertMatrix( TMatrix m );
 
  /*! \fn static double GetSign( double val )
    \brief returns -1 for negative values, +1 for positives, 0 for 0
    \param val the double whose sign is to be retrieved.
  */
  static double GetSign( double val ) {
    if( val < 0 ) return -1.0;
    else if( val > 0 ) return +1.0;
    else return 0;
  }
   
  //! Retrieve histograms matching name from fileselection. Adds them, returns the sum
  static TH1* Get( const char* fileselection, const char* name )
  { return Get( GetFiles( fileselection ), name ); }
  
  //! Retrieve 2D histograms matching name from fileselection. Adds them, returns the sum
  static TH2* Get2D( const char* fileselection, const char* name ) 
  { return Get2D( GetFiles( fileselection ), name ); }
  
  /*! \fn static double HDiv(TH1* h1, TH1* h2, TH1* h3)
    \brief divides 2 histograms bin/bin; stores the result in an third histo; returns ratio of the number of entries in the 2 histo
    \param h1 probe histogram (denominator)
    \param h2 reference histogram (numerator)
    \param h3 output histogram h3[i] = h1[i]/h2[i]; \warning this parameter is changed
  */   
  static double HDiv(TH1* h1, TH1* h2, TH1* h3);
 
  /*! \fn static double HDiv(TH1* h1, TH1* h2, TH1* h3, unsigned int i1, unsigned int i2 )
    \brief divides 2 histograms bin/bin whith specified range; stores the result in an third histo; returns ratio of the number of entries in the 2 histo
    \param h1 probe histogram (denominator)
    \param h2 reference histogram (numerator)
    \param h3 output histogram h3[i] = h1[i]/h2[i]; \warning this parameter is changed
    \param i1 first bin taken into account
    \param i2 last bin taken into account
  */   
  static double HDiv(TH1* h1, TH1* h2, TH1* h3, unsigned int i1, unsigned int i2);
  
  /*! \fn static int GetColorCode( const char* TBName )
    \brief returns color according to detector name
    \param TBName name of the detector. 
  */   
  static int GetColorCode( const char* TBName );
  
    
  //! divide canvas in a rational way to have at least \param nPads pads, returns pointer to the canvas
  static TCanvas* DivideCanvas( TCanvas* cv, int nPads );
 
  #ifndef __CINT__
  // folowing functions have to be hidden from root Cint has it does not like stl vectors
  
  //! retrieve and sum all 1D histograms matching name in files
  static TH1* Get( std::vector<std::string> files, const char* name );

  //! retrieve and sum all 2D histograms matching name in files
  static TH2* Get2D( std::vector<std::string> files, const char* name );

  /*! \fn static string GetTimeStamp( string format="%d_%m_%y.%H_%M_%S" )
    \brief returns a string with date and hour embeded
    \param format the format of the string, based on the 'date' shell command
  */
  static std::string GetTimeStamp( std::string format="%d_%m_%y.%H_%M_%S" );
  
  /*! \fn static string CopyToLocal( string filename )
    \brief copy filename to local directory
    \param return local filename
  */
  static std::string CopyToLocal( std::string filename );

  /*! \fn static vector<string> GetFiles( string fileSelection )
    \brief returns a vector containing all files matching fileSelection using wildcards
    based on 'ls' shell command
    \param fileselection file selection criterion
  */  
  static std::vector<std::string> GetFiles( std::string fileSelection );
  
  /*! \fn static void  MakeBackup( string filename )
     \brief make a backup of a file if it exist, returns the name of the backup file
     the name of the backup file is <filename>.N where N is choosen so that
     the filename does not exist
     \param filename the name of the file for which the backup is done
  */
  static std::string  MakeBackup( std::string filename ); 
  
  /*! \fn static bool  Match( string s1, string s2 );
    \brief compares two strings using 'character wise' wild chards for the second
    \param s1 the first string. Wild cards are not accounted to.
    \param s2 the second string, which may contain wildcards
    example "AuuB" matches "A**B", "***B", "A*" and "A" but not "A*B" nor "*B" 
  */
  static bool  Match( std::string s1, std::string s2 );
  
  /*! \fn static string RemovePath( string filename )
    \brief removes the path name if any from filename. Returns the results
    \param the filename in which the pathname is looked for
  */
 static std::string RemovePath( std::string filename );  
 
 
 #endif  

 ClassDef(Utils,0)
 
};

#endif  
