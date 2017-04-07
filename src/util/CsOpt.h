// $Id: CsOpt.h,v 1.12 2011/02/03 14:20:50 tnagel Exp $

/*!
   \file    CsOpt.h
   \brief   Compass Options Interpreter Class.
   \author  Benigno Gobbo 
   \version $Revision: 1.12 $
   \date    $Date: 2011/02/03 14:20:50 $
*/

#ifndef CsOpt_h
#define CsOpt_h

#include "coral_config.h"
#include <sstream>
#include <string>
#include <list>
#include <vector>


/*! \class CsOpt CsOpt.h
    \brief  Compass Options Interpreter Class. 

    It takes argc and argv from main program and interpretes them.
    It assumes that (at least) an option file is given as argument.
    Polymorphic getOpt methods are available to retrieve information 
    from input file. 
    */

class CsOpt {

 public:

  /*! \fn static CsOpt* Instance( int argc, char **argv );;
    \brief First singleton instantiation.
    At the moment a single argument is accepted and mandatory: the options
    file name. Only the -h option is accepted, if specified Coral prints a 
    short help message and exits.
    \param argc parameter from main function
    \param argv parameter from main function
   */
  static CsOpt* Instance( int argc, char **argv );

  /*! \fn static CsOpt* Instance();
    \brief singleton instantiation (but first).
  */
  static CsOpt* Instance();

  /*! \fn void dump( int dumplevel );
      \brief dumps the list of strings read from input file, and eventally
      the list of unused lines and the list of lines with errors inside.
      \param dumplevel if 0 only the list of read lines is printed;
      if 1 also the list of unused lines is printed; if greather that 1
      also the list of lines with errors is printed.
  */
  void dump( int dumplevel );

  /*! \fn bool getOpt( std::string tag, std::string key );   
      \brief looks for a line with tag+key string. Returns \c TRUE if string
      was found.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
  */
  bool getOpt( std::string tag, std::string key );   


  /*! \fn bool getOpt( std::string tag, std::string key, int& inum ); 
      \brief looks for a line with tag+key string ad get the following
      integer number. 

      Returns \c TRUE if string and integer number were found 
      and the integer number is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param imun reference to the variable where to write the read
      integer number.
  */
  bool getOpt( std::string tag, std::string key, int& inum ); 

  /*! \fn bool getOpt( std::string tag, std::string key, bool& bnum ); 
      \brief looks for a line with tag+key string ad get the following
      boolean number. 

      Returns \c TRUE if string and boolean number were found 
      and the boolean is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param bmun reference to the variable where to write the read
      boolean number.
  */
  bool getOpt( std::string tag, std::string key, bool& bnum ); 

  /*! \fn bool getOpt( std::string tag, std::string key, double& rnum ); 
      \brief looks for a line with tag+key string ad get the following
      real number. 

      Returns \c TRUE if string and double precision real number were found 
      and the number is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param rmun reference to the variable where to write the read
      double precision real number.
  */
  bool getOpt( std::string tag, std::string key, double& rnum ); 

  /*! \fn bool getOpt( std::string tag, std::string key, float& rnum ); 
      \brief looks for a line with tag+key string ad get the following
      real number. 

      Returns \c TRUE if string and single precision real number were found 
      and the number is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param rmun reference to the variable where to write the read
      single precision real number.
  */
  bool getOpt( std::string tag, std::string key, float& rnum ); 

  /*! \fn bool getOpt( std::string tag, std::string key, std::string& str ); 
      \brief looks for a line with tag+key string ad get the following
      character string. 

      Returns \c TRUE if strings and were found 
      and the output string is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param str reference to the variable where to write the read
      character string.
  */
  bool getOpt( std::string tag, std::string key, std::string& str );

  /*! \fn bool getOptRec( std::string tag, std::string key, std::string& str ); 
      \brief looks recursively for a line with tag+key string ad get the 
      following character string. 

      Returns \c TRUE if strings and were found 
      and the output string is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param str reference to the variable where to write the read
      character string.
  */
  bool getOptRec( std::string tag, std::string key, std::string& str );


  /*! \fn bool getOpt( std::string tag, std::string key, std::list<int>& numbers );
      \brief looks for a line with tag+key string ad get the following
      list of integer numbers. 

      Returns \c TRUE if string and at least an 
      integer number were found and the list of integer numbers is read 
      correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param numbers reference to the variable where to write 
      the read list of integer numbers.
  */
  inline bool getOpt( std::string tag, std::string key, std::list<int>& numbers )
    { return _myGetOpt( tag, key, numbers ); }

  /*! \fn bool getOpt( std::string tag, std::string key, std::list<double>& numbers );
      \brief looks for a line with tag+key string ad get the following
      list of real numbers. 

      Returns \c TRUE if string and at least a double precision
      real number were found and the list of numbers is read 
      correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param numbers reference to the variable where 
      to write the read list of numbers.
  */
  inline bool getOpt( std::string tag, std::string key, std::list<double>& numbers )
    { return _myGetOpt( tag, key, numbers ); }

  /*! \fn bool getOpt( std::string tag, std::string key, std::list<float>& numbers );
      \brief looks for a line with tag+key string ad get the following
      list of real numbers. 

      Returns \c TRUE if string and at least a single precision 
      real number were found and the list of numbers is read 
      correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param numbers reference to the variable where 
      to write the read list of numbers.
  */
  inline bool getOpt( std::string tag, std::string key, std::list<float>& numbers )
    { return _myGetOpt( tag, key, numbers ); }

  /*! \fn bool getOpt( std::string tag, std::string key, std::list<std::string>& strlist );
      \brief looks for a line with tag+key string ad get the following
      list of character strings. 

      Returns \c TRUE if tag, key and at least
      one more string were found and the list of strings is read 
      correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param strlist reference to the variable where 
      to write the read list of character strings.
      \param cumulative if true, instead of returning strings only from the
      last matching line, strings from all matching lines are returned
  */
  bool getOpt( std::string tag, std::string key, std::list<std::string>& strlist, bool cumulative=false );

  /*! \fn bool getOptRec( std::string tag, std::string key, std::list<std::string>& strlist );
      \brief looks recursively for a line with tag+key string ad get the 
      following list of character strings. 

      Returns \c TRUE if tag, key and at least
      one more string were found and the list of strings is read 
      correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param strlist reference to the variable where 
      to write the read list of character strings.
  */
  bool getOptRec( std::string tag, std::string key, std::list<std::string>& strlist );


  /*! \fn bool getOpt( std::string tag, std::string key, std::vector<int>& numbers );
      \brief looks for lines with tag+key string. In addition it expects
      a list of vector components and a list of integer values. 
      e.g. "mytag mykey [1,4,6-8]  0 7 10 11 12".
      
      Returns \c TRUE if the tag and key string, the list of vector 
      components and the list of integer values were found and everything 
      is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param numbers reference to the variable where 
      to write the read vector elements.
  */
  inline bool getOpt( std::string tag, std::string key, std::vector<int>& numbers )
    { return _myGetOpt( tag, key, numbers ); }

  /*! \fn bool getOpt( std::string tag, std::string key, std::vector<double>& numbers );
      \brief looks for lines with tag+key string. In addition it expects
      a list of vector components and a list of double precision real values. 
      e.g. "mytag mykey [1,4,6-8]  0 7 10 11 12".
      
      Returns \c TRUE if the tag and key string, the list of vector 
      components and the list of double precision real values were found 
      and everything is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param numbers reference to the variable where 
      to write the read vector elements.
  */
  inline bool getOpt( std::string tag, std::string key, std::vector<double>& numbers )
    { return _myGetOpt( tag, key, numbers ); }

  /*! \fn bool getOpt( std::string tag, std::string key, std::vector<float>& numbers );
      \brief looks for lines with tag+key string. In addition it expects
      a list of vector components and a list of single precision real values. 
      e.g. "mytag mykey [1,4,6-8]  0 7 10 11 12".
      
      Returns \c TRUE if the tag and key string, the list of vector 
      components and the list of single precision real values were 
      found and everything is read correctly.
      \param tag is a (list of) string(s) that should identify the 
      calling package.
      \param key is a (list of) string(s) to be used to identify the
      opportune line in the input file.
      \param numbers reference to the variable where 
      to write the read vector elements.
  */
  inline bool getOpt( std::string tag, std::string key, std::vector<float>& numbers )
    { return _myGetOpt( tag, key, numbers ); }

  // undocumented...
  inline bool getOptRec( std::string tag, std::string key, std::vector<float>& numbers )
    { return _myGetOptRec( tag, key, numbers ); }

  /* \fn void setIntMode( std::string mode );
     \brief allows Oct, Hex and Dec integer number format
     \param mode dec, hex or oct.
  */
  void setIntMode( std::string mode );

  /* \fn static bool expand( std::string filename );
     \brief expands a string (expands environments). if the 
     string contains spaces, it is untouched. returns true if 
     changes are made.
     \param filename the string with the file path
  */
  static bool expand( std::string& filename ) ;
 
  /* \fn static bool regexpMatch(const char* string, const char* regexp);
     \brief match string with a regular expression
  */
  static bool regexpMatch(const char* string, const char* regexp);

 protected:
  CsOpt();
  CsOpt( int argc, char **argv );

 private:
  static CsOpt* _instance; 

  std::list<std::string> _options;
  std::list<std::string> _unused;
  std::list<std::string> _wrong;

  int _intMode; // 10: dec, 16: hex, 8: oct

  void _usage( char* cp );
  void _getOptions( const char* optfile );
  void _clearComments( char* line );
  void _clearFinalBlanks( char* line );
  bool _checkTags( std::string tag, std::string key, std::istringstream& is );
  template <class T> bool _getNumber( std::istringstream& is, T& num ); 
  template <class T> bool _getNumbersList( std::istringstream& is, 
					   std::list<T>& numbers );
  bool _getElementsList( char* text, std::list<int>& elem );
  template <class T> bool _myGetOpt(std::string tag, std::string key, 
				    std::list<T>& numbers);
  template <class T> bool _myGetOpt(std::string tag,std::string key,
				    std::vector<T>& numbers);
  template <class T> bool _myGetOptRec( std::string tag, std::string key,
					std::vector<T>& numbers );
};
#endif // CsOpt_h





