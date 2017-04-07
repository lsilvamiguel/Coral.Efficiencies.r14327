// $Id: Opt.cc,v 1.8 2010/04/26 20:06:45 tnagel Exp $

/*!
   \file    Opt.cc
   \brief   Compass Options Interpreter Class.
   \author  Benigno Gobbo 
   \version $Revision: 1.8 $
   \date    $Date: 2010/04/26 20:06:45 $
*/

#if defined(__GNUG__) && defined(__i386__)
#  include <getopt.h>
#  include <stdio.h>
#endif

#include <wordexp.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <regex.h>

#ifdef __SUNPRO_CC
#  include <strings.h>
#  include <stdio.h>
#endif

#include "Opt.h"

#include <math.h>

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <typeinfo>

using namespace std;

//________________________________________
Opt* Opt::_instance = NULL;

//________________________________________
Opt* Opt::Instance() {
 if( _instance != 0 ) 
   return _instance;
 else {
   cerr << "Opt FATAL: wrong singleton usage." << endl;
   exit(1);
 }
}


//________________________________________
Opt* Opt::Instance( int argc, char **argv ) {
 if( _instance == NULL ) {
   _instance = new Opt( argc, argv );
 }
 return _instance; 
}

Opt::Opt( int argc, char **argv ) {
  
  int rc;
  char* cp;
  char optfile[256] = "";
  setIntMode( "dec" );

  // Program name ...
  if(( cp = (char*) strrchr( argv[0], '/' )) != NULL )
    ++cp;
  else
    cp = argv[0];
   
  // Arguments and Options
  if( argc > 1 ) {
    while(( rc = getopt( argc, argv, "hrq" )) != EOF ) {
      switch( rc ) {
      case 'h':
        _usage(cp);
	exit(0);
	break;
      default:
	break;
      }
    }  
  }

  for( ; optind<argc; optind++ ) {
    if( access( argv[optind], F_OK ) == 0 ) {
      strcpy( optfile, argv[optind] );
    }
  }

  if( access( optfile, F_OK ) == -1 ) {
    perror( cp );
    _usage( cp );
    exit( 1 );
  }

  _options.clear();
  _unused.clear();
  _getOptions( optfile );

}

//________________________________________
void Opt::_usage( char* cp ) {
  cout << "\n usage: " << cp << " [ -h ] <filename> \n\n"
"   Initializer: \n"
"           -h this help \n\n" << endl;
}

//________________________________________
void Opt::_clearComments( char* line ) {
  // It looks for "//" and replaces it with "\0"...
  regex_t re;
  regmatch_t pmatch;
  (void) regcomp( &re, "//", 0 );
  int i = regexec( &re, line, (size_t) 1, &pmatch, 0 ); 
  regfree( &re );
  if( i == 0 ) line[pmatch.rm_so] = '\0';
  // replace tabs with spaces...
  //for( unsigned int j=0; j<strlen(line); j++ ) 
    //if( line[j] == '\t' ) line[j] = ' ';
}

//________________________________________
void Opt::_clearFinalBlanks( char* line ) {
  // From end of string, find last ' ' or '\t' (tab) and replece it with '\0'
  int i;
  for( i=strlen(line)-1; i>=0 && (line[i]==' '||line[i]=='\t'); i-- );
  if( line[i+1]  == ' ' || line[i+1] == '\t' ) line[i+1] = '\0';
}

//________________________________________
void Opt::_getOptions( const char* optfile ) {

  const int lineSize = 512; 
  char line[lineSize];
  string str;
  string key = "include"; 

  // Read all lines in option file. Remove comments and final blanks and tabs

  fstream f( optfile, ios::in );
  if( f.fail() ) {
    cerr << "Opt FATAL: input file " << optfile << " not found." << endl;
    exit(1);
  }

  while( f.getline( line, lineSize ) ) {

    _clearComments( line );
    _clearFinalBlanks( line );

    // look for included files...    
    istringstream is(line);
    is >> str;
    if( str == key ) {
      is >> str; 
      // check existence of the file to be included...
      wordexp_t* exp = (wordexp_t*) malloc(512);
      if( wordexp( str.c_str(), exp, 0 ) == 0 ) {
	char* name = exp->we_wordv[0];
	if( access( name, R_OK ) == 0 ) { 
	  _getOptions( name );
	}
	else {
	  cerr << "Opt FATAL: input file " << name << " not found." << endl;
	  exit(1);
	}
      }
      else {
	cerr << "Opt FATAL: error in include statement syntax" << endl;
	exit(1);
      }
      wordfree( exp );
    }
    else {
      str = line;
      if( !str.empty() ) {
	_options.push_back( str );
	_unused.push_back( str );
      }
    }
  }
}

//________________________________________
void Opt::dump( int dumplevel ) {

  // Obvious, isn't it?
  list<string>::iterator is;
  cout << "----------------------" << endl;
  cout << "   Option File Dump   " << endl;
  cout << "----------------------" << endl;
  for(is=_options.begin(); is!=_options.end(); is++ )
    cout << *is << endl;
  if( dumplevel > 0 ) {
    cout << "----------------------" << endl;
    cout << "   not used options   " << endl;
    cout << "----------------------" << endl;
    for(is=_unused.begin(); is!=_unused.end(); is++ )
      cout << *is << endl;
    if( dumplevel > 1 ) {
      cout << endl;
      cout << "----------------------" << endl;
      cout << "     wrong options    " << endl;
      cout << "----------------------" << endl;
      for(is=_wrong.begin(); is!=_wrong.end(); is++ )
	cout << *is << endl;
      cout << endl;
    }
  }
}


//________________________________________
bool Opt::getOpt( string tag, string key ) {
  // search for a line with tag and key at the beginning
  bool status = false;
  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is((*i).c_str());
    if( _checkTags( tag, key, is ) ) {
      _unused.remove( *i );
      status = true;
    }
  }
  return( status );
}


//________________________________________
bool Opt::getOpt( string tag, string key, int& num ) {
  // read an integer number following tag and key strings in a line
  bool status = false;
  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      status = _getNumber( is, num );
      _unused.remove( *i );
      if( !status ) {
	_wrong.push_back( *i );
      }
    }
  }
  return( status );
}


//________________________________________
bool Opt::getOpt( string tag, string key, double& num ) {
  // read a real number following tag and key strings in a line
  bool status = false;
  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      status = _getNumber( is, num );
      _unused.remove( *i );
      if( !status ) {
	_wrong.push_back( *i );
      }
    }
  }
  return( status );
}

//________________________________________
bool Opt::getOpt( string tag, string key, float& num ) {
  // read a real number following tag and key strings in a line
  bool status = false;
  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      status = _getNumber( is, num );
      _unused.remove( *i );
      if( !status ) {
	_wrong.push_back( *i );
      }
    }
  }
  return( status );
}


//________________________________________
bool Opt::getOpt( string tag, string key, string& str ) {
  // read a string following tag and key strings in a line
  bool status = false;
  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      is >> str;
      if( str.size() != 0 ) status = true;
      _unused.remove( *i );
      if( !status ) {
	_wrong.push_back( *i );
      }
    }
  }
  return( status );
}

//________________________________________
bool Opt::getOptRec( string tag, string key, string& str ) {
  // read a string following tag and key strings in a line
  bool status = false;
  list<string>::iterator i;
  for(i=_unused.begin(); i!=_unused.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      is >> str;
      if( str.size() != 0 ) status = true;
      _unused.remove( *i );
      if( !status ) {
	_wrong.push_back( *i );
      }
      return( status );
    }
  }
  return( status );
}


//________________________________________
bool Opt::getOpt( string tag, string key, list<string>& strlist ) {
  // read a list of strings following tag and key strings in a line
  bool status = false;
  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      strlist.clear();
      status = true;
      while ( !is.eof() && status ) {
	string str;
	is >> str;
	strlist.push_back( str );
      }
      _unused.remove( *i );
      if( strlist.empty() ) {
	_wrong.push_back( *i );
	status = false;
      }
    }
  }
  return( status );
  } 


//________________________________________
bool Opt::getOptRec( string tag, string key, list<string>& strlist ) {
  // read a list of strings following tag and key strings in a line
  bool status = false;
  list<string>::iterator i;
  for(i=_unused.begin(); i!=_unused.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      strlist.clear();
      status = true;
      while ( !is.eof() && status ) {
	string str;
	is >> str;
	strlist.push_back( str );
      }
      _unused.remove( *i );
      if( strlist.empty() ) {
	_wrong.push_back( *i );
	status = false;
      }
      return( status );
    }
  }
  return( status );
  } 


//________________________________________
template <class T> 
bool Opt::_myGetOpt( string tag, string key, list<T>& Tlist ) {

  // read a list of numbers (type T) following tag and key strings in a line
  bool status = false;
  Tlist.clear();
  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      status = _getNumbersList( is, Tlist );
      _unused.remove( *i );
      if( !status ) {
	_wrong.push_back( *i );
      }
    }
  }
  return( status );
}


//________________________________________
template <class T> bool Opt::_myGetOpt( string tag, string key, vector<T>& numbers ) {
  // read a vector of numbers (type T) following tag and key strings in a line
  bool status = false;
  numbers.clear();

  list<int> elem;
  list<T> values;

  list<string>::iterator i;
  for(i=_options.begin(); i!=_options.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      char* text = new char[(*i).size()];
      is.getline( text, 512, '\0' );
      // get the list of vector elements
      status = _getElementsList( text, elem );
      // get the corresponding list of numbers
      istringstream is( text );
      status = _getNumbersList( is, values );
      _unused.remove( *i );
      if( status ) {
	// the number of elements must be equal to the number of numbers
	if( elem.size() != values.size() ) {
	  status = false;
	  _wrong.push_back( *i );
	  numbers.clear();
	}
	else { // store numbers in the right positions
	  list<int>::iterator j;
	  typename list<T>::iterator k;
	  for( j=elem.begin(), k=values.begin(); 
	       j!=elem.end() && k!=values.end(); j++, k++ ) {
	    if( int(numbers.size()) < ((*j)+1) )
	      numbers.resize( (*j)+1 );
	    numbers[(*j)] = (*k);
	  }
	}
      }
      else {
	_wrong.push_back( *i );
      }
      delete []text;
    }
  }
  
  return( status );
}

template bool Opt::_myGetOpt(string,string,vector<int>&);

//________________________________________
template <class T>
bool Opt::_myGetOptRec( string tag, string key, vector<T>& numbers ) {
  // read a vector of numbers (type T) following tag and key strings in a line
  bool status = false;
  numbers.clear();

  list<int> elem;
  list<T> values;

  list<string>::iterator i;
  for(i=_unused.begin(); i!=_unused.end(); i++ ) {
    istringstream is( (*i).c_str() );
    if( _checkTags( tag, key, is ) ) {
      char* text = new char[(*i).size()];
      is.getline( text, 512, '\0' );
      // get the list of vector elements
      status = _getElementsList( text, elem );
      // get the corresponding list of numbers
      istringstream is( text );
      status = _getNumbersList( is, values );
      _unused.remove( *i );
      if( status ) {
	// the number of elements must be equal to the number of numbers
	if( elem.size() != values.size() ) {
	  status = false;
	  _wrong.push_back( *i );
	  numbers.clear();
	}
	else { 	// store numbers in the right positions
	  list<int>::iterator j;
	  typename list<T>::iterator k;
	  for( j=elem.begin(), k=values.begin(); 
	       j!=elem.end() && k!=values.end(); j++, k++ ) {
	    if( int(numbers.size()) < ((*j)+1) ) {
	      numbers.resize( (*j)+1 );
	    }
	    numbers[(*j)] = (*k);
	  }
	}
      }
      else {
	_wrong.push_back( *i );
      }
      delete []text;
      return( status );
    }
  }
  
  return( status );
}


//________________________________________
bool Opt::_checkTags( string tag, string key, istringstream& is ) {

  // verify if line starts with <tag> + <key> strings 
  bool status = false;
  string s1, s2;
  istringstream it(tag.c_str());
  do {
    it >> s2;
    if( !it.fail() ) is >> s1;
  } while( s1 == s2 && !it.fail() );
  if( it.fail() ) {
    istringstream ik(key.c_str());
    do {
      ik >> s2;
      if( !ik.fail() ) is >> s1;
    } while( s1 == s2 && !ik.fail() );
    if( ik.fail() ) {
      status = true;
    }
  }
  return( status );
}

//________________________________________
template <class T> bool Opt::_getNumber( istringstream& is, T& num ) {

  // copy a number from <is> to <num> and check if all is ok
  bool status = true;
  is >> setbase( _intMode ) >> num;
  if( !is.fail() ) {
    string st; is >> st;
    if( st.size() != 0 ) {
      num = 0;
      status = false;
    }
  }
  else {
    num = 0;
    status = false;
  }
  return( status );
} 

//________________________________________
template <class T>
bool Opt::_getNumbersList( istringstream& is, list<T>& numbers ) {

  // copy a list of numbers from <is> to <numbers> and check if all is ok
  bool status;

  if( is.eof() ) {
    numbers.clear();
    status = false;
  }
  else {
    T prev = 0, act = 0;
    bool first = true;
    numbers.clear();
    status = true;
    while( !is.eof() ) {
      string st;
      is >> st;
      if( first ) {
	first = false;
	istringstream isi(st.c_str());
	act = 0; 
	isi >> setbase( _intMode ) >> act;
	if( !isi.fail() ) {
	  st.erase(); isi >> st;
	  if( st.size() == 0 ) {
	    numbers.push_back( act );
	    prev = act;
	  }
	  else {
	    numbers.clear();
	    status = false;
	  }
	}
	else {
	  numbers.clear();
	  status = false;
	}
      }
      else {
	// "N to M", "N - M", "N -> M" can be used to select 
	// numbers from N to M. Spaces are mandatory
	if( st == "to" || st == "-" || st == "->" ) {
	  // allow this only for integer numbers
	  if( typeid(numbers) == typeid(list<int>) ) {
	    is >> st;
	    istringstream isi(st.c_str());
	    act = 0; 
	    isi >> setbase( _intMode ) >> act;
	    st.erase(); isi >> st;
	    if( st.size() == 0 ) {
	      // N must be less than M 
	      if( prev < act ) {
		for( int i=int(prev)+1; i<=int(act); i++ ) {
		  numbers.push_back( i );
		  prev = 0;
		  first = true;
		}
	      }
	      else if( prev > act ) {
		numbers.clear();
		status = false;
	      }
	    }
	    else if( prev > act ) {
	      numbers.clear();
	      status = false;
	    }
	  }
	  else {
	    numbers.clear();
	    status = false;
	  }
	}
	// "N times M" or "N # M" can be used to select
	// N times the M number
	else if( st == "times" || st == "#" ) {
	  is >> st;
	  istringstream isi(st.c_str());
	  act = 0; 
	  isi >> setbase( _intMode ) >> act;
	  st.erase(); isi >> st;
	  // "N" must be integer and positive
	  if( st.size() == 0 && prev > 0 && floor((double) prev) == prev ) {
	    numbers.pop_back();
	    for( int i=0; i<prev; i++ ) {
	      numbers.push_back( act );
	    }
	    prev = 0;
	    first = true;
	  }
	  else {
	    numbers.clear();
	    status = false;
	  }
	}
	else {
	  // At this point next string must be a number
	  istringstream isi(st.c_str());
	  act = 0; 
	  isi >> setbase( _intMode ) >> act;
	  if( !isi.fail() ) {
	    st.erase(); isi >> st;
	    if( st.size() == 0 ) {
	      numbers.push_back( act );
	      prev = act;
	    }
	    else {
	      numbers.clear();
	      status = false;
	    }
	  }
	  else {
	    numbers.clear();
	    status = false;
	  }
	}
      }
    }
  }
  return( status );
}

//________________________________________
bool Opt::_getElementsList( char* text, list<int>& elem ) {

  //Extract numbers inside brackets. Separators can be "," or "-".
  //"N-M" means "from N to M".
  // e.g.: [1,3,4,6-9] = 1 3 4 6 7 8 9
  // spaces and tabs are ignored. Numbers are stored in a list.
  bool status = true;
  regex_t re;
  regmatch_t pm;
  int an, pn;

  char* linebuff = new char[strlen(text)+1];
  char* line = linebuff;

  strcpy( line, text );
  elem.clear();

  // remove initial spaces (in any )
  while( line[0] == ' ' || line[0] == '\t' ) line += 1;
  // At this point "[" must be at the beginning
  (void) regcomp( &re, "[[]", REG_EXTENDED );
  int i = regexec( &re, line, (size_t) 1, &pm, 0 ); 
  regfree( &re );
  if( i == 0 && pm.rm_so == 0 ) {
    line += 1;
    while( line[0] == ' ' || line[0] == '\t' ) line += 1;
    // At this point a number must be at the beginning
    (void) regcomp( &re, "[0-9]+", REG_EXTENDED );
    int i = regexec( &re, line, (size_t) 1, &pm, 0 ); 
    regfree( &re );
    if( i==0 && pm.rm_so==0 ) {
      pn = atoi( line );
      if( pn >= 0 ) {
	elem.push_back( pn ); 
      }
      else {
	status = false;
	elem.clear();
	line[0] = '\0';
      }
      line +=pm.rm_eo;
      while( line[0] == ' ' || line[0] == '\t' ) line += 1;
      while( strlen( line ) > 0 ) {
	// Look for a "]"
	(void) regcomp( &re, "[]]", REG_EXTENDED );
	int i = regexec( &re, line, (size_t) 1, &pm, 0 ); 
	regfree( &re );
	if( i == 0 && pm.rm_so == 0 ) {
	  line += 1;
	  while( line[0] == ' ' || line[0] == '\t' ) line += 1;
	  if( strlen( line ) != 0 ) {
	    strcpy( text, line );
	    line[0] = '\0';
	  }
	}
	else { 
	  // Look for a ","
	  (void) regcomp( &re, "[,]", REG_EXTENDED );
	  int i = regexec( &re, line, (size_t) 1, &pm, 0 ); 
	  regfree( &re );
	  if( i == 0 && pm.rm_so == 0 ) {
	    line += 1;
	    while( line[0] == ' ' || line[0] == '\t' ) line += 1;
	    // At this point a number must be at the beginning
	    (void) regcomp( &re, "[0-9]+", REG_EXTENDED );
	    int i = regexec( &re, line, (size_t) 1, &pm, 0 ); 
	    regfree( &re );
	    if( i==0 && pm.rm_so==0 ) {
	      pn = atoi( line );
	      if( pn >= 0 ) {
		elem.push_back( pn ); 
	      }
	      else {
		status = false;
		elem.clear();
		line[0] = '\0';
	      }
	      line +=pm.rm_eo;
	      while( line[0] == ' ' || line[0] == '\t' ) line += 1;
	    }
	    else {
	      status = false;
	      elem.clear();
	      line[0] = '\0';
	    }
	  }
	  else {
	    // Look for a "-"
	    (void) regcomp( &re, "[-]", REG_EXTENDED );
	    int i = regexec( &re, line, (size_t) 1, &pm, 0 ); 
	    regfree( &re );
	    if( i == 0 && pm.rm_so == 0 ) {
	      line += 1;
	      while( line[0] == ' ' || line[0] == '\t' ) line += 1;
	      // At this point a number must be at the beginning
	      (void) regcomp( &re, "[0-9]+", REG_EXTENDED );
	      int i = regexec( &re, line, (size_t) 1, &pm, 0 ); 
	      regfree( &re );
	      if( i==0 && pm.rm_so==0 ) {
		an = atoi( line );
		if( pn < an ) {
		  for( int i=pn+1; i<=an; i ++ ) {
		    elem.push_back( i );
		  }
		}
		else {
		  status = false;
		  elem.clear();
		  line[0] = '\0';
		}
		line +=pm.rm_eo;
		while( line[0] == ' ' || line[0] == '\t' ) line += 1;
	      }
	      else {
		status = false;
		elem.clear();
		line[0] = '\0';
	      }
	    }
	    else {
	      status = false;
	      elem.clear();
	      line[0] = '\0';
	    }
	  }
	}
      }
    }
  }
  delete []linebuff;
  return( status );
}

//________________________________________
void Opt::setIntMode( string mode ) {
  if( mode == "dec" ) {
    _intMode = 10;
  }
  else if( mode == "hex" ) {
    _intMode = 16;
  }
  else if( mode == "oct" ) {
    _intMode = 8;
  }
  else { // decimal in case of wrong settigs...
    _intMode = 10;
  }
}

//________________________________________
bool Opt::expand( string& filename ) {

  bool status = false;
  wordexp_t* exp = (wordexp_t*) malloc(512);
  if( wordexp( filename.c_str(), exp, 0 ) == 0 ) {
    char* name = exp->we_wordv[0];
    if( access( name, R_OK ) == 0 ) { 
      status = true;
      filename = name;
    }
  }
  wordfree( exp );
  return( status );
}
