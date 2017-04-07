//-*-Mode: C++;-*-
#ifndef _FfrReader_h_
#define _FfrReader_h_
/*!
  \file		FfrReader.h
  \brief	COMGEANT FFR key card reading object definition file.
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:44 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include "CsEndOfJob.h"

//---- FfrReader -----------------------------------------------------------
/*!
	\class	FfrReader
	\brief	FFR card reading class
*/
class FfrReader : public CsEndOfJob {
public:
  static FfrReader* Instance();	//!< Instansiation.
  bool end();	//!< called by CsRegistrySing::callEndMethod()

private:
  FfrReader();	//<! default constructor


public:
  ~FfrReader();	//<! default destructor

  bool read();	//!< read ffr card specified in coral.option file.

//! dump information to ostream
  friend ostream& operator<<(ostream& os, FfrReader& ffr);

private:
  bool getFileName(); //!< get file location from coral option file.
  bool read( istream& is ); //!< read information from stream.
  string checkLine(istream& is);	//!< check one line 

private:

  static FfrReader* singleton_;	//!< singleton pointer
  string ffrDirectory_;	//!< FFR file directory.
  string general_;		//!< General FFR file
  string muon_;			//!< for Muon program
  string hadron_;			//!< for hadron program

  map<string, string> dataLine_;	//! data storage.

  string gVersion_;	//!< geometry version tag

public:
  string ffrDirectory() const {return ffrDirectory_;}	//!< FFR file directory.
  string general() const {return general_;}		//!< General FFR file
  string muon() const {return muon_;}			//!< for Muon program
  string hadron() const {return hadron_;}			//!< for hadron program

  void ffrDirectory(const string& inputString) {ffrDirectory_ = inputString;}	//!< FFR file directory.
  void general(const string& inputString) {general_ = inputString;}		//!< General FFR file
  void muon(const string& inputString) {muon_ = inputString;}			//!< for Muon program
  void hadron(const string& inputString) {hadron_ = inputString;}			//!< for hadron program

  map<string, string>& dataLine() {return dataLine_;}

  string	gVersion(); //!< return Geometry Versioning TAG

  void dump(ostream& os); //!< dump all data line to ostream.

  void strip(string& data, const char& ch);	//!< remove unused char from string
  int stoi(const string& str);	//!< string to int 
  double stod(const string& str);	//!< string to double
  double multi(const string& str, int& n, double& val); //!< analyze *

  static const char nDelim;		//<! delimiter character for string in ffr cards
  static const char space;		//<! space character in ffr cards

	
};


#endif
