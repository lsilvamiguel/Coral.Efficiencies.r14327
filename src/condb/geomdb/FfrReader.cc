#if __GNUG__ >= 2
#  pragma implementation
#endif
/*!
  \file		FfrReader.h
  \brief	COMGEANT FFR key card reading object implementation file.
  \author	$Author: miyachi $
  \version	$Revisopn$
  \date		$Date: 2000/06/07 08:03:17 $
*/

#include "FfrReader.h"
#include "CsRegistry.h"
#include "CsOpt.h"
#include "CsErrLog.h"

# ifdef COMPASS_USE_OSPACE_STD
#  include <ospace/std/algorithm>
# else
#  include <algorithm>
# endif // COMPASS_USE_OSPACE_STD

const char FfrReader::nDelim = '\'';
const char FfrReader::space =	' ';

FfrReader* FfrReader::singleton_ = NULL; // static pointer initialization.

FfrReader* FfrReader::Instance(){
	if( singleton_ == NULL ) {
		singleton_ = new FfrReader();

		CsRegistry csReg_;
		if( csReg_.EOJRegistration((CsEndOfJob*) singleton_) ){
			CsErrLog::Instance()->mes(elDebugging, 
				"FFrReader has been registerd successfully.");;
		} else {
			CsErrLog::Instance()->mes(elError, 
				"FFrReader has not been registerd successfully.");;
		}
	}
	return singleton_;
}

bool FfrReader::end(){
	return true; // if FfrReader has to do someting at the end of its life...
}

FfrReader::FfrReader() : 
	ffrDirectory_(""), 
	general_(""), muon_(""), hadron_(""),
	dataLine_(), gVersion_("") {
	if( this->read() ){
		CsErrLog::Instance()->mes(elDebugging,
			"FFR reader succeeded to read ffr cards.");
	} else {
	}
}

FfrReader::~FfrReader(){
	CsErrLog::Instance()->mes(elDebugging, "FfrReader is going to die.");	
}

bool FfrReader::read(){
	if( ! this->getFileName() ) return false;

	string generalFfr = 
		this->ffrDirectory() + '/' + this->general();

	ifstream gFfr( generalFfr.c_str() );
	if(gFfr.fail()){
		CsErrLog::Instance()->mes(elError,
			"General ffr card(" + generalFfr + ") could not be opend.");
		return false;
	} else {

		if( ! this->read(gFfr) ) {
			CsErrLog::Instance()->mes(elError,
				"General ffr card(" + generalFfr + ") " +
				"reading error occured.");
			return false;
		}
		
	}
	
	string programFfr = this->ffrDirectory() + '/';
	if(this->muon() != ""){
		programFfr += this->muon();
	} else {
		programFfr += this->hadron();
	}

	ifstream pFfr( programFfr.c_str() );
	if(pFfr.fail()){
		CsErrLog::Instance()->mes(elError,
			"Program ffr card(" + programFfr + ") could not be opend.");
		return false;
	} else {
		
		if( ! this->read(pFfr) ) {
			CsErrLog::Instance()->mes(elError,
				"Program ffr card(" + programFfr + ") " +
				"reading error occured.");
			return false;
		}

	}

	return true;
}


bool FfrReader::read(istream& is){
	do {
		string line( this->checkLine( is ) );
    
	    istrstream stringStream(line.c_str());
	    string opt;
	    stringStream >> opt;
    
	    if( opt == "comment") {
			continue;
	    } else {
			char rest[2048];
			char* pRest = rest;
			
//			stringStream.getLine(rest, 2048, '\n');

			while( stringStream.get( *pRest ) ){
				pRest++;
			}
			*pRest = '\0';

			this->dataLine()[opt] = rest ;
//			this->dataLine()[opt] =  stringStream.str() ;

		}
		
	} while( ! is.eof() );
	
	return true;
}

string FfrReader::checkLine(istream& is){
	const int lineSize = 2048;
	char line[lineSize];
  
	is.getline( line, lineSize, '\n' );
  
	if(	is.eof()		||
		line[0] == '\0'	||
		line[0] == '#'	||
		line[0] == 'C'	||
		line[0] == 'c'	||
		line[0] == ' '
	) return "comment";

	return line;
}


bool FfrReader::getFileName(){
	string directory, general, muon, hadron;

	if( CsOpt::Instance()->getOpt("Geom", "ffr_directory",	directory) ){
		this->ffrDirectory(directory);
	}
	
	if( CsOpt::Instance()->getOpt("Geom", "ffr_general",		general) ){
		this->general(general);
	}
	
	if( CsOpt::Instance()->getOpt("Geom", "ffr_muon",		muon) ){
		this->muon(muon);
	}
	
	if( CsOpt::Instance()->getOpt("Geom", "ffr_hadron",		hadron) ){
		this->hadron(hadron);
	}


	if( this->ffrDirectory() == "" ) {
		CsErrLog::Instance()->mes(elError, 
			"FFR directory has not been specified.");
		return false;
	}


	if( this->general() == "" ){
		CsErrLog::Instance()->mes(elError, 
			"General ffr card has not been specified.");
		return false;		
	}


	if( this->muon() == "" && this->hadron() == "" ){
		CsErrLog::Instance()->mes(elError, 
			"Program specified ffr card has not been specified.");
		return false;		
		
	} else if( this->muon() != "" && this->hadron() != "" ) { 
		CsErrLog::Instance()->mes(elError, 
			"Muon and Hadron ffr cards both have been specified.");
		return false;		
	}
 

	ostrstream out;
	
	out	<<	"FFR Directory:\t"	<< this->ffrDirectory()	<< endl;
	out	<<	"    General:\t"	<< this->general()		<< endl;
	
	if(this->muon() != "")
	out	<<	"    Muon:\t"		<< this->muon()			<< endl;
	if(this->hadron() != "" )
	out	<<	"    hadron:\t"		<< this->hadron()		<< endl;

	CsErrLog::Instance()->mes(elDebugging, out.str() );


	return true;	
}


string	FfrReader::gVersion(){
	if( gVersion_ == "" ) {
		const string tag("GVERSION");
		gVersion_ = this->dataLine()[tag];

		this->strip( gVersion_, FfrReader::space);
		this->strip( gVersion_, FfrReader::nDelim);
		
	}
	return gVersion_;
}

void FfrReader::dump(ostream& os){
	os	<<	"C\tFFR Directory:\t"	<< this->ffrDirectory()	<< endl;
	os	<<	"C\t    General:\t"		<< this->general()		<< endl;
	os	<<	"C\t    Muon:\t"		<< this->muon()			<< endl;
	os	<<	"C\t    hadron:\t"		<< this->hadron()			<< endl;
	os	<<	"C\nC\nC\n"	<< endl;
	os	<<	"C------------------------------------------------------" << endl;
	os	<<	"C [TAG] [Data]" << endl ;
	map<string, string>::iterator iData;
	for(	iData =  this->dataLine().begin();
		iData != this->dataLine().end();
		iData++){
		os	<< (*iData).first	<< "\t"	<< (*iData).second << endl;
	}
}

void FfrReader::strip(string& data, const char& ch){
	string::iterator iStr = data.begin();
	while( iStr != data.end() ){	
		if( (*iStr) == ch ) {
			iStr = data.erase( iStr );
		} else {
			iStr++;
		} 
	}	
}

int FfrReader::stoi(const string& str){
	istrstream iStStr( str.c_str() );
	int rval; iStStr >> rval;
	return rval;
}

double FfrReader::stod(const string& str){
	istrstream iStStr( str.c_str() );
	double rval; iStStr >> rval;
	return rval;
}


double FfrReader::multi(const string& str, int& n, double& val){

	istrstream is( str.c_str() );
	char charBuff[16];
	is.getline( charBuff, 16, '*' );

	n = this->stoi( charBuff );

	string strBkg; is >> strBkg ;				 		
	val = this->stod( strBkg );
	
	return n * val;
	
}

ostream& operator<<(ostream& os, FfrReader& ffr){
	os	<<	"CARDS:\t" 
		<< ffr.ffrDirectory() << '/' << 	ffr.general() << endl;
	os	<< ffr.ffrDirectory() << '/' << 	ffr.muon() << endl;
	os	<< ffr.ffrDirectory() << '/' << 	ffr.hadron() << endl;
	os	<<	"GVERSION:\t	("	<< ffr.gVersion()	<< ")" << endl;
	
	return os;
}

