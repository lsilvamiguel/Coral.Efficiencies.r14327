#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "BmsContainerConf.h"
#include "CsErrLog.h"

//---- BmsContainerConf ---------------------------------------------------------

BmsContainerConf::BmsContainerConf() : 
	filename_(""), sourceList_() {
}


BmsContainerConf::BmsContainerConf(const string& filename) : 
	filename_(filename), sourceList_() {
}

BmsContainerConf::BmsContainerConf(const BmsContainerConf& conf) : 
	filename_(conf.filename()), 
	sourceList_(conf.sourceList()) {
}

BmsContainerConf::~BmsContainerConf()
{
}




bool BmsContainerConf::read(){
	if( this->filename() == "" ){
		return false;
	}

	ifstream f( this->filename().c_str(), ios::in );
	if( !f ) {
		CsErrLog::Instance()->mes(elError, 
			this->filename() + " is not found.");
		return( false );
	} else {
		CsErrLog::Instance()->mes(elDebugging, 
			this->filename() + " is successfully opend.");
	}

	return read(f);
}

bool BmsContainerConf::read(istream& is){
	do {
		string line( this->checkLine( is ) );
    
	    istrstream stringStream(line.c_str());
	    string opt;
	    stringStream >> opt;
    
	    if( opt == "comment") {
	      continue;
	    }
	    else 
	    if( opt == "data") {
			Source source;
			stringStream >> source.datafile;

			string strSTime, srtETime;
			stringStream >> strSTime >> srtETime;
			
			
			source.startTime	= this->convertToTime(strSTime);
			source.endTime		= this->convertToTime(srtETime);
			

			sourceList_.push_back(source);

		}
	    else {
			CsErrLog::Instance()->mes(elDebugging, 
				"Unknown comand:\t[" + opt + "] in line:" + line);
		}

	} while( ! is.eof() );

  	return true;
}

string BmsContainerConf::checkLine(istream& is){
  const int lineSize = 2048;
  char line[lineSize];
  
  is.getline( line, lineSize, '\n' );
  
  if(	is.eof()		||
	line[0] == '\0'	||
	line[0] == '#'	||
	line[0] == ' '
	) return "comment";
  
  return line;
}

CsTime BmsContainerConf::convertToTime(const string& strTime){
	istrstream stream(strTime.c_str());
	const int clength(8);

	char 	cYear[clength], cMonth[clength], cDate[clength], 
			cHour[clength], cMinut[clength], cSec[clength], 
			cMicroSec[clength], cNano[clength];

	int	iYear(0), iMonth(0), iDate(0), iHour(0),
		iMinut(0), iSec(0), iMicroSec(0), iNano(0);


	if( stream.getline(cYear, clength, '/') ){
		istrstream buffer( cYear );
		buffer >> iYear;		
	} else {
	}
	
	if( stream.getline(cMonth, clength, '/') ){
		istrstream buffer( cMonth );
		buffer >> iMonth;		
	} else {
	}

	if( stream.getline(cDate, clength, '/') ){
		istrstream buffer( cDate );
		buffer >> iDate;		
	} else {
	}

	if( stream.getline(cHour, clength, '/') ){
		istrstream buffer( cHour );
		buffer >> iHour;		
	} else {
	}

	if( stream.getline(cMinut, clength, '/') ){
		istrstream buffer( cMinut );
		buffer >> iMinut;		
	} else {
	}
	
	if( stream.getline(cSec, clength, '/') ){
		istrstream buffer( cSec );
		buffer >> iSec;		
	} else {
	}

	if( stream.getline(cMicroSec, clength, '/') ){
		istrstream buffer( cMicroSec );
		buffer >> iMicroSec;
	} else {
	}

	if( stream.getline(cNano, clength, '/') ){
		istrstream buffer( cNano );
		buffer >> cNano;
	} else {
	}


//	cout << strTime;
//	cout 	<< "Year:\t"	<< cYear	<< '\t' 	<< iYear
//			<< "Month:\t"	<< cMonth	<< '\t' 	<< iMonth
//			<< "Day:\t"	<< cDate	<< '\t'	<< iDate
//			<< endl;



	return CsTime(	iYear, iMonth, iDate, 
					iHour, iMinut, iSec, iMicroSec);

}


ostream& operator<<(ostream& os, BmsContainerConf& conf){
	list<BmsContainerConf::Source> lSource = conf.sourceList();
	BmsContainerConf::iSource iS;

	os	<< "#------------------------------------------\n"
		<< "#\tBMS container contents configuration file.\n"
		<< endl;
			
	for(	iS = lSource.begin();
		iS != lSource.end();
		iS++){

		os	<< *(iS)	<< endl;
	}

	return os;
}
