//-*-Mode: C++;-*-
#ifndef _BmsContainerConf_h_
#define _BmsContainerConf_h_

#if __GNUG__ >= 2
#  pragma interface
#endif

#include "CsSTD.h"
#include "CsTime.h"

//---- BmsContainerConf -----------------------------------------------------------
/*!
	\class	BmsContainerConf
	\brief	Bms DB container configuration file.
*/
class BmsContainerConf {
public:
	BmsContainerConf();		//!< default constructor
	BmsContainerConf(const string& filename);	//!< constructor with filename
	BmsContainerConf(const BmsContainerConf& conf);	//!< copy constructor

    ~BmsContainerConf();	//!< default destructor

	/*!
		\struct	Source
		\brief	
	*/
	struct Source {
		string datafile;	//!<	source file.
		CsTime startTime;	//!<	varidate period start time
		CsTime endTime;		//!<	varidate period end time
		
		friend ostream& operator<<(
			ostream& os, const Source& source ){
			os	<< "date\t"
				<< source.datafile	<< '\t'
				<< source.startTime	<< '\t'
				<< source.endTime;
			return os;
		}	//!< output operator
		
	};

	typedef list<Source>::iterator iSource;	//!< type definition for iterator

private:
	string filename_;			//!< container configuration file name
	list<Source> sourceList_;	//!<	source list

public:
	inline string filename() const {return filename_;}	//!< return filename
	inline void filename(const string& filename) {
		filename_ = filename;
	} //!< set filename


	inline list<Source> sourceList() const {return sourceList_;} //!< return source list
	inline void sourceList(const list<Source>& sourceList) {
		sourceList_ = sourceList;
	} //!< set source list
	

//----- Published method to world 
	bool read();	//!< read configuration file.


	friend ostream& operator<<(ostream& os, BmsContainerConf& conf); //!< output operator

private:
	bool read(istream& is); //!< read from stream
	string checkLine(istream& is);	//!< check given line
	CsTime convertToTime(const string& strTime); //! convert string to CsTime

};

#endif
