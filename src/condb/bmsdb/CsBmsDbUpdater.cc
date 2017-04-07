#if __GNUG__ >= 2
#  pragma implementation
#endif

#include "CsSTD.h"
#include "CsOpt.h"
/*!
   \file    CsBmsDbUpdater.cc
   \brief   Compass BMS CDB updater class implementation file.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.4 $
   \date    $Date: 2000/06/06 09:48:35 $
*/
#include "CsErrLog.h"
#include "CsBmsDbUpdater.h"
#include "CsBMSpConstants.h"
#include "CsBMSconstants.h"

#include "BmsContainerConf.h"

CsBmsDbUpdater::CsBmsDbUpdater(const string& dbName) :
	CsCondDbUpdater(dbName) {
	this->init();
}


CsBmsDbUpdater::~CsBmsDbUpdater() {
	CsErrLog::Instance()->mes(elDebugging, "CsBmsDbUpdater is going to die.");
}

bool CsBmsDbUpdater::updateAll() {

// get container list by tagging remark data
  list<CsCondDbConf::Container> 
	containerList = configuration_.findByRemark("BMSconstants");

  if(containerList.size() == 0) {
	CsErrLog::Instance()->mes( elInfo, "No container defined in configuration file." );
	ostrstream out;
	out	<< configuration_ << endl;
	CsErrLog::Instance()->mes(elDebugging, out.str() );
	return false;
  }

// open DB
  this->open();

  CsCondDbConf::iContainer iCon;
  for( iCon  = containerList.begin(); iCon != containerList.end(); iCon++){

	BmsContainerConf contConf((*iCon).dataSource );

	if( contConf.read() ) {
	} else {
	  CsErrLog::Instance()->mes( elError, 
		"BMS container information can not be read for " + (*iCon).name );
	  return false;
	}


	list<BmsContainerConf::Source> lSource = contConf.sourceList();
	BmsContainerConf::iSource iS;
	int icheck(1), ncheck = lSource.size();
	
	for(	iS = lSource.begin(); iS != lSource.end(); iS++){

	  CsBMSconstants tConsts;
	  tConsts.readBMSconst( (*iS).datafile.c_str() );

// create persistent detector in memory
	  HepRef(CsBMSpConstants) refconsts =
		new( this->getDatabase()->hint() ) CsBMSpConstants( tConsts );


	  ostrstream out;
	  out << "Persistent object["	<< icheck++ 
		  << "] of total " << ncheck
		  << " has been created."	<< endl;
//			out	<< "Persistents objetcts: \n"
//				<< (*refconsts)		<< '\n'
//				<< " has been created."	<< endl;
	  CsErrLog::Instance()->mes(elDebugging, out.str() );


// convert time
	  HepTime sTime = this->convertTime( (*iS).startTime );
	  HepTime eTime = this->convertTime( (*iS).endTime );

// store object
	  this->getDatabase()->store( 
		refconsts, (*iCon).name.c_str(), sTime, eTime);

	  ostrstream out2;
	  out2	<< "Persistents objetcts has been created in "
			<< container << " of BMS DB (bms) with interval which is from "
			<< sTime << " to " << eTime << "." << endl;

	  CsErrLog::Instance()->mes(elDebugging, out2.str() );
	}
  }

// set status flag to HAVE
  status_=HAVE;

// close DB
  this->close();

  return true;
}

