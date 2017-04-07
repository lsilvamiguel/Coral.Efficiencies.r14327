#include "CsScifiDbUpdater.h"
#include "CsPersGeneralVector.h"
CsScifiDbUpdater::CsScifiDbUpdater(const string &dbName) :
  CsCondDbUpdater(dbName){}
CsScifiDbUpdater::~CsScifiDbUpdater(){}
bool CsScifiDbUpdater::store(CsCalibScifi &scifi,
			 CsTime &startTime,CsTime &endTime){

  string test = "$Revision: 1.3 $"; 

  //get container list by tagging remark data
  list<CsCondDbConf::Container>
    containerList = configuration_.findByRemark("SFD");
  if(containerList.size() == 0) {
    CsErrLog::Instance()->mes( elInfo,
	 	       "No container defined in configuration file." );
    ostrstream out;
    out     << configuration_ << endl;
    CsErrLog::Instance()->mes(elDebugging, out.str() );
    return false;
  }
 
  CsCalibScifiPlane *plane;
  double *dataPointer;
  CsCondDbConf::iContainer iCon;
  vector<int> meanList,sigmaList;
  for(iCon  = containerList.begin();iCon != containerList.end();iCon++){
     plane = scifi.getPlane((*iCon).dataSource);
     cout << (*iCon).name << endl;
 
    //
    //Convert transient object to persistent object
    int length = 0;
    meanList  = plane->getMeanList();
    sigmaList = plane->getSigmaList();
    length+=meanList.size();
    length+=sigmaList.size();

    dataPointer = new double[length];

    int i,j;
    for(i=0;i<meanList.size();i++){
      dataPointer[i] = meanList[i];
    }
    for(j=0;j<meanList.size();j++){
      dataPointer[j+i] = sigmaList[j];
    }

    open();

    HepRef(CsPersGeneralVector)  theContainer 
      = new(getDatabase()->hint()) CsPersGeneralVector(length,dataPointer);

    HepTime sTime = convertTime( startTime );
    HepTime eTime = convertTime( endTime );

    getDatabase()->store(theContainer, (*iCon).name.c_str(), sTime, eTime);
    delete[] dataPointer;
  }

  return true;
}






