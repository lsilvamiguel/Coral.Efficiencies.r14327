#include "CsMagDbUpdater.h"
#include "CsErrLog.h"
#include "CsPersFieldSol.h"
#include "CsPersFieldSM1.h"
#include "CsPersFieldSM2.h"
#include "CsMagFieldSM2.h"

CsMagDbUpdater::CsMagDbUpdater(const string& dbName) :
	CsCondDbUpdater(dbName) {
	this->init();
}


CsMagDbUpdater::~CsMagDbUpdater(){
	CsErrLog::Instance()->mes(elDebugging, 
				"CsMagDbUpdater is going to die.");
}

int CsMagDbUpdater::estimateGrid(string path){
  ifstream f(path.c_str());
  int lineSize = 200;
  char line[lineSize]; 
  int count = 0;
  while( f.getline( line, lineSize ) ) count++;
  f.close();
  return count;
}


bool CsMagDbUpdater::storeFieldCalcSol(string path,string containerName){
  this->open();

  HepRef(CsPersFieldSol) thePersField;
    thePersField = new (getDatabase()->hint()) CsPersFieldSol();
  HepTime theBeginTime = HepTime(2000,1,1,0,0,0);
  HepTime theEndTime   = HepTime(2100,1,1,0,0,0);

  string msg = "Read Calculated Field Maps : "+path;
  CsErrLog::Instance()->mes(elInfo,msg);
    
  ifstream f(path.c_str());
  int lineSize = 200;
  char line[lineSize]; 
  float z,r,Bz,Br,Bt,Brela,Homogeneity;
  int NoOfGrid = estimateGrid(path) - 11;
  thePersField->changeGridNumber(NoOfGrid);
  //Read text header
  for(int i=0;i<11;i++){
    f.getline( line, lineSize );
  }

  for(int i=0;i<NoOfGrid;i++){
    f.getline( line, lineSize );
    istrstream is(line);
    is>>z>>r>>Bz>>Br>>Bt>>Brela>>Homogeneity;
    CsFieldGridSol grid(z,r,Bz,Br);
    thePersField->addGrid(grid);
  }
  getDatabase()->store(thePersField,containerName.c_str(),
  	       theBeginTime,theEndTime);
  msg = "Store Calculated Field Maps to Database as  container "+containerName;
  CsErrLog::Instance()->mes(elInfo,msg);

  this->close();
}

bool CsMagDbUpdater::storeFieldCalcSM1(string path,string containerName){
  this->open();
  HepRef(CsPersFieldSM1) thePersField;
  thePersField = new (getDatabase()->hint()) CsPersFieldSM1();
  //CsPersFieldSM1 *thePersField = new CsPersFieldSM1();
  HepTime theBeginTime = HepTime(2000,1,1,0,0,0);
  HepTime theEndTime   = HepTime(2100,1,1,0,0,0);

  string msg = "Read Calculated Field Maps : "+path;
  CsErrLog::Instance()->mes(elInfo,msg);

  ifstream f(path.c_str());
  int lineSize = 200;
  char line[lineSize]; 
  float x,y,z,Bx,By,Bz;
  int NoOfGrid = estimateGrid(path) - 1;
  thePersField->changeGridNumber(NoOfGrid);
  //Read text header
  f.getline(line,lineSize);
  //

  for(int i=0;i<NoOfGrid;i++){
    f.getline( line, lineSize );
    istrstream is(line);
    is>>x>>y>>z>>Bx>>By>>Bz;
    CsFieldGridSM1 grid(x,y,z,Bx,By,Bz);
    thePersField->addGrid(grid);
    //if(i%1000 == 0) cout << "i = " << i << endl;
  }
  getDatabase()->store(thePersField,containerName.c_str(),
  		       theBeginTime,theEndTime);
  msg = "Store "+path+" to Database as  container "+containerName;
  CsErrLog::Instance()->mes(elInfo,msg);

  this->close();
}

bool CsMagDbUpdater::storeFieldSM2(string path,string containerName){
  this->open();

  HepRef(CsPersFieldSM2) thePersField;
  thePersField = new (getDatabase()->hint()) CsPersFieldSM2();
  HepTime theBeginTime = HepTime(2000,1,1,0,0,0);
  HepTime theEndTime   = HepTime(2100,1,1,0,0,0);

  CsMagFieldSM2 _infoField;

//In CsField

  string msg = "Read Calculated Field Maps : "+path;
  CsErrLog::Instance()->mes(elInfo,msg);

  ifstream f(path.c_str());
  int lineSize = 200,i;
  char line[lineSize]; 
  f.getline( line, lineSize );
  f>>_infoField.FSMA01;
  f>>_infoField.FSMA11;
 
  for( int i=0; i<1770; i++ ) f>>_infoField.FSMAX1[i];
  for( int i=0; i<472 ; i++ ) f>>_infoField.FSMBX1[i];
  for( int i=0; i<2360; i++ ) f>>_infoField.FSMAY1[i];
  for( int i=0; i<300 ; i++ ) f>>_infoField.FSMBY1[i];
  for( int i=0; i<3540; i++ ) f>>_infoField.FSMAZ1[i];
  for( int i=0; i<84  ; i++ ) f>>_infoField.FSMCX1[i];
  for( int i=0; i<126 ; i++ ) f>>_infoField.FSMCY1[i];
  for( int i=0; i<168 ; i++ ) f>>_infoField.FSMCZ1[i];
  for( int i=0; i<312 ; i++ ) f>>_infoField.FSMDX1[i];
  for( int i=0; i<390 ; i++ ) f>>_infoField.FSMDY1[i];
  for( int i=0; i<468 ; i++ ) f>>_infoField.FSMDZ1[i];

  thePersField->setField(_infoField);
  getDatabase()->store(thePersField,containerName.c_str(),
		       theBeginTime,theEndTime);
  msg = "Store "+path+" to Database as container "+containerName;
  CsErrLog::Instance()->mes(elInfo,msg);

  this->close();
}

//bool CsMagDbUpdater::update(const string& Container_name,  const CsTime& startTime,
//  const CsTime& endTime){
//
//return true;
//}

//bool CsMagDbUpdater::checkType(const CsCalibAbstract& TempConst){
//
//return true;
//}
