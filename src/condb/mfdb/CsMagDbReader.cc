#include <stdlib.h>
#include "CsMagDbReader.h"
#include "CsPersFieldSol.h"
#include "CsPersFieldSM1.h"
#include "CsPersFieldSM2.h"


CsMagDbReader::CsMagDbReader(const string& dbName) :
	CsCondDbReader(dbName) {
	this->init();
}

CsMagDbReader::~CsMagDbReader() {
	CsErrLog::Instance()->mes(elDebugging, 
		     "CsMagDbReader is going to die.");
}

CsMagFieldSol CsMagDbReader::getFieldSol(string type){
  CsMagFieldSol infoField;
  string containerName = "SOL_"+type;
  string message = "Read solenoid field map from database";
  if(type==string("CALC")) message+=string(" : caluculated data ");
  this->open();
  HepTime theBeginTime = HepTime(2000,1,1,0,0,0);
  HepRef(calibInterval) theInterval;
  int res = getDatabase()->findInterval(theInterval,
				containerName.c_str(),theBeginTime);
  if(res){
    HepRef(CsPersFieldSol) thePersField 
      = (HepRef(CsPersFieldSol)) theInterval->getObject();

   int NoOfGrid = thePersField->getNoGrid();
   CsFieldGridSol grid;
   for(int i=0;i<NoOfGrid;i++){
     grid = thePersField->getGrid(i);
     infoField.addGrid(grid);
   }
   CsErrLog::Instance()->mes(elInfo,message);
  }else{
    CsErrLog::Instance()->mes(elError,
		 "Solenoid field map dose not exist in database.");
  }
  this->close();
  return infoField;
}

CsMagFieldSM1 CsMagDbReader::getFieldSM1(string type,int poleDistance){
  CsOpt *opt = CsOpt::Instance();
  CsMagFieldSM1 infoField;
  char charPD[10];
  sprintf(charPD,"%d",poleDistance);
  string poleDistance_(charPD);
 
  this->open(); 
  string containerName = string("SM1")+string("_")
    +type+string("_")+poleDistance_;
  string message = "Read SM1 field map from database";
  if(type==string("CALC")) message+=string(" : caluculated data ");
  if(poleDistance== 82  || 
     poleDistance== 132 || 
     poleDistance== 172 ) message+=string(" : pole distance = ")+poleDistance_;

  HepTime theBeginTime = HepTime(2000,1,1,0,0,0);
  HepRef(calibInterval) theInterval;
  int res = getDatabase()->findInterval(theInterval,
				containerName.c_str(),theBeginTime);
  if(res){
    HepRef(CsPersFieldSM1) thePersField 
      = (HepRef(CsPersFieldSM1)) theInterval->getObject();

   int NoOfGrid = thePersField->getNoGrid();
   CsFieldGridSM1 grid;
   for(int i=0;i<NoOfGrid;i++){
     grid = thePersField->getGrid(i);
     infoField.addGrid(grid);
   }
   CsErrLog::Instance()->mes(elInfo,message);
  }else{
    CsErrLog::Instance()->mes(elError,
		      "SM1 field map dose not exist in database.");
  }

  this->close();
  return infoField;
}

CsMagFieldSM2 CsMagDbReader::getFieldSM2(int current){
  CsOpt *opt = CsOpt::Instance();
  char charC[10];
  sprintf(charC,"%d",current);
  string current_(charC);
  this->open();
//   if( !opt->getOpt( "MagDb", "SM2_current",current) ){
//     CsErrLog::Instance()->mes( elError, 
//      "Tag `MagDb` and key `SM2_current` were not found. ");
//   }

  string containerName = string("SM2_")+current_;
  string message = "Read SM2 field map from database";
  if(current == 2000  || 
     current == 4000  || 
     current == 5000 ) message+=string(" : current = ")+current_;
  else{
    CsErrLog::Instance()->mes( elError,"SM2_current is not correct");
  }

  CsMagFieldSM2 infoField; 
  HepTime theBeginTime = HepTime(2000,1,1,0,0,0);
  HepRef(calibInterval) theInterval;
  int res = getDatabase()->findInterval(theInterval,
				containerName.c_str(),theBeginTime);
  if(res){
    HepRef(CsPersFieldSM2) thePersField 
      = (HepRef(CsPersFieldSM2)) theInterval->getObject();
    infoField = thePersField->getField();
    CsErrLog::Instance()->mes(elInfo,message);
  }else{
    CsErrLog::Instance()->mes(elError,
		"SM2 field map dose not exist in database.");
  }
  this->close();
  return infoField;
}






