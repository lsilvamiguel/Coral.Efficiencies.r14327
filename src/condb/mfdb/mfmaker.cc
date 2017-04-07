#include "CsOpt.h"
#include "CsRegistrySing.h"
#include "CsSTD.h"
#include "CsMagDbUpdater.h"

main(int argc, char *argv[]){

  CsOpt* coralOpt = CsOpt::Instance(argc,argv);
  CsRegistrySing *reg_ = CsRegistrySing::Instance();

  const string path = "/afs/cern.ch/exp/compass/mc/comgeant/data/mag_fields";
  const string pathAlx = "/afs/cern.ch/user/k/korzenev/public/mag_fields";
  string fname;
  string containerName;
  CsMagDbUpdater *mag = new CsMagDbUpdater();

  fname = path + "/SM1h/SM1H.map.01.data";
  mag->storeFieldCalcSM1(fname,string("SM1_CALC_82"));

  fname = path + "/SM1m/SM1M.map.01.data";
  mag->storeFieldCalcSM1(fname,string("SM1_CALC_172"));

  fname = pathAlx + "/SM1m/SM1M.map.172.data";
  mag->storeFieldCalcSM1(fname,string("SM1_MEAS_172"));

  fname = path + "/SM2/FSM.map.2000.data";
  mag->storeFieldSM2(fname,containerName = "SM2_2000");

  fname = path + "/SM2/FSM.map.4000.data";
  mag->storeFieldSM2(fname,containerName = "SM2_4000");

  fname = path + "/SM2/FSM.map.5000.data";
  mag->storeFieldSM2(fname,containerName = "SM2_5000");

  fname = path + "/SOL/OXFORD.map.01.data";
  mag->storeFieldCalcSol(fname,containerName = "SOL_CALC");

  mag->commit();

  reg_->callEndMethods();
}













