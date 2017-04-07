#include <string>
#include <vector>
#include <map>

#include "Reco/DataBase.h"
using namespace Reco;
#include "CsGEMDetector.h"

//
//Converter for calibration file of GEM from old format to new format 
//
// Usage: gem_calib_convert <plane name> <Unuveral Time>
// "gem_calib_convert GM01X1__ 1000072800"
// You can use this with "gem_calib_convert.pl".
//

int main(int argc,char *argv[]){

  string folder = argv[1];
  time_t ltime = atoi(argv[2]);
  string DB = argv[3];
  tm *t = localtime(&ltime);

  string tmp = folder;
  while(tmp.find("/")<tmp.size()){
    tmp.replace(tmp.find("/"),1," ");
  }
  istrstream is(tmp.c_str());
  string plane;
  for(int i=0;i<3;i++)is >> plane;

  DataBase *d = new DataBase(DB.c_str());

  string path = "/COMPASS/";
  for(int i=0;i<2;i++) path+=plane[i];
  path+="/";
  for(int i=0;i<8;i++) path+=plane[i];
  path+="/calibration";

  map<int,vector<float> > data;

  int order = 1;
  int chiporder = 1;
  if(plane == "GM01X1__" ||
     plane == "GM02U1__" ||
     plane == "GM02V1__" ||
     plane == "GM04U1__" ||
     plane == "GM04V1__" ||
     plane == "GM05U1__" ||
     plane == "GM05V1__" ||
     plane == "GM06U1__" ||
     plane == "GM06V1__" ||
     plane == "GM08U1__" ||
     plane == "GM08V1__" ||
     plane == "GM09U1__" ||
     plane == "GM09V1__"
     ){
    order = -1;
  }
  if(plane == "GM01V1__" ||
     plane == "GM02V1__" ||
     plane == "GM02X1__" ||
     plane == "GM04V1__" ||
     plane == "GM04X1__" ||
     plane == "GM05V1__" ||
     plane == "GM05X1__" ||
     plane == "GM06V1__" ||
     plane == "GM06X1__" ||
     plane == "GM08V1__" ||
     plane == "GM08X1__" ||
     plane == "GM09V1__" ||
     plane == "GM09X1__"
     ){
    chiporder = -1;
  }
  
  vector<CsGEMDetector::Calib> calib;
  d->Read(plane,calib,*t);

  int ch = 0;
  for(unsigned i=0;i<calib.size();i++){
    int chChip = 0;
    for(unsigned j=0;j<calib[i].channels.size();j++){
      vector<float> calib1V;
      calib1V.push_back(calib[i].channels[j].flag);
      calib1V.push_back(calib[i].channels[j].pedestal_mean);
      calib1V.push_back(calib[i].channels[j].pedestal_sigma);
      calib1V.push_back(calib[i].channels[j].calibration_mean);
      calib1V.push_back(calib[i].channels[j].calibration_sigma);
      
      if(order == 1 && chiporder == 1){
	ch = 128*i + chChip;
      }
      if(order == 1 && chiporder == -1){
	ch = (128*5) - 128*i+chChip;
      }
      if(order == -1 && chiporder == 1){
	ch = 128*i + 127-chChip;
	  }
      if(order == -1 && chiporder == -1){
	    ch = (128*5) - 128*i + 127-chChip;    
      }
      
      data.insert(pair<int,vector<float> >(ch,calib1V));
      chChip++;
    }
  }

  map<int,vector<float> >::iterator it1;
  for(it1=data.begin();it1!=data.end();it1++){
    cout << (*it1).first << " ";
    cout << (*it1).second[0] << " ";
    cout << (*it1).second[1] << " ";
    cout << (*it1).second[2] << endl;
  }

}


