#ifndef CsBeamReconstructionSupport_h 
#define CsBeamReconstructionSupport_h 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdio.h>

#include "TMath.h"

#include "CsTypes.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsTrack.h"
#include "CsBeamBackPropagation.h"

#include <CLHEP/Matrix/Matrix.h>

class CsBeamBackPropagation;

//Function dedicated to reading of a one key from the configuration file 
void read_configuration_key(std::string tag, std::string key, double& value);
void read_configuration_key(std::string tag, std::string key, int& value);
void read_configuration_key(std::string tag, std::string key, bool& value);
void read_configuration_key(std::string tag, std::string key, std::string& value);

//BMS cluster
class CsBMScluster{

   private:
      bool _is_empty;                        //if empty cluster
      double _t, _dt;                        //time
      double _y, _dy;                        //position along y
      double _z;                             //position along z
      CsCluster* _cs_cluster;                //CsCluster object created in case when BMS clusters are stored in mDSTs

   public:
      inline bool isEmpty(){   return _is_empty;  }
      inline double get_t(){   return _t;   }
      inline double get_dt(){   return _dt;   }
      inline double get_y(){   return _y;   }
      inline double get_dy(){   return _dy;   }
      inline double get_z(){   return _z;   }

      void set_cs_cluster(CsDetector* det, std::vector<CsDigit*> digits);   //set CsCluster object
      inline CsCluster* get_cs_cluster(){ return _cs_cluster; }             //get CsCluster object

      CsBMScluster();
      ~CsBMScluster();
      CsBMScluster(double t, double dt, double y, double dy, double z);
};

//BMS track
class CsBMStrack{

   private:
      bool _is_wrong;                        //if is wrong

      double _t, _dt;                        //time
      double _tchi2;                         //time chi2

      double _mom, _dmom;                    //momentum
      double _schi2;                         //spatial chi2

      double _LH;                            //total LH

      int _nPlanes;                          //number of involved planes
      std::vector<CsBMScluster*> _clusters;  //vector with hits/clusters

   public:
      inline void set_is_wrong(){ _is_wrong = true; }
      inline bool is_wrong(){  return _is_wrong; } 

      void calculate_time(double BT_time, double* BTvsBMSResolution, double* BTvsBMSCorrection);
      inline double get_t(){   return _t;  }
      inline double get_dt(){   return _dt;  }
      inline double get_tchi2(){  return _tchi2; }

      inline void set_mom(double mom){ _mom = mom;  }
      inline void set_dmom(double dmom){ _dmom = dmom;  }
      inline double get_mom(){    return _mom;   }
      inline double get_dmom(){    return _dmom;   }

      inline void set_schi2(double schi2){ _schi2 = schi2;  }
      inline double get_schi2(){  return _schi2; }

      inline void set_LH(double LH){ _LH = LH;  }
      inline double get_LH(){  return _LH; }

      inline int get_nPlanes(){   return _nPlanes;  }
      inline const std::vector<CsBMScluster*> get_clusters(){  return _clusters;    }

      void show_clusters();
      int get_cluster_mask();

      CsBMStrack(int nPlanes, CsBMScluster* hit1, CsBMScluster* hit2, CsBMScluster* hit3, 
                              CsBMScluster* hit4, CsBMScluster* hit5, CsBMScluster* hit6);
      ~CsBMStrack();
};

//BMS slice (in time and space)
class CsBMSslice{

   private:
      double _y_min[6];
      double _y_max[6];

      std::vector<CsBMScluster*> _clustersBM01; 
      std::vector<CsBMScluster*> _clustersBM02;
      std::vector<CsBMScluster*> _clustersBM03;
      std::vector<CsBMScluster*> _clustersBM04;
      std::vector<CsBMScluster*> _clustersBM05;
      std::vector<CsBMScluster*> _clustersBM06;

   public:
      inline double get_y_min(int i_plane){  return _y_min[i_plane]; }
      inline double get_y_max(int i_plane){  return _y_max[i_plane]; }

      void add_cluster(int ipl, CsBMScluster* cluster);
      const std::vector<CsBMScluster*> get_clusters(int ipl);

      CsBMSslice(double* y_min, double* y_max);
      ~CsBMSslice();
};

//Beam combination
class CsBeamCombination{
   
   private:
      const CsTrack* _BTtrack;                     //BT track
      CsHelix* _lastBThelix;  

      std::vector<CsBMStrack*> _BMStracks;         //reconstructed BMS tracks 

      std::vector<CsBMSslice*> _BMSslices;         //BMS slices 

      double _time_min;                            //time range
      double _time_max;

      double _slice_mom;                           //momentum interval
      int _slice_n;                                //n slices
      double _slice_overlap;                       //overlap

      CsBeamBackPropagation* BeamBackPropagation;

   public:
      inline const CsTrack* get_BTtrack(){ return _BTtrack; }
      CsBMStrack* get_best_BMStrack(double* best_LH);
      inline CsHelix* get_lastBThelix(){ return _lastBThelix; }

      inline const std::vector<CsBMStrack*> get_BMStracks(){ return _BMStracks; }
      inline const std::vector<CsBMSslice*> get_BMSslices(){ return _BMSslices; }

      bool check_time(double time);
      void add_cluster(int i_plane, CsBMScluster* cluster);
      void add_bms_track(CsBMStrack* bms_track){ _BMStracks.push_back(bms_track); }

      CsBeamCombination(const CsTrack* BTtrack, double time_gate, double int_mom, int n_mom, double slice_overlap, CsBMScluster* empty_cluster);
      ~CsBeamCombination();
};



//container for back propagation coeff
class CsBeamBackPropagationCoeff{
   
   private:
      std::string _plane_name;
      double _correctionFactor;
      double _spaceResolutionBMS; 

   public:
      std::string get_plane_name(){ return _plane_name;  }
      double  get_correctionFactor(){ return _correctionFactor;  }
      double  get_spaceResolutionBMS(){  return _spaceResolutionBMS;   }

      CsBeamBackPropagationCoeff(std::string name, double correction, double resolution){
         _plane_name = name;   _correctionFactor = correction; _spaceResolutionBMS = resolution;
      }
};

#endif
