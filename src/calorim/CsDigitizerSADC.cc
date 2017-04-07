/*!
   \file    CsDigitizerSADC.cc
   \brief   Class to decode time and energy from SADC fit
   \version $Revision: 1.38 $
   \author  Vladimir Kolosov
   \author  Alexander Zvyagin
   \author  Denis Murashev
   \date    $Date: 2010/10/22 08:39:36 $
*/

#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

#include "coral_config.h"
#include <TMath.h>
#include "TGraphAsymmErrors.h"
#include "MyPulse.h"

#include "CsDigitizerSADC.h"
#include "CsCalorimeter.h"
#include "Reco/Exception.h"
#include "TVirtualFitter.h"

#include "CsOpt.h"
#include "CsEvent.h"

using namespace std;
using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

// static member variables

////////////////////////////////////////////////////////////////////////////////

const double                    CsDigitizerSADC::SADC_shape_factor = 0.01;

std::vector< double >           CsDigitizerSADC::SADC_shape[CsDigitizerSADC::N_SADC_shapes];
std::vector< double >           CsDigitizerSADC::D_SADC_shape[CsDigitizerSADC::N_SADC_shapes];
int                             CsDigitizerSADC::SADC_shape_max_position[CsDigitizerSADC::N_SADC_shapes];
int                             CsDigitizerSADC::SADC_shape_max_from_start[CsDigitizerSADC::N_SADC_shapes];
double                          CsDigitizerSADC::SADC_shape_max_from_front[CsDigitizerSADC::N_SADC_shapes];
int                             CsDigitizerSADC::SADC_shape_sadc_front_gate[N_SADC_shapes];
std::vector< double >           CsDigitizerSADC::SADC_shape_value_ower[CsDigitizerSADC::N_SADC_shapes];
std::vector< double >           CsDigitizerSADC::SADC_shape_index_ower[CsDigitizerSADC::N_SADC_shapes];
std::vector< double >           CsDigitizerSADC::SADC_tail[2];

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Constructor of Digitizer
  \param method {Method for digitalisation.}
  \param type {SADC or MSADC}
  \param table {}
*/

CsDigitizerSADC::CsDigitizerSADC ( CsCalorimeter* parent, Method method, DigitizerSADCBase::Type type,  Shape table ) :
                 DigitizerSADCBase ( parent,  type ),
                 method_(method), activated_(false),
                 shape_table_(table), fit_range_to_front_(6.), fit_range_to_back_(3.),
                 fit_ofile_(NULL), mypulse_(NULL), plane_(-1),
                 use_shape_filter_slc_(false), use_shape_filter_line_fit_(true)
{

  // Set overflow acording to SADC type: SADC 10Bit, MSADC 12Bit
  overflow_amplitude_ = type_ == CsDigitizerSADC::MSADC ? 4095 : 1023;
  // common for all methods
  sadc_delta_ = 0;
  sadc_ped_min_ = 0;
  sadc_ped_max_ = sadc_ped_min_+4;
  sadc_signal_min_ = 5;
  sadc_signal_max_ = sadc_signal_min_+20;
  sadc_front_ = 7;


  // Set Methods depending stuff here
  switch( method_ ) {
  case CsDigitizerSADC::MaxSimple:       // Method MaxSimple
    coeff_convert2max_ = 1.00;
    break;
  case CsDigitizerSADC::MaxAdvanced:       // Method MaxAdvanced
    coeff_convert2max_ = 1.00;
    break;
  case CsDigitizerSADC::FitPulse:       // Method PulseFit
    coeff_convert2max_ = 1.00;
    break;
  case CsDigitizerSADC::SummRange:       //  Method SummRange
    coeff_convert2max_ = 0.12;  // empiric constant
    break;
  case CsDigitizerSADC::ShapeFit:       //  Method ShapeFit
    coeff_convert2max_ = 0.12;  // empiric constant
    break;
  default:      //  Other Methods
    throw Reco::Exception("Unknown method '%i' in CsDigitizerSADC::CsDigitizerSADC", method_ );
    break;
  }

  // Init. static memebers
  if( SADC_shape[0].size() <= 0 ) {
    InitStatic();
  }
  option_store_stat_info_ = false;
  StatInfo zero = StatInfo(0.,0.,0.);
  stat_ped_odd_=zero;
  stat_ped_even_=zero;
  stat_ped_=zero;
  stat_dped_=zero;
  stat_time_=zero;

  ped_odd_old_ = 0.;
  ped_even_old_ = 0.;
  dped_old_ = 0.;
  normref_old_ = 1.;
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Initiate the hardcoded shapes
*/

void  CsDigitizerSADC::InitStatic ( void )
{
  bool debug = false;
  // Init. static memebers

  int all_shapes_sofar = N_SADC_shapes;

  if( SADC_shape[0].size() > 0 )
  {
    cerr <<" Wrong usage CsDigitizerSADC::InitStatic was called already!! Please Fix. " << endl;
    exit(1);
  }

  int jshape = 0;
  double shape_ECAL1[] = {   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
                             0.,    0.01,0.046,0.174, 0.526, 0.885,    1.,0.9345, 0.806,0.6585,
                             0.5353,0.4348,0.3556,0.2951,0.2482,0.2109,0.1809,0.1575,0.1395,0.1253,
                             0.1093,0.0919, 0.073, 0.061, 0.051, 0.045, 0.040, 0.035, 0.030, 0.025,
                             0.020, 0.015, 0.010, 0.008, 0.006, 0.004, 0.002, 0.001, 0.0007,0.0005};
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL1[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 6;
  SADC_shape_max_from_front[jshape] = 2.1;
  SADC_shape_sadc_front_gate[jshape] = 3;


  // Shape from MSADC ECAL2
  jshape = 1;
  double shape_ECAL2[] =     {0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,
                              0.,     0.0006, 0.0390, 0.2695, 0.6540, 0.9220, 1.0000, 0.9430, 0.8353, 0.7170,
                              0.6153, 0.5250, 0.4456, 0.3805, 0.3278, 0.2805, 0.2425, 0.2100, 0.1843, 0.1618,
                              0.1435, 0.1307, 0.1185, 0.1038, 0.0900, 0.0800, 0.0751, 0.0727, 0.0670, 0.0610,
                              0.0555, 0.0500, 0.0450, 0.0410, 0.0370, 0.0330, 0.0300, 0.0270, 0.0250, 0.0225 };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL2[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 6;
  SADC_shape_max_from_front[jshape] = 2.5;
  SADC_shape_sadc_front_gate[jshape] = 3;

  // RD Shape SADC ECAL1 GAMS ShapeNM100SADC
  jshape = 2;
  double shape_ECAL1_GAMS_100[] = {0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,
                                   0.,    0., 0.001, 0.054, 0.393, 0.851, 1.000, 0.926, 0.792, 0.658,
                                   0.548, 0.459, 0.388, 0.331, 0.286, 0.248, 0.218, 0.195, 0.178, 0.158,
                                   0.135, 0.112, 0.095, 0.079, 0.068, 0.058, 0.052, 0.045, 0.041, 0.037,
                                   0.033, 0.030, 0.027, 0.025, 0.022, 0.019, 0.017, 0.015, 0.013, 0.011 };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL1_GAMS_100[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 5;
  SADC_shape_max_from_front[jshape] = 2.95;
  SADC_shape_sadc_front_gate[jshape] = 3;

  // RD Shape SADC ECAL1 MT ShapeNM100SADC
  jshape = 3;
  double shape_ECAL1_MAINZ_100[] = {0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,
                                    0.000,  0.000, 0.001, 0.027, 0.414, 0.859, 1.000, 0.924, 0.779, 0.622,
                                    0.485,  0.370, 0.283, 0.221, 0.177, 0.146, 0.120, 0.099, 0.083, 0.074,
                                    0.066,  0.056, 0.044, 0.029, 0.018, 0.012, 0.008, 0.006, 0.004, 0.001,
                                    0.,     0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0. };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL1_MAINZ_100[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 5;
  SADC_shape_max_from_front[jshape] = 1.8;
  SADC_shape_sadc_front_gate[jshape] = 3;

  // RD Shape SADC ECAL1 OS ShapeNM100SADC
  jshape = 4;
  double shape_ECAL1_OLGA_100[] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.001, 0.001, 0.002,
                                    0.003, 0.010, 0.092, 0.350, 0.685, 0.923, 1.000, 0.946, 0.825, 0.688,
                                    0.563, 0.459, 0.378, 0.315, 0.265, 0.227, 0.192, 0.165, 0.146, 0.127,
                                    0.109, 0.090, 0.073, 0.057, 0.046, 0.036, 0.030, 0.033, 0.030, 0.027,
                                    0.025, 0.022, 0.019, 0.017, 0.015, 0.013, 0.011, 0.008, 0.006, 0.004 };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL1_OLGA_100[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 8;
  SADC_shape_max_from_front[jshape] = 2.3;
  SADC_shape_sadc_front_gate[jshape] = 3;

  // LED Shape SADC ECAL2 GAMS ShapeNM1000SADC
  jshape = 5;
  double shape_ECAL2_LED_1000[] = { 0.000, 0.000, 0.000, 0.000, 0.001, 0.012, 0.050, 0.125, 0.227, 0.338,
                                    0.447, 0.544, 0.633, 0.707, 0.774, 0.828, 0.878, 0.918, 0.956, 0.984,
                                    1.000, 0.973, 0.906, 0.809, 0.707, 0.610, 0.526, 0.451, 0.388, 0.333,
                                    0.294, 0.260, 0.219, 0.195, 0.178, 0.158, 0.135, 0.112, 0.095, 0.079,
                                    0.068, 0.058, 0.052, 0.045, 0.041, 0.037, 0.033, 0.030, 0.027, 0.025 };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL2_LED_1000[i]);
  }
  SADC_shape_max_position[jshape] = 20;
  SADC_shape_max_from_start[jshape] = 16;
  SADC_shape_max_from_front[jshape] = 9.2;
  SADC_shape_sadc_front_gate[jshape] = 9;

  // LED Shape SADC ECAL1 GAMS ShapeNM500SADC
  jshape = 6;
  double shape_ECAL1_GAMS_LED_500[] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                        0.000, 0.000, 0.011, 0.119, 0.453, 0.853, 1.000, 0.924, 0.768, 0.611,
                                        0.479, 0.376, 0.297, 0.237, 0.191, 0.157, 0.130, 0.110, 0.099, 0.089,
                                        0.071, 0.049, 0.026, 0.008, 0.001, 0.000, 0.000, 0.000, 0.000, 0.000,
                                        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL1_GAMS_LED_500[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 5;
  SADC_shape_max_from_front[jshape] = 1.8;
  SADC_shape_sadc_front_gate[jshape] = 3;

  // LED Shape SADC ECAL1 MT ShapeNM500SADC
  jshape = 7;
  double shape_ECAL1_MAINZ_LED_500[] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                         0.000, 0.000, 0.001, 0.093, 0.510, 0.888, 1.000, 0.917, 0.760, 0.595,
                                         0.455, 0.345, 0.262, 0.203, 0.162, 0.133, 0.111, 0.094, 0.080, 0.074,
                                         0.066, 0.057, 0.043, 0.028, 0.008, 0.003, 0.001, 0.000, 0.000, 0.000,
                                         0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL1_MAINZ_LED_500[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 5;
  SADC_shape_max_from_front[jshape] = 1.95;
  SADC_shape_sadc_front_gate[jshape] = 3;

  // LED Shape SADC ECAL1 OS ShapeNM500SADC
  jshape = 8;
  double shape_ECAL1_OLGA_LED_500[] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.005,
                                        0.020, 0.044, 0.140, 0.401, 0.719, 0.933, 1.000, 0.954, 0.838, 0.695,
                                        0.559, 0.444, 0.352, 0.280, 0.225, 0.183, 0.151, 0.125, 0.107, 0.093,
                                        0.078, 0.063, 0.057, 0.043, 0.028, 0.008, 0.003, 0.001, 0.000, 0.000,
                                        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };
  for( size_t i=0; i<50; i++ ) {
      SADC_shape[jshape].push_back(shape_ECAL1_OLGA_LED_500[i]);
  }
  SADC_shape_max_position[jshape] = 16;
  SADC_shape_max_from_start[jshape] = 8;
  SADC_shape_max_from_front[jshape] = 2.6;
  SADC_shape_sadc_front_gate[jshape] = 3;

  //  Shapes processing
  for ( int ishape = 0; ishape < all_shapes_sofar; ishape++ ) {
      if( debug ) cout << " New Shape " << ishape << endl;
      for ( int is = 0; is< (int)SADC_shape[ishape].size(); is++ ) {
          if( is == 0 ) {
              D_SADC_shape[ishape].push_back(SADC_shape[ishape][1]-SADC_shape[ishape][0]);
          }
          else if( is == (int)SADC_shape[ishape].size()-1 ) {
              D_SADC_shape[ishape].push_back(SADC_shape[ishape][(int)SADC_shape[ishape].size()-1]-SADC_shape[ishape][(int)SADC_shape[ishape].size()-2]);
          }
          else {
              D_SADC_shape[ishape].push_back((SADC_shape[ishape][is+1]-SADC_shape[ishape][is-1])/2.);
          }
      }

      vector<double> value_l;
      vector<double> index_l;
      vector<double> width;
      for ( int il=0; il< SADC_shape_max_position[ishape]; il++ )
        {
        double vl = SADCshape( il, ishape);
        if( vl <= 0. ) continue;
        int ir = il;
        for ( ir=il+1; ir<50; ir++ )
          {
            double vr = SADCshape( ir, ishape );
            if( vr < vl ) break;
          }
        if( debug ) cout <<" For v[" << il <<"]=" << vl <<" Width = " << ir-il+1 << endl;
        value_l.push_back(vl);
        index_l.push_back(double(il));
        width.push_back(double(ir-il+1));
        }
      if( debug ) cout <<" For v[" << SADC_shape_max_position[ishape] <<"]=" <<
                 SADCshape( SADC_shape_max_position[ishape], ishape ) <<" Width = " << 1 << endl;
      value_l.push_back(SADCshape(SADC_shape_max_position[ishape], ishape ) );
      index_l.push_back(double(SADC_shape_max_position[ishape]));
      width.push_back(1.);

      if( debug ) cout << " Try to extend " << ishape << endl;
      for ( int iw=1; iw < (int)width.size() ; iw++ )
      {
        int inv = width.size() - iw -1;
        double width_min =  width[inv+1];
        double width_max =  width[inv];
        double dw = width_max - width_min;
        double delta =  1./dw;

        double v_min = value_l[inv+1];
        double v_max = value_l[inv];
        double dv = value_l[inv]-value_l[inv+1];

        double index_min = index_l[inv+1];
        double index_max = index_l[inv];
        double dindex = index_l[inv]-index_l[inv+1];


        if( debug ) cout << " Try to find between " << width_min << " -- " << width_max << endl;
        if( debug ) cout <<" For v[" << index_min <<"]=" << v_min <<" Width = " << width_min << endl;
        SADC_shape_value_ower[ishape].push_back(v_min);
        SADC_shape_index_ower[ishape].push_back(index_min);
        for ( int iad = 1; width_min+double(iad) < width_max; iad++ )
          {

            double valuel = v_min + double(iad)*delta*dv;
            double indexl = index_min + double(iad)*delta*dindex;
            double widthl = width_min + double(iad)*delta*dw;
            if( debug ) cout <<" For v[" << indexl <<"]=" << valuel <<" Width = " << widthl << endl;
            SADC_shape_value_ower[ishape].push_back(valuel);
            SADC_shape_index_ower[ishape].push_back(indexl);

          }
        if( iw == (int)width.size() -1 )
          {
            if( debug ) cout <<" For v[" << index_max <<"]=" << v_max <<" Width = " << width_max << endl;
            SADC_shape_value_ower[ishape].push_back(v_max);
            SADC_shape_index_ower[ishape].push_back(index_max);
          }
      }

      if( debug ) cout << " X -Check " << ishape << endl;
      for ( int i = 0; i < (int)SADC_shape_value_ower[ishape].size(); i++ )
        {
          if( debug ) cout << " Width " << i+1 <<" Vl " << SADC_shape_value_ower[ishape][i] <<" index = " << SADC_shape_index_ower[ishape][i] << endl;
        }
      if( debug ) cout << " //////////////// " << endl;
    }

  //     cout << " Vse v poryadke, just for debuging ! " << endl;
  //     exit(0);

    int jtail = 0;
    double tail0[] = { 26689, 24393, 21449, 18233, 15619, 13193, 11299, 9573,
                       8279,  7083,  6169,  5323,  4689,  4073,  3649, 3383,
                       3199,  2883,  2529,  2153,  1889,  1703,  1589, 1493,
                       1399,  1263,  1189,  1113,  1089,  1013,   959,  893 };

  cout << " tail0 " << endl;
  for( size_t i=0; i<32; i++ )
    {
      double t = tail0[i]/tail0[0];
      cout << int(100.*t) <<" ";
      SADC_tail[jtail].push_back(t);
    }
  cout << endl;
  jtail = 1;
  double tail1[] = { 2609, 2451, 2139, 1831, 1579, 1371, 1189, 1051,
                     909,  811,  709,  641,  579,  511,  459,  431,
                     389,  361,  329,  311,  269,  261,  239,  221,
                     209,  181,  169,  161,  169,  161,  149,  131 };

  cout << " tail1 " << endl;
  for( size_t i=0; i<32; i++ )
    {
      double t = tail1[i]/tail1[0];
      cout << int(100.*t) <<" ";
      SADC_tail[jtail].push_back(t);
    }
  cout << endl;

}

////////////////////////////////////////////////////////////////////////////////

double  CsDigitizerSADC::SADCshape ( int is, int table )
{
  if( is < 0 || is >= 50 )
    return 0.;
  else
    return SADC_shape[table][is];
}

////////////////////////////////////////////////////////////////////////////////

double  CsDigitizerSADC::Vshape ( int is, double max, double ped ) const
{
  if( is < 0 || is >= 50 )
    return 0.;
  else
    return max*SADC_shape[shape_table_][is];
}

////////////////////////////////////////////////////////////////////////////////

double  CsDigitizerSADC::DVshape ( int is, double max, double ped ) const
{
  if( is < 0 || is >= 50 )
    return 0.;
  else
    return max*D_SADC_shape[shape_table_][is];
}

////////////////////////////////////////////////////////////////////////////////

bool  CsDigitizerSADC::LoadSettings ( const std::vector< double > &settings )
{
//   bool debug = true;
//   if( debug ) cout << " Nu i gde tut vyletetat ? " << settings.size() << endl;
  if( method_ == CsDigitizerSADC::MaxSimple )       // Method MaxSimple
  {
    coeff_convert2max_ = settings[0];
    sadc_ped_min_ = (int)settings[1];
    sadc_ped_max_ = (int)settings[2];
  }  // End  Method MaxSimple
  else if( method_ == CsDigitizerSADC::SummRange )       //  Method SummRange
  {
    coeff_convert2max_ = settings[0];
    sadc_ped_min_ = (int)settings[1];
    sadc_ped_max_ = (int)settings[2];
    sadc_signal_min_ = (int)settings[3];
    sadc_signal_max_ = (int)settings[4];
    sadc_front_ = (int)settings[5];
//    sadc_front_gate_ = (int)settings[6];
  }  // End  Method SummRange

  activated_=true;
  return activated_;
}

////////////////////////////////////////////////////////////////////////////////

//const std::vector< double > &CsDigitizerSADC::Fit ( const std::vector<CS::uint16> &sample )
bool CsDigitizerSADC::Fit ( const std::vector<uint16> &sample, const unsigned int icell)
{

  switch( method_ ) {
  case CsDigitizerSADC::MaxSimple:       // Method MaxSimple
    FitMaxSimple(sample);
    break;
  case CsDigitizerSADC::MaxAdvanced:       //  Method MaxAdvance
    FitMaxAdvanced(sample);
    break;
  case CsDigitizerSADC::SummRange:       //  Method SummRange
    FitSummRange(sample);
    break;
  case CsDigitizerSADC::ShapeFit:       //  Method ShapeFit
    FitShape(sample);
    break;
  case CsDigitizerSADC::FitPulse:
    PulseFit(sample, icell, this->GetParentName());
    break;
  default:
    // This should never happen
    throw Reco::Exception("Unknown 'method_' = '%s' in CsDigitizerSADC::Fit!", method_);
    break;
  }
  if( option_store_stat_info_) StoreStatInfo( sample);
  return ResultIsReady();
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Determine amplitude and time by searching for maximum amplitude
*/

const CsDigitizerSADC::result_t &CsDigitizerSADC::FitMaxSimple ( const std::vector<uint16> &sample )
{

  // Control debuging output
  bool debug = debug_;
  double signal_cut2print4debug = 2.;

  // Reset memory
  Clear();


  double time(0.), ped(0.);

  // Calculate pedestals from baseline
/// \todo  TODO: Update method for new hardware feature (pedestal subtraction)
  for( int is=sadc_ped_min_; is <= sadc_ped_max_; is++)
    {
      ped += sample[is];
    }
  ped /= (sadc_ped_max_-sadc_ped_min_+1);


  // Get the maximum bin and its content
  double signal = -100000.;
  int max_position = -1;
  for( int is=0; is < (int)sample.size(); is++)
    {
      if( sample[is] > signal )
        {
          signal = sample[is];
          max_position = is;
        }
    }

  // Pedestal subtraction
  signal -= ped;

  // Check if amplitude is above threshold
  if( signal > sadc_delta_ )  // !!!!! SADC Delta
    {

      time = CalcTime( Edge, sample, max_position);

    }   // !!!!! End SADC Delta

  // Some debuging output
  if( debug && signal > signal_cut2print4debug )
    cout << " MaxSimple signal " << signal <<
      " pedestal " << ped << " time " << time <<" max " << max_position <<
      " Xcheck " <<  sample[max_position]-ped <<" Sample size " << sample.size() << endl;
  if( debug )
    cout << " MaxSimple final time " << time << endl;

  // Fill result structure
  result_.ampl = signal;
  result_.time = time;
  result_.base = ped;
  result_.max_pos = max_position;
  result_.ready = 1;

  return result_;
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Generate corrected samples and make basic shape analysis
  \todo Simplify structure
*/

void CsDigitizerSADC::GetCcSample ( const std::vector<uint16> &sample ) {
  bool debug = false;
  DigitizerSADCBase::GetCcSample(sample);
  // Find slope
  double signal_slope = -100000.;
  int front_position = -1;
  for( unsigned int is=1; is < sample.size(); is++) {
      // Get position of maximal slope
      // except for the very first sample, when there is no predecessor
      double w =  csample_[is] - csample_[is-1];
      if( w > signal_slope ) {
          signal_slope = w;
          front_position = is;
      }
  }
  // Get new pedestals
  double cped_odd = 0.;
  double cped_even = 0.;
  double stat = 0.;
  int iped_max = -1;
  if(shape_table_ >= 0) {
      iped_max = result_.max_pos - SADC_shape_max_from_start[shape_table_];

      if( iped_max >= 0 ){
          for( int is=0; is <= iped_max; is+=2 ) {
              cped_odd += csample_[is];
              cped_even += csample_[is+1];
              stat++;
          }
          if( stat > 0 ) {
              cped_odd /= stat;
              cped_even /= stat;
          }
      }
  }
  else {
      cerr << " SADC shape was not provided! Sorry it is a MUST for the moment ! shape_table_ = " << shape_table_ << endl;
      exit(1);
  }

  // Make dynamic pedestals corrections in any case
  //VK ???? Check usefull? And what is recepy?
  if( debug ) {
      if( fabs(float(cped_odd)) > 10. || fabs(float(cped_even)) > 10. ) {
          cerr << " CsDigitizerSADC::GetCcSampl eWARNING!! Deviation of ped_even " << cped_even <<
              " ped_odd " << cped_odd << " was detected ped_even_old = " <<
              ped_even_old_ << " ped_odd_old = " << ped_odd_old_ <<endl;
          for( int is=0; is < (int)sample.size(); is++) {
              cout <<" " << sample[is];
          }
          cout << endl;
      }
  }

  // Correct samples and information by dynamic pedestals
  if( stat > 0 ) {
      for( int is = 0; is < (int)csample_.size(); is+=2) {
          csample_[is] -= cped_odd;
          csample_[is+1] -= cped_even;
      }
      result_.base = GetDynamicPed();
      result_.base_diff = ped_odd_old_ + cped_odd -  result_.base;
  }
  else {
      result_.base = 0.;
      result_.base_diff = 0.;
  }

  // More technical info to output
  front_position_= front_position;
  signal_slope_= signal_slope;

  if( debug ) cout << endl;
}

//////////////////////////////////////////////////////////////////////////////////
bool CsDigitizerSADC::Shape_FilterSLC ( const std::vector<uint16> &smpl ) const
{
  if(smpl.size() != 32 ) {
    cerr <<" Shape_FilterSLC Bad sample size " << smpl.size() << endl;
    return false;
  }
  bool noise = false;
  int jmax = std::max_element(smpl.begin(), smpl.end()) - smpl.begin();
  int jmin = std::min_element(smpl.begin(), smpl.end()) - smpl.begin();
  int max = smpl[jmax];
  int min = smpl[jmin];

  if (jmax <= 24 && jmax >= 5) {
      if ((max - min) > 25) {
          if (((smpl[jmax-5] - min) >= (max - min) / 2) || ((smpl[jmax-5] - min) >= 10)) {
              noise = true;
          }
          else {
              noise = false;
          }
      }
      else { // here amp<=25
          if (((smpl[jmax-4] - min) >= (max - min) /2 ) || ((smpl[jmax-4] - min) >= 3)) {
              noise = true;
          }
          else {
              if (jmax >= 8) {
                  if ( (smpl[jmax-8] - min) <= (max-min) / 3 ) {
                      noise = false;
                  }
                  else {
                      noise=true;
                  }
              }
              else { // jmax<=7
                  if ((smpl[jmax+12]-min)<=5) {
                      noise = false;
                  }
                  else {
                      noise = true;
                  }
              }
          }
      }
  }
  else {  // here jmax>=25 or jmax <=4, only noise
      noise=true;
  }
  return noise;
}

//////////////////////////////////////////////////////////////////////////////////

const CsDigitizerSADC::result_t &CsDigitizerSADC::FitMaxAdvanced ( const std::vector<uint16> &sample )
{

  // Control debuging
  bool debug = debug_;
  double signal_cut2print4debug = 2.;

  if( debug )
    cout << " CsDigitizerSADC::FitMaxAdvanced debug " << endl;

  // Reset
  Clear();

  if( debug )
    cout << " CsDigitizerSADC::FitMaxAdvanced Clear() OK " << endl;
  if( debug )
    PrintSettings();

  bool shape_filter_slc = Shape_FilterSLC(sample);
  if(use_shape_filter_slc_) is_noise_by_shape_anal_ = shape_filter_slc;
  // Production variant
  //   if(use_shape_filter_slc_)
  //   {
  //     is_noise_by_shape_anal_ = Shape_FilterSLC(sample);
  //   }


  // Do basline subtraction and base shape analysis
  GetCcSample( sample );


  if( debug ) cout << " CsDigitizerSADC::FitMaxAdvanced overflow_amplitude_overflow_amplitude_overflow_amplitude_ " << overflow_amplitude_ << endl;


  // Detect overflows
  {

    unsigned int new_over = DetectOverflow(sample);

    // Treat overflows
    if( new_over  )
      {
        if( debug )
          {
            cout << " NOverflow =  " << overflow_samples_.size() << " Overflow? " << endl;
            for( int is=0; is < (int)sample.size(); is++)
              {
                cout <<" " << sample[is];
              }
            cout << endl;

            cout << " Csample " << endl;
            for ( int is=0; is < (int)sample.size(); is++ )
              {
                cout << " " << (int)(10.*csample_[is]);
              }
            cout << endl;
          }

        if( new_over > 1 )
          {
            cerr << " Two samples? And both in overflow ?? Ne byvaet !.. " << endl;
            PrintSample( sample );
          }

        // Calculate time
        int plato_size = overflow_samples_.back() - overflow_samples_.front() + 1;

        if( debug )
          cout << " Overflow detected from sample["<< overflow_samples_.front() <<"]="<<sample[overflow_samples_.front()] <<
            "to sample["<< overflow_samples_.back() <<"]="<<sample[overflow_samples_.back()] <<
            " with plato size = "<< plato_size << endl;

        double value_ower = SADC_shape_value_ower[shape_table_][plato_size];
        double index_ower = SADC_shape_index_ower[shape_table_][plato_size];
        double time = overflow_samples_.front() + SADC_shape_max_position[shape_table_]
          - SADC_shape_max_from_front[shape_table_] - index_ower;

        if( debug )
          cout << " SADC_shape_max_position = " << SADC_shape_max_position[shape_table_] <<
            " SADC_shape_max_from_front = " << SADC_shape_max_from_front[shape_table_] <<
            " from index_ower to front = " << SADC_shape_max_position[shape_table_]-SADC_shape_max_from_front[shape_table_]-index_ower << endl
               << " Time ower in clocks before TCS correction = " << time << endl;

        // Convert to ns
        time *= sadc_clock_;

        if( debug )
          cout << " Over Value_ower " << value_ower << " index_ower " << index_ower << endl;

        double ped = (ped_odd_old_+ped_even_old_)/2.;
        double csignal = overflow_amplitude_ - ped;

        if( debug )
          cout <<" Csample value at first over point(calculated!) " << csignal << endl;

        double signal = csignal/value_ower;

        // Prepare for precise corrections, need more investigations, skip for a while Fri Feb 13 13:34:05 CET 2009
        if( debug )
          cout << " Value " << signal << endl
               << "********************" << endl
               << " Value + ped " << signal + ped <<  endl
               << "--------------------" << endl;

        // Set return structure
        result_.ampl = signal;
        result_.time = time;
        result_.ready = 1;

        // Return after overflow
        return result_;

      }  // End  Method FitOverflow

  }

  double time = 0.;

  if( debug ) cout << " sample " << sample[0] <<" " << sample[1] <<" " << sample[2] <<" " << sample[3] <<" " << sample[4] <<" "<< endl;


  double signal = result_.ampl;
  int max_position = result_.max_pos;
  int front_position = front_position_;


  // calculate chi square
  int ishape0 = (int)(max_position - SADC_shape_max_position[shape_table_]);
  double max = signal;
  double chisq = 0.;
  const double ped_ = result_.base + result_.base_diff / 2.;

  for ( int is=0; is < (int)sample.size(); is++ )
    {
      double f = Vshape( is-ishape0, max, ped_ );
      double r = csample_[is];
      double dvt = DVshape( is-ishape0, max, ped_ );
      dvt = dvt*dvt;
      double d =  r - f;
      chisq += d*d/(1.+dvt);
    }

  if( debug )
    {
      cout << " Csample " << endl;
      for ( int is=0; is < (int)sample.size(); is++ )
        {
          cout << " " << (int)(10.*csample_[is]);
        }
      cout << endl;
      cout << " Vshape " << endl;
      for ( int is=0; is < (int)sample.size(); is++ )
        {
          cout << " " << (int)(10.*Vshape( is-ishape0, max, ped_ ));
        }
      cout << endl;
      cout << " Initial chisq " << chisq << " prob " << TMath::Prob( chisq, 31)  << endl;
      }

  if(use_shape_filter_line_fit_)
  { // Calculate streight line chi-square to compare to shape
    double chisq_line_new = 0.;
    double time_line = 0.;
    double signal_line = 0.;
    {

      double sam_s = 0.;
      double sam_ss = 0.;
      double sam_sf = 0.;
      double sam_f = 0.;
      double sam_ff = 0.;

      double sam_n = (double)csample_.size();
      for ( unsigned int i=0; i < csample_.size(); i++ )
        {
          double s = (double)(csample_[i]);
          double f = (double)(i);
          sam_s += s;
          sam_ss += s*s;
          sam_sf += s*f;
          sam_f += f;
          sam_ff += f*f;
        }
      sam_s /= sam_n;
      sam_ss /= sam_n;
      sam_sf /= sam_n;
      sam_f /= sam_n;
      sam_ff /= sam_n;

      double a_new = ( sam_sf - sam_f*sam_s )/( sam_ff - sam_f*sam_f );
      double ped_new = sam_s -  a_new*sam_f;
      for ( int i=0; i < (int)csample_.size(); i++ )
        {
          double d =  (double)((double)(csample_[i])-ped_new ) - a_new*double(i);
          chisq_line_new += d*d;
        }
      if( a_new >= 0. )
        {
          time_line = (double)(csample_.size());
          signal_line = a_new*(double)(csample_.size());
          time_line *= sadc_clock_;
          //       double tcs_cor = CsEvent::Instance()->getTCSPhaseTime()-TCS_T0;
          //       if( MAKE_TCS_CORRECTIONS ) time_line += tcs_cor;
        }
      else
        {
          time_line = -2.;
          //      signal_line = ped_new + a_new*time_line;
          signal_line = -a_new*(double)(csample_.size());
          time_line *= sadc_clock_;
          //       double tcs_cor = CsEvent::Instance()->getTCSPhaseTime()-TCS_T0;
          //       if( MAKE_TCS_CORRECTIONS ) time_line += tcs_cor;
        }
      if( debug )   cout << " chisq_line_new " << chisq_line_new << endl;
    }

    if( chisq_line_new < chisq )
      {
        result_.time = time_line;
        result_.ampl = signal_line;
        result_.ready = 1;
// ??
        is_noise_by_shape_anal_= true;
        if( debug )
          {
            cout << " Staight line fit looks better then Shape fit " << endl;
            cout << " Considering this amlitude as Garbaje we need to return something Thinking ... " << endl;
            cout << " Signal out of range! MaxAdvanced signal " << GetSignal() <<
              " pedestal " << GetPed() << " time(~ in clocks) " << (int)(GetTime()/sadc_clock_) <<
              " max " << max_position <<
              " Xcheck " <<  sample[max_position]-ped_ <<" Sample size " << sample.size() << endl;
          }
        return result_;
      }
  }

  if( debug ) cout << " MaxAdvanced signal " << signal << " max_position " << max_position <<
                " front_position " << front_position <<" diff " << max_position-front_position-0.5 << endl;

  double expected_front_position = max_position - SADC_shape_max_from_front[shape_table_];

  if( debug ) cout << " front_position " << front_position  <<
                " expected front_position " << max_position - SADC_shape_max_from_front[shape_table_] <<
                " sadc_front_gate " << SADC_shape_sadc_front_gate[shape_table_] << endl;

  // Start of signal fits expected start of signal
  if( fabs(float(expected_front_position - front_position)) < SADC_shape_sadc_front_gate[shape_table_] )
    {

/// \todo  TODO: remove old time calculation
      // Not used at the moment
      int front_min = front_position - SADC_shape_sadc_front_gate[shape_table_] - 1;
      int front_max = front_position + SADC_shape_sadc_front_gate[shape_table_] + 1;

    }
  // Start of signal does not match expectations
/// \todo  TODO calculate time for mismatch between expected and real signal time
  time = CalcTime( HalfMax, sample, (unsigned int)max_position);

  if( !isfinite(time) )
    {
      time = expected_front_position;
      // GoTo ns
      time *= sadc_clock_;
    }


  double sxx=0.;
  double sy=0.;
  double fmax = csample_[max_position];
  for( int is=max_position-2; is <= max_position+2; is++)
    {
      if( is <= 0 ) continue;
      if( is >= (int)sample.size() ) break;
      double dx = is-max_position;
      sxx += dx*dx;
      double dy = csample_[is]-fmax;
      sy += dy;
    }

  double a = 0.;
  if( sxx > 0 ) a = sy/sxx;

  if( debug ) cout << " MaxAdvanced a in max " << a << " max[" <<max_position<<"] = "<< csample_[max_position] << endl;
  if( a < 0 )
    {
      double ss = 0.;
      double sx = 0.;
      double sy = 0.;
      double syx = 0.;
      double sxx = 0.;
      double sxxx = 0.;

      // Recalculate signal
      for( int is=max_position-2; is <= max_position+2; is++)
        {
          if( is <= 0 ) continue;
          if( is >= (int)sample.size() ) break;
          double dx = is-max_position;
          double dx2 = dx*dx;
          double dy = csample_[is]-fmax;
          if( debug ) cout << " dx " << dx <<" dy " << dy << endl;
          ss += 1.;
          sx += dx;
          sy += dy;
          syx += dy*dx;
          sxx += dx2;
          sxxx += dx2*dx;
        }

      if( ss <= 0 )
        {
          cerr << "  EROOR EROOR EROOR EROOR MaxAdvanced internal error ss = " << ss << endl;
        }
      else
        {
          sx /= ss;
          sy /= ss;
          syx /= ss;
          sxx /= ss;
          sxxx /= ss;
          double b =  ( syx - a*sxxx - sy*sx + a*sxx*sx )/(sxx-sx*sx);
          double c =  sy - a*sxx - sx*b;
          double xmax = -b/(2*a);
          double ymax = c- a*xmax*xmax;
          if( debug ) cout << " MaxAdvanced xmax = " << xmax << " ymax = " << ymax << endl;
          if( debug )
            {
              cout << " X check max_position = " << max_position << " fmax " << fmax <<
                " X check1 " << csample_[max_position] << " X check2 signal " << signal << endl;
            }

          if( xmax < -2 ||  xmax > 2 || ymax < -1. || ymax > 0.3*signal || a > - 0.19 )
            {
              if ( debug ) {
                cout << " Not expected: MaxAdvanced xmax = " << xmax << " ymax = " << ymax <<" a= " << a << endl;
                cout << " X check max_position = " << max_position << " fmax " << fmax <<
                  " X check1 " << csample_[max_position] << " X check2 signal " << signal << endl;
              }
            }
          else
            {
              if( debug )
                {
                  cout << " Expected: MaxAdvanced xmax = " << xmax << " ymax = " << ymax <<" a= " << a <<  endl;
                  cout << " X check max_position = " << max_position << " fmax " << fmax <<
                    " X check1 " << csample_[max_position] << " X check2 signal " << signal << endl;
                }
              signal += ymax;
            }
        }
    }

    // Prepare result
  result_.ampl = signal;
  result_.time = time;
  result_.ready = 1;
  if( debug && signal > signal_cut2print4debug )
    cout << " MaxAdvanced signal " << GetSignal() <<
      " pedestal " << GetPed() << " time(~ in clocks) " << (int)(time/sadc_clock_) <<" max " << max_position <<
      " Xcheck " <<  sample[max_position]-ped_ <<" Sample size " << sample.size() << endl;
  if( debug ) cout << " MaxAdvanced time " << GetTime() << endl;

  return result_;
}

////////////////////////////////////////////////////////////////////////////////

const CsDigitizerSADC::result_t &CsDigitizerSADC::FitSummRange ( const std::vector<uint16> &sample )
{

  // Control debug output
  bool debug = false;
  double signal_cut2print4debug = 2.;

  // Reset module
  Clear();


  double time = 0.;

  // Calculate pedestal
  double ped =0.;
  for( int is=sadc_ped_min_; is <= sadc_ped_max_; is++)
    {
      ped += sample[is];
    }
  ped /= (sadc_ped_max_-sadc_ped_min_+1);

  // Calculate signal amplitude as sum over range from sadc_signal_min_ to sadc_signal_max_
  double signal = 0.;
  for( int is=sadc_signal_min_; is <= sadc_signal_max_; is++)
    {
      signal += sample[is] - ped;
    }

  // Check threshold
  if( signal*coeff_convert2max_ > sadc_delta_ )  // !!!!! SADC Delta
    {

      // Calculate time edge finder
      double ww = 0.;
      for( int is=sadc_front_ - SADC_shape_sadc_front_gate[shape_table_]; is <= sadc_front_ + SADC_shape_sadc_front_gate[shape_table_]; is++)
        {
          double w =  sample[is] - sample[is-1];
        if( w > 0 )
          {
            time += w*(is - sadc_front_);
            ww += w;
          }
        }
      if( ww > 0 ) time /= ww;

      if( debug && signal > signal_cut2print4debug )
        cout << " Relative to front signal " << signal << " time " << time << " ww " << ww << endl;
      time += sadc_front_; // switch back to absolute SADC time ( counting from 0 sample
      time *= sadc_clock_;

      if( debug && signal > signal_cut2print4debug ) cout << " time " << time << endl;

      signal = coeff_convert2max_*signal;

    }   // !!!!! End SADC Delta


  if( debug && signal > signal_cut2print4debug )
    cout << " SummRange signal " << signal << " pedestal " << ped << " time " << time << " sadc_clock_ " << sadc_clock_ << endl;

/// \todo  signal maybe unconverted, is this correct?
  result_.ampl = signal;
  result_.time = time;
  result_.base = ped;
  result_.ready = 1;

  return result_;
}

////////////////////////////////////////////////////////////////////////////////

const CsDigitizerSADC::result_t &CsDigitizerSADC::FitShape ( const std::vector<uint16> &sample )
{

  // Control debugging output
  bool debug = debug_;
  double signal_cut2print4debug = 2;

  // Reset module
  Clear();


  double time = 0.;

  if( debug ) cout << " sample " << sample[0] <<" " << sample[1] <<" " << sample[2] <<" " << sample[3] <<" " << sample[4] <<" "<< endl;

  // Get pedestals
  double ped_odd =0.;
  double stat_odd =0.;
  for( int is=0; is <= sadc_ped_max_; is+=2 )
    {
      ped_odd += sample[is];
      stat_odd++;
    }
  if( stat_odd > 0 ) ped_odd /= stat_odd;

  double ped_even =0.;
  double stat_even =0.;
  for( int is=1; is <= sadc_ped_max_; is+=2 )
    {
      ped_even += sample[is];
      stat_even++;
    }
  if( stat_even > 0 ) ped_even /= stat_even++;

  double ped( (ped_odd+ped_even)/2. );
  double dpedd = (ped_even-ped_odd)/2.;

  if( debug ) cout << " ped " << ped << " ped_odd " << ped_odd << " ped_even " << ped_even <<" dpedd " << dpedd  << endl;

  // correct sample
  vector<double> csample;
  if( debug ) cout << " Csample ";
  for( int is=0; is < (int)sample.size(); is++)
    {
      if( (is%2) > 0 )
        csample.push_back( double(sample[is]) - dpedd - ped);
      else
        csample.push_back( double(sample[is]) + dpedd - ped );
      if( debug ) cout << (int)(10.*csample.back() ) <<" ";
    }
  if( debug ) cout << endl;


  // Do basic pulse shape analysis
  double signal = -100000.;
  int max_position = -1;
  double signal_slope = -100000.;
  int front_position = -1;
  for( int is=0; is < (int)sample.size(); is++)
    {
      if( csample[is] > signal )
        {
          signal = csample[is];
          max_position = is;
        }
      if( is > 0 )
        {
          double w =  csample[is] - csample[is-1];
          if( w > signal_slope )
            {
              signal_slope = w;
              front_position = is;
            }
        }
    }


  if( debug ) cout << " FitShape signal " << signal << " max_position " << max_position <<
                " front_position " << front_position <<" diff " << max_position-front_position-0.5 << endl;

  // Calculate chisq
  int ishape0 = (int)(max_position - SADC_shape_max_position[shape_table_]);
  double max = signal;
  double chisq = 0.;
  for ( int is=0; is < (int)sample.size(); is++ )
    {
      double f = Vshape( is-ishape0, max, ped );
      double r = csample[is];
      double dvt = SADC_shape_factor*r;
      dvt = dvt*dvt;
      double d =  r - f;
      chisq += d*d/(1.+dvt);
    }

  if( debug )
  {
    for ( int is=0; is < (int)sample.size(); is++ )
    {
      if( debug ) cout << " " << (int)(csample[is]);
    }
    cout << endl;
    for ( int is=0; is < (int)sample.size(); is++ )
    {
      cout << " " << (int)(Vshape( is-ishape0, max, ped ));
    }
    cout << endl;
    cout << " Initial chisq " << chisq << endl;
  }


  if( debug && signal > signal_cut2print4debug ) cout << " FitShape signal " << signal << " pedestal " << ped << " time " << time << " sadc_clock_ " << sadc_clock_ << endl;

  // Prepare result
/// \todo  Not time calculated in CsDigitizerSADC::FitShape!
  result_.ampl = signal;
  result_.time = time;
  result_.base = ped;
  result_.dbase = dpedd;
  result_.max_pos = max_position;
  result_.ready = 1;

  return result_;
}

////////////////////////////////////////////////////////////////////////////////

const CsDigitizerSADC::result_t &CsDigitizerSADC::FitOverflow ( const std::vector<uint16> &sample )
{
  bool debug = true;
  Clear();
  if( true )       //  Method FitOverflow
  {

    int max_amplitude_ = overflow_amplitude_;
    vector <int> over;
    vector <int> nover;

    double maxx = 0.;
    for( int is=0; is < (int)sample.size(); is++)
    {
      if( sample[is] > maxx )
        maxx = sample[is];
      if( sample[is] >= max_amplitude_ )
      {
        if( (over.size() == 0) or (over.back() < is-1) )
        {
          nover.push_back(is);
        }
        over.push_back(is);
      }
    }
    if( over.size() <= 0 )
      return result_;
    if( nover.size() > 1 )
    {
      cerr << " Two samples? And both in overflow ?? Ne byvaet !.. " << endl;
      PrintSample( sample );
//      exit(1);
    }


    if( debug )
    {
      cout << " NOverflow =  " << over.size() << " Overflow? " << endl;
      PrintSample(sample);
    }

    int plato_size = over.back() - over[0] + 1;

    double value_ower = SADC_shape_value_ower[shape_table_][plato_size];
    double index_ower = SADC_shape_index_ower[shape_table_][plato_size];

    cout << " Over Value_ower " << value_ower << " index_ower " << index_ower << endl;
    PrintSample( sample );

//    GetCcSample( sample );

    const double ped_ = result_.base+result_.base_diff/2.;
    double value = (overflow_amplitude_-ped_)/value_ower;

    cout << " Value " << value << endl;
    cout << "********************" << endl;
    cout << " Value + ped " << value + ped_ << " Maxx " << maxx <<  endl;
    cout << "--------------------" << endl;

//     for ( int il=0; il< SADC_shape_max_position[shape_table_]; il++ )
//     {
//       double vl = Vshape( il, 1., 0.);
//       if( vl <= 0. ) continue;
//       int ir = il;
//       double vr = vl;
//       for ( ir=il+1; ir<50; ir++ )
//       {
//         double vr = Vshape( ir, 1., 0.);
//         if( vr < vl ) break;
//       }
//       cout <<" For v[" << il <<"]=" << vl <<" Width = " << ir-il << endl;
//     }
    return result_;
  }  // End  Method FitOverflow

  return result_;
}

////////////////////////////////////////////////////////////////////////////////

std::vector< double > CsDigitizerSADC::GetCShape ( int size ) const
{
  std::vector< double > cshape;
  if( method_ == CsDigitizerSADC::ShapeFit && ResultIsReady() )       //  Method FitShape
  {
    int max_position = (int)GetMaxPosition();
    int ishape0 = (int)(max_position - SADC_shape_max_position[shape_table_]);
    double ped = GetPed();
    double max = GetSignal();
    for ( int is=0; is < size; is++ )
    {
      double f = Vshape( is-ishape0, max, ped );
      cshape.push_back(f);
    }

  }
  return cshape;
}

////////////////////////////////////////////////////////////////////////////////

void CsDigitizerSADC::PrintSettings ( void ) const
{
  DigitizerSADCBase::PrintSettings();
//   cout << " sadc_format_version_ " << sadc_format_version_;
//   cout << " sadc_decode_version_ " << sadc_decode_version_;
//  cout << " SADC::  Type " << type_;
  cout << " Method " << method_;
  cout << " sadc_front_ " << sadc_front_;
//  cout << " sadc_front_gate_ " << sadc_front_gate_;
  cout << " sadc_ped_min_ " << sadc_ped_min_;
//  cout << " sadc_ped_max_ " << sadc_ped_max_;
  cout << " sadc_signal_min_ " << sadc_signal_min_;
  cout << " sadc_signal_max_ " << sadc_signal_max_;
  cout << " sadc_delta_ " << sadc_delta_;
//  cout << " sadc_clock_ " << sadc_clock_;
  cout << " coeff_convert2max_ " << coeff_convert2max_;
//   cout << " sadc_max_position_ " << sadc_max_position_;
  cout << endl;
}

////////////////////////////////////////////////////////////////////////////////

void CsDigitizerSADC::StoreStatInfo (  const std::vector<uint16> &sample )
{
  if( !option_store_stat_info_ ) return;
  // Subtract pedestals
  for( int is=0; is < (int)sample.size(); is++)
  {
    if( is >5 ) break;
    double s = double(sample[is]);
    if( (is%2) > 0 ) // even
    {
      stat_ped_odd_.Add(s);
      stat_dped_.Add(s-double(sample[is-1]));
    }
    else   // odd
    {
      stat_ped_even_.Add(s);
    }
//    if( debug ) cout << (int)(10.*sample.back() ) <<" ";
  }

  if( !ResultIsReady() ) return;
  stat_time_.Add( GetTime() );
  bool print_sample = false;
// There is clearly something to watch and improve!!! Wed Oct 26 12:27:01 CEST 2011 VK
//  if(csample_[0] > 10. ) print_sample = true;
  if( print_sample )
  {
    cout <<" Ped deviation in LED sample of Calorimeter " << parent_->GetName() << " Method " << method_ << endl;
  }
  for( int is=0; is < (int)csample_.size(); is++)
  {
    if( is >5 ) break;
    double s = double(csample_[is]);
    stat_ped_.Add(s);
  }
  if( print_sample )
  {
    for( int is=0; is < (int)csample_.size(); is++)
    {
//      if( is >5 ) break;
      cout << (int)(10.*csample_[is] ) <<" ";
    }
    cout << endl;
    cout <<" And here is original sample " << endl;
    for( int is=0; is < (int)csample_.size(); is++)
    {
//      if( is >5 ) break;
      cout << sample[is] <<" ";
    }
    cout << endl;
  }

}

////////////////////////////////////////////////////////////////////////////////

/*!
  \brief Simple function to extract signal using cfd algorithm
  \param sample {A baseline subtracted list of samples}
  \return {Pointer to the result. Null on no signal. (Has to be delted after use.)}
*/

CsDigitizerSADC::cfd_result_t*
CsDigitizerSADC::CFD(const CsDigitizerSADC::cfd_option_t &options,
                   const std::vector<unsigned int> &sample ){

  // Point to be returned initialized at first founfd signal
  cfd_result_t *res(NULL);

  // CFD implementation
  {

    // To remeber last slice
    float last_result(0);

    // Loop over valid samples defined by delay and sample.size()
    // Can be further restricrted by user seting smin and smax
    for( unsigned int i = max(options.delay, options.smin);
         i < min(sample.size(),(size_t)options.smax);
         i++) {

      const float this_result = sample[i] - sample[i - options.delay] * options.ampl;
      // CFD condition
      if( last_result > 0 && this_result <= 0 ) {
        // Threshold check
        if ( sample[i] > options.thr ) {

          // We have a signal so lets check if res exists
          if( !res ) {
            res = new cfd_result_t;
            res->ampl=0;
          }

          // Check for biggest signal
          if ( res->ampl < sample[i] ) {
            res->ampl = sample[i];
            res->time = i + ( this_result / (last_result - this_result) );
          }

        }
      }

      last_result = this_result;

    }

  }

  return res;


};

//////////////////////////////////////////////////////////////////////////////////////////////////////

double CsDigitizerSADC::CalcTime(CsDigitizerSADC::TimeExtractionMethod method,
                               const vector<short unsigned int> &sample,
                               const unsigned int max_position) {

  switch( method ) {
  case Edge:
    return CalcTime_edge(sample, max_position);
    break;
  case HalfMax:
    return CalcTime_halfMax(sample, max_position);
    break;
  default:
    return NAN;
    break;
  }

  return 0.;

}

double CsDigitizerSADC::CalcTime_edge(const vector<short unsigned int> &sample,
                                   const unsigned int max_position) {

  double time(0.);

  // Time calclculation (Edge finder)
  if( sadc_ped_max_ < (int)max_position )
    {
      // Calculate time edge finder

      double ww(0.);
      for( int is=sadc_front_ - SADC_shape_sadc_front_gate[shape_table_];
           is <= sadc_front_ + SADC_shape_sadc_front_gate[shape_table_]; is++)
        {
          double w =  sample[is] - sample[is-1];
          if( w > 0 )
            {
              time += w*(is - sadc_front_);
              ww += w;
            }
        }
      if( ww > 0 )
        time /= ww;
      else return NAN;

      time += sadc_front_; // switch back to absolute SADC time ( counting from 0 sample

    }
  else
    time = max_position - SADC_shape_sadc_front_gate[shape_table_];

  time *= sadc_clock_;

  return time;

}

double CsDigitizerSADC::CalcTime_halfMax(const vector<short unsigned int> &sample, const unsigned int max_position) {

  double time(0.);

  double fmax = csample_[max_position];

  int front_min = front_position_ - SADC_shape_sadc_front_gate[shape_table_] - 1;
  int front_max = front_position_ + SADC_shape_sadc_front_gate[shape_table_] + 1;

  for( unsigned int is = max( 1, front_min);
       (int)is < min( front_max, (int)sample.size());
       is++)
    {

      // Get time of half max crossing
      if( csample_[is-1] <= fmax/2. && csample_[is] > fmax/2. )
      {
          time = is-1 + (fmax/2.- csample_[is-1])/(csample_[is]- csample_[is-1]);
          break;
      }
    }

  time *= sadc_clock_;

  return time;

}

vector<short unsigned int> Correct_Background_Exp(const vector<short unsigned int> &_sample) {

/// \todo  This function only works for hardware based Baseline suptraction (at the moment)
/// \todo  TODO make it configurable in case of advantages

  // Time constant
  double tau[2];
  vector<short unsigned int> sample = _sample;


  // Baseline set by frontends
  short unsigned int baseline = 50;

  // Subtract baseline
  for( vector<short unsigned int>::iterator it = sample.begin();
       it != sample.end(); it++)
    *it -= baseline;

/// \todo  TODO replace tau extraction by fit or by parameter
/// \todo  TODO include error handling

  // Get tau for odd and even (should be identical except for errors)
  tau[0] = -2. / log( sample[2] / sample[0]);
  tau[1] = -2. / log( sample[3] / sample[1]);
  double Tau = ( tau[0] + tau[1] ) / 2.;

  // Correct samples again
  for( unsigned int i = 1;
       i < sample.size(); i++) {
    int corr = (int)round( sample[0] * exp( - i / Tau ) );
    sample[i] -= corr ;
  }
  return sample;

}


//////////////////////////////////////////////////////////

/// \todo  CsDigitizerSADC::PulseFit is not yet in production status!!!!!

const CsDigitizerSADC::result_t &
CsDigitizerSADC::PulseFit(const std::vector<uint16> &sample, const unsigned int icell, const string &DAQDetName ) {

  assert( parent_->GetCalID() != 0 );


  const int id(parent_->GetCalID()),
    x(parent_->GetColumnOfCell( icell )),
    y(parent_->GetRowOfCell( icell ));

  ClearResult();

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,00)

  // Initilize spline
  if( !PulseManager::Instance()->GetPulse(id,0) ) {

    {
      list<string> spline_name;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "SPLINE_FILE", spline_name ) ) {
        unsigned int i = 0;
        for( list<string>::iterator it = spline_name.begin(); it != spline_name.end(); ++it) {
          try {
            cout << "CsDigitizerSADC::PulseFit: Register " << *it
                 << " as " <<  parent_->GetName() << ":" << i << endl;
            MyPulse *mypulse = new MyPulse(it->c_str());
            if( mypulse )
              PulseManager::Instance()->Insert(id, i, mypulse );
            else
              throw -66;
            i++;
          } catch ( int &err ) {
            /*
            for( vector< MyPulse*>::iterator it=parent_->mypulse_.begin();
                 it!=parent_->mypulse_.end(); ++it) {
              if (*it) {
                delete *it;
                *it = NULL;
              }
            }
            parent_->mypulse_.clear();
            */
            throw Reco::Exception("%s in %s at %i : Catched %i while init mypulse_\n",
                                  __FILE__, __FUNCTION__, __LINE__, err);
          }
        }
      } else {
        throw Reco::Exception("%s in %s at %i : Uncomplete configuration : No SPLINE_FILE_NAME for %s\n",
                              __FILE__, __FUNCTION__, __LINE__, parent_->GetName().c_str());
      }

    }

  }


  // Get other configurations
  //cout << DAQDetName << endl;
  if( plane_ == -1 ) {
    if( id == 1 ) // ECAL01
      plane_ = atoi(DAQDetName.substr(5,2).c_str());
    else if ( id == 2 ) { //ECAL02
      float xcell = x - 30.7, ycell = y - 23.4;
      if ( ( xcell > -23 && xcell < -15 && ycell >- 5  && ycell < 5 )
           ||( xcell > -15 && xcell <-11 && ycell > -8  && ycell < 8 )
           ||( xcell > -11 && xcell < 13 && ycell > -11 && ycell < 12 )
           ||( xcell > 13 && xcell < 17 && ycell > -8  && ycell < 8 ) )
        plane_ = 0;
      else if ( xcell > 17 && xcell < 25 && ycell > -5 && ycell < 5 )
        plane_ = 1;
      else
        plane_ = 2;
    } else
      plane_ = 0;
  }
  if( !mypulse_ ) {
    mypulse_ = PulseManager::Instance()->GetPulse(id, plane_);
    {
      string value;
      if( CsOpt::Instance()->getOpt( "CALDEBUG", "FIT_OFILE", value ) )
        {
          asprintf( &fit_ofile_, "%s", value.c_str());
        } else {
          fit_ofile_ = NULL;
        }
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_BASE", value ) )
        {
          fit_base_ = atof(value.c_str());
        } else {
          fit_base_ = 0.;
        }
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_RANGE_TO_FRONT", value ) )
        {
          fit_range_to_front_ = -atof(value.c_str());
        }
      else
        fit_range_to_front_ = -10.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_RANGE_TO_BACK", value ) )
        {
          fit_range_to_back_ = atof(value.c_str());
        }
      else
        fit_range_to_back_ = 10.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_CUT_SLOPE", value ) )
        {
          fit_cut_slope_ = atof(value.c_str());
        }
      else
        fit_cut_slope_ = 1000.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_CUT_OFFSET", value ) )
        {
          fit_cut_offset_ = atof(value.c_str());
        }
      else
        fit_cut_offset_ = 1000.;
       if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_CUT_DEBUG", value ) )
        {
          fit_cut_debug_ = atof(value.c_str());
        }
      else
        fit_cut_debug_ = 0.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_SAMPLE_CERR", value ) )
        {
          fit_sample_cerr_ = atof(value.c_str());
        }
      else
        fit_sample_cerr_ = 7.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_SAMPLE_DERR", value ) )
        {
          fit_sample_derr_ = atof(value.c_str());
        }
      else
        fit_sample_derr_ = 0.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "AC_CUT_MIN", value ) )
        {
          ac_cut_min_ = atof(value.c_str());
        }
      else
        ac_cut_min_ = 0.;

    }


  }


  if( !mypulse_  )
    throw Reco::Exception("%s in %s at %i : Uncomplete configuration : %s: No SPLINE for plane %i (%s)\n",
                          __FILE__, __FUNCTION__, __LINE__, parent_->GetName().c_str(), plane_, DAQDetName.c_str());

  TF1 mypulse_func_("f",mypulse_,0.,32.,3,"MyPulse");


  // Only continue with non empty sample
  if( !sample.size() )
    return result_;

  unsigned short Sample[sample.size()];

  double base[2] = {fit_base_,fit_base_};
  if ( fit_base_ == 0. ) {
    unsigned int use_base_diff = 6; // has to be even
    {
      unsigned int nbase = 0;
      for ( unsigned int i = 0; i < sample.size() && i < use_base_diff; i+=2 ) {
        base[0] += sample[i];
        base[1] += sample[i+1];
        nbase++;
      }
      base[0] /= nbase;
      base[1] /= nbase;
    }
  }

  // Generate histogram
  TGraphAsymmErrors g( sample.size() );

  // Get starting parameters for fit and fill histogram
  list<int> localMax;
  {
    int max = 0;
    for ( unsigned int i = 0; i < sample.size(); i++ ) {
      Sample[i] = sample[i];
      if ( max < 0 && sample[i] >= sample[i] - 1 )
        max = i;
      else if ( sample[max] < sample[i] )
        max = i;
      else if ( ( (sample[max]-base[max&1]) * 0.8 > (sample[i]-base[i&1]) ) ||
                ( i + 1 ==  sample.size() ) ) {
        localMax.push_back(max);
        max = -1;
      }

      g.SetPoint(i, i, sample[i] - base[i&1] );
      double ampl_err[2] = {
        (fit_sample_cerr_+fit_sample_derr_*sqrt(sample[i] - base[i&1])),
        (fit_sample_cerr_+fit_sample_derr_*sqrt(sample[i] - base[i&1])),
      };
      if( sample[i] == overflow_amplitude_ )
        ampl_err[1] = 4. * overflow_amplitude_;
      g.SetPointError(i, 0.3 / 12.6,  0.3 / 12.6, ampl_err[0], ampl_err[1]);

    }
    if( max > 0 && (unsigned int)max < sample.size() -1 )
      localMax.push_back(max);
  }


  if ( !localMax.size() )
    return result_;

  unsigned int Max( 1 << 31);
  for ( list<int>::iterator it = localMax.begin();
        it !=  localMax.end(); ++it ) {

    int use(0), max(*it);
    if ( ac_cut_min_ == -1 ||
         sample[max] - base[max&1] <= ac_cut_min_ ) {
      use = 1;
    } else {
      int pre(0), post(0);
      float sighalf =  (sample[max]-base[max&1])* 0.5;
      while ( max - pre >= 0 &&
              sample[max - pre]-base[(max - pre)&1]  > sighalf )
        pre++;
      while ( (unsigned int)max + post < sample.size() &&
              sample[max + post]-base[(max + post)&1]  > sighalf ) post++;
      if (  pre + post - 1 > 5  ) {
        use = 1;
      }
    }
    if ( use ) {


      if ( Max == (unsigned int)(1 << 31) )
        Max = *it;

      else if ( sample[Max] < sample[*it] )
        Max = *it;

    }


  }

  if ( Max == (unsigned int)(1 << 31) )
    return result_;


  // Prepare function
  mypulse_func_.SetParameters(Max , sample[Max] - base[Max&1], 0. );

  // Execute fit
  g.Fit(&mypulse_func_, "Q", "",
        max(0. , Max - fit_range_to_front_) ,
        min( (double)sample.size(), Max + fit_range_to_back_));

  TF1 *fpulse = g.GetFunction("f");
  double time = fpulse->GetParameter(0);
  double timeErr = fpulse->GetParError(0);
  double ampl = fpulse->GetParameter(1);
  double amplErr = fpulse->GetParError(1);
  double fbase =  fpulse->GetParameter(2);
  double fbaseErr = fpulse->GetParError(2);
  double chi2 = fpulse->GetChisquare();
  double ndf = fpulse->GetNDF();
  double prob = fpulse->GetProb();
  unsigned int run = CsEvent::Instance()->getRunNumber();
  unsigned int eventNb = CsEvent::Instance()->getEventNumberInRun();


  if (time < 9 || time > 22 ) return result_;

  if ( fit_cut_debug_ && fit_cut_debug_ <= fabs( ampl / (sample[Max] - (base[Max&1] + fbase) ) - 1 ) ) return result_;

/// \todo  TODO: Add checks to assert success of fit

  time *= sadc_clock_;
  timeErr *= sadc_clock_;

  time -= 27.8;


  //Analysis output to file

  if ( fit_ofile_ ) {
    static FILE *file = NULL;
    if (!file) {
      char* fname;
      asprintf( &fname, "%s.%f.%f",
                fit_ofile_, -fit_range_to_front_, fit_range_to_back_);
      cout << "Opening file '" << fname << "'\t";
      cout.flush();
      file = fopen( fname , "w");
      free( fname );
      cout << "Done"<< endl;
      assert( file );
    }

    double diff = ampl / (sample[Max] - (base[Max%2] + fbase) ) - 1;

    unsigned int size = sample.size();

    fwrite(&run, sizeof(unsigned int), 1, file);
    fwrite(&eventNb, sizeof(unsigned int), 1, file);
    fwrite(&id, sizeof(int), 1, file);
    fwrite(&plane_, sizeof(int), 1, file);
    fwrite(&x, sizeof(int), 1, file);
    fwrite(&y, sizeof(int), 1, file);
    fwrite(&diff, sizeof(double), 1 , file);
    fwrite(&ampl, sizeof(double), 1 , file);
    fwrite(&chi2, sizeof(double), 1 , file);
    fwrite(&ndf, sizeof(double), 1 , file);
    fwrite(&prob, sizeof(double), 1 , file);
    fwrite(&time, sizeof(double),1 , file);
    fwrite(base, sizeof(double), 2, file);
    fwrite(&fbase, sizeof(double), 1, file);
    fwrite(&size, sizeof(unsigned int), 1, file);
    fwrite(Sample, sizeof(unsigned short), size, file);

  }

  if( ampl && fit_cut_slope_ != 1000. &&  fit_cut_offset_ != 1000. ) {
    if( (chi2 / ndf / ampl * 200.) > ( fit_cut_slope_ * ampl + fit_cut_offset_ ) )
      return result_;
  }

  // Fill result structure


  result_.ampl = ampl;
  result_.time = time;
  result_.base = fbase+base[0];
  result_.base_diff = base[0] - base[1];
  result_.dampl = amplErr;
  result_.dtime = timeErr;
  result_.dbase = fbaseErr;
  result_.chi2 = chi2;
  result_.ndf = ndf;
  result_.ready = 1;

  return result_;

#else
  throw Reco::Exception("CsDigitizerSADC::FitPulse Requires ROOT-Version >= 5.16/00! \n");
  return result_;
#endif
};

// Destructor

CsDigitizerSADC::~CsDigitizerSADC  (void) {
  if ( fit_ofile_ )
    free(fit_ofile_);
};

// Clear results

void CsDigitizerSADC::ClearResult() {

  memset( &result_, 0, sizeof(result_t) );

};
