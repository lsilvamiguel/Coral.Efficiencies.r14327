/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CalorimeterDB.cc,v $
   $Date: 2011/03/01 14:14:00 $
   $Revision: 1.54 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

   Copyright(C): 1999-2000  V.Kolosov,A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

// --- Standard C/C++ library ---
#include <cmath>
#include <sstream>
#include <ctime>

// --- Internal files ----
#include "Exception.h"
#include "Calorimeter.h"
#include "Cell.h"
#include "DataBase.h"

using namespace std;

namespace Reco {

#define CDB_LINE_MAX 132

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SetOutputInfoName ( const string &tag )
{
  for( int i=0; i< (int)options.output_info_names.size(); i++ )
    if( options.output_info_names[i] == tag) return;

  options.output_info_names.push_back( tag );
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SetInputInfoName ( const string &tag )
{
  for( int i=0; i< (int)options.input_info_names.size(); i++ )
    if( options.input_info_names[i] == tag) return;

  options.input_info_names.push_back( tag );
}

////////////////////////////////////////////////////////////////////////////////

const std::vector< std::string >& Calorimeter::GetOutputInfoNames(void) const
{
  return options.output_info_names;
}

////////////////////////////////////////////////////////////////////////////////

const std::vector< std::string >& Calorimeter::GetInputInfoNames(void) const
{
  return options.input_info_names;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputAnyCalibInfo(const string &tag, string &s,const std::string &comment) const
{
  if( tag == "_CALIB" )
  {
    if( !XYRegularGrid() )
      return OutputCalibInfo(s);
    else
      return OutputCalibInfoXY(s);
  }
  else if( tag == "_LED" )
  {
    if( !XYRegularGrid() )
      return OutputLEDInfo(s,comment);
    else
      return OutputLEDInfoXY(s,comment);
  }
//   else if( tag == "_LEDinSpills" )
//   {
//     return OutputLEDInfoInSpills(s,comment);
//   }
//   else if( tag == "_LEDinSpillsXY" )
//   {
//     return OutputLEDInfoInSpillsXY(s,comment);
//   }
//   else if( tag == "_TimeInSpillLED" )
//   {
//     return OutputTimeInSpillLED(s);
//   }
  else if( tag == "_PED" )
  {
    if( !XYRegularGrid() )
      return OutputPEDInfo(s);
    else
      return OutputPEDInfoXY(s);
  }
  else if( tag == "_TimeCALIB" )
  {
    if( !XYRegularGrid() )
      return OutputTimeCalibInfo(s);
    else
      return OutputTimeCalibInfoXY(s);
  }
  else if( tag == "_NOISE" )
  {
    return OutputCellsInfo(s,NOISE,NEW);
  }
  else if( tag == "_EcutCells" )
  {
    if( !XYRegularGrid() )
      return OutputEcutCellsInfo(s);
    else
      return OutputEcutCellsInfoXY(s);
  }
  else if( tag == "_BadCells" )
  {
//     if( !XYRegularGrid() )
      return OutputBadCellsInfo(s);
//     else
//       return OutputEcutCellsInfoXY(s);
  }
  else if( tag == "_ProbEnT" )
  {
    return OutputProbEnT(s);
  }
  else if( tag == "_PositionInBurst" )
  {
    return OutputPositionInBurst( s );
  }
  else
  {
    throw Exception( (string(__func__) + " Bad tag: " + tag).c_str() );
  }
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputAnyCalibInfo(const string &tag,const string &s)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Try Calorimeter::InputAnyCalibInfo " << GetName() << " tag " << tag << endl;

  if( tag == "_CALIB" )
  {
    if( !XYInitialized() )
      return InputCalibInfo(OLD, s);
    else
      return InputCalibInfoXY(OLD, s);
  }
  else if( tag == "_CALIB_PRIM" )
  {
    if( !XYInitialized() )
      return InputCalibInfo(PRIM, s);
    else
      return InputCalibInfoXY(PRIM, s);
  }
  else if( tag == "_LED" )
  {
    if( !XYInitialized() )
      return InputLEDInfo(OLD, s);
    else
      return InputLEDInfoXY(OLD, s);
  }
  else if( tag == "_LEDinSpills" )
  {
    return InputLEDInfoInSpills(OLD, s);
  }
  else if( tag == "_LEDinSpillsXY" )
  {
    return InputLEDInfoInSpillsXY(OLD, s);
  }
//   else if( tag == "_TimeInSpillLED" )
//   {
//     return InputTimeInSpillLED(s);
//   }
  else if( tag == "_LED_PRIM" )
  {
    if( !XYInitialized() )
      return InputLEDInfo(PRIM, s);
    else
      return InputLEDInfoXY(PRIM, s);
  }
  else if( tag == "_TimeCALIB" )
  {
    if( !XYInitialized() )
      return InputTimeCalibInfo( OLD, s);
    else
      return InputTimeCalibInfoXY( OLD, s);
  }
  else if( tag == "_TimeCALIB_PRIM" )
  {
    if( !XYInitialized() )
      return InputTimeCalibInfo( PRIM, s);
    else
      return InputTimeCalibInfoXY( PRIM, s);
  }
  else if( tag == "_TrigGroupCALIB" )
  {
    cerr <<" WARNING!!!  Sorry!! No more trigger groups calibrations " << endl;
    return -1;
//      return InputTrGrCalibInfo(s);
  }
  else if( tag == "_TrigGroupTimeCALIB" )
  {
    cerr <<"  WARNING!!! Sorry!! No more trigger groups time calibrations " << endl;
    return -1;
//      return InputTrGrTimeCalibInfo(s);
  }
  else if( tag == "_NOISE" )
  {
    return InputCellsInfo(s,Reco::Calorimeter::NOISE,Reco::Calorimeter::NEW);
  }
  else if( tag == "_EcutCells" )
  {
    if( !XYInitialized() )
      return InputEcutCellsInfo(s);
    else
      return InputEcutCellsInfoXY(s);
  }
  else if( tag == "_BadCells" )
  {
    return InputBadCellsInfo(OLD, s);
//     if( !XYInitialized() )
//       return InputEcutCellsInfo(s);
//     else
//       return InputEcutCellsInfoXY(s);
  }
  else if( tag == "_MC_BadCells" )
  {
    return InputBadCellsInfo(MC, s);
//     if( !XYInitialized() )
//       return InputEcutCellsInfo(s);
//     else
//       return InputEcutCellsInfoXY(s);
  }
  else if( tag == "_ProbEnT" )
  {
    return InputProbEnT(s);
  }
  else if( tag == "_BaseProbNoise" )
  {
    cerr <<" WARNING!!! Still using BaseProbNoise ?? Seems it is obsolet please revise your code or options !!! " << endl;
    return -1;
//     return InputStatInfo(s,prob_ebase_old_);
  }
  else if( tag == "_PositionInBurst" )
  {
    return InputPositionInBurst( s );
  }
  else if( tag == "_Navigator" )
  {
    return InputPositionInTime( s );
  }
  else if( tag == "_NavigatorR" )
  {
    return InputPositionInTimeRedundant( s );
  }
  else if( tag == "_PositionInTime" )
  {
    return InputPositionInTimeGeneric( s );
  }
  else if( tag == "_JOUJOU" )
  {
    return InputJOUJOUInfo( s );
  }
  else
  {
    return -1;
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::GetInfoFromDataBase(const DataBase &db)
{
  for( int i=0; i<3; i++ )
  {
    DataBase::Element< vector<StatInfo> > b(CellInfoNames[i]);
    db.Read(b);
    cells_info[i][OLD] = cells_info[i][NEW] = *b.GetData();
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::PutInfoToDataBase(DataBase &db) const
{
  for( int i=0; i<3; i++ )
  {
    DataBase::Element< vector<StatInfo> > b(CellInfoNames[i],cells_info[i][NEW]);
    db.Write(b);
  }
}

////////////////////////////////////////////////////////////////////////////////
namespace {

class ForFileReading
{
  public:
    enum {N=500000};
    friend istream & operator >> (istream &,ForFileReading &c);
    char s[N];
};

istream & operator >> (istream &in,ForFileReading &c)
{
  in.read(c.s,ForFileReading::N);
  if( in.gcount()>=ForFileReading::N )
    throw Exception("CsCalorimeter::ReadCalib(): too small internal buffer!");
  c.s[in.gcount()]=0;
  in.clear();
//  cout << "Is it OK? " << in.good() << "\n";
  return in;
}

}

void Calorimeter::ReadFromDataBase(DataBase &db,tm &t)
{
   ForFileReading r;
   db.Read(string(GetName())+"_CALIB",r,t);
   string s(r.s);
//   cout << "String: " << s << "\n";
//#warning " Default Calibration setting should be implemented in case some chanels are missing "
   InputCalibInfo(OLD, s);

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::WriteToDataBase(DataBase &db,tm &time_start,tm &time_finish)
{
   string s("");
   int ok = OutputCalibInfo(s);
   db.Write(string(GetName())+"_CALIB",s,time_start,time_finish);
   s.clear();
   ok = OutputLEDInfo(s);
   db.Write(string(GetName())+"_LED"  ,s  ,time_start,time_finish);
   s.clear();
   ok = OutputPEDInfo(s);
   db.Write(string(GetName())+"_PED"  ,s  ,time_start,time_finish);
   return;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputCalibInfo(size_t when, const string &s )
{
  assert( when == OLD || when == PRIM );

//   bool debug = true;
  bool debug = false;
  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;
  int default_format_level = 1;

  if( XYRegularGrid() && default_format_level == 1 )
    return InputCalibInfoXY(when, s);

  //  Read line of comments
  getline(is,str);
  getline(is,str);
  if( debug ) cout << str << endl;
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);

  if( debug )
      cout << " Calorimeter " << calorim_name <<" with " << dummy  << endl;

  //  Read calibration
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
    {
      if( debug ) cout << str << endl;
      int icell;
      float w,c,sc;
      assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
      if( sscanf(str.c_str(), " %d %g %g %g %s", &icell, &c, &sc, &w, cellname) != 5 )
        throw Exception("FORMAT ERROR for InputCalibInfo in line %s", str.c_str());
      int jcell = FindCellByName( cellname );
      if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
      if( jcell < 0 || (int)NCells() <= jcell )
	{
	  stringstream ss;
	  ss << "Unexpected Cell name " << cellname << " in Calorimeter::InputCalibInfo "
	     << GetName() << " " << when << endl;
	  ss << "Input string: " << endl;
	  ss << str << endl;
	  ss << "It well might be that you use wrong calibration file for this detector "
	     << GetName() << endl;
	  throw Exception( ss.str().c_str() );
	}

      if( icell != jcell ) {
	icell = jcell;
	not_matched_cells_cnt++;
	if (not_matched_cells_cnt==1) {
	  cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	       << "cell id and cell name do not match.  (The calibration file "
	       << "was produced with a different geometry description, which "
	       << "might be an indication that you're using the wrong file.)"
	       << endl;
	}
      }

      if( debug )
	  cout << " icell " << icell <<" " << GetCellName(icell) << " c " << c/1000.
	       << " sc " << sc/1000. << " stat " << w << endl;

      if( icell >=0 && icell < (int)NCells() )
        SetCellInfo(CALIB, when,icell,StatInfo(w,c/1000.,sc/1000.));
    }

  return true;
}


////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ScaleCalibrationByCellCorrections( size_t when )
{
  for( size_t icell=0; icell != NCells(); icell++ )
  {
    double facell = individual_calib_factor_[icell];
    double stat = cells_info[CALIB][when][icell].GetEntries();
    double mean = cells_info[CALIB][when][icell].GetMean()*facell;
    double sigma = cells_info[CALIB][when][icell].GetSigma()*facell;
    SetCellInfo(CALIB, when,icell,StatInfo(stat,mean,sigma));
  }
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputCalibInfoXY(const string &s)
{
  return InputCalibInfoXY(OLD, s);
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputCalibInfoXY(size_t when, const string &s )
{
  assert( XYRegularGrid() );

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );

  getline(is,str);
  float energy_calib=0.;
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  ret = sscanf(str.c_str(), " %s %g", dummy, &energy_calib);
  assert( ret == 2 );
  if(GetName()=="EC00P1__") cout << " Calorimeter " << calorim_name << " energy_calib = " << energy_calib << endl;

//  Read 2 lines of comments
  getline(is,str);
  getline(is,str);

//  Read calibration
  while( getline(is,str) )
  {
    int x,y;
    float w,m,s,c,sc;
 // Need to check if line is a valid line in many places
    if ( sscanf(str.c_str(), " %d %d %g %g %g %g %g", &x, &y, &m, &s, &w, &c, &sc) != 7 )
      throw Exception("InputCalibInfoXY() %s bad line: %s", GetName().c_str(), str.c_str());

    int icell = GetCellOfColumnRow( x, y);
    if ( icell >= 0 )
      SetCellInfo(CALIB, when,icell,StatInfo(1,c/1000.,sc/1000.));
//      if(GetName() == "EC00P1__")cout << "My: ReadCalibration x= " << x << " y= " << y << " c= " << c << " sc= " << sc <<
//      " CALIB= " << CALIB << " icell= " << icell << endl;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputCalibInfo( string &s ) const
{
//   bool debug = true;
  bool debug = false;
//  if(GetName() != "EC00P1__")debug=false;
  if( debug ) cout << " Calorimeter::OutputCalibInfo for " <<
                   (this)->GetName() << " this = " << this << endl;

  char o[500000];
  o[0]=0;
  if( debug ) cout << " At start strlen(o) = " << strlen(o) << endl;

  sprintf(o,"### \n");
  sprintf(o+strlen(o),"Calibration Calorimeter:  %s  Ncells=  %zu \n",GetName().c_str(),NCells());
  sprintf(o+strlen(o),"### \n");

  for( size_t icell=0; icell!=NCells(); icell++ )
  {
     double entries = cells_info[CALIB][NEW][icell].GetEntries();
     double fac=cells_info[CALIB][NEW][icell].GetMean();
     double sfac=cells_info[CALIB][NEW][icell].GetSigma();

     bool tolerance_check = true;
     if(fac>0.) tolerance_check = (fabs( 1.-1./fac ) > options.tolerate2keepold_calib );

     if( entries > options.new_calib_stat_min && fac > 0 && tolerance_check )
     {
       if( debug )
       {
         cout << " Stat OK in cell " << GetCellName(icell) <<
	         " cells_info[CALIB][NEW][icell].GetMean() " << cells_info[CALIB][NEW][icell].GetMean() <<
	         " cells_info[CALIB][OLD][icell].GetMean() " << cells_info[CALIB][OLD][icell].GetMean() << endl;

       }


       double cfnew= cells_info[CALIB][OLD][icell].GetMean()/fac;
       double scfnew= sfac/fac * cfnew;

       if( debug )
       {
         printf(" %zu %7.2f %7.2f %7.2f %s \n",
                    icell,1000.*cfnew,1000.*scfnew,cells_info[CALIB][NEW][icell].GetEntries(), GetCellName(icell).c_str() );
       }
       sprintf(o+strlen(o)," %zu %7.2f %7.2f %7.2f %s \n",
                    icell,1000.*cfnew,1000.*scfnew,cells_info[CALIB][NEW][icell].GetEntries(), GetCellName(icell).c_str() );
     }
     else
     {

       if( debug )
       {
         printf(" %zu %7.2f %7.2f %7.2f %s \n",icell,
          1000.*cells_info[CALIB][OLD][icell].GetMean(),
          1000.*cells_info[CALIB][OLD][icell].GetSigma(),
            cells_info[CALIB][OLD][icell].GetEntries(), GetCellName(icell).c_str());
       }
       sprintf(o+strlen(o)," %zu %7.2f %7.2f %7.2f %s \n",icell,
          1000.*cells_info[CALIB][OLD][icell].GetMean(),
          1000.*cells_info[CALIB][OLD][icell].GetSigma(),
            cells_info[CALIB][OLD][icell].GetEntries(), GetCellName(icell).c_str());
     }
  }
  if( debug ) cout << " At finish strlen(o) = " << strlen(o) << endl;
  string ss(o);
//  s += ss;
  s = ss;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputCalibInfoXY( string& s ) const
{
  assert( XYRegularGrid() );

//   cout << " Calorimeter::OutputCalibInfoXY " << GetName() << " opt.tol2OLD CALIB " << options.tolerate2keepold_calib << endl;
  char o[500000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter Energy Calibration File \n",GetName().c_str());
  sprintf(o+strlen(o)," Calibration_Energy(GeV)    %7.2f \n",options.calibration_energy );
  sprintf(o+strlen(o)," X     Y      MeanPeak   SigmaPeak  Statistic     Coeff    SigmaCoeff \n");
  sprintf(o+strlen(o),"            (ADCchanels)                     (MeV/ADCchanel)          \n");
  for( int x=0; x != GetNColumns(); x++ )
    for( int y=0; y!= GetNRows(); y++ )
    {
      int icell = GetCellOfColumnRow(x,y);
      if( icell >= 0 )
      {
        double entries = cells_info[CALIB][NEW][icell].GetEntries();
        double fac=cells_info[CALIB][NEW][icell].GetMean();
        double sfac=cells_info[CALIB][NEW][icell].GetSigma();
        bool tolerance_check = true;
        if(fac>0.) tolerance_check = (fabs( 1.-1./fac ) > options.tolerate2keepold_calib );

//         if( entries > options.new_calib_stat_min && fac > 0  )
// 	  cout << " cell " << icell << " ftol " << fabs( 1.-1./fac) << " tol " << tolerance_check << endl;
       if( entries > options.new_calib_stat_min && fac > 0 && tolerance_check )
        {

          double cfnew= cells_info[CALIB][OLD][icell].GetMean()/fac;
          double scfnew= sfac/fac * cfnew;

          double peak_adc=0;
          double sigma_peak_adc=0;
          if(options.calibration_energy > 0)
          {
            if(cfnew > 0 )
            {
              peak_adc = options.calibration_energy/cfnew;
              sigma_peak_adc = scfnew/cfnew*peak_adc;
            }
          }
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f %11.2f %11.2f \n",x,y,
                                peak_adc,sigma_peak_adc,entries,1000.*cfnew,1000.*scfnew);
        }
        else
        {
          double old_entries=cells_info[CALIB][OLD][icell].GetEntries();
          double cfnew=cells_info[CALIB][OLD][icell].GetMean();
          double scfnew=cells_info[CALIB][OLD][icell].GetSigma();
          double peak_adc=0;
          double sigma_peak_adc=0;
          if(options.calibration_energy > 0)
          {
            if(cfnew > 0 )
            {
              peak_adc = options.calibration_energy/cfnew;
              sigma_peak_adc = options.calibration_energy*scfnew/cfnew;
            }
          }
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f %11.2f %11.2f \n",x,y,
                                peak_adc,sigma_peak_adc,old_entries,1000.*cfnew,1000.*scfnew);
        }
      }
      else
      {
        sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f %11.2f %11.2f \n",x,y,-1.,-1.,-1.,-1.,-1.);
      }
    }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputLEDInfo(size_t when, const string &s)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);
//  Read calibration
  unsigned int cells_cnt             = 0;
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell,normflag;
    float w,c,sc;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %d %s", &icell, &c, &sc, &w, &normflag, cellname);
    assert( ret == 6 );
    if(debug) printf(" LED for Cell %d Mean=%7.2f Sigma=%7.2f Stat=%7.2f \n",icell,c,sc,w);
    int jcell = FindCellByName( string(cellname) );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      stringstream ss;
      ss << " Unexpected Cell name " << cellname << " in Calorimeter::InputLEDInfo "
	 << " jcell = " <<    jcell    <<"   " <<  GetName() << " " << when  << endl;
      ss << " Expected Cellname " << GetCellName(icell) << endl;
      ss << " Input string: " << endl;
      ss << str << endl;
      ss << " It well might be that you use wrong calibration file for this detector "
	 << GetName() << endl;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if (icell >=0 && icell < (int)NCells() )
      SetCellInfo(LED, when, icell, StatInfo(w, c, sc) );
    else
      throw Exception( (string(__func__)+" "+GetName()+" Error wrong cell number=%i").c_str(), icell);

    cells_cnt++;
  }

  if ( cells_cnt != NCells() )
    throw Exception( (string(__func__)+" "+GetName()+" Mismatch in number of cells, "
                      "expected: %i, read: %i").c_str(), NCells(), cells_cnt);

  return 0;
}


////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputLEDInfoXY( size_t when, const string &s)
{
  assert( XYRegularGrid() );

  bool debug = false;
//   bool debug = true;
  if( when == OLD )
  {
    if( debug ) cout << " InputLEDInfoXY OLD " << GetName() << endl;
  }
  else if( when == PRIM )
  {
    if( debug ) cout << " InputLEDInfoXY PRIM " << GetName() << endl;
  }
  else
  {
    cerr << " Not expected InputLEDInfoXY " << when << " in " << GetName() << endl;
    exit(1);
  }

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
//  Read one more line of comments
  getline(is,str);


//  Read calibration
  unsigned int cells_cnt = 0;
  while( getline(is,str) )
  {
    int x,y;
    float w,m,s;
    ret = sscanf(str.c_str(), " %d %d %g %g %g", &x, &y, &m, &s, &w);
    assert( ret == 5 );
//    cout << " x=" << x << " y=" << y << " m=" << m << " s=" << s << " w=" << w << " c=" << c << endl;
    int icell = GetCellOfColumnRow( x, y );
    if ( icell >= 0 ) {
      SetCellInfo( LED, when, icell, StatInfo(w, m, s) );
      cells_cnt++;
    }
  }

  if ( cells_cnt != NCells() )
    throw Exception( (string(__func__)+" "+GetName()+" Mismatch in number of cells, "
                      "expected: %i, read: %i").c_str(), NCells(), cells_cnt);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputLEDInfoInSpills(size_t when, const string &s)
{
//   bool debug = true;
  bool debug = false;

  if( debug )
     cout << " Calorimeter::InputLEDInfoInSpills " << GetName() << " debug " << endl;

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);

//  Read calibration
  const string run_tag("run");
  const string spill_tag("spills");
  const string nan_tag("nan");
  unsigned int cells_cnt             = 0;
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell,normflag;
    float w,c,sc;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %d %s", &icell, &c, &sc, &w, &normflag, cellname);
    assert( ret == 6 );
    if (debug)
      printf(" LED for Cell %d Mean=%7.2f Sigma=%7.2f Stat=%7.2f Normflag %d \n",icell,c,sc,w,normflag);
    int jcell = FindCellByName( cellname );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      stringstream ss;
      ss << "Unexpected Cell name " << cellname << " in Calorimeter::InputLEDInfoInSpills "
	 << GetName() << " " << when << endl;
      ss << "Input string: " << endl;
      ss << str << endl;
      ss << "It well might be that you use wrong calibration file for this detector "
	 << GetName() << endl;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
      SetCellInfo( LED, when, icell, StatInfo(w, c, sc) );
    else
      throw Exception( (string(__func__)+" "+GetName()+" Error wrong cell number=%i").c_str(), icell);

    cells_cnt++;

    getline(is,str);
    char runtag[CDB_LINE_MAX],spilltag[CDB_LINE_MAX];
    int run, nspills;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %s %d %s %d", runtag, &run, spilltag, &nspills);
    assert( ret == 4 );

    string run_check(runtag);
    string spill_check(spilltag);
    if ( debug )
      cerr << "Check Format InputLEDInfoInSpills " << GetName()
	   << " runtag " << run_check << " spills " << spill_check << endl;

    if ( run_check != run_tag || spill_check != spill_tag )
    {
      stringstream ss;
      ss << "Format ERROR InputLEDInfoInSpills " << GetName() << " runtag " << run_check << " spills " << spill_check;
      throw Exception( ss.str().c_str() );
    }
    if( debug ) cerr << "Format OK InputLEDInfoInSpills " << GetName() << " runtag " << runtag << " spills " << nspills << endl;

    int nlines = nspills/20;
    int morespills =   nspills%20;
    if( debug ) cout << " nlines " << nlines << " morespills " << morespills << endl;

    if( debug ) cout << " Base stat " << w << " mean " << c << " sigma " << sc << endl;
    int ispill = 0;
    for( int iline = 0; iline < nlines; iline++ )
    {
      if( !getline(is,str) ) break;
      istringstream s(str.c_str());

      if( debug ) cout << " permils : ";
      string opt;
      for( int in = 0; in < 20; in++ )
      {
        ispill++;
        s >> opt;
        if( opt == nan_tag )
        {
//          if ( debug ) cout << " No measurements in spill " << ispill << " What to do? " << endl;
          SetInSpillMeanLED(icell, run, ispill, c);
        }
        else
        {
          int permils;
          ret = sscanf(opt.c_str(), " %d", &permils);
          assert( ret == 1 );
          double v =  c*(1. + double(permils)/1000.);
          SetInSpillMeanLED(icell, run, ispill, v);
          if ( debug )
            cout << " " << permils;
        }
//        if( debug ) cout << " " << opt;
      }
      if( debug ) cout << endl;
    }
    if( morespills > 0 )
    {
      if( !getline(is,str) ) break;
      istringstream s(str.c_str());

      if( debug ) cout << " permils : ";
      string opt;
      for( int in = 0; in < 20; in++ )
      {
        ispill++;
        s >> opt;
        if( opt == nan_tag )
        {
//          if ( debug ) cout << " No measurements in spill " << ispill << " What to do? " << endl;
          SetInSpillMeanLED(icell, run, ispill, c);
        }
        else
        {
          int permils;
          ret = sscanf(opt.c_str(), " %d", &permils);
          assert( ret == 1 );
          if( debug ) cout << " " << permils;
          double v =  c*(1. + double(permils)/1000.);
          SetInSpillMeanLED(icell, run, ispill, v);
          if ( debug )
            cout << " " << permils;
        }
//        if( debug ) cout << " " << opt;
      }
      if( debug ) cout << endl;
    }
  }

  if ( cells_cnt != NCells() )
    throw Exception( (string(__func__)+" "+GetName()+" Mismatch in number of cells, "
                      "expected: %i, read: %i").c_str(), NCells(), cells_cnt);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputLEDInfoInSpillsXY(size_t when, const string &s)
{
//   bool debug = true;
  bool debug = false;

  if( debug )
     cout << " Calorimeter::InputLEDInfoInSpillsXY " << GetName() << " debug " << endl;

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  assert( sscanf(str.c_str(), "%s %s", calorim_name, dummy) == 2 );
  getline(is,str);

//  Read calibration
  const string run_tag("run");
  const string spill_tag("spills");
  const string nan_tag("nan");
  unsigned int cells_cnt             = 0;
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell,normflag,x,y;
    float w,c,sc;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    assert( sscanf(str.c_str(), " %d %d %d %g %g %g %d %s", &x, &y, &icell, &c, &sc, &w, &normflag, cellname) == 8 );
    if (debug)
      printf(" LED for Cell %d Mean=%7.2f Sigma=%7.2f Stat=%7.2f Normflag %d \n",icell,c,sc,w,normflag);
    int jcell = FindCellByName( cellname );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      stringstream ss;
      ss << "Unexpected Cell name " << cellname << " in Calorimeter::InputLEDInfoInSpillsXY "
	 << GetName() << " " << when << endl;
      ss << "Input string: " << endl;
      ss << str << endl;
      ss << "It well might be that you use wrong calibration file for this detector "
	 << GetName() << endl;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
      SetCellInfo( LED, when, icell, StatInfo(w, c, sc) );
    else
      throw Exception( (string(__func__)+" "+GetName()+" Error wrong cell number=%i").c_str(), icell);

    cells_cnt++;

    getline(is,str);
    char runtag[CDB_LINE_MAX],spilltag[CDB_LINE_MAX];
    int run, nspills;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    assert( sscanf(str.c_str(), " %s %d %s %d", runtag, &run, spilltag, &nspills) == 4 );

    string run_check(runtag);
    string spill_check(spilltag);
    if ( debug )
      cerr << "Check Format InputLEDInfoInSpillsXY " << GetName()
	   << " runtag " << run_check << " spills " << spill_check << endl;

    if ( run_check != run_tag || spill_check != spill_tag )
    {
      stringstream ss;
      ss << "Format ERROR InputLEDInfoInSpillsXY " << GetName() << " runtag " << run_check << " spills " << spill_check;
      throw Exception( ss.str().c_str() );
    }
    if( debug ) cerr << "Format OK InputLEDInfoInSpillsXY " << GetName() << " runtag " << runtag << " spills " << nspills << endl;

    int nlines = nspills/20;
    int morespills =   nspills%20;
    if( debug ) cout << " nlines " << nlines << " morespills " << morespills << endl;

    if( debug ) cout << " Base stat " << w << " mean " << c << " sigma " << sc << endl;
    int ispill = 0;
    for( int iline = 0; iline < nlines; iline++ )
    {
      if( !getline(is,str) ) break;
      istringstream s(str.c_str());

      if( debug ) cout << " permils : ";
      string opt;
      for( int in = 0; in < 20; in++ )
      {
        ispill++;
        s >> opt;
        if( opt == nan_tag )
        {
//          if ( debug ) cout << " No measurements in spill " << ispill << " What to do? " << endl;
          SetInSpillMeanLED(icell, run, ispill, c);
        }
        else
        {
          int permils;
          assert( sscanf(opt.c_str(), " %d", &permils) == 1 );
          double v =  c*(1. + double(permils)/1000.);
          SetInSpillMeanLED(icell, run, ispill, v);
          if ( debug )
            cout << " " << permils;
        }
//        if( debug ) cout << " " << opt;
      }
      if( debug ) cout << endl;
    }
    if( morespills > 0 )
    {
      if( !getline(is,str) ) break;
      istringstream s(str.c_str());

      if( debug ) cout << " permils : ";
      string opt;
      for( int in = 0; in < 20; in++ )
      {
        ispill++;
        s >> opt;
        if( opt == nan_tag )
        {
//          if ( debug ) cout << " No measurements in spill " << ispill << " What to do? " << endl;
          SetInSpillMeanLED(icell, run, ispill, c);
        }
        else
        {
          int permils;
          assert( sscanf(opt.c_str(), " %d", &permils) == 1 );
          if( debug ) cout << " " << permils;
          double v =  c*(1. + double(permils)/1000.);
          SetInSpillMeanLED(icell, run, ispill, v);
          if ( debug )
            cout << " " << permils;
        }
//        if( debug ) cout << " " << opt;
      }
      if( debug ) cout << endl;
    }
  }

  if ( cells_cnt != NCells() )
    throw Exception( (string(__func__)+" "+GetName()+" Mismatch in number of cells, "
                      "expected: %i, read: %i").c_str(), NCells(), cells_cnt);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputLEDInfo( string &s, const string &comment ) const
{
  char o[500000];
  o[0]=0;
  if( comment == "") {
    sprintf(o,"### \n");
  } else {
    sprintf(o,"%s\n",comment.c_str() );
  }
  sprintf(o+strlen(o),"LEDs Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  sprintf(o+strlen(o),"  Cell      MeanPeak   SigmaPeak  Statistic  FitFlag  CellName \n");

  for( size_t icell=0; icell!=NCells(); icell++ )
  {
     int fitflag = fit_info[LED][icell].fit_ok;
       sprintf(o+strlen(o)," %zu %7.2f %7.2f %7.2f %5d   %s \n",icell,
             cells_info[LED][NEW][icell].GetMean(),
             cells_info[LED][NEW][icell].GetSigma(),
             cells_info[LED][NEW][icell].GetEntries(),
             fitflag,
	     GetCellName(icell).c_str() );
  }
  string ss(o);
//  s += ss;
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputLEDInfoXY( string& s, const string &comment ) const
{
  assert( XYRegularGrid() );

  char o[500000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter LEDs  %s\n",GetName().c_str(),comment.c_str());
  sprintf(o+strlen(o)," X     Y      MeanPeak   SigmaPeak  Statistic  FitFlag  CellName Cell \n");

  for( int x=0; x != GetNColumns(); x++ )
    for( int y=0; y != GetNRows(); y++ )
    {
      int icell = GetCellOfColumnRow(x,y);

      if( icell >= 0 )
      {
        int fitflag = fit_info[LED][icell].fit_ok;
        double entries = cells_info[LED][NEW][icell].GetEntries();
        double peaknew=cells_info[LED][NEW][icell].GetMean();
        double speaknew=cells_info[LED][NEW][icell].GetSigma();
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f %5d   %s  %d \n",x,y,
                                       peaknew,speaknew,entries,fitflag,
                                	     GetCellName(icell).c_str(),
                                             icell );
      }
      else
      {
        sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,-1.,-1.,-1.);
      }
    }

  string ss(o);
  s = ss;
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InputEdepCorr(const string &s) {

  const bool debug = false;

  unsigned nmods               = 0;  // number of cells read
  unsigned nmods_with_E_depend = 0;
  unsigned ndata               = 0;  // number of data points read
  unsigned int not_matched_cells_cnt = 0;

  istringstream is(s);
  string str;
  unsigned cell;     // cell id
  string   cellname; // cell name
  unsigned Ncoeffs;        // number of data points to be read

  while( is >> cell >> cellname >> Ncoeffs) {

    if ( cell >= NCells() )
      throw Exception( (string(__func__)+" "+GetName()+" invalid cell: %i").c_str(), cell);

    if (debug) cout<<cell<<" N = "<<Ncoeffs<<endl;

    int jcell = FindCellByName( cellname );
    if( !(0 <= jcell && jcell < (int)NCells()) )
    {
      stringstream ss;
      ss << "Unexpected Cell name " << cellname << " in " << __func__
	 << GetName() << endl
	 << "Input string: " << endl << str << endl
	 << "It well might be that you use wrong calibration file for this detector "
	 << GetName() << endl;
      throw Exception( ss.str().c_str() );
    }

    if( cell != (unsigned int)jcell ) {
      cell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if (Ncoeffs == 0) continue;

    nmods++;
    if (Ncoeffs > 1) nmods_with_E_depend++;

    if ( cells[cell].GetEdepCorr().NPoints() != 0 )
      throw Exception( (string(__func__)+" "+GetName()+" duplicate cell: %i").c_str(), cell);

    for (unsigned i=0; i < Ncoeffs; i++) {
      double e, a;

      if ( !(is>>e>>a) )
        throw Exception( (string(__func__)+" "+GetName()+" bad line: "+str).c_str() );

      if (debug) cout<<"    "<<e<<"  "<<a<<endl;

      cells[cell].GetEdepCorrNonConst().SetDataPoint(e, a);

      ndata++;
    }
  }

  if (debug) {
    cout << "Correction maps for " << nmods << "cells were read-in." << endl;
    cout << nmods_with_E_depend << " maps contain more than 1 data point." << endl;
    cout << ndata << " data points in total." << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InputTiSdepCorr(const string &s) {

  const bool debug = false;

  unsigned nmods                 = 0;  // number of cells read
  unsigned nmods_with_TiS_depend = 0;
  unsigned ndata                 = 0;  // number of data points read

  istringstream is(s);
  string str;
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is, str) ) {

    istringstream line(str);
    unsigned cell;     // cell id
    string   cellname; // cell name
    unsigned N;        // number of data points to be read

    if ( !( line >> cell >> cellname >> N) )
      throw Exception( (string(__func__)+" "+GetName()+" bad line start: "+str).c_str() );

    if ( cell >= NCells() )
      throw Exception( (string(__func__)+" "+GetName()+" invalid cell: %i").c_str(), cell);

    if (debug) cout<<cell<<" N = "<<N<<endl;

    int jcell = FindCellByName( cellname );
    if( !(0 <= jcell && jcell < (int)NCells()) )
    {
      stringstream ss;
      ss << "Unexpected Cell name " << cellname << " in " << __func__
	 << GetName() << endl
	 << "Input string: " << endl << str << endl
	 << "It well might be that you use wrong calibration file for this detector "
	 << GetName() << endl;
      throw Exception( ss.str().c_str() );
    }

    if( cell != (unsigned int)jcell ) {
      cell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if (N == 0) continue;

    nmods++;
    if (N > 1) nmods_with_TiS_depend++;

    if ( cells[cell].GetTiSdepCorr().NPoints() != 0 )
      throw Exception( (string(__func__)+" "+GetName()+" duplicate cell: %i").c_str(), cell);

    for (unsigned i=0; i < N; i++) {
      double e, a;

      if ( !(line>>e>>a) )
        throw Exception( (string(__func__)+" "+GetName()+" bad line: "+str).c_str() );

      if (debug) cout<<"    "<<e<<"  "<<a<<endl;

      cells[cell].GetTiSdepCorrNonConst().SetDataPoint(e, a);

      ndata++;
    }

    if ( line >> str )
      throw Exception( (string(__func__)+" "+GetName()
                        +" unexpected token at end of line: "+str).c_str() );
  }

  if (debug) {
    cout << "Correction maps for " << nmods << "cells were read-in." << endl;
    cout << nmods_with_TiS_depend << " maps contain more than 1 data point." << endl;
    cout << ndata << " data points in total." << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputPEDInfo(const string &s)
{
  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);

//  Read calibration
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell;
    float w,c,sc;

    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %s", &icell, &c, &sc, &w, cellname);
    assert( ret == 5 );
    int jcell = FindCellByName( cellname );
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cellname << " in Calorimeter::InputPEDInfo " <<
                                                        GetName() << " " << endl;
      cerr << " Input string: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      ss << " Unexpected Cell name " << GetName() << " Calorimeter::InputPEDInfo " << str;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
      SetCellInfo(PED,OLD,icell,StatInfo(w,c,sc));
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputPEDInfoXY(const string &s)
{
  assert( XYRegularGrid() );

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
//  Read one more line of comments
  getline(is,str);

//  Read calibration
  while( getline(is,str) )
  {
    int x,y;
    float w,m,s;
    ret = sscanf(str.c_str(), " %d %d %g %g %g", &x, &y, &m, &s, &w);
    assert( ret == 5 );
//    cout << " x=" << x << " y=" << y << " m=" << m << " s=" << s << " w=" << w << " c=" << c << endl;
    int icell = GetCellOfColumnRow( x, y);
    if ( icell >= 0 )
      SetCellInfo(PED,OLD,icell,StatInfo(w,m,s));
  }
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputPEDInfo( string& s ) const
{
  char o[500000];
  o[0]=0;
  sprintf(o,"### \n");
  sprintf(o+strlen(o)," Pedestals Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  sprintf(o+strlen(o),"### \n");

  for( size_t icell=0; icell!=NCells(); icell++ )
  {
     double entries = cells_info[PED][NEW][icell].GetEntries();
     if( entries > 0 )
     {
       sprintf(o+strlen(o)," %zu %7.2f %7.2f %7.2f %s \n",icell,
             cells_info[PED][NEW][icell].GetMean(),
             cells_info[PED][NEW][icell].GetSigma(),
             cells_info[PED][NEW][icell].GetEntries(),
	     GetCellName(icell).c_str() );
     }
     else
     {
       sprintf(o+strlen(o)," %zu %7.2f %7.2f %7.2f %s \n",icell,
             cells_info[PED][OLD][icell].GetMean(),
             cells_info[PED][OLD][icell].GetSigma(),
             cells_info[PED][OLD][icell].GetEntries(),
	     GetCellName(icell).c_str() );
     }
  }
  string ss(o);
//  s += ss;
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputPEDInfoXY( string& s ) const
{
  assert( XYRegularGrid() );

  char o[500000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter PEDS for SADC readout \n",GetName().c_str());
  sprintf(o+strlen(o)," X     Y      MeanPeak   SigmaPeak  Statistic \n");

  for( int x=0; x != GetNColumns(); x++ )
    for( int y=0; y != GetNRows(); y++ )
    {
      int icell = GetCellOfColumnRow(x,y);
      if( icell >= 0 )
      {
        double entries = cells_info[PED][NEW][icell].GetEntries();
        double peaknew=cells_info[PED][NEW][icell].GetMean();
        double speaknew=cells_info[PED][NEW][icell].GetSigma();
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,
                                peaknew,speaknew,entries);
      }
      else
      {
        sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,-1.,-1.,-1.);
      }
    }

  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputTimeCalibInfo(size_t when, const string &s)
{
//   bool debug = true;
  bool debug = false;
//  if( GetName() == "EC01P1__" ) debug = true;

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );

//  Read 2 lines of comments
  getline(is,str);
  getline(is,str);
//  Read 1 line of comment and time offset
  getline(is,str);
  getline(is,str);
  float t_offset;
  ret = sscanf(str.c_str(), " %g", &t_offset);
  assert( ret == 1 );
  if( debug ) cout << " Time offset " << t_offset << endl;
//  Read calibration
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell;
    float t,st,w;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %s", &icell, &w, &t, &st, cellname);
    assert( ret == 5 );
    int jcell = FindCellByName( cellname );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cellname << " in Calorimeter::InputTimeCalibInfo " <<
                                                        GetName() << " " << when  << endl;
      cerr << " Input string: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      ss << " Unexpected Cell name " << GetName() << "Calorimeter::InputTimeCalibInfo " << str;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
    {
      if( debug )
      {
        cout << "  cell " << icell << " w " << w << " t0 " << t+t_offset <<" st " << st << endl;
      }

      if( w == 0. ) w=1.;
      if( debug )
      {
        cout << "  cell " << icell << " t0 " << t+t_offset << endl;
      }
      SetCellInfo(TIME, when,icell,StatInfo((double)w,(double)(t+t_offset),(double)st));
      if( debug )
      {
        double t0=cells_info[TIME][when][icell].GetMean();
        cout << GetName() <<" InputTimeCalibInfo cell " << icell << " t0 " << t0 << endl;
      }
    }
  }

  if( debug )
  {
    cout << " Calorimeter::InputTimeCalibInfo stop for debug " << GetName() << endl;
    exit(0);
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputTimeCalibInfoXY( size_t when, const string &s)
{
  assert( XYRegularGrid() );

//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::InputTimeCalibInfoXY " << GetName() << " debug " << endl;
  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],offset_id[CDB_LINE_MAX],try_to_read_offset[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );

//  Read 2 lines of comments
  getline(is,str);
  float t_offset = 0.;
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  ret = sscanf(str.c_str(), "%s %s", offset_id, try_to_read_offset);
  assert( ret == 2 );
  if( std::string(offset_id) == "Time_Offset" )
  {
    if( debug )
    {
      cout << " Time_Offset = " << try_to_read_offset << " was identified " << endl;
      cout <<" Try to recognize once more " << str << endl;
    }
    ret = sscanf(str.c_str(), " %s %g", offset_id, &t_offset);
    assert( ret == 2 );
    if( debug ) cout << " And indeed we read Time_Offset = " << t_offset << " Ura !!! " << endl;
  }

//  Read calibration
  while( getline(is,str) )
  {
    int x,y;
    float t,st,w;
    ret = sscanf(str.c_str()," %d %d %g %g %g", &x, &y, &t, &st, &w);
    assert( ret == 5 );
    t += t_offset;
    if( w <= 0.) w = 1.;
    int icell = GetCellOfColumnRow( x, y);
    if ( icell >= 0 )
      SetCellInfo(TIME, when,icell,StatInfo(w,t,st));
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputTimeCalibInfo( string &s ) const
{
  char o[500000];
  o[0]=0;
  sprintf(o," %s Calorimeter Time Calibration File \n",GetName().c_str());
  sprintf(o+strlen(o)," Cell#    Statistic   Time     SigmaTime  CellName    \n");
  sprintf(o+strlen(o),"                      (ns)       (ns)                 \n");

  float t_offset=0.;
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    double entries = cells_info[TIME][NEW][icell].GetEntries();
    if( entries > 10 )
    {
      double t=cells_info[TIME][NEW][icell].GetMean()+
                      cells_info[TIME][OLD][icell].GetMean();
      t_offset += t;
    }
    else
    {
//      double old_entries=cells_info[TIME][OLD][icell].GetEntries();
      double t=cells_info[TIME][OLD][icell].GetMean();
      t_offset += t;
    }
  }

  t_offset = t_offset/float(NCells());
  sprintf(o+strlen(o)," Time  Offset  \n");
  sprintf(o+strlen(o)," %11.2f \n",t_offset );

  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    double entries = cells_info[TIME][NEW][icell].GetEntries();
    if( entries > 10 )
    {
      double t=cells_info[TIME][NEW][icell].GetMean()+
                      cells_info[TIME][OLD][icell].GetMean();
      double st=cells_info[TIME][NEW][icell].GetSigma();
      sprintf(o+strlen(o)," %4zu %11.2f %11.2f %11.2f %s \n",
                           icell,(float)entries,(float)t-t_offset,(float)st, GetCellName(icell).c_str() );
    }
    else
    {
//      double old_entries=cells_info[TIME][OLD][icell].GetEntries();
      double t=cells_info[TIME][OLD][icell].GetMean();
      double st=cells_info[TIME][OLD][icell].GetSigma();
      sprintf(o+strlen(o)," %4zu %11.2f %11.2f %11.2f %s \n",
                         icell,1.,(float)t-t_offset,(float)st, GetCellName(icell).c_str());
    }
  }
  string ss(o);
//  s += ss;
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputTimeCalibInfoXY( string& s ) const
{
  assert( XYRegularGrid() );

  char o[500000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter Time Calibration File \n",GetName().c_str());
  sprintf(o+strlen(o),"   X     Y        Time     SigmaTime  Statistic    \n");
//  sprintf(o+strlen(o),"                  (ns)       (ns)                  \n");
  float t_offset=0.;
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    double entries = cells_info[TIME][NEW][icell].GetEntries();
    if( entries > 10 )
    {
      double t=cells_info[TIME][NEW][icell].GetMean()+
                      cells_info[TIME][OLD][icell].GetMean();
      t_offset += t;
    }
    else
    {
//      double old_entries=cells_info[TIME][OLD][icell].GetEntries();
      double t=cells_info[TIME][OLD][icell].GetMean();
      t_offset += t;
    }
  }

  t_offset = t_offset/float(NCells());
  sprintf(o+strlen(o)," Time_Offset %11.2f  \n",t_offset);

  for( int x=0; x != GetNColumns(); x++ )
    for( int y=0; y != GetNRows(); y++ )
    {
      int icell = GetCellOfColumnRow(x,y);
      if( icell >= 0 )
      {
        double entries = cells_info[TIME][NEW][icell].GetEntries();
        if( entries > 10 )
        {
          double t=cells_info[TIME][NEW][icell].GetMean()+
                          cells_info[TIME][OLD][icell].GetMean();
          double st=cells_info[TIME][NEW][icell].GetSigma();
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,(float)t-t_offset,(float)st,entries);
        }
        else
        {
//          double old_entries=cells_info[TIME][OLD][icell].GetEntries();
          double t=cells_info[TIME][OLD][icell].GetMean();
          double st=cells_info[TIME][OLD][icell].GetSigma();
          sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,(float)t-t_offset,(float)st,1.);
        }
      }
      else
      {
        sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,-1.,-1.,-1.);
      }
    }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputBadCellsInfo(size_t when, const string &s)
{
//   bool debug = true;
  bool debug = false;
  if( debug )
  {
    cout << " Calorimeter::InputBadCellsInfo " << GetName() << " debug " << endl;
    cout << s;
  }
  int nbadcells=0;
  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;
  double scale = options.scale_cells_threshols;
  if( scale != 1. )
  {
    cerr << " WARNING! Cells threshols for " << GetName() << " Scaled by a factor of " << scale << endl;
  }

//  Read line of comments
// ??? Gams-4pi format does not suppose this line of comments
// COMPASS format does, commenting the next line will crash 2008 production (and presumably also 2004)
  getline(is,str);

  // workaround for the difference to Gams-4pi format, only read the second time, if the first line is indeed a comment
  if (str[0]=='#')
    getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);

  getline(is,str);
  float particle_energy_threshold=0.;
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  ret = sscanf(str.c_str(), " %s %g", dummy, &particle_energy_threshold);
  assert( ret == 2 );
  if(fabs(options.particle_energy_threshold-particle_energy_threshold) >0.01)
  {
    cout << " WARNING! In Calorimeter " << calorim_name <<
            " InputBadCellsInfo:: particle_energy_threshold has been changed " << endl;
  }

  for( size_t i=0; i<NCells(); i++ )
  {
    cell_is_bad_[when][i]=false;
    status_bad_cells_[when][i]=0;
  }
  bad_cells_[when].clear();
//  Read Bad Cells Info
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell,status;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %d %s", &icell, &status, cellname);
    assert( ret == 3 );
    int jcell = FindCellByName( cellname );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cellname << " in Calorimeter::InputBadCellsInfo " <<
                                                      GetName() << " " << when  << endl;
      cerr << " Input string: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      ss << " Unexpected Cell name " << GetName() << "Calorimeter::InputBadCellsInfo " << str;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
    {
      bad_cells_[when].push_back(icell);
      status_bad_cells_[when][icell]=status;
      cell_is_bad_[when][icell]=true;
      nbadcells++;
      if( debug )
      {
        cout << " Set cell " << icell << " bad status " << status << " Cell " << GetCellName(icell) << endl;
      }
    }
  }
  if( debug ) cout << " Calorimeter::InputBadCellsInfo " << GetName() << " debug  OK !!!! Nbadcells = " << nbadcells << endl;
  if( !s.empty() &&  nbadcells == 0 )
  {
    cerr << " Suspicious! Calorimeter::InputBadCellsInfo No bad cell in not empty input!!! " << endl;
    exit(1);
  }


  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputBadCellsInfo( string &s ) const
{
  char o[500000];
  o[0]=0;

  size_t nbad_cells = bad_cells_[NEW].size();
  sprintf(o,"### \n");
  sprintf(o+strlen(o)," %s Bad Cells List  Ncells=%zu  NBadCells=%zu \n",GetName().c_str(),NCells(),nbad_cells);
  sprintf(o+strlen(o),"### \n");
  sprintf(o+strlen(o)," Min_Gamma_Energy(MeV)    %7.2f \n",options.particle_energy_threshold );

  for( list<size_t>::const_iterator it=bad_cells_[NEW].begin();it!=bad_cells_[NEW].end(); it++ )
  {
    sprintf(o+strlen(o)," %6zu %4d %s \n",*it, status_bad_cells_[NEW][*it],GetCellName(*it).c_str());
  }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputEcutCellsInfo(const string &s)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::InputEcutCellsInfo " << GetName() << endl;

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX],cellname[CDB_LINE_MAX];
  string str;
  double scale = options.scale_cells_threshols;
  if( scale != 1. )
  {
    cerr << " WARNING! Cells threshols for " << GetName() << " Scaled by a fctor of " << scale << endl;
  }

//  Read line of comments
  getline(is,str);
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);

  getline(is,str);
  float particle_energy_threshold=0.;
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  ret = sscanf(str.c_str(), " %s %g", dummy, &particle_energy_threshold);
  assert( ret == 2 );
  if(fabs(options.particle_energy_threshold-particle_energy_threshold) >0.01)
  {
    cout << " WARNING! In Calorimeter " << calorim_name <<
            " InputEcutCellsInfo:: particle_energy_threshold has been changed " << endl;
  }

  if( debug ) cout << " Calorimeter::InputEcutCellsInfo Clear  " << endl;
  for( size_t i=0; i<NCells(); i++ )
  {
    energy_cut_bad_cells_old[i]=0;
  }
//  Read Bad Cells Info
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    if( debug ) cout << str << endl;
    int icell;
    float eadccut,eadcgamcut;
    float ecut,egamcut;
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %g %s",
                 &icell, &eadccut, &eadcgamcut, &ecut, &egamcut, cellname);
    assert( ret == 6 );
    int jcell = FindCellByName( cellname );
    if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cellname << " in Calorimeter::InputEcutCellsInfo " <<
                                                        GetName() << " " << endl;
      cerr << " Input string: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      ss << " Unexpected Cell name " << GetName() << "Calorimeter::InputEcutCellsInfo " << str;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
    {
      double cf = cells_info[CALIB][OLD][icell].GetMean();
      energy_cut_bad_cells_old[icell]=eadccut*cf*scale;
      energy_gamma_cut_bad_cells_old[icell]=eadcgamcut*cf*scale;

      if(energy_cut_bad_cells_old[icell]==0) energy_cut_bad_cells_old[icell]=0.005;
      if(energy_gamma_cut_bad_cells_old[icell]==0)
          energy_gamma_cut_bad_cells_old[icell]=options.particle_energy_threshold;
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputEcutCellsInfo( string &s ) const
{
  char o[500000];
  o[0]=0;
  sprintf(o,"### \n");
  sprintf(o+strlen(o)," %s Energy threshods  Ncells=%zu \n",GetName().c_str(),NCells());
  sprintf(o+strlen(o),"### \n");
  sprintf(o+strlen(o)," Min_Gamma_Energy(MeV)    %7.2f \n",options.particle_energy_threshold );

  for( size_t i=0; i<NCells(); i++ )
  {
    if( energy_cut_bad_cells_new[i] > 0 )
    {
      double cf = cells_info[CALIB][OLD][i].GetMean();
      sprintf(o+strlen(o)," %6zu %7.2f %7.2f %7.2f %7.2f %s \n",i,
                energy_cut_bad_cells_new[i]/cf, energy_gamma_cut_bad_cells_new[i]/cf,
                1000*energy_cut_bad_cells_new[i],1000*energy_gamma_cut_bad_cells_new[i],
                                                    GetCellName(i).c_str());
    }
  }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputEcutCellsInfoXY(const string &s)
{
  assert( XYRegularGrid() );

  bool debug = false;
  double scale = options.scale_cells_threshols;
  if( scale != 1. )
  {
    cerr << " WARNING! Cells threshols for " << GetName() << " Scaled by a fctor of " << scale << endl;
  }

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
//  Read one more line of comments
  getline(is,str);


//  Read calibration
 for( size_t i=0; i<NCells(); i++ )
  {
    energy_cut_bad_cells_old[i]=0;
    energy_gamma_cut_bad_cells_old[i]=0;
  }
  bad_cells_[OLD].clear();

//  Read Bad Cells Info
  while( getline(is,str) )
  {
    int x,y;
    float eadccut,eadcgamcut;
    float ecut,egamcut;
    ret = sscanf(str.c_str(), " %d %d %g %g %g %g",
                 &x, &y, &eadccut, &eadcgamcut, &ecut, &egamcut);
    assert( ret == 6 );
    int icell = GetCellOfColumnRow( x, y );
    if(icell >=0 && icell < (int)NCells() )
    {
      double cf = cells_info[CALIB][OLD][icell].GetMean();
      energy_cut_bad_cells_old[icell]=eadccut*cf*scale;
      energy_gamma_cut_bad_cells_old[icell]=eadcgamcut*cf*scale;
      if( debug ) cout << GetCellName(icell) << " EcellTh " << energy_cut_bad_cells_old[icell] <<
                                    " EgammaTh " << energy_cut_bad_cells_old[icell] << endl;
//  Should be the same if the calibration didn't change
      if(energy_cut_bad_cells_old[icell]==0) energy_cut_bad_cells_old[icell]=0.005;
      if(energy_gamma_cut_bad_cells_old[icell]==0)
          energy_gamma_cut_bad_cells_old[icell]=options.particle_energy_threshold;
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputEcutCellsInfoXY( string& s ) const
{
  assert( XYRegularGrid() );

  char o[500000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter Bad Cells List   \n",GetName().c_str());
  sprintf(o+strlen(o)," LowThrProb  HighThrProb  \n %9.6f %9.6f \n ",options.prob_cell_noise, options.prob_gamma_noise);
//   sprintf(o+strlen(o)," X     Y     ElowCUT(ADC) EhighCUT(ADC) ElowCUT(MeV) EhighCUT(MeV) Status \n");
  sprintf(o+strlen(o)," X     Y     ElowCUT(ADC) EhighCUT(ADC) ElowCUT(MeV) EhighCUT(MeV) \n");

  for( int x=0; x != GetNColumns(); x++ )
    for( int y=0; y != GetNRows(); y++ )
    {
      int icell = GetCellOfColumnRow(x,y);
      if( icell >= 0 )
      {
          double cf = cells_info[CALIB][OLD][icell].GetMean();
          sprintf(o+strlen(o)," %4d %4d  %9.2f %9.2f %9.2f %9.2f \n",x,y,
                   energy_cut_bad_cells_new[icell]/cf,
                   energy_gamma_cut_bad_cells_new[icell]/cf,
                   1000*energy_cut_bad_cells_new[icell],
                   1000*energy_gamma_cut_bad_cells_new[icell]);
      }
      else
      {
        sprintf(o+strlen(o)," %4d %4d  %9.2f %9.2f \n",x,y,-1.,-1.);
      }
    }

  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputProbEnT(const string &s)
{
//   bool debug = true;
  bool debug = false;

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX], cellname[CDB_LINE_MAX];
  string str;

  //  Read line of comments
  if( debug ) cout << " Calorimeter::InputProbEnT " << GetName() << endl;
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s", calorim_name);
  assert( ret == 1 );
  //  Read 3 more line of comments
  if( debug ) cout << " skip line " << endl;
  getline(is,str);
  if( debug ) cout << " skip line " << endl;
  getline(is,str);
  if( debug ) cout << " skip line " << endl;
  getline(is,str);

  if( debug ) cout << " Now prepare to read and decode lines " << endl;
  bool even = false;
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    if( debug ) cout << " Got new line ";
    int icell;
    float adc[11],ecut[11];
    if( even )        // Read even lines
    {
      if( debug ) cout << " Read even line " << endl;
      even = false;
      assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
      if ( sscanf(str.c_str(), " %d %g %g %g %g %g %g %g %g %g %g %g %s", &icell,
		  &ecut[0], &ecut[1], &ecut[2], &ecut[3], &ecut[4], &ecut[5],
		  &ecut[6], &ecut[7], &ecut[8], &ecut[9], &ecut[10], cellname) != 13 )
	throw Exception("InputProbEnT() %s bad line: %s", GetName().c_str(), str.c_str());

      int jcell = FindCellByName( cellname );
      if(debug) cout << " Cell name " << cellname << " jcell " << jcell << endl;
      if( !(jcell >=0 && jcell < (int)NCells()) )
      {
        cerr << " Unexpected Cell name " << cellname << " in Calorimeter::ProbEnT " <<
                                                        GetName() << " " <<  endl;
        cerr << " Input string: " << endl;
        cerr << str << endl;
        cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
        stringstream ss;
        ss << " Unexpected Cell name " << GetName() << "Calorimeter::ProbEnT " << str;
        throw Exception( ss.str().c_str() );
      }

      if( icell != jcell ) {
	icell = jcell;
	not_matched_cells_cnt++;
	if (not_matched_cells_cnt==1) {
	  cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	       << "cell id and cell name do not match.  (The calibration file "
	       << "was produced with a different geometry description, which "
	       << "might be an indication that you're using the wrong file.)"
	       << endl;
	}
      }

      if(icell >=0 && icell < (int)NCells() )
      {
        double cf = cells_info[CALIB][OLD][icell].GetMean();
        if( adc[0] > 0.)
        {
          energy_cut_[0][icell] = adc[0]*cf;
        }
        else
        {
          energy_cut_[0][icell] = 0.;
        }
        for( int ip=1; ip<11; ip++ )
        {
          energy_cut_[ip][icell] = adc[ip]*cf;
	  if( energy_cut_[ip][icell] < energy_cut_[ip-1][icell] )
          {
            energy_cut_[ip][icell] = energy_cut_[ip-1][icell];
          }
        }
        if( debug )
	  cout << " cell " << icell << " cf " << cf << " energy_cut "
	       << energy_cut_[0][icell] << " " << energy_cut_[1][icell] << " "
	       << energy_cut_[2][icell] << " " << energy_cut_[3][icell] << " "
	       << energy_cut_[4][icell] << " " << energy_cut_[5][icell] << " "
	       << energy_cut_[6][icell] << " " << energy_cut_[7][icell] << " "
	       << energy_cut_[8][icell] << " " << endl;
      }
      continue;
    }
    else
    {
      if( debug ) cout << " and try to decode it " << endl;
    }
    even = true;
    ret = sscanf(str.c_str(), " %d %g %g %g %g %g %g %g %g %g %g %g", &icell,
                 &adc[0], &adc[1], &adc[2], &adc[3], &adc[4], &adc[5],
                 &adc[6], &adc[7], &adc[8], &adc[9], &adc[10]);
    assert( ret == 12 );
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputProbEnT( string& s ) const
{
//  bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter:: " << GetName() << " OutputProbEnT " << endl;
  char o[1000000];
  o[0] = 0;
  sprintf(o," %s COMPASS Calorimeter Probability Energy Thresholds   \n",GetName().c_str());
  sprintf(o+strlen(o)," comment: two line for each cell in first line thresholds in ADC in second in MeV   \n");
  sprintf(o+strlen(o),"                         -1         -2       -2         -3       -3         -4       -4         -5       -5         -6 \n");
  sprintf(o+strlen(o),"  icell     Emin       10       5*10       10       5*10       10       5*10       10       5*10       10       5*10   \n");
  if( debug )
  {
    printf(" %s COMPASS Calorimeter Probability Energy Thresholds   \n",GetName().c_str());
    printf(" comment: two line for each cell in first line thresholds in ADC in second in MeV   \n ");
    printf("                    -1      -2    -2      -3    -3      -4    -4      -5    -5     -6 \n ");
    printf("  icell   Emin    10    5*10    10    5*10    10    5*10    10    5*10    10    5*10  \n ");
  }
  for( int i=0; i < (int)NCells(); i++ )
  {
    if( debug ) cout << " New cell " << i << GetCellName(i) << endl;
    double cf = cells_info[CALIB][OLD][i].GetMean();
    double adc[11];
    double ecut[11];
    for( int ip=0; ip<11; ip++ )
    {
      if( energy_cut_[ip][i] > 0. )
      {
        adc[ip]= energy_cut_[ip][i]/cf;
        ecut[ip] =energy_cut_[ip][i]*1000.;
      }
      else
      {
        adc[ip]= energy_cut_[ip][i];
        ecut[ip] =energy_cut_[ip][i];
      }
    }

    if( debug )
    {
      printf(" %5d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n ",i,
                                adc[0],adc[1],adc[2],adc[3],adc[4],adc[5],adc[6],adc[7],adc[8],adc[9],adc[10]);

    }
    sprintf(o+strlen(o)," %5d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f \n",i,
                               adc[0],adc[1],adc[2],adc[3],adc[4],adc[5],adc[6],adc[7],adc[8],adc[9],adc[10]);
    if( debug )
    {
      printf("       %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %s\n",
          ecut[0],ecut[1],ecut[2],ecut[3],ecut[4],ecut[5],ecut[6],ecut[7],ecut[8],ecut[9],ecut[10],GetCellName(i).c_str());

    }
    sprintf(o+strlen(o),"       %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %s\n",
          ecut[0],ecut[1],ecut[2],ecut[3],ecut[4],ecut[5],ecut[6],ecut[7],ecut[8],ecut[9],ecut[10],GetCellName(i).c_str());

  }

  if( debug ) cout << " Calorimeter:: " << GetName() << " OutputProbEnT seems OK " << endl;
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputStatInfo(string& s, size_t what,size_t when) const
{
//#warning " implement size check"
  bool debug = false;

  if( cells_stat_info[what][when] == NULL )
  {
    cerr << " ERROR Calorimeter::OutputStatInfo what=" << what << " when=" <<
                                         when << " not initialized " << endl;
    return -1;
  }

  if( debug ) cout << " Calorimeter::OutputStatInfo What=" << what << " When=" << when << endl;
  bool recognised_info = true;
  char o[500000];
  o[0]=0;
  sprintf(o,"### \n");
  if( what == CALIB )
    sprintf(o+strlen(o),"Monitoring of CALIB Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == LED )
    sprintf(o+strlen(o),"Monitoring of LED Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == PED )
    sprintf(o+strlen(o),"Monitoring of PED Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == TIME )
    sprintf(o+strlen(o),"Monitoring of TIME Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == NOISE )
    sprintf(o+strlen(o),"Monitoring of NOISE Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  {
    sprintf(o+strlen(o),"Monitoring of UNKNOWN_INFO Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
    recognised_info = false;
  }

  if( recognised_info )
  {
    for( size_t icell=0; icell!=NCells(); icell++ )
    {
      double entries = cells_stat_info[what][when][icell].GetEntries();
      double cfnew=cells_stat_info[what][when][icell].GetMean();
      double scfnew=cells_stat_info[what][when][icell].GetSigma();
      sprintf(o+strlen(o)," %zu %12.6f %12.6f %12.6f %s \n",icell,cfnew,scfnew,entries,GetCellName(icell).c_str());
    }
  }

  string ss(o);
//  s += ss;
  s = ss;
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputStatInfo(const string &s,size_t what,size_t when)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell;
    float mean,sigma,stat;
    char cell_name[CDB_LINE_MAX];
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %s", &icell, &mean, &sigma, &stat, cell_name);
    assert( ret == 5 );
    int jcell = FindCellByName( cell_name );
    if(debug) cout << " Cell name " << cell_name << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cell_name << " in Calorimeter::InputStatInfo " <<
                                                        GetName() << " " << when  << endl;
      cerr << " Input string: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      ss << " Unexpected Cell name " << GetName() << "Calorimeter::InputStatInfo " << str;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
    {
      if(debug) printf(" Cell %d Mean=%g Sigma=%g Stat=%g cell_name=%s \n",icell,mean,sigma,stat,cell_name);
      SetCellInfo(what,when,icell,StatInfo(stat,mean,sigma));
      if(debug) cout << " check mean " << cells_info[what][when][icell].GetMean() <<
              " check sigma " << cells_info[what][when][icell].GetSigma() << endl;
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputStatInfo(string& s,const std::vector<StatInfo> &stat_info) const
{
  char o[500000];
  o[0]=0;
  sprintf(o,"### \n");
  sprintf(o+strlen(o)," Stat Info Calorimeter: %s \n",GetName().c_str());
  sprintf(o+strlen(o)," %zu  vector size \n",stat_info.size());

  for( size_t i=0; i!=stat_info.size(); i++ )
  {
    double entries = stat_info[i].GetEntries();
    double cfnew=stat_info[i].GetMean();
    double scfnew=stat_info[i].GetSigma();
    sprintf(o+strlen(o)," %zu %12.6f %12.6f %12.2f \n",i,cfnew,scfnew,entries);
  }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputStatInfo(const string &s,std::vector<StatInfo> &stat_info)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s);
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  getline(is,str);
  int vector_size;
  int ret = sscanf(str.c_str(), " %d", &vector_size);
  assert( ret == 1 );
  stat_info.clear();

  while( getline(is,str) )
  {
    int i;
    float mean,sigma,stat;
    ret = sscanf(str.c_str(), " %d %g %g %g", &i, &mean, &sigma, &stat);
    assert( ret == 4 );
    if(i>=0 )
    {
      if(debug) printf(" Index %d Mean=%g Sigma=%g Stat=%g \n",i,mean,sigma,stat);
      stat_info.push_back( StatInfo(stat,mean,sigma) );
    }
  }

  if( (int)stat_info.size() != vector_size )
  {
    cerr << " Calorimeter::InputStatInfo " << GetName() << " vector size not consistent " << endl;
    return 1;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputCellsInfo(string& s, size_t what,size_t when) const
{
//  bool debug = true;
  bool debug = false;

  if( debug ) cout << " Calorimeter::OutputStatInfo What=" << what << " When=" << when << endl;
  bool recognised_info = true;
  char o[500000];
  o[0]=0;
  sprintf(o,"### \n");
  if( what == CALIB )
    sprintf(o+strlen(o),"Monitoring of CALIB Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == LED )
    sprintf(o+strlen(o),"Monitoring of LED Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == PED )
    sprintf(o+strlen(o),"Monitoring of PED Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == TIME )
    sprintf(o+strlen(o),"Monitoring of TIME Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  if( what == NOISE )
    sprintf(o+strlen(o),"Monitoring of NOISE Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
  else
  {
    sprintf(o+strlen(o),"Monitoring of UNKNOWN_INFO Calorimeter: %s Ncells=%zu \n",GetName().c_str(),NCells());
    recognised_info = false;
  }
  sprintf(o+strlen(o),"### \n");

  if( recognised_info )
  {
    for( size_t icell=0; icell!=NCells(); icell++ )
    {
      double entries = cells_info[what][when][icell].GetEntries();
      double cfnew=cells_info[what][when][icell].GetMean();
      double scfnew=cells_info[what][when][icell].GetSigma();
      sprintf(o+strlen(o)," %zu %12.6f %12.6f %12.6f %s \n",icell,cfnew,scfnew,entries,GetCellName(icell).c_str());
    }
  }
  string ss(o);
//  s += ss;
  s = ss;
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputCellsInfoXY(string& s, size_t what, size_t when) const
{
  assert( XYRegularGrid() );

  char o[500000];
  o[0]=0;
//  sprintf(o," %s COMPASS Calorimeter  \n",GetName().c_str());
  if( what == CALIB )
    sprintf(o+strlen(o),"Monitoring of CALIB Calorimeter: %s Ncells= %zu \n",GetName().c_str(),NCells());
  else
  if( what == LED )
    sprintf(o+strlen(o),"Monitoring of LED Calorimeter: %s Ncells= %zu \n",GetName().c_str(),NCells());
  else
  if( what == PED )
    sprintf(o+strlen(o),"Monitoring of PED Calorimeter: %s Ncells= %zu \n",GetName().c_str(),NCells());
  else
  if( what == TIME )
    sprintf(o+strlen(o),"Monitoring of TIME Calorimeter: %s Ncells= %zu \n",GetName().c_str(),NCells());
  else
  if( what == NOISE )
    sprintf(o+strlen(o),"Monitoring of NOISE Calorimeter: %s Ncells= %zu \n",GetName().c_str(),NCells());
  else
  {
    sprintf(o+strlen(o),"Monitoring of UNKNOWN_INFO Calorimeter: %s Ncells= %zu \n",GetName().c_str(),NCells());
    return -2;
  }
  sprintf(o+strlen(o)," X     Y      MeanPeak   SigmaPeak  Statistic \n");

  for( int x=0; x != GetNColumns(); x++ )
    for( int y=0; y!= GetNRows(); y++ )
    {
      int icell = GetCellOfColumnRow(x,y);
      if( icell >= 0 )
      {
        double entries = cells_info[what][when][icell].GetEntries();
        double peaknew=cells_info[what][when][icell].GetMean();
        double speaknew=cells_info[what][when][icell].GetSigma();
        sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,
                                peaknew,speaknew,entries);
      }
      else
      {
        sprintf(o+strlen(o)," %4d %4d %11.2f %11.2f %11.2f \n",x,y,-1.,-1.,-1.);
      }
    }

  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputCellsInfo(const string &s,size_t what,size_t when)
{
//  bool debug = true;
  bool debug = false;
  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
  getline(is,str);
  unsigned int not_matched_cells_cnt = 0;
  while( getline(is,str) )
  {
    int icell;
    float mean,sigma,stat;
    char cell_name[CDB_LINE_MAX];
    assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
    ret = sscanf(str.c_str(), " %d %g %g %g %s", &icell, &mean, &sigma, &stat, cell_name);
    assert( ret == 5 );
    int jcell = FindCellByName( cell_name );
    if(debug) cout << " Cell name " << cell_name << " jcell " << jcell << endl;
    if( !(jcell >=0 && jcell < (int)NCells()) )
    {
      cerr << " Unexpected Cell name " << cell_name << " in Calorimeter::InputCellsInfo "
	   << GetName() << "  " << when  << endl;
      cerr << " Input string: " << endl;
      cerr << str << endl;
      cerr << " It well might be that you use wrong calibration file for this detector " << GetName() << endl;
      stringstream ss;
      ss << " Unexpected Cell name " << GetName() << "Calorimeter::InputCellsInfo " << str;
      throw Exception( ss.str().c_str() );
    }

    if( icell != jcell ) {
      icell = jcell;
      not_matched_cells_cnt++;
      if (not_matched_cells_cnt==1) {
	cerr << "Notice: " <<  __func__ << " " << GetName() << " "
	     << "cell id and cell name do not match.  (The calibration file "
	     << "was produced with a different geometry description, which "
	     << "might be an indication that you're using the wrong file.)"
	     << endl;
      }
    }

    if(icell >=0 && icell < (int)NCells() )
    {
      if(debug) printf(" Cell %d Mean=%g Sigma=%g Stat=%g cell_name=%s \n",icell,mean,sigma,stat,cell_name);
      SetCellInfo(what,when,icell,StatInfo(stat,mean,sigma));
      if(debug) cout << " check mean " << cells_info[what][when][icell].GetMean() <<
              " check sigma " << cells_info[what][when][icell].GetSigma() << endl;
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputCellsInfoXY(const string &s,size_t what, size_t when)
{
  assert( XYRegularGrid() );

  istringstream is(s);
  char calorim_name[CDB_LINE_MAX],dummy[CDB_LINE_MAX];
  string str;

//  Read line of comments
  getline(is,str);
  assert( str.length() < CDB_LINE_MAX );  // avoid buffer overflow
  int ret = sscanf(str.c_str(), "%s %s", calorim_name, dummy);
  assert( ret == 2 );
//  Read one more line of comments
  getline(is,str);


//  Read calibration
  while( getline(is,str) )
  {
    int x,y;
    float w,m,s;
    ret = sscanf(str.c_str(), " %d %d %g %g %g", &x, &y, &m, &s, &w);
    assert( ret == 5 );
//    cout << " x=" << x << " y=" << y << " m=" << m << " s=" << s << " w=" << w << " c=" << c << endl;
    int icell = GetCellOfColumnRow( x, y);
    if ( icell >= 0 )
      SetCellInfo(what,when,icell,StatInfo(w,m,s));
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::UpdateFrontEndDependentSettings( void )
{

  bool debug = false;
  if( GetName() == "EC02P1__" ) debug = true;
//  bool debug = false;
  if( debug )
  {
    cout << " Calorimeter::UpdateFrontEndDependentSettings " << GetName() <<  endl;
    cout << "  CALIB ready ? " << cells_info[CALIB][OLD].size() << endl;
    cout << "  options.delta_spars_mode_ = " << options.delta_spars_mode_ <<  "cf[0]= " << cells_info[CALIB][OLD][0] << endl;
    cout << " options.default_calibration  = " << options.default_calibration  << endl;

  }

  for ( size_t i=0; i<NCells(); i++ )
  {
    energy_cut_sparse_mode [i] = options.delta_spars_mode_*cells_info[CALIB][OLD][i].GetMean();
  }

  if( debug )
  {
    cout << " Calorimeter::UpdateFrontEndDependentSettings finished OK " << endl;
//     if( GetName() == "EC02P1__" ) exit(0);
  }

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::UpdateInternalAfterNewSettings( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug )
  {
    cout << " Calorimeter::UpdateInternalAfterNewSettings " << GetName() << endl;
    cout << " update_ecut_using_prob_cell_noise = " << options.update_ecut_using_prob_cell_noise << endl;
  }

  UpdateFrontEndDependentSettings();

  if( !options.update_ecut_using_prob_cell_noise ) return;

  if( debug )
  {
    cout << " prob_cell_noise = " << options.prob_cell_noise << endl;
    cout << " prob_gamma_noise = " << options.prob_gamma_noise << endl;
  }

  int indx_pp = 1;
  double pp =options.prob_cell_noise;
  if( pp >= 1.e-01 )
    indx_pp = 1;
  else if( pp < 1.e-01 && pp >= 5.e-02 )
    indx_pp = 2;
  else if( pp < 5.e-02 && pp >= 1.e-02 )
    indx_pp = 3;
  else if( pp < 1.e-02 && pp >= 5.e-03 )
    indx_pp = 4;
  else if( pp < 5.e-03 && pp >= 1.e-03 )
    indx_pp = 5;
  else if( pp < 1.e-03 && pp >= 5.e-04 )
    indx_pp = 6;
  else if( pp < 5.e-04 && pp >= 1.e-04 )
    indx_pp = 7;
  else if( pp < 1.e-04 && pp >= 5.e-05 )
    indx_pp = 8;
  else if( pp < 5.e-05 && pp >= 1.e-05 )
    indx_pp = 9;
  else if( pp < 1.e-05 && pp >= 5.e-06 )
    indx_pp = 10;
  else
    indx_pp = 10;

  int indx_prob_low =indx_pp;

  pp =options.prob_gamma_noise;
  if( pp >= 1.e-01 )
    indx_pp = 1;
  else if( pp < 1.e-01 && pp >= 5.e-02 )
    indx_pp = 2;
  else if( pp < 5.e-02 && pp >= 1.e-02 )
    indx_pp = 3;
  else if( pp < 1.e-02 && pp >= 5.e-03 )
    indx_pp = 4;
  else if( pp < 5.e-03 && pp >= 1.e-03 )
    indx_pp = 5;
  else if( pp < 1.e-03 && pp >= 5.e-04 )
    indx_pp = 6;
  else if( pp < 5.e-04 && pp >= 1.e-04 )
    indx_pp = 7;
  else if( pp < 1.e-04 && pp >= 5.e-05 )
    indx_pp = 8;
  else if( pp < 5.e-05 && pp >= 1.e-05 )
    indx_pp = 9;
  else if( pp < 1.e-05 && pp >= 5.e-06 )
    indx_pp = 10;
  else
    indx_pp = 10;

  int indx_prob_high =indx_pp;

  if( debug )
  {
    cout << " index prob_cell_noise = " << indx_prob_low << endl;
    cout << " index prob_gamma_noise = " << indx_prob_high << endl;
  }

  for( int i=0; i < (int)NCells(); i++ )
  {
    if( energy_cut_[indx_prob_low][i] > 0. )
    {
      energy_cut_bad_cells_old[i] = energy_cut_[indx_prob_low][i];
    }
    else
    {
      int ii=0;
      for( ii=0; ii<indx_prob_low; ii++ )
      {
        if( energy_cut_[ii][i] > 0 ) break;
      }
      if( ii == indx_prob_low )
         energy_cut_bad_cells_old[i] = 0.005;      // set at some "absolute" minimium
      else
         energy_cut_bad_cells_old[i] = energy_cut_[ii][i];

    }

    if( energy_cut_[indx_prob_high][i] > 0. )
    {
      energy_gamma_cut_bad_cells_old[i] = energy_cut_[indx_prob_high][i];
    }
    else
    {
      int ii=0;
      for( ii=0; ii<indx_prob_high; ii++ )
      {
        if( energy_cut_[ii][i] > 0 ) break;
      }
      if( ii == indx_prob_high )
         energy_gamma_cut_bad_cells_old[i] = 0.010;      // set at some "absolute" minimium
      else
         energy_gamma_cut_bad_cells_old[i] = energy_cut_[ii][i];
    }
    if( debug ) cout << " cell " << i << " " << GetCellName(i) << " ecut " << energy_cut_bad_cells_old[i] << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputPositionInBurst(const string &s)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s);
  string str;
//  map_burst2position_.clear();
  if( debug ) cout <<  " InputPositionInBurst " << GetName() << endl;

//  Read line of comments
  getline(is,str);
//  Read one more line of comments
  getline(is,str);


//  Read calibration
  while( getline(is,str) )
  {
    int run,from_burst,to_burst,eventmin, eventmax;
    float x,y;
    int ret = sscanf(str.c_str(), " %d %d %d %d %d %g %g",
                     &run, &from_burst, &to_burst, &eventmin, &eventmax, &x, &y);
    assert( ret == 5 );

    if( debug ) cout << " From burst " << from_burst <<" to burst " << to_burst << " x= " << x << " y="<< y << endl;
    vector <double> v;
    v.push_back(x);
    v.push_back(y);
    for (int ib=from_burst; ib<= to_burst; ib++ )
    {
      map_burst2position_.insert(std::pair< int, vector<double> > (ib,v) );
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::OutputPositionInBurst ( string& s) const
{
  bool debug = true;
  if( debug ) cout << " OutputPositionInBurst " << GetName() << endl;
  char o[500000];
  o[0]=0;
  sprintf(o," %s COMPASS Calorimeter  Position in burst \n",GetName().c_str());
  sprintf(o+strlen(o),"## Run Spills Range              X       Y     Z(not used) \n");
  for( size_t n=0; n< positionX_in_burst_new_.size(); n++ )
  {
    double x = positionX_in_burst_new_[n].GetMean();
    double y = positionY_in_burst_new_[n].GetMean();
    double z=0.;
    int evmin = -1;
    int evmax = -1;
    sprintf(o+strlen(o)," %d %zu %zu %d %d %11.2f %11.2f %11.2f \n", -1, n+1, n+1, evmin, evmax, x, y, z);
  }
  string ss(o);
  s = ss;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputJOUJOUInfo(const string &s)
{
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputPositionInTime(const string &s)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s);
  string str;
  if( debug ) cout <<  " InputPositionInTime " << GetName() << endl;

  double xoffset = options.position_offset[0];
  double yoffset = options.position_offset[1];

//  Read calibration
  while( getline(is,str) )
  {
    char dummy[CDB_LINE_MAX];
    float x_rel,y_rel,x_abs,y_abs;
    assert( 24 < CDB_LINE_MAX );  // avoid buffer overflow
    int ret = sscanf(str.c_str(), " %g %g %g %g %24c", &x_abs, &x_rel, &y_abs, &y_rel, dummy);
    assert( ret == 5 );
    dummy[24] = '\0';  // need to set trailing zero manually with %c
    if( debug ) cout << " " << str << " " << endl;

    vector < double > pos;
    pos.push_back(-x_rel*10.+xoffset);
    pos.push_back(-y_rel*10.+yoffset);
    pos.push_back(-x_abs*10.+xoffset);
    pos.push_back(-y_abs*10.+yoffset);
    tm ttm;
    if( strptime(dummy,"%a %b %d %T %Y", &ttm) == NULL )
    {
      cerr << " InputPositionInTime Format problems " << GetName() << endl;
      cerr << " " << str << " " << endl;
      exit(1);
      continue;
    }
    else
    {
//      if( debug ) cout << " No Format problems time " << endl;
      time_t timet = mktime( &ttm );
      if( debug ) cout << " No Format problems time " << timet << endl;
      alignment_in_time_.map_time2position_.push_back( std::pair< time_t, vector<double> > (timet,pos) );
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputPositionInTimeRedundant(const string &s)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s.c_str());
  string str;
  if( debug ) cout <<  " InputPositionInTimeRedundant " << GetName() << endl;

  double xoffset = options.position_offset[0];
  double yoffset = options.position_offset[1];
  if( debug ) cout <<  " xoffset = " << xoffset << " yoffset = " << yoffset << endl;

//  Read calibration
  while( getline(is,str) )
  {
    char dummy[132];
    float x_rel,y_rel,x_abs,y_abs,x_xabs,y_xabs;
    sscanf(str.c_str()," %g %g %g %g %g %g %24c",&x_abs,&x_rel,&y_abs,&y_rel,&x_xabs,&y_xabs,dummy);
    if( debug ) cout << " " << str << " " << endl;

    vector < double > pos;
    pos.push_back(-x_rel*10.+xoffset);
    pos.push_back(-y_rel*10.+yoffset);
    pos.push_back(-x_abs*10.+xoffset);
    pos.push_back(-y_abs*10.+yoffset);
//    string time_string(dummy);
//    if( debug ) cout << " Try to decode time stamp: " << time_string << " " << endl;
    if( debug ) cout << " Try to decode time stamp: " << dummy << " " << endl;
    tm ttm;
    if( strptime(dummy,"%a %b %d %T %Y", &ttm) == NULL )
    {
      cerr << " InputPositionInBurst Format problems " << GetName() <<  endl;
      cerr << " " << str << " " << endl;
      exit(1);
      continue;
    }
    else
    {
//      if( debug ) cout << " No Format problems time " << endl;
      time_t timet = mktime( &ttm );
// Stupid fix for ??????????????????
//      timet += 3600*2;
      if( debug ) cout << " No Format problems time " << timet << endl;
      alignment_in_time_.map_time2position_.push_back( std::pair< time_t, vector<double> > (timet,pos) );
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InputPositionInTimeGeneric(const string &s)
{
//   bool debug = true;
  bool debug = false;
  istringstream is(s.c_str());
  string str;
  if( debug ) cout <<  " InputPositionInTimeGeneric " << GetName() << endl;

  double xoffset = options.position_offset[0];
  double yoffset = options.position_offset[1];
  if( debug ) cout <<  " xoffset = " << xoffset << " yoffset = " << yoffset << endl;

//  Read calibration
  while( getline(is,str) )
  {
    char dummy[132];
    float x_abs,y_abs;
    sscanf(str.c_str()," %g %g %24c",&x_abs,&y_abs,dummy);
    if( debug ) cout << " " << str << " " << endl;

    vector < double > pos;
    pos.push_back(x_abs+xoffset);
    pos.push_back(y_abs+yoffset);
    pos.push_back(x_abs+xoffset);
    pos.push_back(y_abs+yoffset);
    if( debug ) cout << " Try to decode time stamp: " << dummy << " " << endl;
    tm ttm;
    if( strptime(dummy,"%a %b %d %T %Y", &ttm) == NULL )
    {
      cerr << " InputPositionInTimeGeneric Format problems " << GetName() <<  endl;
      cerr << " " << str << " " << endl;
      exit(1);
      continue;
    }
    else
    {
//      if( debug ) cout << " No Format problems time " << endl;
      time_t timet = mktime( &ttm );
// Stupid fix for ??????????????????
      timet += 3600*2;
      if( debug ) cout << " No Format problems time " << timet << endl;
      alignment_in_time_.map_time2position_.push_back( std::pair< time_t, vector<double> > (timet,pos) );
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::OperationsWithCalibInfo( void )
{
  if( options.scale_old_calib_info ) ScaleCalibration(OLD);
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SpecialOperationsWithCalibInfo( void )
{
  bool debug = false;
//   bool debug = true;

  std::vector<StatInfo> cf0;
  std::vector<StatInfo> cf1;
  std::vector<StatInfo> cf2;
  std::vector<StatInfo> cf3;

  cout << " Start SpecialOperationsWithCalibInfo for " << GetName() << endl;
  {
    const std::string file1_name("EC02P1___CALIB~~start-Run_57260~~finish-Run_99999");
    std::cout << " Try to open file1 " <<  file1_name << std::endl;
    std::ifstream f1(file1_name.c_str());
    string oo;
    char cc[1000000];

    if( !f1.is_open() )
    {
      std::cerr << " ERROR in SpecialOperationsWithCalibInfo " << GetName() <<
                                       " Can not open file1 " << file1_name << std::endl;
      return;
    }
    else
    {

      for( size_t line_n=1; true; line_n++ )
      {
        string line;
        getline(f1,line);
        if( f1.eof() )
          break;
        sprintf(cc+strlen(cc),"%s\n",line.c_str());
      }

      string s(cc);
      if( debug ) cout << s;
      InputCalibInfo(OLD, s);

      for( size_t ic=0; ic<NCells(); ic++ )
      {
        cf1.push_back( cells_info[CALIB][OLD][ic] );
      }
    }
    cout << " cf1 " << cf1.size() << endl;

  }
  {

    const std::string file1_name("EC02P1___CALIB_PRIM~~start-Run_57260~~finish-Run_99999");
    std::cout << " Try to open file1 " <<  file1_name << std::endl;
    std::ifstream f1(file1_name.c_str());
    string oo;
    char cc[1000000];

    if( !f1.is_open() )
    {
      std::cerr << " ERROR in SpecialOperationsWithCalibInfo " << GetName() <<
                                       " Can not open file1 " << file1_name << std::endl;
      return;
    }
    else
    {

      for( size_t line_n=1; true; line_n++ )
      {
        string line;
        getline(f1,line);
        if( f1.eof() )
          break;
        sprintf(cc+strlen(cc),"%s\n",line.c_str());
      }

      string s(cc);
      if( debug ) cout << s;
      InputCalibInfo(OLD, s);

      for( size_t ic=0; ic<NCells(); ic++ )
      {
        cf2.push_back( cells_info[CALIB][OLD][ic] );
      }
    }
    cout << " cf2 " << cf2.size() << endl;

  }
  {

    const std::string file1_name("EC02P1___CALIB~~start-Run_56100~~finish-Run_57259");
    std::cout << " Try to open file1 " <<  file1_name << std::endl;
    std::ifstream f1(file1_name.c_str());
    string oo;
    char cc[1000000];

    if( !f1.is_open() )
    {
      std::cerr << " ERROR in SpecialOperationsWithCalibInfo " << GetName() <<
                                       " Can not open file1 " << file1_name << std::endl;
//     exit(1);
      return;
    }
    else
    {

      for( size_t line_n=1; true; line_n++ )
      {
        string line;
        getline(f1,line);
        if( f1.eof() )
          break;
        sprintf(cc+strlen(cc),"%s\n",line.c_str());
      }

      string s(cc);
      if( debug ) cout << s;
      InputCalibInfo(OLD, s);

      for( size_t ic=0; ic<NCells(); ic++ )
      {
        cf3.push_back( cells_info[CALIB][OLD][ic] );
      }
    }
    cout << " cf3 " << cf3.size() << endl;

  }

  for( size_t ic=0; ic<NCells(); ic++ )
  {
    double newcf = cf1[ic].GetMean()/cf2[ic].GetMean()*cf3[ic].GetMean();

    StatInfo newcff(1.,newcf,1.);
    cells_info[CALIB][OLD][ic]=newcff;
  }



  cout << " Finish SpecialOperationsWithCalibInfo for " << GetName() << endl;


}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
