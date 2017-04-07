/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/Calorimeter.cc,v $
   $Date: 2011/02/04 17:54:37 $
   $Revision: 1.120 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch, Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )
     Denis     Murashev  ( Denis.Mourachev@cern.ch, Murashev@sirius.ihep.su )

   Copyright(C): 1999-2001  V.Kolosov, A.Zvyagin, D.Murashev

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
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

// --- Internal files ----
#include "Calorimeter.h"

#include "CalorimeterParticle.h"
#include "CellDataRaw.h"
#include "CellType.h"
#include "DataBase.h"
#include "Exception.h"
#include "GUICalorimeter.h"
#include "MCConstruction.h"
#include "Reconstruction.h"
#include "ReconstructionCombined.h"
#include "ReconstructionKolosov.h"
#include "ReconstructionLednev.h"

using namespace std;

namespace Reco {

template <class T> T sqr(const T &x) {return x*x;}

const char *Calorimeter::CellInfoNames[]={"CALIB","LED","PED","TIME","NOISE","CLUSTER"};

////////////////////////////////////////////////////////////////////////////////

Calorimeter::~Calorimeter(void)
{
  delete gui;
  if ( options.blackbox_chk_file )
    delete options.blackbox_chk_file;

  if (energy_cut_sparse_mode!=NULL)
      delete [] energy_cut_sparse_mode;
  if (energy_cut_bad_cells_old!=NULL)
      delete [] energy_cut_bad_cells_old;
  if (energy_cut_bad_cells_new!=NULL)
      delete [] energy_cut_bad_cells_new;
  if (energy_gamma_cut_bad_cells_old!=NULL)
      delete [] energy_gamma_cut_bad_cells_old;
  if (energy_gamma_cut_bad_cells_new!=NULL)
      delete [] energy_gamma_cut_bad_cells_new;

  for (size_t i=0; i<CellInfoTypeSize; i++)
    for (size_t j=0; j<CalibTimeTypeSize; j++) {
      if (cells_stat_info[i][j]!=NULL)
        delete [] cells_stat_info[i][j];
      if (cells_store[i][j]!=NULL)
        delete [] cells_store[i][j];
    }

  if (reconstruction_!=NULL)
    delete reconstruction_;
}

////////////////////////////////////////////////////////////////////////////////

Calorimeter::Calorimeter(const string &the_name, const string &geom_file) :
  name(the_name),
  gui(NULL),
  calo_hist(NULL),
  internal_hist_(NULL),
  profiles_hist_(NULL),
  internal_correlations_hist_(NULL),
  external_correlations_hist_(NULL),
  trig_group_hist_(NULL),
  make_particles_histo_(NULL),
  fit_info_histo_(NULL),
  store_led_(NULL),
  fem_(NULL),
  energy_cut_sparse_mode(NULL),
  energy_cut_bad_cells_old(NULL),
  energy_cut_bad_cells_new(NULL),
  energy_gamma_cut_bad_cells_old(NULL),
  energy_gamma_cut_bad_cells_new(NULL),
  final_initialization_(false),
  evid_updated_(false),
  monteconstruction_(NULL),
  monteconstruction_cleared_(false),
  reconstruction_(NULL),
  histograms_basedir_(NULL)
{

  bool debug = false;
  if( debug ) cout << " Calorimeter constructor from " << geom_file << endl;

  position_atstartofrun_[0] = 0;
  position_atstartofrun_[1] = 0;
  position_atstartofrun_[2] = 0;
  position_correction_[0]   = 0;
  position_correction_[1]   = 0;
  position_correction_[2]   = 0;
  SetPosition(0, 0, 0);
  SetVertexPosition(0, 0, 0);
  InitOptions();

  if ( geom_file.empty() )
    return;

  ifstream f(geom_file.c_str());

  if( !f.is_open() )
    throw Exception("Calorimeter::Calorimeter():  can not open file \"%s\"",
		    geom_file.c_str());

  bool calo_found=false;

  for( size_t line_n=1; true; line_n++ )
    {
      char line[555];

      f.getline( line, sizeof(line) );
      if( f.eof() )
        break;

      if( line[0]==0 || line[0]=='#' )
        continue;

      istringstream s(line);
      if( debug )
      {
        cout << " new line " << string(line) << endl;
      }
      string opt;
      if( !(s >> opt) )
        throw Exception("Calorimeter::Calorimeter(): bad format: \"%s\"",line);

      if( opt == "cell" ) {
        if( debug ) cout << " new  CellType definition " << endl;
        // CellType definition
        string cell_name;
        double RadiationLength, NuclearLength, sizeXYZ[3];
	double ConstantTerm, StochasticTerm, ReadoutTerm, active_material, density;

        if( !( s >> cell_name
	       >> sizeXYZ[0] >> sizeXYZ[1] >> sizeXYZ[2]
	       >> RadiationLength
	       >> NuclearLength
	       >> ConstantTerm
	       >> StochasticTerm
	       >> ReadoutTerm
	       >> active_material
	       >> density) )
          throw Exception("Calorimeter::Calorimeter(): bad format in cell setting: \"%s\"",line);

	// critical energy is only contained in new versions of detectors.dat,
	// therefore we default to a reasonable (average) value
	double crit_energy = 11.;
	s >> crit_energy;
	crit_energy /= 1000.;  // convert to GeV

	string hw_type;
	s >> hw_type;

        if( debug ) cout << " CellType  " << cell_name << endl;

        for( list<CellType>::const_iterator it=cells_type.begin(); it!=cells_type.end(); it++ )
	  if( it->GetName()==cell_name )
	    throw Exception("Calorimeter::Calorimeter():  second cell type "
			    "\"%s\" definition on line %d of file \"%s\"",
			    cell_name.c_str(),line_n,geom_file.c_str());

        if( debug )
	  cout << " insert  CellType  " << cell_name <<" RadL " << RadiationLength
	       << " NuclL " << NuclearLength << endl;

        cells_type.insert( cells_type.end(), CellType(cell_name,
						      sizeXYZ[0], sizeXYZ[1], sizeXYZ[2],
						      sizeXYZ[0], sizeXYZ[1],
                                                      RadiationLength, NuclearLength,
						      crit_energy,
                                                      StochasticTerm, ConstantTerm, ReadoutTerm,
                                                      active_material, density, hw_type) );
      }
      else if( opt == "calo" ) {
	if( debug ) cout << " new  Calorimeter definition " << endl;
	    // Calorimeter definition

        string name;

        if( !(s >> name) )
          throw Exception("Calorimeter::Calorimeter(): bad format: \"%s\"",line);
        if( debug ) cout << " Calorimeter " << name << endl;

        if( name!=GetName() )
          continue;       // This is another calorimeter

        if( calo_found )
          throw Exception("Calorimeter::Calorimeter():  second calorimeter "
			  "definition on line %d of file \"%s\"",
			  line_n, geom_file.c_str());
        calo_found=true;

        double pos[3];
        if( !(s >> pos[0] >> pos[1] >> pos[2]) )
          throw Exception("Calorimeter::Calorimeter(): bad format: \"%s\"",line);

        SetPosition(pos[0],pos[1],pos[2]);
      }
      else if( opt == "calm" ) {
	if( debug ) cout << " new  Calorimeter matrix definition " << endl;
	    // Calorimeter matrixdefinition

        string name;

        if( !(s >> name) )
          throw Exception("Calorimeter::Calorimeter(): bad format: \"%s\"",line);
        if( debug ) cout << " Calorimeter " << name << endl;

        if( name!=GetName() )
          continue;       // This is another calorimeter

        std::string matrixName, cellTypeName;
        unsigned int nH, nV, nHoles;
        double pos[3];
        if ( !(s >> matrixName >> cellTypeName >> nH >> nV >> pos[0] >> pos[1] >> pos[2] >> nHoles) )
          throw Exception("Calorimeter::Calorimeter(): bad format: \"%s\"",line);

        std::list<CellType>::const_iterator it=cells_type.begin();
        for (; it!=cells_type.end(); it++)
          if (it->GetName()==cellTypeName)
              break;
        if (it==cells_type.end())
          throw Exception("Calorimeter::Calorimeter(): unknown cell type: \"%s\"", cellTypeName.c_str());
        const CellType& cellType = *it;

        matrixes.push_back(CellsMatrix(matrixName, cellType, nH, nV, pos[0], pos[1], pos[2]));
        CellsMatrix& matrix = matrixes.back();

        for (unsigned int i=0; i<nHoles; i++) {
          unsigned int x1, x2, y1, y2;
          if ( !(s >> x1 >> x2 >> y1 >> y2) )
            throw Exception("Calorimeter::Calorimeter(): bad format: \"%s\"",line);

          matrix.AddHole(CellsMatrix::Hole(x1, x2, y1, y2));
        }

        const double xStart = pos[0] - (((double)nH)/2.-0.5)*cellType.GetSizeX();
        const double yStart = pos[1] - (((double)nV)/2.-0.5)*cellType.GetSizeY();

        for (unsigned int v=0; v<nV; v++) {
          for (unsigned int h=0; h<nH; h++) {
            if (!matrix.IsInHole(h, v))
              cells.push_back(Cell(cellType, true, xStart + h*cellType.GetSizeX(), yStart + v*cellType.GetSizeY(), pos[2]));
          }
        }
      }
      else
        cerr << "Calorimeter::Calorimeter():   Unknown keyword " << opt << endl;
    }

  if( !calo_found )
    throw Exception("Calorimeter::Calorimeter():  calorimeter \"%s\" was not found in file \"%s\"",
		    GetName().c_str(),geom_file.c_str());
}

////////////////////////////////////////////////////////////////////////////////
void Calorimeter::InitOptions(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InitXY(void) {
//  bool debug = true;
  bool debug = false;
  if( debug )  cout << " Calorimeter " << GetName() << " InitXY started" << endl;

  assert( NCells() > 0 );

  const double epsilon = .01;  // epsilon used for arithmetic comparison [mm]

  // Make sure that the z position of the front surface is identical.
  // (This must be checked for ECAL1, too, therefore we have to put it before
  // the regular grid test which aborts in the ECAL1 case.)

  double tolerate_in_xy_structure_check = GetOptions().tolerance_for_nearby_cells;

  if( options.check_front_surface_position )
  {
    double zmin =  1e9;
    double zmax = -1e9;
    for (vector<Cell>::const_iterator cell=cells.begin(); cell!=cells.end(); cell++) {
      // z position of front edge
      const double zfront = cell->GetZ() - cell->GetCellType().GetSizeZ()/2.;
      if ( zfront < zmin )
        zmin = zfront;
      if ( zfront > zmax )
        zmax = zfront;
    }

    if ( zmax-zmin > epsilon )
      throw Exception("Calorimeter %s front surface z delta %gmm!  Fix detectors.dat!",
              GetName().c_str(), zmax-zmin);

    if ( debug )
      cout << "Calorimeter " << GetName() << " front surface z delta: " << zmax-zmin << endl;
  }

  // Make sure that there is no overlap between the active area of the cells.
  // (Might fail to detect overlap if the shape of the cells isn't quadratic.)
  // Padding must be excluded from this test because it is missing at some
  // boundaries between cell matrices (cmtx'es), like eg. LG10/LG18 of ECAL 1.

  // Definition of the model of the geometry:  The center of the cell as
  // defined in detectors.dat is assumed to be the same as the center of the
  // active part of the cell.  In other words:  Padding is assumed to be
  // distributed symmetrically around the active part (but possibly
  // differently in X and Y directions).
  {
    int i=0;
    for (vector<Cell>::const_iterator cell1=cells.begin(); cell1!=cells.end(); cell1++, i++) {
      const double left   = cell1->GetLeftActive()   + epsilon;
      const double right  = cell1->GetRightActive()  - epsilon;
      const double bottom = cell1->GetBottomActive() + epsilon;
      const double top    = cell1->GetTopActive()    - epsilon;

      int j=0;
      for (vector<Cell>::const_iterator cell2=cells.begin(); cell2!=cells.end(); cell2++, j++) {

	// we have to check all combinations (not just cell1 < cell2) except cell1==cell2
	if ( cell1 == cell2 ) continue;

	if (    cell2->InTrueActiveArea(left,  bottom, cell1->GetZ())
	     || cell2->InTrueActiveArea(right, bottom, cell1->GetZ())
	     || cell2->InTrueActiveArea(left,     top, cell1->GetZ())
	     || cell2->InTrueActiveArea(right,    top, cell1->GetZ()) )
	  throw Exception( "Calorimeter::InitXY(): %s has overlap between cells "
			   "#%i (%.3f,%.3f)-(%.3f,%.3f) and "
			   "#%i (%.3f,%.3f)-(%.3f,%.3f)!",
			   GetName().c_str(),
			   i, cell1->GetLeftActive(), cell1->GetBottomActive(),
			   cell1->GetRightActive(), cell1->GetTopActive(),
			   j, cell2->GetLeftActive(), cell2->GetBottomActive(),
			   cell2->GetRightActive(), cell2->GetTopActive() );
      }
    }
  }

  // Stop here for inhomogenous calorimeters.
  if ( options.mixed_blocks ) {
    fNcols_ = -1;
    fNrows_ = -1;
    fXStep_ = 0.;
    fYStep_ = 0.;
    return;
  }

  // Make sure that cells lie on a regular grid.
  fXStep_            = cells[0].GetCellType().GetStepX();
  fYStep_            = cells[0].GetCellType().GetStepY();
  const double XSize = cells[0].GetCellType().GetTrueSizeX();
  const double YSize = cells[0].GetCellType().GetTrueSizeY();
  fNcols_            = (int) rint( (GetXmax()-GetXmin()) / fXStep_ );
  fNrows_            = (int) rint( (GetYmax()-GetYmin()) / fYStep_ );
  const double xmin  = GetTrueXmin();
  const double ymin  = GetTrueYmin();
  size_t i=0;
  for (vector<Cell>::const_iterator cell=cells.begin(); cell!=cells.end(); cell++, i++) {
    // center of cell in calorimeter coordinates [mm]
    const double xc     = cell->GetX() - xmin;
    const double yc     = cell->GetY() - ymin;
    // cell position in integer coordinates
    const int     x     = (int) rint( (xc-XSize/2.) / fXStep_ );
    const int     y     = (int) rint( (yc-YSize/2.) / fYStep_ );
    // nominal (expected) center of cell in calo coords [mm]
    const double xc_nom = x*fXStep_ + XSize/2.;
    const double yc_nom = y*fYStep_ + YSize/2.;
    if( debug )
      cout << "InitXY:: Cell " << i
           << " dx =" << xc - xc_nom << " dy =" << yc - yc_nom << endl;

    if( fabs( xc - xc_nom ) > tolerate_in_xy_structure_check || fabs( yc - yc_nom ) > tolerate_in_xy_structure_check )
      throw Exception("%s: InitXY Check failed at cell index %i, X=%i, Y=%i, "
		      "xc=%f, yc=%f, nominal xc=%f, nominal yc=%f, "
		      "nominal xstep=%f, nominal ystep=%f",
		      GetName().c_str(), i, x, y, xc, yc,
		      xc_nom, yc_nom, fXStep_, fYStep_);
  }

  if ( debug )
    cout << "Calorimeter " << GetName() << " InitXY Check OK " << endl;

  map_cell_xy_.clear();
  map_xy_cell_[0].clear();
  map_xy_cell_[1].clear();
  for( int it=0; it!=fNcols_*fNrows_; it++ ) {
    map_cell_xy_.push_back(-1);
  }
  for( size_t it=0; it<NCells(); it++ ) {
    map_xy_cell_[0].push_back(-1);
    map_xy_cell_[1].push_back(-1);
  }

  i=0;
  for( vector<Cell>::const_iterator cell=cells.begin(); cell!=cells.end(); cell++, i++ ) {
    const int x = (int) rint( (cell->GetX() - xmin - fXStep_/2.) / fXStep_ );
    const int y = (int) rint( (cell->GetY() - ymin - fYStep_/2.) / fYStep_ );
    if( x < 0 || x >= fNcols_ || y < 0 || y >= fNrows_ )
      throw Exception("%s: ERROR xy setting: Cell outside calorimeter: x=%i, y=%i",
		      GetName().c_str(), x, y);
    if ( map_cell_xy_[x + fNcols_*y] != -1 )
      throw Exception("%s: ERROR xy setting: Duplicate cell: x=%i, y=%i",
		      GetName().c_str(), x, y);
    map_cell_xy_[x + fNcols_*y] = i;
    map_xy_cell_[0][i] = x;
    map_xy_cell_[1][i] = y;
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::UpdateAfterOptionsSettings ( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug) cout << "  Calorimeter::UpdateAfterOptionsSettings debug " << GetName() << endl;
  if( debug)
  {
    cout <<" cells_info[CALIB][OLD].size() = " << cells_info[CALIB][OLD].size() << endl;
    cout <<" cells_info[CALIB][OLD].size() = " << cells_info[CALIB][OLD].size() << endl;
    cout <<" cells_info[CALIB][PRIM].size() = " << cells_info[CALIB][PRIM].size() << endl;

    cout <<" cells_info[CALIB][MC].size() = " << cells_info[CALIB][MC].size() << endl;

    cout <<" cells_info[TIME ][OLD].size() = " << cells_info[TIME ][OLD].size() << endl;
    cout <<" cells_info[TIME ][NEW].size() = " << cells_info[TIME ][NEW].size() << endl;
    cout <<" cells_info[TIME ][PRIM].size() = " << cells_info[TIME ][PRIM].size() << endl;
  }

  for ( int i=0; i<3; i++) {
    vertex_position[i]      = options.vertex_position[i];
    position_correction_[i] = options.position_correction[i];
  }

  for( std::list<CellType>::iterator it=cells_type.begin(); it!=cells_type.end(); it++ )
    it->SetReadoutTerm( options.readout_term );

  for ( size_t i=0; i<NCells(); i++ )
  {
    cells_info[CALIB][OLD][i] = StatInfo(1,options.default_calibration,1.);
    cells_info[CALIB][PRIM][i] = StatInfo(1,options.default_calibration,1.);

    cells_info[CALIB][MC][i] = StatInfo(1,options.mc_default_calibration,1.);

    cells_info[TIME ][OLD][i] = StatInfo(1,options.default_time0_calibration,1.);
    cells_info[TIME ][NEW][i] = StatInfo(1,options.default_time0_calibration,1.);
    cells_info[TIME ][PRIM][i] = StatInfo(1,options.default_time0_calibration,1.);

  }

  if( debug) cout << " Implement electronis delta energy cuts settings in RD and MC for " <<
                           GetName() << "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << prob_cut_.size() <<  endl;
  for( size_t i=0; i<prob_cut_.size()+1; i++ )
    for ( size_t icell=0; icell<NCells(); icell++ )
    {
      energy_cut_[i][icell] = 0;
    }
  leds_in_spills_old_ = new std::map<long unsigned int,double>[NCells()];

  if( debug) cout << "  Calorimeter::UpdateAfterOptionsSettings debug OK " << GetName() << endl;

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InitializeWithOptions( const std::list<std::string>& reco_opt_list )
{

//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::InitializeWithOption " << GetName() << " reco_opt_list size = " << reco_opt_list.size() << endl;
  InitOptions();
  if( debug ) cout << " SetOptions " << endl;
  SetOptions( reco_opt_list );
  SetOptionsOK();
  if( !options.Ok() )
  {
    cerr << " Internal error in Calorimeter::Initialize " << GetName() <<
              " you must initialize and read options first " << endl;
    exit(1);
  }
  if( debug ) cout << " Reco::Calorimeter::Init() " << endl;
  Reco::Calorimeter::Init();
  if( debug ) cout << " Reco::Calorimeter::InitXY() ?? default " << XYRegularGrid() <<  endl;
  if( options.init_xy_structure) Reco::Calorimeter::InitXY();
// Here we assume that Options initialization and reding is implemented externally
  if( debug ) cout << " XYRegularGrid " << XYRegularGrid() << endl;
  if( debug ) cout << " UpdateAfterOptionsSettings " << endl;
  UpdateAfterOptionsSettings();

  if( debug ) cout << " PrintGeneralInfo ? " << options.print_general_info << endl;
  if( options.print_general_info ) PrintGeneralInfo();
  if( debug ) cout << " Calorimeter::InitializeWithOption " << GetName() << " OK!!! " << endl;
  final_initialization_ = true;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Initialize(void)
{
//   InitOptions();
//   ReadOptions();
//   SetOptionsOK();
  if( !options.Ok() )
  {
    cerr << " Internal error in Calorimeter::Initialize " << GetName() <<
              " you must initialize and read options first " << endl;
    exit(1);
  }

  // Here we assume that Options initialization and reading is implemented externally

  UpdateAfterOptionsSettings();

  Reco::Calorimeter::Init();
  Reco::Calorimeter::InitXY();
  if( options.print_general_info ) PrintGeneralInfo();
  final_initialization_ = true;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InsertCellType  ( CellType cell_type)
{
  list<CellType>::const_iterator it_cell_type;
  for( it_cell_type=cells_type.begin(); it_cell_type!=cells_type.end(); it_cell_type++ )
  {
    if( it_cell_type->GetName() == cell_type.GetName() )
      {
	cerr << " Now Fatal in " << GetName() << " ! New CellType with not unic name: "
	     << cell_type.GetName() << endl;
      }
  }
  cells_type.insert( cells_type.end(),cell_type );
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InsertCell  ( Cell cell )
{
  list<CellType>::const_iterator cell_type;
  for( cell_type=cells_type.begin(); cell_type!=cells_type.end(); cell_type++ )
    if( cell_type->GetName() == cell.GetCellType().GetName() )
            break;

   if( cell_type==cells_type.end() )
   {
     cerr << "Calorimeter::InsertCell(): in calorimeter " << GetName() <<
      " Can not find cell type. And dont want to insert new.  " << endl;
     exit(1);
   }

   cells.push_back( Cell(*cell_type, true, cell.GetX(), cell.GetY(), cell.GetZ()) );
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::InitMemory(void)
{
  bool debug = false;
//   bool debug = true;
  if( debug ) cout << " Calorimeter::InitMemory debug " << GetName() << endl;
  map_cell_xy_.clear();
  map_xy_cell_[0].clear();
  map_xy_cell_[1].clear();
  fNcols_            = -1;
  fNrows_            = -1;
  fXStep_            = 0.;
  fYStep_            = 0.;

  prob_cut_.push_back(1.e-01);
  prob_cut_.push_back(5.e-02);
  prob_cut_.push_back(1.e-02);
  prob_cut_.push_back(5.e-03);
  prob_cut_.push_back(1.e-03);
  prob_cut_.push_back(5.e-04);
  prob_cut_.push_back(1.e-04);
  prob_cut_.push_back(5.e-05);
  prob_cut_.push_back(1.e-05);
  prob_cut_.push_back(5.e-06);

  ebin_stat_.push_back(0.01);
  ebin_stat_.push_back(0.02);
  ebin_stat_.push_back(0.05);
  ebin_stat_.push_back(0.10);
  ebin_stat_.push_back(0.20);
  ebin_stat_.push_back(0.50);
  ebin_stat_.push_back(1.00);
  ebin_stat_.push_back(2.50);
  ebin_stat_.push_back(5.00);
  ebin_stat_.push_back(10.00);

  StatInfo zero = StatInfo(0.,0.,0.);
  FitInfo  fit_empty = FitInfo();

  // allocate memory for calorimeter data
  mc_input_.clear();
  signals_.clear();
  if( NCells() <= 0 ) assert(false);
  for ( size_t i=0; i<NCells(); i++ )
  {
    cells_info[CALIB][OLD].push_back(StatInfo(1,options.default_calibration,1.));
    cells_info[CALIB][PRIM].push_back(StatInfo(1,options.default_calibration,1.));
    cells_info[CALIB][NEW].push_back(zero);
    cells_info[CALIB][MONITOR].push_back(zero);
    cells_info[CALIB][MC].push_back(StatInfo(1,options.default_calibration,1.));

    cells_info[LED  ][OLD].push_back(zero);
    cells_info[LED  ][NEW].push_back(zero);
    cells_info[LED  ][PRIM].push_back(zero);
    cells_info[LED  ][MONITOR].push_back(zero);

    cells_info[PED  ][OLD].push_back(zero);
    cells_info[PED  ][NEW].push_back(zero);
    cells_info[PED  ][PRIM].push_back(zero);
    cells_info[PED  ][MONITOR].push_back(zero);

    cells_info[TIME ][OLD].push_back(StatInfo(1,options.default_time0_calibration,1.));
    cells_info[TIME ][NEW].push_back(StatInfo(1,options.default_time0_calibration,1.));
    cells_info[TIME ][PRIM].push_back(StatInfo(1,options.default_time0_calibration,1.));
    cells_info[TIME ][MONITOR].push_back(zero);

    cells_info[NOISE][OLD].push_back(zero);
    cells_info[NOISE][NEW].push_back(zero);
    cells_info[NOISE][MONITOR].push_back(zero);

    cells_info[CLUSTER][OLD].push_back(zero);
    cells_info[CLUSTER][NEW].push_back(zero);
    cells_info[CLUSTER][MONITOR].push_back(zero);

    cells_info[RAW][OLD].push_back(zero);
    cells_info[RAW][NEW].push_back(zero);
    cells_info[RAW][MONITOR].push_back(zero);

    fit_info  [CALIB].push_back(fit_empty);
    fit_info  [LED  ].push_back(fit_empty);
    fit_info  [PED  ].push_back(fit_empty);
    fit_info  [TIME ].push_back(fit_empty);
    fit_info  [NOISE].push_back(fit_empty);
    fit_info  [CLUSTER].push_back(fit_empty);
    fit_info  [RAW].push_back(fit_empty);

    individual_calib_factor_.push_back(1.);
  }
  // Noise cells info
  cells_noise_statistic=0.;
  cells_gamma_noise_statistic=0.;
  led_statistic_=0.;
  noise_histo_booked = false;

//   bad_cells_old.clear();
//   bad_cells_new.clear();

  bad_cells_[NEW].clear();
  bad_cells_[OLD].clear();
  bad_cells_[MC].clear();

//   cell_is_bad_old  = new bool[NCells()];
//   cell_is_bad_new  = new bool[NCells()];
//   status_bad_cells_old = new int[NCells()];
//   status_bad_cells_new = new int[NCells()];

  energy_cut_sparse_mode = new double[NCells()];
  energy_cut_bad_cells_old = new double[NCells()];
  energy_cut_bad_cells_new = new double[NCells()];
  energy_gamma_cut_bad_cells_old = new double[NCells()];
  energy_gamma_cut_bad_cells_new = new double[NCells()];
  bad_ped_ = new int[NCells()];
  bad_ped_cells_.clear();

  for ( size_t i=0; i<NCells(); i++ )
  {
//     cell_is_bad_new          [i]=false;
//     cell_is_bad_old          [i]=false;

    cell_is_bad_[NEW].push_back(false);
    cell_is_bad_[OLD].push_back(false);
    cell_is_bad_[MC].push_back(false);

//     status_bad_cells_old     [i]=0;
//     status_bad_cells_new     [i]=0;

    status_bad_cells_[NEW].push_back(0);
    status_bad_cells_[OLD].push_back(0);
    status_bad_cells_[MC].push_back(0);

    energy_cut_sparse_mode [i]=0;
    energy_cut_bad_cells_old [i]=0;
    energy_cut_bad_cells_new [i]=0;
    energy_gamma_cut_bad_cells_old [i]=0;
    energy_gamma_cut_bad_cells_new [i]=0;
    bad_ped_[i] = 0;
  }

  for( size_t i=0; i < CellInfoTypeSize; i++ )
    for( size_t j=0; j < CalibTimeTypeSize; j++ )
    {
      cells_stat_info[i][j] = NULL;
      cells_store[i][j] = NULL;
    }

  StatInfo *p = new StatInfo[ NCells() ];
  cells_stat_info[CALIB  ][MONITOR]= p;

  cells_store[LED][OLD] = new StatInfoStore[NCells()];
  cells_store[CALIB][MONITOR] = new StatInfoStore[NCells()];

  if( prob_cut_.size() <= 0 )
  {
    cerr << " prob arry not init: change design " << endl;
    exit(1);
  }

  if( debug) cout << " prob arry " << GetName() << " prob_cut_.size() = "  << prob_cut_.size() << endl;

  for( size_t i=0; i<prob_cut_.size()+1; i++ )
    energy_cut_[i] = new double[NCells()];

  for( size_t i=0; i<prob_cut_.size()+1; i++ )
    for ( size_t icell=0; icell<NCells(); icell++ )
      energy_cut_[i][icell] = 0.;

  if( debug ) cout << " Calorimeter::InitMemory debug OK " << GetName() << endl;

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Init(void)
{
//  bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::Init " << GetName() << endl;
  if( ! options.Ok() )
  {
    cerr << " Calorimeter::Init FATAL internal error in " << GetName() << " Options not initialized !! " << endl;
    exit(1);
  }

  DetectCellNeighbors();
  DetectCellsOnBoundary();

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CreateGUI()
{
  if( gui!=NULL )
  {
    cout << "Calorimeter::CreateGUI(): there is a GUI for calorimeter \""
         << GetName() << "\"\n";
    return;
  }

  gui = new GUICalorimeter(*this,
			   NULL,
			   this->GetName().c_str(),
			   0,
			   this->XYRegularGrid() );
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::UpdateGUI( void )
{
  if( gui!=NULL )
  {
    gui->UpdateForPreviousSpill();
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Clear(void)
{
    #if calorim_debug>1
       printf(" Clear: Calorimeter %s \n",GetName().c_str());
    #endif
  mc_input_.clear();
  signals_.clear();

  hint_particles.clear();
// ???
//   reco_particles.clear();

  for( size_t it=0; it< calorimeter_sub_sets.size(); it++ )
  {
    calorimeter_sub_sets[it].Clear();
  }
  evid_updated_ = false;
  monteconstruction_cleared_ = false;
}

////////////////////////////////////////////////////////////////////////////////

const EventID & Calorimeter::GetEventID(void ) const
{
  if( !evid_updated_ ) cerr << " Warning: Calorimeter " << GetName() << " usage of not initialized EventID " << endl;
  return evid_;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::UpdateEventID(int ev_type, int run_num, int spill_num, int ev_in_spill, int ev_in_run, const unsigned int& trigger_mask,
                                const time_t &time, const double time_in_spill, bool is_led_event )
{
  evid_.Update( ev_type, run_num, spill_num, ev_in_spill, ev_in_run, trigger_mask, time, time_in_spill, is_led_event );
  evid_updated_ = true;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::UpdateEventID( const EventID &evid )
{
  evid_ = evid;
  evid_updated_ = true;
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::SetHintParticles (const vector<CalorimeterParticle> &particles)
{
// # warning " SetHintParticles particles in DRS !!! "
   hint_particles = particles;
// As now we are somewhere inside Calorimeter we may find hint_particles hited_cell
   for( vector<CalorimeterParticle>::iterator it=hint_particles.begin(); it!=hint_particles.end(); it++ )
   {
      vector <double> x;
      x.push_back(it->GetX());
      x.push_back(it->GetY());
      x.push_back(it->GetZ());

      int m = FindCell( x );
//       int m = FindCell( it->GetX(), it->GetY(), it->GetZ() );
      if( m >= 0 ) it->SetHitedCell( m );
   }
}

/////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::AddHintParticle (const CalorimeterParticle &p)
{
// # warning " Now it is internal and importat that: SetHintParticles particles now in MRS !!! "
   CalorimeterParticle p_in = ConvertParticleMRS2DRS(p);
   vector <double> x;
   x.push_back(p_in.GetX());
   x.push_back(p_in.GetY());
   x.push_back(p_in.GetZ());

   int m = FindCell( x );

//    int m = FindCell( p_in.GetX(), p_in.GetY(), p_in.GetZ() );

   if( m >= 0 ) p_in.SetHitedCell( m );
   hint_particles.push_back(p_in);
}


void  Calorimeter::SetHitedCells ( vector<CalorimeterParticle> &p_in)
{
// Here we assume particles in MRS!
  bool debug = true;
//   bool debug = false;
  if( debug ) cout << GetName() << " Calorimeter::SetHitedCells for N particles=" <<
                                                                 p_in.size() << endl;
  for( vector<CalorimeterParticle>::iterator it=p_in.begin(); it!=p_in.end(); it++ )
  {
    SetHitedCells( *it );
  }
}

////////////////////////////////////////////////////////////////////////////////

bool  Calorimeter::SetHitedCells ( CalorimeterParticle &p_in)
{
// // Here we assume particles in MRS!
  bool debug = true;
// //   bool debug = false;
//   if( debug ) cout << GetName() << " Calorimeter::SetHitedCells for N particles=" <<
//                                                                  p_in.size() << endl;
//   for( vector<CalorimeterParticle>::iterator it=p_in.begin(); it!=p_in.end(); it++ )
//   {
//     Particle p = ConvertParticleMRS2DRS(Particle((Particle::ParticleID)it->GetID(),it->GetProb(),it->GetE(),
//                                         it->GetX(),it->GetY(),it->GetZ(),
//                                         it->GetEerr(),it->GetXerr(),it->GetYerr(),it->GetZerr(),0.,0.) );
    CalorimeterParticle p =
      ConvertParticleMRS2DRS( CalorimeterParticle( (CalorimeterParticle::ParticleID)p_in.GetID(),
						   p_in.GetProb(), p_in.GetE(),
						   p_in.GetX(), p_in.GetY(), p_in.GetZ(),
                                                   this,
						   p_in.GetEerr(), p_in.GetXerr(),
						   p_in.GetYerr(),p_in.GetZerr(),
						   0., 0.) );
//    int m = FindCell( it->GetX()-GetPositionX(), it->GetY()-GetPositionY(), it->GetZ()-GetPositionZ() );


   vector <double> x;
   x.push_back(p.GetX());
   x.push_back(p.GetY());
   x.push_back(p.GetZ());

   int m = FindCell( x );
//     int m = FindCell( p.GetX(), p.GetY(), p.GetZ() );
    if( debug ) cout << " m= " << m  << endl;
    if( m >= 0 )
    {
      p_in.SetHitedCell( m );
      return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::SetHitedCellsUltimate ( CalorimeterParticle &p_in)
{
// // Here we assume particles in MRS!
  bool debug = false;
//   bool debug = true;
    CalorimeterParticle p = ConvertParticleMRS2DRS(CalorimeterParticle((CalorimeterParticle::ParticleID)p_in.GetID(),p_in.GetProb(),p_in.GetE(),
                                        p_in.GetX(),p_in.GetY(),p_in.GetZ(),
                                        this,
                                        p_in.GetEerr(),p_in.GetXerr(),p_in.GetYerr(),p_in.GetZerr(),0.,0.) );
   vector <double> x;
   x.push_back(p.GetX());
   x.push_back(p.GetY());
   x.push_back(p.GetZ());

   int m = FindCell( x );
   if( debug ) cout << " m= " << m  << endl;
   if( m >= 0 )
   {
     p_in.SetHitedCell( m );
     return;
   }
   m= FindCellUltimate( x );
   if( m < 0 )
   {
     if( debug ) cout << " Calorimeter::SetHitedCellsUltimate This is expected m = " << m << endl;
     p_in.SetHitedCell( -m );
     return;
   }
   else
   {
     cerr << " Calorimeter::SetHitedCellsUltimate This is not expected m = " << m << endl;
     p_in.SetHitedCell( m );
     return;
   }

}

////////////////////////////////////////////////////////////////////////////////

CalorimeterParticle Calorimeter::ConvertParticleDRS2MRS (const CalorimeterParticle &particle) const
{
  CalorimeterParticle p=particle;
  p.SetX( p.GetX()+GetPositionX(), p.GetXerr() );
  p.SetY( p.GetY()+GetPositionY(), p.GetYerr() );
  p.SetZ( p.GetZ()+GetPositionZ(), p.GetZerr() );

  std::pair<bool, double> miscData;
  miscData = p.GetMiscInfo(CalorimeterParticle::UNCORR_X);
  if (miscData.first) p.SetMiscInfo(CalorimeterParticle::UNCORR_X, miscData.second+GetPositionX());
  miscData = p.GetMiscInfo(CalorimeterParticle::UNCORR_Y);
  if (miscData.first) p.SetMiscInfo(CalorimeterParticle::UNCORR_Y, miscData.second+GetPositionY());

  return p;
}

////////////////////////////////////////////////////////////////////////////////

CalorimeterParticle Calorimeter::ConvertParticleMRS2DRS (const CalorimeterParticle &particle) const
{
  CalorimeterParticle p=particle;
  p.SetX( p.GetX()-GetPositionX(), p.GetXerr() );
  p.SetY( p.GetY()-GetPositionY(), p.GetYerr() );
  p.SetZ( p.GetZ()-GetPositionZ(), p.GetZerr() );

  std::pair<bool, double> miscData;
  miscData = p.GetMiscInfo(CalorimeterParticle::UNCORR_X);
  if (miscData.first) p.SetMiscInfo(CalorimeterParticle::UNCORR_X, miscData.second-GetPositionX());
  miscData = p.GetMiscInfo(CalorimeterParticle::UNCORR_Y);
  if (miscData.first) p.SetMiscInfo(CalorimeterParticle::UNCORR_Y, miscData.second-GetPositionY());

  return p;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ConvertDRS2MRS ( const vector <double> &x_in, vector<double> &x_out ) const
{
  assert( x_in.size() == 3);
  x_out = x_in;
  x_out[0] += GetPositionX();
  x_out[1] += GetPositionY();
  x_out[2] += GetPositionZ();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ConvertMRS2DRS ( const vector <double> &x_in, vector<double> &x_out ) const
{
  assert( x_in.size() == 3);
  x_out = x_in;
  x_out[0] -= GetPositionX();
  x_out[1] -= GetPositionY();
  x_out[2] -= GetPositionZ();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CellsInfoClear(int what,int when)
{
  if( what<0 )
  {
    for( size_t i=0; i < CellInfoTypeSize; i++ )
      CellsInfoClear(i,when);

    return;
  }

  if( when<0 )
  {
    CellsInfoClear(what,OLD);
    CellsInfoClear(what,NEW);
    return;
  }

  vector<StatInfo> &v=cells_info[what][when];
  for( vector<StatInfo>::iterator it=v.begin(); it!=v.end(); it++ )
    it->Clear();
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SaveFit2CellInfo(int what)
{
//  bool debug = true;
  bool debug = false;
//   if( what == LED && GetName() == "HC02P1__" ) debug = true;
  if( debug ) cout << " Calorimeter " << GetName() << " SaveFit2CellInfo " << what << endl;
  if( what < 0 )
  {
    SaveFit2CellInfo(CALIB);
    SaveFit2CellInfo(LED  );
    SaveFit2CellInfo(PED  );
    SaveFit2CellInfo(TIME );
    SaveFit2CellInfo(NOISE);
    SaveFit2CellInfo(CLUSTER);
  }
  else if ( what < (int)CellInfoTypeSize )
  {
    bool use_fit_info4calibration = true;
    bool use_raw_fit_info4calibration = true;
//    CellsInfoClear(what, NEW);
    if( debug ) cout << " Calorimeter " << GetName() << " SaveFit2CellInfo " <<
                      what << " size " << fit_info[what].size() <<
		      " options.calib_histo_units_in_gev " << options.calib_histo_units_in_gev << endl;
//     if( what == LED ) cout << " Leds fit !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   " << endl;
    for( size_t it=0; it != NCells(); it++ )
    {
      if(fit_info[what][it].fit_ok)
      {
        if( debug ) cout << " Fit OK! cell " << GetCellName(it) <<
                           " Stat " <<fit_info[what][it].fit_parameters[0].first <<
                           " Mean " << fit_info[what][it].fit_parameters[1].first <<
                           " Sigma " <<fit_info[what][it].fit_parameters[2].first << endl;


	if( what == CALIB && options.calib_histo_units_in_gev && options.calibration_energy > 0 )
	{
	  double facell =1.;
	  if( options.use_cell_calib_correction )
	  {
            facell = individual_calib_factor_[it];
	  }
          cells_info[what][NEW][it].Set(fit_info[what][it].fit_parameters[0].first,
                                      fit_info[what][it].fit_parameters[1].first/(options.calibration_energy*facell),
                                      fit_info[what][it].fit_parameters[2].first/(options.calibration_energy*facell));
    	}
	else
	{
          cells_info[what][NEW][it].Set(fit_info[what][it].fit_parameters[0].first,
                                      fit_info[what][it].fit_parameters[1].first,
                                      fit_info[what][it].fit_parameters[2].first);
    	}
      }
      else
      {
        if( what == CALIB && use_fit_info4calibration && use_raw_fit_info4calibration )
	{

          if( debug == true ) cout << " Fit info is not OK ! Check for special fit for RAW " << endl;
	  if( options.calib_histo_units_in_gev && options.calibration_energy > 0 )
	  {
	    double facell =1.;
	    if( options.use_cell_calib_correction )
	    {
              facell = individual_calib_factor_[it];
	    }
            if(fit_info[RAW][it].fit_ok)
            {
              cells_info[what][NEW][it].Set(fit_info[RAW][it].fit_parameters[0].first,
                                      fit_info[RAW][it].fit_parameters[1].first/(options.calibration_energy*facell),
                                      fit_info[RAW][it].fit_parameters[2].first/(options.calibration_energy*facell));
	    }
	  }
	}
      }
    }
  }
  else cout << " ERROR in SaveFit2CellInfo::" <<
         " Sorry only CALIB,LED,PED,TIME,NOISE and CLUSTER for the moment, but? " << what << endl;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::PrintGeneralInfo( void ) const
{
  cout << " General Info For Calorimeter " << GetName() << endl;
  cout << " Position X = " << GetPositionX() << " Y = " << GetPositionX() << " Z = " << GetPositionZ() << endl;
  cout << " NCells " << NCells() << endl;
  for ( size_t ic=0; ic < NCells(); ic++ )
  {
    cout << " " << ic << " " << GetCellName(ic) <<
            " x = " << cells[ic].GetX() << " y = " << cells[ic].GetY() << " x = " << cells[ic].GetZ() << endl;
  }

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CellsInfoPrint(int what,int when)
{
  if( what<0 )
  {
    CellsInfoPrint(CALIB,when);
    CellsInfoPrint(LED  ,when);
    CellsInfoPrint(PED  ,when);
    CellsInfoPrint(TIME ,when);
    CellsInfoPrint(NOISE,when);
    CellsInfoPrint(CLUSTER,when);
    return;
  }

  if( when<0 )
  {
    CellsInfoPrint(what,OLD);
    CellsInfoPrint(what,NEW);
    return;
  }

  static const char
    *str_when[]={"OLD","NEW"};

  printf("Cells %s %s information for calorimeter %s\n",
             CellInfoNames[what],str_when[when],GetName().c_str());
  cout << cells_info[what][when];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CellsInfoUpdate(CellInfoType what)
{
  if( what!=CALIB && what!=LED && what!=PED  && what!=TIME && what!=NOISE)
    throw Exception("Calorimeter::CellsInfoUpdate():  bad request");

  cells_info[what][OLD] = cells_info[what][NEW];
}

////////////////////////////////////////////////////////////////////////////////

const StatInfo &Calorimeter::GetCellInfo(CellInfoType what, CalibTimeType when, size_t cell) const
{
  if( ( (unsigned)(what) >= CellInfoTypeSize) || ( (unsigned)(when) >= CalibTimeTypeSize ) )
    throw Exception("Calorimeter::GetCellInfo():  unknown requested type of information: what=%d when=%d",what,when);
  if( cell>NCells() )
    throw Exception("Calorimeter::GetCellInfo(): bad cell number %d. NCells=%d",cell,NCells());
  if(cells_info[what][when].size() != NCells() )
    throw Exception("Calorimeter::GetCellInfo():what=%d when=%d cell number %d NCells=%d not initialised ?? ",what,when,cell );

  return cells_info[what][when][cell];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CopyCellsInfo(size_t what, size_t when, size_t what_target, size_t when_target)
{
  if( ( (unsigned)(what) >= CellInfoTypeSize) || ( (unsigned)(when) >= CalibTimeTypeSize ) )
    throw Exception("Calorimeter::GetCellInfo():  unknown requested type of information: what=%d when=%d",what,when);

  if( (what_target >= CellInfoTypeSize) || (when_target >= CalibTimeTypeSize ) )
    throw Exception("Calorimeter::GetCellInfo():  unknown requested type of information: what=%d when=%d",what_target,when_target);

  for( size_t i=0; i<NCells(); i++ )
  {
    cells_info[what_target][when_target][i] = cells_info[what][when][i];
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SetCellInfo(size_t what, size_t when, size_t cell, const StatInfo &ci)
{
  if( ( (unsigned)(what) >= CellInfoTypeSize) || ( (unsigned)(when) >= CalibTimeTypeSize ) )
    throw Exception("Calorimeter::GetCellInfo():  unknown requested type of information: what=%d when=%d",what,when);

  if( cell==size_t(-1) )
  {
    for( size_t i=0; i<NCells(); i++ )
      Calorimeter::SetCellInfo(what,when,i,ci);
    return;
  }

  if( cell >= NCells() )
    throw Exception("Calorimeter::SetCellInfo(): bad cell number %d. NCells=%d",cell,NCells());

  cells_info[what][when].at(cell) = ci;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StoreData(const vector<CellDataRaw>& data)
{
  signals_.clear();
  for (vector<CellDataRaw>::const_iterator it=data.begin(); it!=data.end(); it++)
  {
    size_t address = it->GetCellIdx();
    if( address>=NCells() )
      throw Exception("Reco::Calorimeter::StoreData() Cell Address %d is out of calorimeter size %d",
                       address, NCells() );

    signals_.push_back(*it);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StorebackData(const vector<CalorimeterParticle>& p)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::StorebackData " << GetName() << " option store_back_data " << GetOptions().store_back_data <<  endl;

  if( !GetOptions().store_back_data ) return;
  signals_.clear();
  if( debug ) cout << " Calorimeter::StorebackData debug Try to extract raw data from " << p.size() << " CalorimeterParticles in " << GetName() <<  endl;
  if( p.size() <= 0 ) return;

  double r[NCells()];
  double t[NCells()];
  for (size_t ic=0; ic<NCells(); ic++)
  {
    r[ic]= 0.;
  }

  for( vector<CalorimeterParticle>::const_iterator p_in = p.begin(); p_in != p.end(); p_in++ )
  {
//    if(p_in->HasTime())  p.SetTime( p_in->GetTime(),p_in->GetTimeErr() );
    for( vector< pair< size_t,double> >::const_iterator d = p_in->GetClusterData().begin(); d != p_in->GetClusterData().end(); d++ )
    {
      size_t icell = d->first;
      double ecell = d->second;
      r[icell] += ecell;
      if(p_in->HasTime()) t[icell] += ecell*p_in->GetTime();
    }
  }

  for (size_t ic=0; ic < NCells(); ic++)
  {
    if(r[ic] > 0.)
    {
      CellDataRaw data(ic);
      data.SetEnergy(r[ic]);
      if(t[ic] > 0.) data.SetTime(t[ic]/r[ic]);
      signals_.push_back( data );
    }
  }

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::RecalculateCalibData( void )
{
//   bool debug = true;
  bool debug = false;
//   if(GetName() !="EC00P1__") debug=false;
  if( debug ) cout <<" Calorimeter::RecalculateCalibData " << GetName() << endl;
//  Set new shema since SignalToEnergy function was introduced
  double scale_e = GetScaleCalibration();
  if( scale_e == 1. && ! options.recalc_from_primary_calib ) return;
  for (vector<CellDataRaw>::iterator it=signals_.begin(); it!=signals_.end(); it++) {
    size_t icell = it->GetCellIdx();
    if( icell>=NCells() )
    {
      cerr << " FATAL!! in Reco::Calorimeter::StoreCalibData() Cell Address " << icell <<
                           " is out of calorimeter size " << NCells() << endl;
      exit(1);
    }
    if( icell>=NCells() )
      throw Exception("Reco::Calorimeter::StoreCalibData() Cell Address %d is out of calorimeter size %d",
                       icell, NCells() );
    double data = it->GetEnergy();
    if( options.recalc_from_primary_calib )
    {
      double cf_old = cells_info[CALIB][PRIM][icell].GetMean();
      double amp = it->GetEnergy()/cf_old;   // It might be wrong in case of LED corrections implementred during DST production
      double energy = SignalToEnergy( icell, amp, -1. );
      it->SetEnergy( energy );
    }

//     double corr = 1.;
//     if( options.recalc_from_primary_calib )
//     {
//       double cf_old = cells_info[CALIB][PRIM][icell].GetMean();
//       double cf_new = cells_info[CALIB][OLD][icell].GetMean();
//       corr = cf_new/cf_old;
//     }
//     it->SetEnergy(data*corr*scale_e);
  }
  if( debug ) cout <<" Calorimeter::RecalculateCalibData was indeed recalculated in " << GetName() << endl;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::RecalculateTimeCalibData( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout <<" Calorimeter::RecalculateTimeCalibData " << GetName() <<
                    " reqwest = "  << options.recalc_from_primary_time_calib <<
                                                 " for " <<signals_.size() << " signals " << endl;
  if( !options.recalc_from_primary_time_calib ) return;
  for (vector<CellDataRaw>::iterator it=signals_.begin(); it!=signals_.end(); it++) {
    if( !it->HasTime() ) continue;
    size_t icell = it->GetCellIdx();
    if( icell>=NCells() )
    {
      cerr << " FATAL!! in Reco::Calorimeter::StoreCalibData() Cell Address " << icell <<
                           " is out of calorimeter size " << NCells() << endl;
      exit(1);
    }
    double time = it->GetTime();
    double corr = 1.;
    double time_old = cells_info[TIME][PRIM][icell].GetMean();
    double time_new = cells_info[TIME][OLD][icell].GetMean();
    corr = time_new-time_old;
    if( debug ) cout << " Cell " << icell << " Time old " << time_old <<" Time new " << time_new << " Time corr " << corr << endl;
    it->SetTime( time-corr );
  }
  if( debug ) cout <<" Calorimeter::RecalculateTimeCalibData was indeed recalculated in " << GetName() << endl;
}

////////////////////////////////////////////////////////////////////////////////

std::vector <Reco::CalorimeterParticle> Calorimeter::RecalculateParticles(const std::vector <Reco::CalorimeterParticle> &particles_in  ) const
{
//   bool debug = true;
  bool debug = false;
  std::vector <Reco::CalorimeterParticle> particle_out;
  if( debug ) cout <<" Calorimeter::RecalculateParticles " << GetName() <<
                    " correct energy = "  << GetOptions().recalc_from_primary_calib <<
                      " correct time = "  << GetOptions().recalc_from_primary_time_calib <<
                                         " for " << particles_in.size() << " particles " << endl;
  for( int p=0; p < (int)particles_in.size(); p++)
  {
    Reco::CalorimeterParticle p_out = particles_in[p];
    const std::vector<size_t> &mcells = p_out.GetHitedCells();
    if( mcells.size() <= 0 )
    {
      cout <<" Not expected! in Calorimeter::RecalculateParticles " << GetName() <<" main cell not provided for " << p <<" particle " << endl;
      continue;
    }
    size_t icell = mcells[0];

    if( GetOptions().recalc_from_primary_calib )
    {
      double cf_old = cells_info[CALIB][PRIM][icell].GetMean();
//      double cf_new = GetActualCalibration(icell);
      double cf_new = cells_info[CALIB][OLD][icell].GetMean();
      double corr = cf_new/cf_old;

      double e = p_out.GetE();
      double se = p_out.GetEerr();

      if( debug ) cout << " # " << p << " E= " << p_out.GetE() <<
	                                   " T= " << p_out.GetTime() <<
	                                   " X= " << p_out.GetX() <<
	                                   " Y= " << p_out.GetY() <<
	                                   " Z= " << p_out.GetZ() << endl;
      if( debug ) cout <<  " ecorr= " << corr << endl;
      p_out.SetE( e*corr, se);
//       double sum=0.;
//       const vector < pair< size_t,double> > &clust = p_out.GetClusterData();
//       for( unsigned pc=0; pc<clust.size(); pc++)
//         sum += clust[pc].second;
//
//       cout << " Cluster Sum =" << sum << endl;
    }

    if( GetOptions().recalc_from_primary_time_calib )
    {
      double time_old = cells_info[TIME][PRIM][icell].GetMean();
      double time_new = cells_info[TIME][OLD][icell].GetMean();
      double tcorr = time_new-time_old;
      double time = p_out.GetTime();
      double stime = p_out.GetTimeErr();
      if( debug ) cout << " Cell " << icell << " Time old " << time_old <<" Time new " << time_new << " Time corr " << tcorr << endl;
      p_out.SetTime( time-tcorr, stime );
    }

    if( p_out.GetE() >= GetOptions().particle_energy_threshold ) particle_out.push_back( p_out );

  }
  if( debug ) cout <<" Calorimeter::Recalculateparticles  " << GetName() << " Np_out = " << particle_out.size() << endl;
  return particle_out;

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StoreRawFromCalibData( void )
{
  for (vector<CellDataRaw>::iterator it=signals_.begin(); it!=signals_.end(); it++) {
    size_t address = it->GetCellIdx();
    if( address>=NCells() )
      throw Exception("Reco::Calorimeter::StoreCalibData() Cell Address %d is out of calorimeter size %d",
                       address, NCells() );
    double data = it->GetEnergy();
//    cout << " Calorimeter " << GetName() << " Adress " << address << " data " << data <<
//                                           " Calib " << cells_info[CALIB][OLD][address].GetMean() << endl;
    it->SetAmplitude( data/cells_info[CALIB][OLD][address].GetMean() );
  }
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::EnergyDispersionInCellReadout(double energy_in_cell,size_t jcell) const
{
  double ro = cells[jcell].GetCellType().GetReadoutTerm();
  return ro*ro;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::ScaleDispersionInSparceMode(double e, double dispe, size_t jcell) const
{
  double de = dispe;
//  it seems we have options doubling
//  if( options.rd_fiadc_digitization &&  dispe > 0. )
  if( GetOptions().readout_sparsified_ &&  dispe > 0. )
  {
//  it seems we have options doubling
//    double e_sparce_cut = options.rd_fiadc_sparce_delta*cells_info[CALIB][OLD][jcell].GetMean();
    double e_sparce_cut = GetEnergyCutSparseMode()[jcell];
    if( e > e_sparce_cut )
    {
      double edev = fabs( e - e_sparce_cut );
      edev = edev*edev/dispe;
      if( edev < 1. )
      {
        de = dispe/2.5;
//        de = 2.*de; // 50
      }
      else if( edev < 4. )
      {
        de = dispe/1.8;
      }
      else if( edev < 6. )
      {
        de = dispe/1.4;
      }
    }
  }
  return de;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::EnergyDispersionInCellDigitization(double energy_in_cell,size_t jcell) const
{
  double rs = 0.;
  double e_sparce_cut = cells_info[CALIB][OLD][jcell].GetMean()/0.5;
//  it seems we have options doubling
//  if( options.rd_fiadc_digitization )
  if( GetOptions().readout_sparsified_ )
  {
//  it seems we have options doubling
//    e_sparce_cut = (0.5+options.rd_fiadc_sparce_delta)*cells_info[CALIB][OLD][jcell].GetMean();
    e_sparce_cut = GetEnergyCutSparseMode()[jcell];
    if( energy_in_cell < e_sparce_cut ) rs = e_sparce_cut;
  }
  double ro =  cells_info[CALIB][OLD][jcell].GetMean()/3.;
  double factor = 2.00/(e_sparce_cut/0.06);  //   150 GeV
// 100 GeV
// factor = factor/4.; // e_sparce_cut = 6.*0.6  at 50. GeV
  factor = factor/4.; // e_sparce_cut = 6.*0.6  at 5. GeV


  return (ro*ro+rs*rs)*factor;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::EnergyDispersionInCellBasic(double energy_in_cell,size_t jcell) const
{
  double StochasticTerm = cells[jcell].GetCellType().GetStochasticTerm();
  double ConstantTerm = cells[jcell].GetCellType().GetConstantTerm();
  return energy_in_cell*StochasticTerm*StochasticTerm +
           energy_in_cell*energy_in_cell*(ConstantTerm*ConstantTerm + options.calibration_uncertainty*options.calibration_uncertainty);

}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::EnergyDispersionInCell(double energy_in_cell,size_t jcell) const
{
  return EnergyDispersionInCellBasic(energy_in_cell, jcell) + EnergyDispersionInCellReadout(energy_in_cell, jcell);
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::EnergySigmaInCellReadout(double energy_in_cell,size_t jcell) const
{
  return sqrt( EnergyDispersionInCellReadout(energy_in_cell, jcell) );
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::EnergySigmaInCellBasic(double energy_in_cell,size_t jcell) const
{
  return sqrt( EnergyDispersionInCellBasic(energy_in_cell, jcell) );
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::EnergySigmaInCell(double energy_in_cell,size_t jcell) const
{
  return sqrt( EnergyDispersionInCell(energy_in_cell, jcell) );
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::StoreRawDataStatistic(void)
{
  // Check consistency
  if( cells_info[RAW][NEW].size() != NCells() )
  {
    cerr << " StoreRawDataStatistic " << GetName() << " Bad size of storage " <<  cells_info[RAW][NEW].size() << endl;
    exit(1);
  }

  // Add raw amplitude to cells info
  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++) {
    // Get cell info
    size_t icell = it->GetCellIdx();
    //Get amplitude
    double ecell = it->GetEnergy();
    // Check range
    if( icell>=NCells() )
      throw Exception("Reco::Calorimeter::InsertData() Cell Address %d is out of calorimeter size %d",
                       icell, NCells() );
    // Add amplitude to raw data
    cells_info[RAW][NEW][icell].Add(ecell);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ResetStatistic(void)
{
  StatInfo zero = StatInfo(0.,0.,0.);
  for( size_t icell=0; icell < NCells(); icell++ )
  {
    if( cells_info[RAW][NEW].size() == NCells() ) cells_info[RAW][NEW][icell] = zero;
    if( cells_info[CALIB][NEW].size() == NCells() )cells_info[CALIB][NEW][icell] = zero;
    if( cells_info[LED][NEW].size() == NCells() )cells_info[LED][NEW][icell] = zero;
    if( cells_info[PED][NEW].size() == NCells() )cells_info[PED][NEW][icell] = zero;
  }
}

////////////////////////////////////////////////////////////////////////////////

vector<size_t> Calorimeter::BadCellsAround( size_t icell) const
{
  vector<size_t> cnt;
  for( vector<size_t>::const_iterator it=cells[icell].GetNeighbors().begin(); it!=cells[icell].GetNeighbors().end(); it++ )
  {
//     if( cell_is_bad_old[*it] > 0 ) cnt.push_back(*it);
    if( cell_is_bad_[OLD][*it] > 0 ) cnt.push_back(*it);
  }
  return cnt;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::CellIsBad(  size_t icell, size_t when) const
{
  for( list <size_t>::const_iterator it = bad_cells_[when].begin(); it!=bad_cells_[when].end(); it++ )
  {
    if( icell == *it ) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::GetBadCellStatus(  size_t icell, CalibTimeType when) const
{
// TODO: checks on validity of input

  return status_bad_cells_[when][icell];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FilterOutBadCells(const vector<CellDataRaw>& data_in, vector<CellDataRaw>& data_out, size_t when) const
{
//   bool debug = true;
  bool debug = false;

  if( &data_in == &data_out )
  {
    cerr << " Error Calorimeter::FilterOutBadCells can not filter at the same container " << endl;
    exit(1);
  }
  if( debug ) cout << " Calorimeter::FilterOutBadCells debug NbadCells = " << bad_cells_[when].size() << endl;
  if( debug ) cout << " Calorimeter::FilterOutBadCells debug NbadCells OLD = " << bad_cells_[OLD].size() << endl;
  if( debug ) cout << " Calorimeter::FilterOutBadCells debug NbadCells MC = " << bad_cells_[MC].size() << endl;

  data_out.clear();
  for (vector<CellDataRaw>::const_iterator it=data_in.begin(); it!=data_in.end(); it++)
  {
    if (!CellIsBad(it->GetCellIdx(), when))
      data_out.push_back(*it);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ResetMonitorLEDInBurstHisto ( void )
{
  if( p1_led_in_spills_.size() == 0 ) return;
  if( p1_led_in_spills_.size() != NCells() ) assert(false);
  for( unsigned icell=0; icell < NCells(); icell++ )
  {
    p1_led_in_spills_[icell]->Reset();
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::MonitorLEDInBurst( int spill )
{
//   bool debug = true;
  bool debug = false;

  if( debug ) cout << " Calorimeter::MonitorLEDInBurst " << GetName() << endl;
  if( p1_led_in_spills_.size() == 0 )
  {
    bool ok=false;
    TDirectory *dir_save=gDirectory; // Save current ROOT directory.
    char dir_name[111];
    ok = gDirectory->cd("/");
    assert(ok);
    sprintf(dir_name,"Calorimeter_%s",GetName().c_str());
    TDirectory *root_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok = root_dir->cd();
    assert(ok);
    sprintf(dir_name,"MonitorLED");
    TDirectory *mon_dir = myROOT_utils::TDirectory_checked(dir_name);
    ok=mon_dir->cd();
    for( unsigned icell=0; icell < NCells(); icell++ )
    {
      char hname[132];
      char hist_name[132];
      sprintf(hname,"MonitorLED_%s", GetCellName(icell).c_str());
      sprintf(hist_name, "MonitorLEDinSpills in Cell %s ", GetCellName(icell).c_str() );
      TProfile *p1 = myROOT_utils::TProfile_checked(hname, hist_name, 250, -0.5, -0.5+250. );
      p1_led_in_spills_.push_back(p1);
    }
    c_led_in_spills_=NULL;
    ok=dir_save->cd();
    assert(ok);
  }

  if( p1_led_in_spills_.size() != NCells() ) assert(false);

  for( unsigned icell=0; icell < NCells(); icell++ )
  {
    p1_led_in_spills_[icell]->Fill( spill, cells_info[LED][MONITOR][icell].GetMean() );
//    cells_info[LED][MONITOR][icell]=cells_info[LED][NEW][icell];
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ProcLED(void)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout <<" Calorimeter::ProcLED debug " << GetName() <<" DataSize " << signals_.size() <<
                                                " options.store_leds_all " << options.store_leds_all << endl;

// ??? Wild Monitoring activity
  for( unsigned icell=0; icell < NCells(); icell++ )
  {
    cells_info[LED][MONITOR][icell]=StatInfo(0.,0.,0.);
  }

  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++ )
  {
    size_t icell = it->GetCellIdx();
    double amp_led = it->GetAmplitude();

//    cout << " Adress " << icell << " Amplitude " << amp_led << endl;
    cells_info[LED][NEW][icell].Add(amp_led,1.);

    cells_info[LED][MONITOR][icell] = StatInfo(1.,amp_led, 1.);
  }

  FillHistoLED();
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::MakeInspections( void )
{
  int status = 0;
  int status_max = 0;
  status = InspectCalibration();
  if(status > status_max )  status_max = status;
  status = InspectCellsNoise();
  if(status > status_max )  status_max = status;
  status = InspectCellsGamNoise();
  if(status > status_max )  status_max = status;
  status = InspectCellsNoise4RandomTrigger();
  if(status > status_max )  status_max = status;
  status = InspectInternalCorrelations();
  if(status > status_max )  status_max = status;
  status = InspectExternalCorrelations();
  if(status > status_max )  status_max = status;
  status = InspectLED();
  if(status > status_max )  status_max = status;
//   InspectPedestal() is obsolet ???
//   status = InspectPedestal();
//   if(status > status_max )  status_max = status;
  return status_max;                            // Yee, Just Fun ! Don't care this has no sense
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectLED( void )
{
  bool debug=false;
//   bool debug=true;
  if( options.led_statistic_inspect < 0.1 ) return 0;
  SummaryInspectLED summary;
  led_statistic_++;
  if( options.led_statistic_inspect > led_statistic_ ) return 0;
  double led_statistic = led_statistic_;
  led_statistic_ = 0;

  bool print_every_channel=false;
  bool print_summary_statistic=true;
  double led_min_value=25;

  double stat_max=0;
  for( size_t icell=0; icell < NCells(); icell++ )
  {
    double stat  = cells_info[LED][NEW][icell].GetEntries();
    if(stat > stat_max )  stat_max = stat;
  }
  if(debug) cout << " **** CALORIMETER " << GetName() << " InspectLED STATISTIC Total=" << stat_max << endl;


// ??   if( stat_max < led_statistic )
// ??   {
// ??     cout << " InspectLED()::ERROR Terminate due to low total statistic " << endl;
// ??   }

  int no_ref = 0;
  int no_entries = 0;
  int low_statistic = 0;
  int small_value = 0;
  for( size_t icell=0; icell < NCells(); icell++ )
  {
    bool cell_ok = true;
    double stat  = cells_info[LED][NEW][icell].GetEntries();
    double value = cells_info[LED][NEW][icell].GetMean();
    double ref_value = cells_info[LED][OLD][icell].GetMean();
    double tolerance = 0.2;

    if( stat == 0 )  // This cell has no led signals !
    {
      if(print_every_channel) cout << GetCellName(icell) << " no LED statistic " << endl;
      no_entries++;
      continue;
    }

    if( stat < led_statistic/10 )
    {
      low_statistic++;
      continue;
    }

    if( value < led_min_value )
    {
      small_value++;
    }


    if( stat < led_statistic )
    {
      cell_ok = false;
    }

    if( cell_ok ) summary.cells_ok_++;
    if( ref_value > 0 )
    {
       double delta = value/ref_value - 1.;
       if( fabs(delta) > tolerance ) cell_ok = false;

    }
    else
    {
       no_ref++;
       cell_ok = false;
    }

    if( !cell_ok ) if(print_every_channel) cout << " Cell " << GetCellName(icell) << " LED " << (int)value << " REF " << (int)ref_value << endl;

  }


  summary.stat_max_ = (int)stat_max;
  summary.no_entries_ = no_entries;
  summary.low_statistic_ = low_statistic;
  summary.small_value_ = small_value;
  summary.init_ = true;

  if(print_summary_statistic)
  {
//     cout << " **** CALORIMETER " << GetName() << " InspectLED STATISTIC Total=" << stat_max <<
//             " Dead cells " << no_entries << " Low statistic " << low_statistic << " Small LEDs " << small_value <<endl;
    cout << " **** CALORIMETER " << GetName() << endl;
    summary.Print();

  }

  prev_summary_led_ = summary_led_;
  summary_led_ = summary;

  if(options.reset_after_led_inspect)
  {
    for ( size_t i=0; i<NCells(); i++ )
    {
      cells_info[LED  ][NEW][i] = StatInfo(0.,0.,0.);
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::MakeInspectionsDB( void ) const
{
  int status = 0;
  int status_max = 0;
  status = InspectCalibrationDB();
  if(status > status_max )  status_max = status;
//   status = InspectCellsNoiseDB();
//   if(status > status_max )  status_max = status;
  status = InspectLEDDB();
  if(status > status_max )  status_max = status;
//   status = InspectPedestalDB();
//   if(status > status_max )  status_max = status;
  status = InspectPorogiDB();
  if(status > status_max )  status_max = status;
  return status_max;                            // Yee, Just Fun ! Don't care this has no sense
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectCalibrationDB(void) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " InspectCalibrationDB " << GetName() << " debug " << endl;
  int printout_level = options.monitor_db_calib_level;
  int n_low_mean = 0;
  int n_high_mean = 0;
//   double mean_ref = 0.005;
//   double mean_min = mean_ref/5.;
//   double mean_max = mean_ref*5.;
  double mean_min = options.monitor_db_calib_min_value;
  double mean_max = options.monitor_db_calib_max_value;
  if( printout_level >= 0 )
  {
    cout << " **** CALORIMETER " << GetName() << " InspectCalibrationDB " << endl;
  }

  for( size_t icell=0; icell < NCells(); icell++ )
  {
    double mean  = cells_info[CALIB][OLD][icell].GetMean();
//     double stat  = cells_info[CALIB][OLD][icell].GetEntries();
    if(mean>0)
    {
      if(  mean < mean_min )
      {
        if( printout_level >= 2 )
        {
          cout << " Small cf in cell " << icell << " " << GetCellName(icell) << " cf = " << mean << endl;
        }
        n_low_mean++;
      }
      if(  mean > mean_max )
      {
        if( printout_level >= 2 )
        {
          cout << " Big cf in cell " << icell << " " << GetCellName(icell) << " cf = " << mean << endl;
        }
        n_high_mean++;
      }
      if( debug ) cout << " cell " << icell << " " << GetCellName(icell) << " cf = " << mean << endl;
    }
    else
    {
      if( debug ) cerr << " BAD cf in cell " << icell << " " << GetCellName(icell) << " cf = " << mean << endl;
    }
  }

  if( printout_level >= 0 )
  {
    cout << " Summary: Ncells (cf < " << mean_min << " ) = " << n_low_mean <<
                      " Ncells (cf > " << mean_max << " ) = " << n_high_mean << endl;
    if( printout_level > 0  )
    {
    }
    cout << endl;
  }
  if( debug ) cout << " InspectCalibrationDB " << GetName() << " debug  OK " << endl;
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectPorogiDB(void) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " InspectPorogiDB " << GetName() << " debug " << endl;
  int printout_level = options.monitor_db_porog_level;
  int n_low_mean = 0;
  int n_high_mean = 0;
//   double mean_ref = 0.005;
//   double mean_min = mean_ref/5.;
//   double mean_max = mean_ref*5.;
  double mean_min = options.monitor_db_porog_min_value;
  double mean_max = options.monitor_db_porog_max_value;
  if( printout_level >= 0 )
  {
    cout << " **** CALORIMETER " << GetName() << " InspectPorogiDB LOW thresholds inspection " << endl;
  }

  for( size_t icell=0; icell < NCells(); icell++ )
  {
    double mean  = energy_cut_bad_cells_old[icell];
//    double mean  = energy_gamma_cut_bad_cells_old[icell];
    if(mean>0)
    {
      if(  mean < mean_min )
      {
        if( printout_level >= 2 )
        {
          cout << " Small threshol in cell " << icell << " " << GetCellName(icell) << " thr = " << mean << endl;
        }
        n_low_mean++;
      }
      if(  mean > mean_max )
      {
        if( printout_level >= 2 )
        {
          cout << " Big threshol in cell " << icell << " " << GetCellName(icell) << " thr = " << mean << endl;
        }
        n_high_mean++;
      }
      if( debug ) cout << " cell " << icell << " " << GetCellName(icell) << " thr = " << mean << endl;
    }
    else
    {
      if( debug ) cerr << " BAD  threshol in cell " << icell << " " << GetCellName(icell) << " cf = " << mean << endl;
    }
  }

  if( printout_level >= 0 )
  {
    cout << " Summary: Ncells (cell thr < " << mean_min << " ) = " << n_low_mean <<
                      " Ncells (cell thr > " << mean_max << " ) = " << n_high_mean << endl;
    if( printout_level > 0  )
    {
    }
    cout << endl;
  }


  n_low_mean = 0;
  n_high_mean = 0;
  mean_min = options.monitor_db_porog_min_value;
  mean_max = options.monitor_db_porog_max_value;
  if( printout_level >= 0 )
  {
    cout << " **** CALORIMETER " << GetName() << " InspectPorogiDB HIGH thresholds inspection " << endl;
  }

  for( size_t icell=0; icell < NCells(); icell++ )
  {
//     double mean  = energy_cut_bad_cells_old[icell];
    double mean  = energy_gamma_cut_bad_cells_old[icell];
    if(mean>0)
    {
      if(  mean < mean_min )
      {
        if( printout_level >= 2 )
        {
          cout << " Small threshol in cell " << icell << " " << GetCellName(icell) << " thr = " << mean << endl;
        }
        n_low_mean++;
      }
      if(  mean > mean_max )
      {
        if( printout_level >= 2 )
        {
          cout << " Big threshol in cell " << icell << " " << GetCellName(icell) << " thr = " << mean << endl;
        }
        n_high_mean++;
      }
      if( debug ) cout << " cell " << icell << " " << GetCellName(icell) << " thr = " << mean << endl;
    }
    else
    {
      if( debug ) cerr << " BAD  threshol in cell " << icell << " " << GetCellName(icell) << " cf = " << mean << endl;
    }
  }

  if( printout_level >= 0 )
  {
    cout << " Summary: Ncells (gamma thr < " << mean_min << " ) = " << n_low_mean <<
                      " Ncells (gamma thr > " << mean_max << " ) = " << n_high_mean << endl;
    if( printout_level > 0  )
    {
    }
    cout << endl;
  }
  if( debug ) cout << " InspectCalibrationDB " << GetName() << " debug  OK " << endl;
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::InspectLEDDB(void) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " InspectLEDDB " << GetName() << " debug " << endl;

  int printout_level = options.monitor_db_led_level;
  int n_low_mean = 0;
  int n_high_mean = 0;
  double mean_min = options.monitor_db_led_min_value;
  double mean_max = options.monitor_db_led_max_value;
  if( printout_level >= 0 )
  {
    cout << " **** CALORIMETER " << GetName() << " InspectLEDDB " << endl;
  }

  for( size_t icell=0; icell < NCells(); icell++ )
  {
    double mean  = cells_info[LED][OLD][icell].GetMean();
    if(mean >0)
    {
      if( debug ) cout << " cell " << icell << " " << GetCellName(icell) << " led = " << mean << endl;
      if(  mean < mean_min )
      {
        if( printout_level >= 2 )
        {
          cout << " Small led in cell " << icell << " " << GetCellName(icell) << " led = " << mean << endl;
        }
        n_low_mean++;
      }

      if(  mean > mean_max )
      {
        if( printout_level >= 2 )
        {
          cout << " Big led in cell " << icell << " " << GetCellName(icell) << " led = " << mean << endl;
        }
        n_high_mean++;
      }

    }
    else
    {
      if( debug ) cerr << " NO led in cell " << icell << " " << GetCellName(icell) << " led = " << mean << endl;
    }
  }

  if( printout_level >= 0 )
  {
    cout << " Summary: Ncells (led < " << mean_min << " ) = " << n_low_mean <<
                      " Ncells (led > " << mean_max << " ) = " << n_high_mean << endl;
    if( printout_level > 0)
    {
    }
    cout << endl;
  }
  if( debug ) cout << " InspectLEDDB " << GetName() << " debug  OK " << endl;
  return 0;

}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::CheckCalibData(const std::vector<CellDataRaw>& data, double threshold, double sum_energy_cut, int nmax_cells_cut) const
{
  double sum = 0.;
  int n=0;
  for( size_t it=0; it < data.size(); it++ )
  {
    if( data[it].GetEnergy() > threshold )
    {
      sum += data[it].GetEnergy();
      n += 1;
    }
  }
  bool ok = true;
  if( sum > sum_energy_cut ) ok = false;
  if( n > nmax_cells_cut ) ok = false;
  if( !ok && options.print_bad_data )
  {
    cerr << GetName() << " Bad Data !!!! ESum = " << sum << " N cells " << n << endl;
  }
  return ok;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::CheckPedestalIsOK( const CellDataRaw& data ) const
{
  double amp_ped = data.GetAmplitude();
  if( amp_ped > 511 ) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::RecoverBadPed (void)
{
  for (vector<CellDataRaw>::iterator it=signals_.begin(); it!=signals_.end(); it++) {
    size_t cell = it->GetCellIdx();
    if(bad_ped_[cell] != 0)
    {
      double data = it->GetAmplitude()-511;
      it->SetAmplitude(data);
      it->SetEnergy(data*cells_info[CALIB][OLD][cell].GetMean());
//       cout << " Bad pedestal in cell " << cell << " Ampl " << data << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ProcPedestal(void)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Debug ProcPedestal for calorimeter " << GetName() << endl;
  if( debug ) cout << " bad_ped_cells_.size() = " << bad_ped_cells_.size() << endl;
  for ( list <size_t>::iterator it=bad_ped_cells_.begin(); it!=bad_ped_cells_.end(); it++ )
  {
    if( debug ) cout << " cleaning bad_ped_[" << *it <<"] " << endl;
    bad_ped_[*it]=0;
  }
  bad_ped_cells_.clear();

  if( debug ) cout << " signals_.size() = " << signals_.size() << endl;
  for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++) {
    size_t icell = it->GetCellIdx();
    double amp_ped = it->GetAmplitude();
    if( debug ) cout << " Adress " << icell << " Amplitude " << amp_ped << endl;
    if( icell>=NCells() )
      throw Exception("Reco::Calorimeter::ProcPedestal() Cell Address %d is out of calorimeter size %d",
                       icell, NCells() );
    cells_info[PED][NEW][icell].Add(amp_ped,1.);
    if( debug ) cout << " Try to CheckPedestalIsOK " << endl;
    if( !CheckPedestalIsOK(*it) )
    {
      cout << " BIG !!! Pedestal in " << GetName() << " cell= " << GetCellName(icell)
                                           << " ped= " << (int)amp_ped  << endl;
      if(bad_ped_[icell] != 0)
      {
        cerr << " Calorimeter::ProcPedestal Bad DATA the same cell " << icell << endl;
        continue;
      }
      bad_ped_cells_.push_back(icell);
      bad_ped_[icell]=(int)amp_ped;
    }
  }
  if( debug ) cout << " Try to  FillHistoPED " << endl;

  FillHistoPED();
  if( debug ) cout << " Debug ProcPedestal for calorimeter " << GetName() << " is OK !!!!!!! " << endl;

}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SubstractPED(vector<CellDataRaw>& data) const
{
  for (vector<CellDataRaw>::iterator it=data.begin(); it!=data.end(); it++)
  {
    size_t icell = it->GetCellIdx();
    if( icell>=NCells() )
      throw Exception("Reco::Calorimeter::SubstractPED() Cell Address %d is out of calorimeter size %d",
                       icell, NCells() );
    it->SetAmplitude(it->GetAmplitude() - cells_info[PED][OLD][icell].GetMean());
  }
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCalibrationFactorBase(size_t cell ) const
{
  return cells_info[CALIB][OLD][cell].GetMean() * options.scale_calibration;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCalibrationFactorLED(size_t cell ) const
{
//   bool debug = true;
  bool debug = false;
  double cf=1.;
  if( options.use_led_ref_calib_correction )
  {
    if( debug ) cout << " start " << GetName() << " cell " << cell
		     << " LED PRIM " << cells_info[LED][PRIM][cell].GetMean()
		     << " LED OLD " <<  cells_info[LED][OLD][cell].GetMean() << endl;

    if( cells_info[LED][PRIM][cell].GetMean() > 10. )
      {
        double current_led = cells_info[LED][OLD][cell].GetMean();
        if( options.use_led_ref_calib_correction_in_spills )
	  {
	    double current_led_in_spill = GetInCurrentSpillMeanLED( cell );
	    if( current_led_in_spill > 0. ) current_led = current_led_in_spill;
	  }

        if( current_led > 10. )
	  {
	    double cfled = cells_info[LED][PRIM][cell].GetMean() / current_led;
	    if( cfled > 0.05 && cfled < 20. )
	      {
		if( debug ) cout << GetName() << " cell " << cell << " LED correction " << cfled << endl;
		cf *= cfled;
	      }
	  }
      }
  }
  return cf;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCalibrationFactorEdep(size_t cell, double energy ) const
{
  if ( options.edep_corr != "" )
    return cells[cell].GetEdepCorr().Interpolate(energy);
  else
    return 1.;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetCalibrationFactorSdep(size_t cell, double tis ) const
{
  if ( options.tisdep_corr != "" && tis != -1. )
    return cells[cell].GetTiSdepCorr().Interpolate(tis);
  else
    return 1.;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetECal0_nonlinear_calib_correction(size_t cell, double signal) const
{
  if ( options.ecal0_nonlinear_calib_correction) {

  double p0;

  if(cell<54) { p0=4200;
  } else if(cell>=54 && cell<63){ p0=4100; }
  else { p0=6300;}

  double ttt=(1-signal/p0);
  signal = (-p0)*log(ttt);
}
return signal;
}


double Calorimeter::SignalToEnergy(size_t cell, double signal, double tis) const
{
//   bool debug = true;
  bool debug = false;
  if(GetName() != "EC00P1__")debug=false;
  double cf = GetCalibrationFactorBase(cell);

  if( options.ecal0_nonlinear_calib_correction) {
     signal = GetECal0_nonlinear_calib_correction(cell, signal );
     }

  double energy = cf*signal;
  if( debug ) cout << " Calorimeter::GetActualCalibration " << GetName()
		   << " options.use_led_ref_calib_correction "
		   << options.use_led_ref_calib_correction << endl;

  if( options.use_led_ref_calib_correction )
    energy *= GetCalibrationFactorLED(cell);

  if ( options.edep_corr != "" )
    energy *= GetCalibrationFactorEdep(cell, energy);

  if ( options.tisdep_corr != "" && tis != -1. )
    energy *= GetCalibrationFactorSdep(cell, tis);

  return energy;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::DetectCellNeighbors(void)
{
//  bool debug = true;
  bool debug = false;

  if( debug ) cout << " DetectCellNeighbors " << GetName() << " eps= " << options.tolerance_for_nearby_cells << endl;
  for( size_t it=0; it<NCells(); it++ )
    cells[it].GetNeighbors().clear();

  for( size_t it=0; it<NCells(); it++ )
    cells[it].GetNeighbors().push_back(it);

  if( NCells()>1 )
    for( size_t it1=0; it1<NCells()-1; it1++ )
      for( size_t it2=it1+1; it2<NCells(); it2++ )
      {
        if( cells[it1].IsNeighbor(cells[it2], options.tolerance_for_nearby_cells) )
        {
          cells[it1].GetNeighbors().push_back(it2);
          cells[it2].GetNeighbors().push_back(it1);
        }
      }


  for ( size_t i=0; i<NCells(); i++ )
  {
    if( cells[i].GetNeighbors().size() <= 1 )
    {
      cout << "Sorry the code was not tested in case the cell has no neighbors "
	   << "it might be FATAL! The program terminated! " << endl;
      exit(1);
    }
  }

  if( debug )
  {
    cout << " DetectCellNeighbors " << GetName() << endl;
    for ( size_t i=0; i<NCells(); i++ )
    {
      cout << " cell " << i;
      cout << " size " << cells[i].GetCellType().GetSizeX();
      cout << " x " << cells[i].GetX();
      cout << " y " << cells[i].GetY();
      cout << " Neighbors " << cells[i].GetNeighbors().size() << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::DetectCellsOnBoundary(void)
{
  for( vector<Cell>::iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double xc = it->GetX();
    double yc = it->GetY();
    double zc = it->GetZ();
    double sxc = it->GetCellType().GetSizeX();
    double syc = it->GetCellType().GetSizeY();
    for( int ix=-1; ix < 2; ix++)
      for( int iy=-1; iy < 2; iy++)
      {
        if( !(ix == 0 && iy == 0) )
        {
          double x = xc + (double)(ix)*sxc;
          double y = yc + (double)(iy)*syc;
	  vector <double> xx;
	  xx.push_back(x);
	  xx.push_back(y);
	  xx.push_back(zc);
	  vector <double> x_out;
	  ConvertDRS2MRS(xx, x_out);
//          if( !InActiveArea(x+GetPositionX(),y+GetPositionY(),zc+GetPositionZ()) ) // This is a mess! InActiveArea(x,y,zc) method require coordinates in MRS !
          if( !InActiveArea(x_out) ) // This is a mess! InActiveArea(x,y,zc) method require coordinates in MRS !
          {
//            it->SetBoundaryRegion(x,y,zc);
            it->SetBoundaryRegion(x_out);
          }
        }
      }
  }
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetXmin(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetXmin():  there are no cells "
		    "in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetX() - cells[0].GetCellType().GetSizeX()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetX() - it->GetCellType().GetSizeX()/2;
    if( v < m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetTrueXmin(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetXmin():  there are no cells "
		    "in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetX() - cells[0].GetCellType().GetTrueSizeX()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetX() - it->GetCellType().GetTrueSizeX()/2;
    if( v < m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetXmax(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetXmax():  there are no cells "
                    "in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetX() + cells[0].GetCellType().GetSizeX()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetX() + it->GetCellType().GetSizeX()/2;
    if( v > m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetTrueYmin(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetYmin():  there are no cells "
		    "in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetY() - cells[0].GetCellType().GetTrueSizeY()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetY() - it->GetCellType().GetTrueSizeY()/2;
    if( v < m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetYmin(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetYmin():  there are no cells "
		    "in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetY() - cells[0].GetCellType().GetSizeY()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetY() - it->GetCellType().GetSizeY()/2;
    if( v < m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetYmax(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetYmax():  there are no cells "
		    "in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetY() + cells[0].GetCellType().GetSizeY()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetY() + it->GetCellType().GetSizeY()/2;
    if( v > m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetZmin(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetZmin():  there are no cells "
		    "in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetZ() - cells[0].GetCellType().GetSizeZ()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetZ() - it->GetCellType().GetSizeZ()/2;
    if( v < m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetZmax(void) const
{
  if( cells.empty() )
    throw Exception("Reco::Calorimeter::GetZmax():  there are no "
		    "cells in calorimeter \"%s\"", GetName().c_str());

  double m =  cells[0].GetZ() + cells[0].GetCellType().GetSizeZ()/2;
  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
    double v = it->GetZ() + it->GetCellType().GetSizeZ()/2;
    if( v > m )
      m = v;
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetMinCellSizeX(void) const
{
  double m = 1.e+10;
  list<CellType>::const_iterator cell_type;
  for( cell_type=cells_type.begin(); cell_type!=cells_type.end(); cell_type++ )
  {
    if( cell_type->GetSizeX() < m  ) m = cell_type->GetSizeX();
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetMaxCellSizeX(void) const
{
  double m = 1.e-10;
  list<CellType>::const_iterator cell_type;
  for( cell_type=cells_type.begin(); cell_type!=cells_type.end(); cell_type++ )
  {
    if( cell_type->GetSizeX() > m  ) m = cell_type->GetSizeX();
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetMinCellSizeY(void) const
{
  double m = 1.e+10;
  list<CellType>::const_iterator cell_type;
  for( cell_type=cells_type.begin(); cell_type!=cells_type.end(); cell_type++ )
  {
    if( cell_type->GetSizeY() < m  ) m = cell_type->GetSizeX();
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetMaxCellSizeY(void) const
{
  double m = 1.e-10;
  list<CellType>::const_iterator cell_type;
  for( cell_type=cells_type.begin(); cell_type!=cells_type.end(); cell_type++ )
  {
    if( cell_type->GetSizeY() > m  ) m = cell_type->GetSizeX();
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetMinCellSizeZ(void) const
{
  double m = 1.e+10;
  list<CellType>::const_iterator cell_type;
  for( cell_type=cells_type.begin(); cell_type!=cells_type.end(); cell_type++ )
  {
    if( cell_type->GetSizeZ() < m  ) m = cell_type->GetSizeX();
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetMaxCellSizeZ(void) const
{
  double m = 1.e-10;
  list<CellType>::const_iterator cell_type;
  for( cell_type=cells_type.begin(); cell_type!=cells_type.end(); cell_type++ )
  {
    if( cell_type->GetSizeZ() > m  ) m = cell_type->GetSizeX();
  }
  return m;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::FindCell(const CalorimeterParticle &p) const
{
  vector <double> vp;
  vp.push_back(p.GetX());
  vp.push_back(p.GetY());
  vp.push_back(p.GetZ());
  return FindCell(vp);
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::FindCellInternal(double x,double y,double z) const
{
  vector <double> vp;
  vp.push_back(x);
  vp.push_back(y);
  vp.push_back(z);
  return FindCell(vp);
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::FindCell(const vector <double> &x) const
{
//   bool debug=true;

  for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
  {
//     if( it->InActiveArea(x,y,z) ) return it-cells.begin();
    if( it->InActiveArea(x) ) return it-cells.begin();
  }
//   if( debug )
//   {
//     if( GetName() == "GS" )
//     {
//
//       cout << " MMissed in GS x=" << x[0] << " y=" << x[1] << " z=" << x[2] << endl;
// //       for( vector<Cell>::const_iterator it=cells.begin(); it!=cells.end(); it++ )
// //       {
// //         double xc = it->GetX();
// //         double yc = it->GetY();
// //         double zc = it->GetZ();
// // //        cout << " dx=" << x-xc << " dy=" << y-yc << " dz=" << z-zc << endl;
// //       }
//     }
//  }
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::FindCellUltimate(const vector <double> &x) const
{
  double distmin = 1e+10;
  int nearest_cell = -1;
  for( size_t it=0; it<cells.size(); it++ )
  {
    double dist = cells[it].Distance(x);
    if( dist < 0. ) return it;
    if( dist < distmin )
    {
      nearest_cell = it;
      distmin = dist;
    }
  }
  if( nearest_cell < 0 )
  {
    cerr << " Calorimeter::FindCellUltimate internal error in " << GetName() <<
              " At x= " << x[0] <<" y= " << x[1] << " z= " << x[2] << endl;
  }
  return -nearest_cell;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::Associated(const CalorimeterParticle &pe, const std::vector<double>& vpar, int syst, double gate ) const
{
  if( pe.GetMainCells().size() <= 0 ) return false; // Particle missed the Calorimeter???

  double Nsigma = gate;

  assert(vpar.size() >= 5);
  double zcal0 = GetPositionZ();

//  double eg = pe.GetE();
  double xcalc = pe.GetX();
  double ycalc = pe.GetY();
  double zcalc = pe.GetZ();

  double excal = pe.GetXerr();
  double eycal = pe.GetYerr();

  double ztrk = vpar[0]; // Z goes firts!!!! drrrr!!!! brrr!!!
  double xtrk = vpar[1];
  double ytrk = vpar[2];

  const static double delta = 0.01;

  if( fabs(ztrk-zcal0) > delta )
  {
    cerr << " Calorimeter::Associated " << GetName() << " oshibochka ztrk = " << ztrk <<
                                   " zcal0 = " << zcal0 << endl;
  }

  double dxdztrk = vpar[3];
  double dydztrk = vpar[4];

  double gate_x = Nsigma*excal;
  double gate_y = Nsigma*eycal;

  bool relaxed_gate = true;
  if( relaxed_gate  )
  {
    int icell = pe.GetMainCells()[0];
    double sx = GetCells()[icell].GetCellType().GetSizeX();
    double sy = GetCells()[icell].GetCellType().GetSizeY();
    gate_x = 2.8*sx/4.;
    gate_y = 2.8*sy/4.;
  }

  if(fabs(zcal0-zcalc) > 0.1 ) // OK precision first!
  {
//     Extrapolate(trackpar_at_destination,zcalc, trackpar_precise, false);
//     xtrk = trackpar_precise.GetX();
//     ytrk = trackpar_precise.GetY();

    xtrk = xtrk+dxdztrk*(zcalc-ztrk);
    ytrk = ytrk+dydztrk*(zcalc-ztrk);
  }
// It seems no need to search for nearest cluster to track. Try to associate now.
  bool associate = false;
  if( fabs( xcalc-xtrk) < gate_x && fabs( ycalc-ytrk) < gate_y ) // accepted
  {
    associate = true;
  }

  return associate;
}

////////////////////////////////////////////////////////////////////////////////

pair<bool,double> Calorimeter::InActiveAreaExternal(const std::vector<double>& vpar, int syst ) const
{
//   bool debug = true;
  assert(vpar.size() >= 5);
//  if( debug ) cout << " Kak by snova ne lopuxnut'sya !!!! " << vpar.size() <<" Z = " << vpar[0] << endl;
//  return InActiveAreaExternal(vpar[0],vpar[1],vpar[2],vpar[3],vpar[4]); !!! lopuxnulsya-to kak !!!!
  return InActiveAreaExternal(vpar[1],vpar[2],vpar[0],vpar[3],vpar[4]);
}

////////////////////////////////////////////////////////////////////////////////

pair<bool,double> Calorimeter::InActiveAreaExternal(double x_mrs,double y_mrs,double z_mrs, double dxdz, double dydz ) const
{
  vector <double> x_in;
  x_in.push_back(x_mrs);
  x_in.push_back(y_mrs);
  x_in.push_back(z_mrs);
  return pair<bool,double>(InActiveArea(x_in), 0.);
}

////////////////////////////////////////////////////////////////////////////////

// bool Calorimeter::InActiveArea(double x_mrs,double y_mrs,double z_mrs) const
bool Calorimeter::InActiveArea(const vector <double> &x_in) const
{
//   bool debug = true;
//   if( debug ) cout << "  Calorimeter::InActiveArea MRS x = " << x_in[0] << " y="<< x_in[1] << " z="<< x_in[2] << endl;
  vector <double> x_out;
  ConvertMRS2DRS( x_in, x_out);
//  if( FindCell( x_mrs-GetPositionX(),y_mrs-GetPositionY(),z_mrs-GetPositionZ() ) >= 0 ) return true;
//   if( debug ) cout << "  Calorimeter::InActiveArea local x = " << x_out[0] << " y="<< x_out[1] << " z="<< x_out[2] << endl;
  if( FindCell( x_out ) >= 0 )
  {
//     if( debug ) cout << "  InActiveArea " << endl;
    return true;
  }
  else
  {
//     if( debug ) cout << " NOT InActiveArea " << endl;
    return false;
  }
}

////////////////////////////////////////////////////////////////////////////////

//bool Calorimeter::OnBoundary(double x_mrs,double y_mrs,double z_mrs) const
bool Calorimeter::OnBoundary(const vector<double> &x_mrs) const
{
  vector <double> x_out;
  ConvertMRS2DRS( x_mrs, x_out);
//   int m = FindCell( x_mrs-GetPositionX(),y_mrs-GetPositionY(),z_mrs-GetPositionZ() );
  int m = FindCell( x_out );
  if( m >= 0 )
  {
//    if(  cells[m].InBoundaryRegion(x_mrs-GetPositionX(),y_mrs-GetPositionY(),z_mrs-GetPositionZ()) ) return true;
    if(  cells[m].InBoundaryRegion( x_out ) ) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////

// bool Calorimeter::WellInActiveArea(double x_mrs,double y_mrs,double z_mrs) const
bool Calorimeter::WellInActiveArea(const vector<double> &x_mrs) const
{
  vector <double> x_out;
  ConvertMRS2DRS( x_mrs, x_out);
//   int m = FindCell( x_mrs-GetPositionX(),y_mrs-GetPositionY(),z_mrs-GetPositionZ() );
  int m = FindCell( x_out );
  if( m >= 0 )
  {
//     if(  !cells[m].InBoundaryRegion(x_mrs-GetPositionX(),y_mrs-GetPositionY(),z_mrs-GetPositionZ()) ) return true;
    if(  !cells[m].InBoundaryRegion( x_out ) ) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////

vector<size_t> Calorimeter::FindCells(double x,double y,double z,double r) const
{
  vector<size_t> c;
  for( size_t i=0; i<cells.size(); i++ )
  {
    const Cell &cell = cells[i];
    if( sqrt(sqr(cell.GetX()-x) + sqr(cell.GetY()-y) + sqr(cell.GetZ()-z)) <= r )
      c.push_back(i);
  }
  return c;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetEtotalCalib(void) const {
    double ETotal(0.);

    for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++)
        ETotal += it->GetEnergy();

    return ETotal;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetEtotalRaw(void) const {
    double ESumADC(0.);

    for (vector<CellDataRaw>::const_iterator it=signals_.begin(); it!=signals_.end(); it++)
        ESumADC += it->GetAmplitude();

    return ESumADC;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SetOptions ( const list<string>& reco_opt_list)
{
  bool debug = false;
  for( list<string>::const_iterator it=reco_opt_list.begin(); it!=reco_opt_list.end(); it++ )
  {
    if( debug ) cout << Reco::Calorimeter::GetName() <<" RecoOptions " << *it << endl;
    SetOptions(*it);
  }
}

////////////////////////////////////////////////////////////////////////////////

int   Calorimeter::FindCellByName         (const std::string &cell_name)
{
  // if this method is called for the first time
  // fill the map of cell names to id associations
  if (name_to_id_.empty())
    for( int i=0; i<(int)NCells(); i++ )
      name_to_id_[GetCellName(i)] = i;

  std::map<std::string, int>::const_iterator it;
  if ( (it=name_to_id_.find(cell_name))==name_to_id_.end() )
    return -1;

  return it->second;
}

////////////////////////////////////////////////////////////////////////////////

string   Calorimeter::GetCellName          (int icell) const
{
  char cell_name[132];
  if( XYRegularGrid() )
  {
    int x = GetColumnOfCell(icell);
    int y = GetRowOfCell(icell);
    int res = snprintf(cell_name, 132, "X_%d_Y_%d", x, y);
    assert ( 0 <= res && res < 132 );
  }
  else
  {
    int res = snprintf(cell_name, 132, "%s_Cell_%d", GetName().c_str(), icell);
    assert ( 0 <= res && res < 132 );
  }
  return cell_name;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::GetCellOfColumnRow(int x, int y) const
{
  assert( XYRegularGrid() );

  if ( 0 <= x && x < fNcols_ && 0 <= y && y < fNrows_ )
    return map_cell_xy_[x+y*fNcols_];

  return -1;
}

////////////////////////////////////////////////////////////////////////////////

int  Calorimeter::GetColumnOfCell(int icell) const
{
  // provide backwards compatibility for PHAST, should be replaced by
  // the assertion sooner or later
  //  assert( XYRegularGrid() );
  if ( !XYRegularGrid() )
    return -1;

  if ( icell >= 0 && icell < (int)NCells() )
    return map_xy_cell_[0][icell];

  return -1;
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::GetRowOfCell(int icell) const
{
  // provide backwards compatibility for PHAST, should be replaced by
  // the assertion sooner or later
  //  assert( XYRegularGrid() );
  if ( !XYRegularGrid() )
    return -1;

  if ( icell >= 0 && icell < (int)NCells() )
    return map_xy_cell_[1][icell];

  return -1;
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::ScaleCalibration( size_t when )
{
  if( options.scale_calibration == 1. ) return;
  cout << " WARNING ! Actual calibrations in " << GetName() <<
                " were altered by factor " << options.scale_calibration << endl;
  for( size_t i=0; i<NCells(); i++ )
  {
    double stat  = cells_info[CALIB][when][i].GetEntries();
    double value = cells_info[CALIB][when][i].GetMean()*options.scale_calibration;
    double sigma = cells_info[CALIB][when][i].GetSigma()*options.scale_calibration;
    cells_info[CALIB][when][i].Set(stat, value, sigma);
  }
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::RecalculateParticlesPosition( std::vector<CalorimeterParticle> &particles ) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " RecalculateParticlesPosition " << GetName() << " x= " << GetPositionCorrectionX() <<
                                           " y= " << GetPositionCorrectionY() <<" z= " << GetPositionCorrectionZ() << endl;
  if( particles.size() == 0 ) return;
  double correct_x = GetPositionCorrectionX();
  double correct_y = GetPositionCorrectionY();
  double correct_z = GetPositionCorrectionZ();
  if( correct_x == 0. &&  correct_y == 0. &&  correct_z == 0.) return;
  if( debug ) cout << " RecalculateParticlesPosition " << GetName() << " x= " << correct_x <<
                                           " y= " << correct_y <<" z= " << correct_z << endl;
  for( size_t it=0; it< particles.size(); it++)
  {
    CalorimeterParticle &pp = particles[it];
    double x = pp.GetX()+correct_x;
    double sx = pp.GetXerr();
    pp.SetX(x,sx);

    double y = pp.GetY()+correct_y;
    double sy = pp.GetYerr();
    pp.SetY(y,sy);

    double z = pp.GetZ()+correct_z;
    double sz = pp.GetZerr();
    pp.SetZ(z,sz);
  }
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::RecalculateActualToStartOfRunPosition( std::vector<CalorimeterParticle> &particles ) const
{
  bool debug = true;
  if( particles.size() == 0 ) return;
  double correct_x = position[0] - position_atstartofrun_[0];
  double correct_y = position[1] - position_atstartofrun_[1];
  double correct_z = position[2] - position_atstartofrun_[2];
  if( debug ) cout << " RecalculateActualToStartOfRunPosition for " << GetName() <<
                 " start_x " << position_atstartofrun_[0] <<
                   " start_y " << position_atstartofrun_[1] <<
                    " start_z " << position_atstartofrun_[2] <<  endl;
  if( correct_x == 0. &&  correct_y == 0. &&  correct_z == 0.)
  {
    if( debug ) cout << " No need to RecalculateActualToStartOfRunPosition for " << GetName() << endl;
    return;
  }
  if( debug ) cout << " Start RecalculateActualToStartOfRunPosition for " << GetName() <<
                 " correct_x " << correct_x <<  " correct_y " << correct_y << " correct_z " << correct_z <<  endl;
  for( size_t it=0; it< particles.size(); it++)
  {
    CalorimeterParticle &pp = particles[it];
    double x = pp.GetX()+correct_x;
    double sx = pp.GetXerr();
    pp.SetX(x,sx);

    double y = pp.GetY()+correct_y;
    double sy = pp.GetYerr();
    pp.SetY(y,sy);

    double z = pp.GetZ()+correct_z;
    double sz = pp.GetZerr();
    pp.SetZ(z,sz);
  }
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::RecalculateParticlesEnergy( std::vector<CalorimeterParticle> &particles ) const
{
  if( particles.size() == 0 ) return;
  double scale_e = GetScaleCalibration();
  if( scale_e == 1. && !options.recalc_from_primary_calib ) return;
  for( size_t it=0; it< particles.size(); it++)
  {
    CalorimeterParticle &pp = particles[it];
    const std::vector < std::pair< size_t,double> >  &cluster_data = pp.GetClusterData();
    std::vector < std::pair< size_t,double> >  new_cluster_data;
    double e_old = 0.;
    double e_new = 0.;
    for( int it_data = 0;  it_data < (int)cluster_data.size();  it_data++ )
    {
      size_t icell = cluster_data[it_data].first;
      double ecell = cluster_data[it_data].second;
      double corr = 1.;
      if( options.recalc_from_primary_calib )
      {
        double cf_old = cells_info[CALIB][PRIM][icell].GetMean();
        double cf_new = cells_info[CALIB][OLD][icell].GetMean();
        corr = cf_new/cf_old;
//        if( debug ) cout << " For cell " << icell << " cf_new = " << cf_new << " cf_old = " << cf_old <<" corr = " << cf_new/cf_old << endl;
        e_old +=  ecell;
        e_new +=  ecell*corr;
      }
      new_cluster_data.push_back( pair< size_t,double>( icell, ecell*corr*scale_e) );
    }
    pp.ReplaceClusterData(new_cluster_data);
    double factor = scale_e;
    if( options.recalc_from_primary_calib )
    {
      if( e_old > 0. ) factor = factor*e_new/e_old;
    }
    double e = pp.GetE()*factor;
    double se = pp.GetEerr()*factor;
    pp.SetE(e,se);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SetPosition( double x, double y, double z )
{
  position[0]=x;
  position[1]=y;
  position[2]=z;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::SetVertexPosition( double x, double y, double z )
{
  vertex_position[0]=x;
  vertex_position[1]=y;
  vertex_position[2]=z;
}

////////////////////////////////////////////////////////////////////////////////

void   Calorimeter::ScaleCalibrationByLED( void )
{
  bool debug = true;
//   bool debug = false;

  double min_led_we_trust = 20.; // could be an option

  if( debug ) cout << " Calorimeter::ScaleCalibrationByLED options.scale_calibration_by_led = " <<
                                                            options.scale_calibration_by_led  << endl;
  if( !options.scale_calibration_by_led ) return;
  for( size_t i=0; i<NCells(); i++ )
  {
//    double ref_led_stat  = cells_info[LED][OLD][i].GetEntries();
    double ref_led_value = cells_info[LED][OLD][i].GetMean();
//    double ref_led_sigma = cells_info[LED][OLD][i].GetSigma();

//    double current_led_stat  = cells_info[LED][MONITOR][i].GetEntries();
    double current_led_value = cells_info[LED][MONITOR][i].GetMean();
//    double current_led_sigma = cells_info[LED][MONITOR][i].GetSigma();

    if( debug ) cout << " Cell " << i << " refLED= " << ref_led_value <<
                                   " CurrentLED " << current_led_value << endl;
    double scale = 1;
    if( current_led_value > min_led_we_trust && ref_led_value > min_led_we_trust ) scale = ref_led_value/current_led_value;

    double stat  = cells_info[CALIB][OLD][i].GetEntries();
    double value = cells_info[CALIB][OLD][i].GetMean()*scale;
    double sigma = cells_info[CALIB][OLD][i].GetSigma()*scale;
    if( debug && scale != 1 ) cout << " Cell " << i <<  " Old Cf=" << cells_info[CALIB][OLD][i].GetMean() <<
                        " New Cf=" <<  value <<" refLED= " << ref_led_value <<
                        " CurrentLED " << current_led_value << endl;
    if( current_led_value <= min_led_we_trust )
    {
        cerr << " Calorimeter::ScaleCalibrationByLED " << GetName() <<" Cell " << i <<" "<<GetCellName(i) <<
	                                                                  " Current led value is too low = " << current_led_value << endl;
    }
    if( ref_led_value <= min_led_we_trust )
    {
        cerr << " Calorimeter::ScaleCalibrationByLED " << GetName() <<" Cell " << i <<" "<<GetCellName(i) <<
	                                                                  " Reference led value is too low = " << ref_led_value << endl;
    }

    cells_info[CALIB][OLD][i].Set(stat, value, sigma);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::ProcessPedestalsAtEndOfJob( void )
{
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::FilterBadCellsOutOfStatisticLED( void )
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout <<" FilterBadCellsOutOfStatisticLED " << GetName() << endl;
  double led_min = 20.;
  double led_max = 1000.;
  if( GetName() == "EC02P1__" ) led_max = 4000.;
  double stat_led_max = 0.;
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    double stat_led =  cells_info[LED][NEW][icell].GetEntries();
    if( stat_led_max < stat_led ) stat_led_max = stat_led;
  }

  if( stat_led_max < 10. )
  {
    if( debug ) cout <<" No LED statistic FilterBadCells  canceled " << GetName() << endl;
    return;
  }
  double stat_led_max2tolerate = 0.8*stat_led_max;
  double stat_led_min2tolerate = 0.2*stat_led_max;

  int bad_cnt = 0;
  for( size_t icell=0; icell!=NCells(); icell++ )
  {
    int fitflag = fit_info[LED][icell].fit_ok;
    double led =  cells_info[LED][NEW][icell].GetMean();
    double sigma_led =  cells_info[LED][NEW][icell].GetMean();
    double stat_led =  cells_info[LED][NEW][icell].GetEntries();
    if( stat_led > stat_led_max2tolerate ) // statisic OK
    {
      if( led < led_min )
      {
        if( debug ) cout <<" Problems in cell " << icell <<" Stat OK, but low led  " << GetCellName(icell) <<" led " << led << endl;
        bad_cnt++;
        AddNewBadCell(icell, 1);
      }
    }
    else if( stat_led > stat_led_min2tolerate ) // suspiciously low statisic cell
    {
      if( led < led_min )
      {
        if( debug ) cout <<" Problems in cell " << icell <<" Low stat, no led  " << GetCellName(icell) <<" led " << led << endl;
        bad_cnt++;
        AddNewBadCell(icell, 1);
      }
    }
    else // Clear problem with LED
    {
      if( led < led_min )
      {
        if( debug ) cout <<" Problems in cell " << icell <<" No stat, no led " << GetCellName(icell) <<" led " << led << endl;
        bad_cnt++;
        AddNewBadCell(icell, 1);
      }
    }
  }
  if( debug ) cout <<" FilterBadCellsOutOfStatisticLED " << GetName() <<" BadCells = " << bad_cnt << " found " << endl;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::AddNewBadCell( int icell, unsigned status )
{
  if(icell >=0 && icell < (int)NCells() )
  {
    if( cell_is_bad_[NEW][icell] )
    {
      status_bad_cells_[NEW][icell] = status_bad_cells_[NEW][icell]|status;
    }
    else
    {
      bad_cells_[NEW].push_back(icell);
      status_bad_cells_[NEW][icell] = status_bad_cells_[NEW][icell]|status;
      cell_is_bad_[NEW][icell]=true;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::EndOfJob( void )
{
  bool debug = false;
//   bool debug = true;

  if( debug ) cout << " Reco::Calorimeter " << GetName() << " EndOfJob " << endl;

  ProcessPedestalsAtEndOfJob();
  FilterBadCellsOutOfStatisticLED();

  if( debug ) cout << " FitAllCellsHisto " << endl;
  if( options.fit_led_histo_at_end_of_job) FitAllCellsHisto(LED,"");

  if( debug ) cout << " SaveFit2CellInfo " << endl;
  SaveFit2CellInfo(-1); // Should be optional !
  FillFitInfoHist();
}

////////////////////////////////////////////////////////////////////////////////

int Calorimeter::GetSubSetIndex(const std::string &name) const
{
  for( int i=0; i < (int)calorimeter_sub_sets.size(); i++)
    if( name == calorimeter_sub_sets[i].GetName() ) return i;
  return -1;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::UpdateOnNewRun( int previous_run, int new_run )
{
  if( positionX_in_burst_new_.size() > 0 )
  {
    cerr << " !!!!!!!!! WARNING !!! For the moment we suppose to store spill by spill position only for one RUN !!! " << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
#include "Shower.h"

double Calorimeter::MakePreshowerCorrections( double eg, double xpart, double ypart, CalorimeterParticle::ParticleID pid ) const
{
  if( options.preshower_length > 0. )
  {
    size_t ind = (size_t)(fabs(eg)/0.1);
    if( preshower_corrections_.size() == 0 )
    {
      cerr << " Fatal you  must run InitPreshowerCorrections before using them  " << GetName() << "  !!! " << endl;
      exit(1);
    }
    if( preshower_corrections_.size() <= ind )
    {
      cerr << " Failed to MakePreshowerCorrections in " << GetName() << " for Egamma " << eg << " !!! " << endl;
      return eg;
    }
    double corr = preshower_corrections_[ind];
    return eg*(1.+corr);
  }
  else
  {
    return eg;
  }
}

////////////////////////////////////////////////////////////////////////////////

double  Calorimeter::GetInCurrentSpillMeanLED(int icell) const
{
//   bool debug = true;
  bool debug = false;
  const EventID &evid = GetEventID();
  int run    = evid.GetRunNumber();
  int ispill = evid.GetBurstNumber();
  double v   = GetInSpillMeanLED( icell, run, ispill );

  if ( debug )
    cout << " In Calorimeter " << GetName() << "  GetInCurrentSpillMeanLED run " << run
	 <<" spill " << ispill <<" cell " << icell << " led " << v  << endl;

  return v;
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::SetInSpillMeanLED(int icell, int run, int ispill, double v)
{
  if( v >= 0. )
  {
    double old_value = GetInSpillMeanLED(icell, run, ispill);
    if ( old_value != 0. )
      throw Exception("%s: SetInSpillMeanLED(cell %s, run %i, spill %i, value %f) "
		      "trying to overwrite existing value %f",
		      GetName().c_str(), GetCellName(icell).c_str(), run, ispill, v, old_value);
    if ( ispill >= 10000 )
      throw Exception("%s: SetInSpillMeanLED(cell %s, run %i, spill %i, value %f) "
		      "exceeding maximum spill number 9999",
		      GetName().c_str(), GetCellName(icell).c_str(), run, ispill, v);

    long unsigned int key = 10000ul * run + ispill;
    leds_in_spills_old_[icell].insert( pair<long unsigned int,double>(key, v) );
  }
}

////////////////////////////////////////////////////////////////////////////////

double  Calorimeter::GetInSpillMeanLED( int icell, int run, int ispill ) const
{
  long unsigned int key = 10000ul * run + ispill;
  std::map<long unsigned int,double>::const_iterator end = leds_in_spills_old_[icell].end();
  std::map<long unsigned int,double>::const_iterator e   = leds_in_spills_old_[icell].find(key);
  if ( e != end )
    return e->second;
  else
    return 0.;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::Reconstruction( void ) {
  if ( reconstruction_==NULL ) {
    if ( options.fortran_reconstruction ) {
      if ( reconstruction_==NULL ) {
        reconstruction_ = new ReconstructionLednev(this);
      } else {
        delete reconstruction_;
        reconstruction_ = NULL;

        throw Exception("Calorimeter::Reconstruction: Multiple reconstruction algorithms "
			"selected in options file for %s, it will not be reconstructed!", GetName().c_str());
      }
    }

    if ( options.kolosov_reconstruction ) {
      if ( reconstruction_==NULL ) {
        reconstruction_ = new ReconstructionKolosov(this);
      } else {
        delete reconstruction_;
        reconstruction_ = NULL;

        throw Exception("Calorimeter::Reconstruction: Multiple reconstruction algorithms selected in options file for %s, it will not be reconstructed!", GetName().c_str());
      }
    }

    if ( options.combined_reconstruction ) {
      if ( reconstruction_==NULL ) {
        reconstruction_ = new ReconstructionCombined(this);
      } else {
        delete reconstruction_;
        reconstruction_ = NULL;

        throw Exception("Calorimeter::Reconstruction: Multiple reconstruction algorithms selected in options file for %s, it will not be reconstructed!", GetName().c_str());
      }
    }

    if ( reconstruction_ == NULL )
      throw Exception("Calorimeter::Reconstruction: No reconstruction algorithm selected in options file for %s, it will not be reconstructed!", GetName().c_str());

    // create histograms
    reconstruction_->BookHistograms();

    // read calibrations
    reconstruction_->ReadCalibrations();
  }

  if ( options.calib_event_reconstruction || !GetEventID().IsLEDEvent())
    reco_particles = reconstruction_->DoReconstruction(signals_);
  else
    reco_particles.clear();

  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::RepeatReconstruction() {
  bool debug = false;
//   bool debug = true;
  if( debug ) cout <<" RepeatReconstruction debug  repeat_option = " << GetOptions().repeat_reconstruction << endl;
  if ( reconstruction_ == NULL ) {
      reconstruction_ = new ReconstructionKolosov(this);
      cout <<" Warning !!! For the moment only ReconstructionKolosov is available in RepeatReconstruction function " << endl;
  }
//       throw Exception("Calorimeter::RepeatReconstruction: No reconstruction object existing for %s, so no reconstruction was executed. Will not re-reconstruct!", GetName().c_str());

  if( GetOptions().repeat_reconstruction  < 0 ) return false;

  if ( options.calib_event_reconstruction || !GetEventID().IsLEDEvent())
    reco_particles = reconstruction_->DoRepeatReconstruction(signals_);
  else
    reco_particles.clear();

  return true;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetRandomFlat(void) const {
    const double r=drand48();

    return r;
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::GetRandomGaus(void) const {
    const double u1=2*M_PI*drand48();
    const double u2=sqrt(-2*log(drand48()));

    const double r1=sin(u1)*u2;
    const double r2=cos(u1)*u2;

    return r1;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::MonteConstruction(const vector<CalorimeterParticle>& particles, vector<CellDataRaw>& data) {
  if ( monteconstruction_==NULL ) {
    monteconstruction_ = new MCConstruction(this);

    if ( monteconstruction_ == NULL )
      throw Exception("Calorimeter::MonteConstruction: MCConstruction object could not be created for %s, MC will fail!", GetName().c_str());
  }

  // get ideal energy distribution of the particles in the calorimeter from a shower profile
  mc_input_ = monteconstruction_->DoDigitization(particles);

  // now simulate the detector response
  MakeMCDigitization();

  data = signals_;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::MakeMCDigitization() {
  if ( monteconstruction_==NULL ) {
    monteconstruction_ = new MCConstruction(this);

    if ( monteconstruction_ == NULL )
      throw Exception("Calorimeter::MakeMCDigitization: MCConstruction object could not be created for %s, MC will fail!", GetName().c_str());
  }

  signals_ = monteconstruction_->DoConstruction(mc_input_);
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::MakeMCHits(const CalorimeterParticle& particle_in, CalorimeterParticle& particle_out,
                             std::vector<CalorimeterParticle> &other_particles_out) {
    if ( monteconstruction_==NULL ) {
        monteconstruction_ = new MCConstruction(this);

        if ( monteconstruction_ == NULL )
            throw Exception("Calorimeter::MakeMCHits: MCConstruction object could not be created for %s, MC will fail!", GetName().c_str());
    }

    if ( !monteconstruction_cleared_ ) {
        monteconstruction_->ClearExpected();
        monteconstruction_cleared_ = true;
    }

    double depositedEnergy = monteconstruction_->AddParticle(particle_in);

    // The following is a quick and dirty patch
    particle_out=particle_in;

    double einit = particle_in.GetE();
    if ( einit-depositedEnergy < 0.001 ) // stop tracking, this is not our particle any more
        particle_out.SetE(0.001);
    else
        particle_out.SetE(einit-depositedEnergy);
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::MakeMCResponse() {
    MakeMCDigitization();
}

////////////////////////////////////////////////////////////////////////////////

/// Output function:  Called from CsEvent to collected reconstructed clusters.
vector<CalorimeterParticle>  Calorimeter::GetCalorimeterParticles(void) const
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " Calorimeter::GetCalorimeterParticles " << GetName() <<   endl;
  if( debug ) cout << "  reco_particles.size() = " <<  reco_particles.size() << endl;
  vector<CalorimeterParticle> calorimeter_particles;
  for( size_t ip = 0; ip < reco_particles.size(); ip++ )
  {
    CalorimeterParticle p = ConvertParticleDRS2MRS(reco_particles[ip]);
// #warning "Calorimeter::GetCalorimeterParticles temporary setting AngleX() and AngleY() to zero "
//     calorimeter_particles.push_back( CalorimeterParticle( (CalorimeterParticle::ParticleID)p.GetID(),p.GetProb(),p.GetE(),
//                                         p.GetX(),p.GetY(),p.GetZ(),
//                                         p.GetEerr(),p.GetXerr(),p.GetYerr(),p.GetZerr(),p.GetAngleX(),p.GetAngleY()) );
    calorimeter_particles.push_back( CalorimeterParticle( (CalorimeterParticle::ParticleID)p.GetID(),p.GetProb(),p.GetE(),
							  p.GetX(), p.GetY(), p.GetZ(),
                                                          this,
							  p.GetEerr(), p.GetXerr(), p.GetYerr(), p.GetZerr(),
							  p.GetAngleX(), p.GetAngleY() ) );
    if(p.HasTime())  calorimeter_particles.back().SetTime( p.GetTime(),p.GetTimeErr() );
    calorimeter_particles.back().SetMainCells( p.GetMainCells() );

    const vector< pair< size_t,double> >& data = p.GetClusterData();
    calorimeter_particles.back().SetClusterData(data);

    if( p.GetMainCells().size() <= 0 )
    {
      cerr << " ERROR! in Calorimeter::GetCalorimeterParticles: Main cells vector is empty " << endl;
      exit(1);
    }

    int icell = p.GetMainCells()[0];
    double gatex = GetCells()[icell].GetCellType().GetSizeX();
    double gatey = GetCells()[icell].GetCellType().GetSizeY();
    double gatez = GetCells()[icell].GetCellType().GetSizeZ();
    calorimeter_particles.back().SetGateX(gatex);
    calorimeter_particles.back().SetGateY(gatey);
    calorimeter_particles.back().SetGateZ(gatez);
    if( debug )
    {
      cout << " gatex " << gatex << " gatey " << gatey << " gatez " << gatez << endl;
      cout << " Check gatex " << calorimeter_particles.back().GetGateX() <<
              " gatey " << calorimeter_particles.back().GetGateY() <<
	      " gatez " << calorimeter_particles.back().GetGateZ() << endl;
    }

    std::map< CalorimeterParticle::MiscDataID, double > &mapin = p.GetMiscData();
    std::map< CalorimeterParticle::MiscDataID, double > &mapout = calorimeter_particles.back().GetMiscData();
//     if( debug )
//     {
//       cout << GetName() << " Calorimeter::GetCalorimeterParticles mapin.size " << mapin.size() <<
//                                                                " mapout.size " << mapout.size() << endl;
//     }
    for( std::map<CalorimeterParticle::MiscDataID, double>::iterator mit = mapin.begin(); mit != mapin.end(); mit++ )
    {
      mapout.insert( pair< CalorimeterParticle::MiscDataID, double > ( CalorimeterParticle::MiscDataID(mit->first), mit->second) );
    }
//     if( debug )
//     {
//       cout << GetName() << " Calorimeter::GetCalorimeterParticles mapout.size after insertion " <<  mapout.size() << endl;
//     }


//     cout << " Calorimeter " << GetName() << p.HasTime() << " Time " << p.GetTime() << endl;
//     cout << " CalorimeterParticle Has Time " << calorimeter_particles.back().HasTime() << " Time " <<
//                                                                     calorimeter_particles.back().GetTime() << endl;
  }
  return  calorimeter_particles;
}

////////////////////////////////////////////////////////////////////////////////

void  Calorimeter::SetCalorimeterParticles(const vector<CalorimeterParticle> &particles)
{
  reco_particles.clear();
  for( vector<CalorimeterParticle>::const_iterator p_in = particles.begin(); p_in != particles.end(); p_in++ )
  {
    CalorimeterParticle p = ConvertParticleMRS2DRS(CalorimeterParticle((CalorimeterParticle::ParticleID)p_in->GetID(),p_in->GetProb(),p_in->GetE(),
                                        p_in->GetX(),p_in->GetY(),p_in->GetZ(),
                                        this,
                                        p_in->GetEerr(),p_in->GetXerr(),p_in->GetYerr(),p_in->GetZerr(),0.,0.) );
    if(p_in->HasTime())  p.SetTime( p_in->GetTime(),p_in->GetTimeErr() );
    p.SetMainCells( p_in->GetMainCells() );
    reco_particles.push_back(p);
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Print(ostream &o, const string &prefix) const
{
  o << prefix << GetName() << "  Total matrixes amount is " << matrixes.size() << endl;
  for( list<CellsMatrix>::const_iterator it=matrixes.begin(); it!=matrixes.end(); it++ )
    it->Print(o,prefix+"matrix: ");
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::AddMCHit(int id_hit, const void *data)
{
  bool debug = false;
  if( debug ) cout << " Nunu poprobuem AddMCHit id_hit=" << id_hit << endl;
  bool result = false;

  if( id_hit!=4 )
    return false;                     // No, this is not a calorimter hit...

  struct CalorimeterMCData
  {
    unsigned int cell_id;
    double dE;
    double dT;
  };
  const CalorimeterMCData &d = *(CalorimeterMCData*)(data);

//cout << "cell_id = " << d.cell_id << "  E=" << d.dE << "\n";

  list<CellsMatrix>::iterator mtx = matrixes.begin();

  for( size_t i=0, cell_n=0; i<comgeant_first_cell_n_.size(); i++,mtx++ )
  {
    assert(mtx!=matrixes.end());
    size_t first_cell = comgeant_first_cell_n_[i];

//printf("~~~~  i=%d  first_cell=%d  size=%d\n",i,first_cell,mtx->Size());

    if( d.cell_id>first_cell && (d.cell_id-first_cell)<=mtx->Size() )
    {
      // OK, this is hit from THIS calorimeter and THIS matrix.

      cell_n += d.cell_id-first_cell-1;         // This is cell number in the range [0,NCells())
      if( cell_n>=NCells() )
        throw Exception("CsCalorimeter::AddMCHit(): d.cell_id=%d  first_cell=%d  cell_n=%d   NCells=%d",
			d.cell_id,first_cell,cell_n,NCells());

      if( debug )  cout << "HIT: Calorimeter " << GetName() << " matrix " << i << " cell_n=" <<
                       cell_n << " dE=" << d.dE << endl;

      mc_input_.push_back( CellDataRaw(cell_n, d.dE, d.dT) );
//       cout << " Hit Energy " << d.dE << " and Time " << d.dT << endl;

      result = true;
      break;
    }
    else
      cell_n += mtx->Size();
  }

//   // Fill histogram.
//   if( hist_MC_hitE!=NULL )
//     hist_MC_hitE->Fill(GetAmplitudesSum());

  return result;
}

////////////////////////////////////////////////////////////////////////////////

} /// using namespace Reco
