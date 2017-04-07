#include <iostream>
#include "Exception.h"
#include "CalorimeterGAMS.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

CalorimeterGAMS::CalorimeterGAMS(const string &the_name,int type) :
  Calorimeter(the_name)
{
  cells_type.push_back(Reco::CellType("lead glass"   ,3.8,3.8,45,2,3,0.09*0.09));
  cells_type.push_back(Reco::CellType("PbWO4 crystal",1.9,1.9,20,2,3,0.06*0.06)); // ?? 

  vector<size_t> holesXmin, holesXmax, holesYmin, holesYmax;

  switch( type )
  {
    case 0:

      GAMS_Nx = GAMS_Ny = 48;
      GAMS_hole = Hole(GAMS_Nx/2-1,GAMS_Nx/2,GAMS_Ny/2-1,GAMS_Ny/2);

      SAD_Nx = SAD_Ny = 0;
      SAD_hole = Hole(0,0,0,0);
      
      break;

    case 2000:

      GAMS_Nx = GAMS_Ny = 48;
      GAMS_hole = Hole(GAMS_Nx/2-2,GAMS_Nx/2+1,GAMS_Ny/2-2,GAMS_Ny/2+1);

      SAD_Nx = SAD_Ny = 8;
      SAD_hole = Hole(3,4,3,4);
      
      break;

    default:
      throw Exception("CalorimeterGAMS::CalorimeterGAMS()   Unknown GAMS calorimeter type %d",type);
  }

  for( size_t y=0; y<GAMS_Ny; y++ )
    for( size_t x=0; x<GAMS_Nx; x++ )
      if( x>=GAMS_hole.Xmin && x<=GAMS_hole.Xmax &&
          y>=GAMS_hole.Ymin && y<=GAMS_hole.Ymax  )
        continue;   // hole
      else
      {
        double
          pos_x = (x-(GAMS_Nx-1)/2.)*cells_type[0].GetSizeX(),
          pos_y = (y-(GAMS_Ny-1)/2.)*cells_type[0].GetSizeY(),
          pos_z = 0;//cells_type[0].GetSizeZ()/2;
        cells.push_back( Cell(cells_type[0],true,pos_x,pos_y,pos_z) );
      }

  SAD_cells = cells.size();

  for( size_t y=0; y<SAD_Ny; y++ )
    for( size_t x=0; x<SAD_Nx; x++ )
      if( x>=SAD_hole.Xmin && x<=SAD_hole.Xmax &&
          y>=SAD_hole.Ymin && y<=SAD_hole.Ymax  )
        continue;   // hole
      else
      {
        double
          pos_x = (x-(SAD_Nx-1)/2.)*cells_type[1].GetSizeX(),
          pos_y = (y-(SAD_Ny-1)/2.)*cells_type[1].GetSizeY(),
          pos_z = 10+cells_type[1].GetSizeZ()/2;
        cells.push_back( Cell(cells_type[1],true,pos_x,pos_y,pos_z) );
      }

  Init();
}

////////////////////////////////////////////////////////////////////////////////

int CalorimeterGAMS::FindCell(double x,double y,double z) const
{
//   double
//     cell_size = cells_type[0].GetSizeX(),
//     x0        = GetCells()[0].GetY()-cell_size/2,
//     y0        = GetCells()[0].GetX()-cell_size/2;
//   // (x0,y0) - left down corner.
// 
//   size_t
//     i = (size_t) ((x-x0)/cell_size),
//     j = (size_t) ((y-y0)/cell_size);
// 
//   // i,j always >=0
//   if( i>=Nx || j>=Ny )
//     return -1;  // Outside the calorimeter.
// 
//   if( i>=hole.Xmin && i<=hole.Xmax && j>=hole.Ymin && j<=hole.Ymax )
//     return -1;  // Inside the hole.
// 
//   if( j<hole.Ymin )
//     return i+Ny*j;
// 
//   if( j>hole.Ymax )
//     return i+Ny*j - hole.size();

  // We are around the hole. Call general method.
  return Calorimeter::FindCell(x,y,z);
}

////////////////////////////////////////////////////////////////////////////////

int CalorimeterGAMS::GetGAMSCell(size_t x,size_t y) const
{
  if( x>=GAMS_Nx || y>=GAMS_Ny )
    //throw Exception("CalorimterGAMS::GetGAMSCell()  outside the calorimeter x=%d y=%d",x,y);
    return -1;

  if( x>=GAMS_hole.Xmin && x<=GAMS_hole.Xmax &&
      y>=GAMS_hole.Ymin && y<=GAMS_hole.Ymax )
    //throw Exception("CalorimterGAMS::GetGAMSCell()  inside the hole:  x=%d y=%d",x,y);
    return -1;

  if( y<GAMS_hole.Ymin )
    return x+GAMS_Nx*y;

  if( y>GAMS_hole.Ymax )
    return x+GAMS_Nx*y - GAMS_hole.Size();
  
  return GAMS_Nx*GAMS_hole.Ymin +
         (y-GAMS_hole.Ymin)*(GAMS_Nx-GAMS_hole.SizeX()) +
          x - ( x<GAMS_hole.Xmin ? 0 : GAMS_hole.SizeX() );
}

////////////////////////////////////////////////////////////////////////////////

int CalorimeterGAMS::GetSADCell(size_t x,size_t y) const
{
  if( x>=SAD_Nx || y>=SAD_Ny )
    //throw Exception("CalorimterSAD::GetSADCell()  outside the calorimeter x=%d y=%d",x,y);
    return -1;

  if( x>=SAD_hole.Xmin && x<=SAD_hole.Xmax &&
      y>=SAD_hole.Ymin && y<=SAD_hole.Ymax )
    //throw Exception("CalorimterSAD::GetSADCell()  inside the hole:  x=%d y=%d",x,y);
    return -1;

  if( y<SAD_hole.Ymin )
    return x+SAD_Nx*y;

  if( y>SAD_hole.Ymax )
    return x+SAD_Nx*y - SAD_hole.Size();
  
  return SAD_Nx*SAD_hole.Ymin +
         (y-SAD_hole.Ymin)*(SAD_Nx-SAD_hole.SizeX()) +
          x - ( x<SAD_hole.Xmin ? 0 : SAD_hole.SizeX() );
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
