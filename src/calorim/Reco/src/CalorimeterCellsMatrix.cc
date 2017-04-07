#include "Calorimeter.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

Calorimeter::CellsMatrix & Calorimeter::CellsMatrix::operator = (const Calorimeter::CellsMatrix &m)
{
  if( &m!=this )
    {
      matrix_name = m.matrix_name;
      cell_type   = m.cell_type;
      size[0]     = m.size[0];
      size[1]     = m.size[1];
      pos[0]      = m.pos[0];
      pos[1]      = m.pos[1];
      pos[2]      = m.pos[2];
      holes       = m.holes;
      first_cell_ = m.first_cell_;
      last_cell_  = m.last_cell_;
    }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::CellsMatrix::IsInHole(size_t x, size_t y) const
{
  for( std::vector<Hole>::const_iterator it=holes.begin(); it!=holes.end(); it++ )
    if( it->IsInside(x,y) )
      return true;
  return false;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CellsMatrix::Print(std::ostream &o,const std::string &prefix) const
{
  o << prefix << matrix_name << " size: " << size[0] << " x " << size[1] <<
          " pos=( " << pos[0] << "," << pos[1] << "," << pos[2] << ")   CellType=";
  cell_type.Print(o);

  for( std::vector<Hole>::const_iterator it=holes.begin(); it!=holes.end(); it++ )
  {
    it->Print(o,prefix+"hole: ");
    o<<"\n";
  }
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::CellsMatrix::Hole::Print(std::ostream &o,const std::string &prefix) const
{
  o << prefix << "[Xmin,Xmax]x[Ymin,Ymax]=[" <<
       Xmin() << "," << Xmax() << "]x["<< Ymin() << "," << Ymax() <<"]";
}

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco
