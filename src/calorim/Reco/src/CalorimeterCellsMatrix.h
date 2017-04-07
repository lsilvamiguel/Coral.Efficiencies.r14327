#ifndef inside_Calorimeter_h
#  error Only for inclusion from inside Calorimeter.h!
#endif


/*! \brief CellsMatrix is a rectangular matrix of identical cells, may be
  with holes - such structure is quite usual.
*/
class CellsMatrix
{
 public:

  /*! \brief Matrix hole with size [Xmin,Xmax]x[Ymin,Ymax]
   */
  class Hole
    {
    public:

      /// Destructor
      ~Hole                    (void) {}

      /// Constructor
      Hole                    (void) : x_min(1), x_max(0), y_min(1), y_max(0) {}
      Hole                    (size_t x_min_, size_t x_max_, size_t y_min_, size_t y_max_)
	: x_min(x_min_), x_max(x_max_), y_min(y_min_), y_max(y_max_) {
	if( x_min>x_max || y_min>y_max )
	  throw "Hole::Hole()  bad hole size (x_min>x_max or y_min>y_max)";
      }

      size_t    Xmin                    (void) const {return x_min;}
      size_t    Xmax                    (void) const {return x_max;}
      size_t    Ymin                    (void) const {return y_min;}
      size_t    Ymax                    (void) const {return y_max;}
      size_t    SizeX                   (void) const {return (x_max-x_min+1);}
      size_t    SizeY                   (void) const {return (y_max-y_min+1);}
      size_t    Size                    (void) const {return SizeX()*SizeY();}
      bool      IsInside                (size_t x, size_t y) const {
	return x>=x_min && x<=x_max && y>=y_min && y<=y_max; }
      void      Print                   (std::ostream &o=std::cout,
					 const std::string &prefix="") const;

    private:
      size_t    x_min, x_max, y_min, y_max;

    };  // end class Hole



  /// Constructor in case no holes
  CellsMatrix (const std::string &name,
	       const CellType &ct,
	       int Size_x, int Size_y,
	       double Pos_x, double Pos_y, double Pos_z,
	       const std::vector<Hole> &h = std::vector<Hole>() )
    : matrix_name  (name),
    cell_type    (ct),
    holes        (h),
    first_cell_(0),
    last_cell_(0)
    { pos[0]=Pos_x; pos[1]=Pos_y; pos[2]=Pos_z; size[0]=Size_x; size[1]=Size_y; }

  /// Assignment operator
  CellsMatrix      &operator =         (const CellsMatrix &m);

  /// Copy constructor
  CellsMatrix (const CellsMatrix &m) { *this = m; }

  const std::string &GetName( void ) const { return matrix_name; }

  void AddHole( const Hole &h ) {holes.push_back(h);}

  /// \return Total number of cells in the Matrix.
  size_t Size( void ) const { return size[0]*size[1]; }

  const std::vector<Hole> &GetHoles (void) const { return holes; }

  /// \return true if cell (x,y) is in Hole
  bool IsInHole( size_t x,size_t y ) const;

  const CellType &GetCellType( void ) const { return cell_type; }

  /// Print properties
  void            Print(std::ostream &o=std::cout,
			const std::string &prefix="") const;

  size_t           GetFirstCellIndex ( void ) const {return first_cell_;}
  size_t           GetLastCellIndex  ( void ) const {return last_cell_;}
  void             SetFirstCellIndex ( size_t first_cell ) {first_cell_ = first_cell;}
  void             SetLastCellIndex  ( size_t last_cell ) {last_cell_ = last_cell;}

 private:

  /// The matrix's name.
  std::string          matrix_name;

  /// Matrix Cells step.
  CellType        cell_type;

  /// Matrix size.
  int             size[2];

  /// Matrix position in mother Calorimeter reference system.
  double          pos[3];

  /// Vector of matrix's holes.
  std::vector<Hole>                        holes;

  size_t           first_cell_;
  size_t           last_cell_;

};



