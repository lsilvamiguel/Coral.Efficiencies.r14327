#ifndef inside_Calorimeter_h
#  error Only for inclusion from inside Calorimeter.h!
#endif


/*! \brief List of calorimeter cells with some special properties or functions.
  In some cases may be considered as sub-detector description.
*/
class SubSet
{
  // =========================================================================
  // Constructors and destructor
  // =========================================================================
 public:

  /// Destructor
  virtual            ~SubSet                (void) {}

  SubSet                ( Calorimeter *c, const std::string &the_name,
			  const std::vector <size_t> &sub_set) :
    calorimeter(c),
    name(the_name),sub_set_cells(sub_set),
    xy_structure_init_(false) {
  }

  // =========================================================================
  // Operators
  // =========================================================================

 public:

  /// Test for equality
  bool                operator ==             (SubSet &s);

  /// Print SubSet info to output stream
  friend std::ostream     &operator <<             (std::ostream &o,SubSet &s);

  // =========================================================================
  // Methods
  // =========================================================================

 public:

  /// \return SubSet name.
  const std::string         &GetName           (void) const {return name;}

  /// Clear
  void          Clear           (void) {info_.clear(); finfo_.clear();}

  /// \return vector of cells numbers
  const std::vector<size_t> &GetCells          (void) const {return sub_set_cells;}

  /// \return  amount of cells in sub_set
  size_t                Size                    (void) const {return sub_set_cells.size();}

  /// \return  Vector of associated integer info
  std::vector < int >   &GetIntData             (void) {return info_;}

  /// \return  Vector of associated double info
  std::vector < double >   &GetData             (void) {return finfo_;}

  /// \return calibrated  energy as default_info
  double                    GetEnergy            (void) {if(finfo_.size() >0 )
    return finfo_[0];
  else
    return 0.;
  }
  /// \return raw  energy as default_info
  double                    GetAmp               (void) {if(info_.size() >0 )
    return info_[0];
  else
    return 0;
  }

  /// Write sub_set configuration to the file
  bool                  WriteConfig             (const std::string &file_name, Calorimeter &c) const;

  /// Read sub_set configuration from the file
  bool                  ReadConfig              (const std::string &file_name, Calorimeter &c);

  /// Init XY structure
  bool                InitXY             ( double tolerance );
  int 	              GetCellOfColumnRow (int x, int y) const;
  int 	              GetColumnOfCell    (int icell) const;
  int 	              GetRowOfCell       (int icell) const;

  int                 GetNColumns        ( void ) const { return fNcols_; }
  int                 GetNRows           ( void ) const { return fNrows_; }

  double              GetXStep           ( void ) const { return fXStep_; }
  double              GetYStep           ( void ) const { return fYStep_; }

  /// Left-most position occupied by a block (not including padding)
  double              GetXmin  ( void ) const;
  /// Right-most position occupied by a block (not including padding)
  double              GetXmax  ( void ) const;
  /// Lower-most position occupied by a block (not including padding)
  double              GetYmin  ( void ) const;
  /// Top-most position occupied by a block (not including padding)
  double              GetYmax  ( void ) const;
  bool                XYInitialized ( void ) const { return xy_structure_init_;}
  bool	              IsMyCell ( int icell) const;

 protected:
  int 	            GetCellIndex ( int icell) const;


  // =========================================================================
  // Data Members
  // =========================================================================

 protected:

  /// This is our Calorimeter
  Calorimeter                       *calorimeter;

  /// The sub_set's name.
  std::string                       name;

  /// Vector of cells numbers
  std::vector < size_t >            sub_set_cells;

  /// Vector associated integer info
  std::vector < int >               info_;

  /// Vector associated double info
  std::vector < double >            finfo_;

  /// SubSet XY structure
  bool                                xy_structure_init_;
  std::vector<int>                    map_cell_xy_;
  std::vector<int>                    map_xy_cell_[2];
  int                                 fNcols_;
  int                                 fNrows_;
  double                              fXStep_;
  double                              fYStep_;
};

