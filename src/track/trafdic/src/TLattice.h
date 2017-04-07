// $Id: TLattice.h 13163 2012-01-09 03:41:19Z ybedfer $

/*!
  \file    Tlattice.h
  \brief   Lattice == Traffic class for fitting TTrack's using QNewton fit
  \author  Yann.Bedfer@cern.ch
  \version $Revision: 13163 $
  \date    $Date: 2012-01-09 04:41:19 +0100 (Mon, 09 Jan 2012) $
  */

#ifndef TLattice_h
#define TLattice_h

#include <cstdlib>   // for exit()
#include "THit.h"

/*! \class TLattice
  \brief Class devoted to the providing a look-up mechanism for
  determining hit positions given track parameters and then fitting.
  A priori using a file, viz. Dico, via a table. But could be using a
  track propagator.
  */

// ***** STRUCTURE FOR GRID DEFINITION
typedef struct{ float min; float max; int num; } parbound;

// ***** STRUCTURE FOR TRACK == PARAMETERS + HIT LIST
#define NROWS 142

#define NRW32 NROWS/32
typedef float Hits[NROWS];
typedef struct {
  float h[NROWS];
  int   n;
  unsigned int hp[NRW32+1];  /* hitpattern */
} Track_Hits;
typedef struct {
  float x, y, pxzi, tx, ty;
} Parameters;
typedef struct {
  Track_Hits hits;
  Parameters par;
  unsigned int par_pattern;  /* Which track-parameters should be fitted */
  float      z;     /* z-vertex position to use in tracking */
  float      chi2;  /* chi-square for the given set of parameters */
} Track;

class TLattice {

  unsigned int sieve[6];    //<! Bit pattern sieve (Dimension=6 is supposed to cover all detectors entitled to enter some dico)

  Hits       *grid;
  int        global_grid_point;

  parbound   bound[5];
  int        jump[5], multi_jump[32];

  Parameters il_diff;

  float al_v[5],ar_v[5];
  
  float max_par[5];   //! Maximum parameter values to avoid floating overflows

public:

  float *sigmas;

  static int IplLast;

  float dZStart;      //! Distance of dico ref. plane w.r.t. its origin. That distance is crossed by straight line extrapolation, cf. "TLattice::GenInit".  

  /*! \fn void TLattice(int);
    \brief Default Constructor. Reading is done here.
    \param read True for reading an existing dico. False for when creating a new one.
    */
  TLattice(int read);

  // copy constructor
  TLattice(const TLattice& t){
    std::cerr<<"Copy constructor is not available!"<<std::endl;
    exit(1);
  }

  /*! \fn TLattice Ptr();
    \brief Accessor: Returns pointer to TLattice class
    */
  static TLattice* Ptr();

  // assignment operator
  TLattice& operator = (const TLattice& t){
    std::cerr<<"= operator is not available\n"
	<<"Are you really needing 2 lattices at once!?"<<std::endl;
    exit(1);
  } 

  std::vector<int> dico_idx;

  // ********** PUBLIC METHODS **********

  // ***** LOOKUP FUNCTIONS
  void get_coords(Track *);
  int get_point(Parameters *);
  int get_point_par(Parameters *,Parameters *);
  void get_coords_auto(Track *);
  bool OutDico(Parameters *);
  bool FarOutDico(Parameters *);

  // ***** INTERPOLATION FUNCTIONS
  void interpolate_lin0(Hits,Hits*);
  void interpolate_lin1(Hits,Hits*);
  void interpolate_lin2(Hits,Hits*);
  void interpolate_lin3(Hits,Hits*);
  void interpolate_lin4(Hits,Hits*);

  // ***** FIT FUNCTION
  int Chi2Fit(Track *);
  int Chi2Fit_(Track *);
  int FastChi2Fit(Track *,int);
  int FastFit(Track *,float **,int,float *,float *, int,int);
  void derive_lattice(float *,float *,int *,int,Parameters *);

  // ***** GENERATION FUNCTIONS

  /*! \fn void GenInit();
    \brief Init the generation of the Dico.
    */
  void GenInit();
  /*! \fn void GenServer(float*,float*,int*);
    \brief Track server for the generation of the Dico.
    */
  void GenServer(float*,float*,int*);
  /*! \fn void GenOutput(float *);
    \brief Output the Dico to file.
    */
  void GenOutput(float *);
  /*! \fn void GenClose();
    \brief Close the generation of the Dico.
    */
  void GenClose();
  /*! \fn unsigned int getSize();
    \brief Return the size of the lattice grid.
    */
  unsigned int getSize();

private:
  static TLattice* address;  //!< pointer to itself  
};

inline bool TLattice::OutDico(Parameters *p) {
  float *pi; parbound *bi;
  for (pi = &p->x, bi = &bound[0]; pi<&p->x+5; pi++, bi++) {
    if (*pi<bi->min-1.5*(bi->max-bi->min)/bi->num ||
        *pi>bi->max+1.5*(bi->max-bi->min)/bi->num ) return true;
  }
  return false;
}
inline bool TLattice::FarOutDico(Parameters *p) {
  float *pi; parbound *bi;
  for (pi = &p->x, bi = &bound[0]; pi<&p->x+2; pi++, bi++) {
    if (*pi<bi->min-3.0*(bi->max-bi->min)/bi->num ||
        *pi>bi->max+3.0*(bi->max-bi->min)/bi->num ) return true;
  }
  for (                          ; pi<&p->x+5; pi++, bi++) {
    if (*pi<bi->min-1.5*(bi->max-bi->min)/bi->num ||
        *pi>bi->max+1.5*(bi->max-bi->min)/bi->num ) return true;
  }
  return false;
}
#endif
