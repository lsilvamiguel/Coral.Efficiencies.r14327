#ifndef MYPULSE_H_
#define MYPULSE_H_
#include <cassert>
#include <math.h>
#include <map>
/*!
  \brief Pulse form functional description
*/

#include "APV_timing.cc"

class MyPulse {

 private:

  double *fDots; /*! Array of points */
  unsigned int fNDots; /*! size of the array */

  double fXmin; /*! xaxis describtion */
  double fXmax; /*! xaxis describtion */
  double fMax; /*! Maximum/scaling factor */
  unsigned int fMaxBin; /*! Bin number of the maximum */

  double GetRPos(const double &x) {

    // Check if x in range
    if ( x < fXmin || x > fXmax )
      return fNDots;

    // find position
    return (fNDots - 1) * ( x - fXmin ) / ( fXmax - fXmin ) ;

  };

  unsigned int GetBin(const double &x) {

    return lround( GetRPos(x) );
  };

 public:


  double Evaluate(double *x, double *par) {

    if ( fXmin == fXmax && fNDots == 8 ) {

      double t = *x + fXmax - par[0];

      return par[1]*APV_SI_ampl_f(&t, fDots)+par[2];


    } else {

      // Get relative position
      double pos = GetRPos(*x - par[0]);

      unsigned int bin = (unsigned int) pos;

      // only the position between bins is now intresting
      pos -= bin;

      if ( bin < fNDots ) {
        double correction(0.);
        // Linear interpolation between bins
        if ( bin + 1 < fNDots )
          correction = pos * ( fDots[bin+1] - fDots[bin] );

        return par[1] * ( fDots[bin] + correction )/ fMax + par[2];
      } else
        // Elswhere 0;
        return par[2];
    }

  };

  double EvaluatePileUp(double *x, double *par){

    double* par1 = par, *par2 = par + 3;

    return Evaluate(x, par1) + Evaluate(x, par2);

  };


  double operator() (double *x, double *par) {

    return Evaluate(x, par);

  };

  MyPulse(double *_dots, unsigned int _ndots,
          const double &xmin, const double &xmax) :
    fDots(_dots), fNDots(_ndots), fXmin(xmin), fXmax(xmax) {

    // Search max
    fMax = fDots[0];
    fMaxBin = 0;
    for ( unsigned int i = 1; i < fNDots; i++)
      if ( fDots[i] > fMax ) {
        fMax = fDots[i];
        fMaxBin = i;
      }

  };

  MyPulse(double *_dots, unsigned int _ndots,
          const double &range) {
     Init(_dots, _ndots, range);
  }

  void Init(double *_dots, unsigned int _ndots,
       const double &range) {
    fDots = _dots;
    fNDots = _ndots;

    // APV  parametristion
    if( _ndots == 8 && range == 0. ) {

      double t = 0.;

      fMax = APV_SI_ampl_f(&t, _dots);

      fMaxBin = 0;

      fXmin = t;
      fXmax = t;

    } else {

      fMax = fDots[0];
      fMaxBin = 0;
      // Search max
      for ( unsigned int i = 1; i < fNDots; i++) {
        if ( fDots[i] > fMax ) {
          fMax = fDots[i];
          fMaxBin = i;
        }
      }

      // Set boarders
      fXmin = - range * fMaxBin  / fNDots;
      fXmax = range + fXmin;
    }

  };

  MyPulse(const char * filename) {

    // Open file
    FILE *f = fopen(filename, "r");

    // Not oppened file
    if ( !f ) {
      fprintf(stderr, "MyPulse::MyPulse(char*): Cannot open file '%s'!\n" , filename );
      throw 1;

    }


    // Get file size
    if (fseek(f, 0, SEEK_END))
      throw 2;
    long size = ftell(f);
    if (size == -1)
      throw 3;
    // Back to start of file
    if ( fseek(f, 0, SEEK_SET) )
      throw 2;

    // Get number of doubles in array
    size /= sizeof(double);

    double range(0.);
    int status = 0;

    if ( size != 8 ) {

      // First argument is total length
      size -= 1;

      // Get range
      status = fread( &range, sizeof(double), 1, f );
      if ( status != 1 )
      throw 4;
    }

    // Get Array
    double *a = new double[size];
    status = fread( a, sizeof(double), size, f );
    if ( status != size ) {
      if ( ferror(f) )
        throw 6;
      throw 5;
    }

    // Use constructor
    Init(a, size, range);

    fclose(f);

  };

  ~MyPulse() {

    // Delete array
    for ( unsigned int i = 0; i < fNDots; i++) {
      double *ptr = fDots;
      fDots++;
      //free(ptr);
    }

  };
};


class PulseManager {

 private:

  std::map< unsigned int, MyPulse* >  fMapPulse;
  std::map< unsigned int, void* > fMapVoid;

  static PulseManager *instance;

  PulseManager() {

    assert( !instance );

  };

 public:

  static PulseManager *Instance() {

    if( !instance )
      instance = new PulseManager();

    assert( instance );

    return instance;

  };

  MyPulse *GetPulse(unsigned int CalID, unsigned int Plane) {

    unsigned int Id = ( CalID << 16 ) | (Plane & 0xffff);

    std::map< unsigned int, MyPulse* >::iterator it = fMapPulse.find(Id);

    if( it == fMapPulse.end() )
      return NULL;

    return it->second;

  };

  void *GetVoid(unsigned int CalID, unsigned int Plane) {

    unsigned int Id = ( CalID << 16 ) | (Plane & 0xffff);

    std::map< unsigned int, void* >::iterator it = fMapVoid.find(Id);

    if( it == fMapVoid.end() )
      return NULL;

    return it->second;

  };

  int Insert(const unsigned int CalID, const unsigned int Plane,
             MyPulse* PULSE, void* VOID= NULL) {

    unsigned int Id = ( CalID << 16 ) | (Plane & 0xffff);

    if( PULSE ) {
      std::map< unsigned int, MyPulse* >::iterator it = fMapPulse.find(Id);

      if( it != fMapPulse.end() )
        return 1;
    }

    if( VOID ){
      std::map< unsigned int, void* >::iterator it = fMapVoid.find(Id);

      if( it != fMapVoid.end() )
        return 2;
    }

    if( PULSE )
      fMapPulse.insert(std::make_pair(Id,PULSE));

    if( VOID )
      fMapVoid.insert(std::make_pair(Id,VOID));

    return 0;

  };

};

PulseManager *PulseManager::instance = NULL;

#endif
