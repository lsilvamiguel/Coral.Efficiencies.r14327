// $Id: CsRTRelation.h,v 1.6 2006/06/01 11:57:33 corallib Exp $

/*!
   \file    CsRTRelation.h 
   \brief   RT relation for drift-like detectors
   \author  Hugo Pereira
   \version $Revision: 1.6 $
   \date    $Date: 2006/06/01 11:57:33 $
*/

#ifndef CsRTRelation_h
#define CsRTRelation_h


/*! \class CsSpacePoint 
    \brief Compass SpacePoint Class.
*/

#include "Coral.h"
#include "CsSTD.h"
#include "CsDetector.h"
#include "CsRTGridPoint.h"


class CsRTRelation {
 public:
  CsRTRelation( CsDetector &det );                     //!< default constructor 
  CsRTRelation( const CsRTRelation& rt ) {*this = rt;} //!< copy constructor
  CsRTRelation& operator = ( const CsRTRelation& );    //!< assignment operator  
  virtual ~CsRTRelation() {}                           //!< destructor.
   
  inline CsDetector* getDetector( void ) const { return d_; }
  inline bool hasGrid( void ) const { return hasGrid_; }
  inline bool hasRegularGrid( void ) const { return hasRegularGrid_; }
  inline bool hasRTParameters( void ) const { return hasRTPars_; }
  
  inline unsigned int gridSize( void ) const { return ( hasRegularGrid_ ) ? rtGReg_.size() : rtGrid_.size(); }
  inline unsigned int parSize( void ) const { return rtPars_.size(); }

  bool readFromOpt( void );
  inline void clear( void ){ rtGrid_.clear(); }
  double getTfromR( const double r, bool &error );  
  double getRfromT( const double t, bool &error );  
  double getRfromT( const double t, const int wire, bool &error );  
  double getResAtT( const double t, bool &error );  
  double getDRDT  ( const double t, bool &error );  

  void dump( FILE* out );
  void dump( void ) { dump( stdout ); }
  
 private:
  void _makeRegular( void );                   //!< modifies grid to have regular points and speed up getRfromT (or ResAtT)
  void _sortGrid( std::vector< RTGridPoint > &rt ); //!< used to sort grid vector (tricky and slightly dirty)
  
  double dt_;                    //!< time interval for regular grid
  CsDetector* d_;                //!< associate detector
  bool hasGrid_;                 //!< rtGrid_ make sense
  bool hasRegularGrid_;          //!< rtGReg_ make sense
  bool hasRTPars_;               //!< rtPars_ make sense

  std::vector< RTGridPoint > rtGrid_;   //!< rt relation
  std::vector< RTGridPoint > rtGReg_;   //!< regular rt relation
  std::vector< double > rtPars_;        //!< rtRelation parametrisation (_modified_ polynom)
  bool _FIT1_;                     //!< if true, fit function is polynom
  bool _FIT2_;                     //!< if true, fit function is hyperbolic tangeant.
  
  struct sortRTGridPoints_;      //!< structure used to sort grid points

  friend std::istream& operator>>(std::istream& in, CsRTRelation& rt); //!< used to read RT relations from DB
};

#endif 
