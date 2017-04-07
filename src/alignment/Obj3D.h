// $Id: Obj3D.h,v 1.10 2010/02/04 15:16:07 tnagel Exp $
#ifndef Obj3D_h
#define Obj3D_h
/*!
   \file    Obj3D.h
   \brief   base 3D object
   \author  Hugo Pereira
   \version $Revision: 1.10 $
   \date    $Date: 2010/02/04 15:16:07 $
*/

#include "Point.h"

#include <TMatrix.h>

#include <string>
#include <list>
#include <vector>
#include <iostream>

class TRotMatrix;
class TNode;
class TShape;
class Point;
//__________________________________________________________________________________________
/*!
   \class   Obj3D
   \brief   base 3D object
*/
class Obj3D {
  public:
  
  //! shapes definition
  enum Shape3D { pave, circle, polyline, tunnel };
  
  /*! \fn Obj3D( std::string TBName, Shape3D shape, Point center )
    \brief constructor with unit rotation matrix
    \param TBName name of the object
    \param shape one of the shape id
    \param center center position point object
  */
  Obj3D( std::string TBName, Shape3D shape, Point center );
  
  /*! \fn Obj3D( std::string TBName, Shape3D shape, Point center, TMatrix m )
    \brief constructor with specified rotation matrix
    \param shape one of the shape id
    \param center center position point object
    \param m rotation matrix wrt main reference system
  */
  Obj3D( std::string TBName, Shape3D shape, Point center, TMatrix m );
  
  virtual ~Obj3D( void ); //!< destructor
    
  //! create nodes associated to object, with specified color
  virtual void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1 ) { return; }
  
  void DeleteTObjects( void );   //!< to delete TNode/TShapes  pointed to by nodes_/shapes_
  bool Drawn( void ) const {return drawn_; } //!< tells if nodes/shapes have been built
  std::string TBName_;    //!< obj3d name
  
  protected:  
  Shape3D shape_;    //!< shape of current object
  TRotMatrix *rotM_; //!< rotation matrix wrt main reference system
  Point  center_ ;   //!< center position 
  bool drawn_;       //!< True if nodes are built

  std::vector< TNode* >  nodes_;  //!< vector of pointer to booked TNodes, via MakeNodes
  std::vector< TShape* > shapes_; //!< vector of pointer to booked TShapes, via MakeNodes

  friend std::ostream &operator << ( std::ostream &, const Obj3D & );
};

//__________________________________________________________________________________________
/*!
   \class   Pave3D
   \brief   rotated pave object (for rectangle dead zone/active area)
*/
class Pave3D: public Obj3D {
  public:
  
  /*! \fn Pave3D( std::string TBName, Point center, TMatrix m, double du, double dv, double dz )
    \brief constructor with specified rotation matrix
    \param TBName name of the object
    \param center center position point object
    \param m rotation matrix wrt main reference system
    \param du half width  (x) in local reference system
    \param dv half height (y) in local reference system
    \param dz half depth  (z) in local reference system
  */
  Pave3D( std::string TBName, Point center, TMatrix rotM,
    double du, double dv, double dz ): 
    Obj3D( TBName, pave, center, rotM ),
    du_(du), dv_(dv), dz_(dz)
    {}
    
  /*! \fn void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1)
    \brief create/draw nodes corresponding to object with given color and line width
    \param parent the parent node in wich new ones are to be created
    \param color color to be used for the object
    \param width line width
  */
  void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1 );
     
  double du_;    //!< half width DRS [mm]
  double dv_;    //!< half height DRS [mm]
  double dz_;    //!< half depth [mm]

};

//__________________________________________________________________________________________
/*!
   \class   Tunnel3D
   \brief   rotated tunnel object: tube with rectangular section. used for magnets
*/
class Tunnel3D: public Obj3D {
  public:
  
  /*! \fn Tunnel3D( std::string TBName, Point center, TMatrix m, 
    double duIn,  double dvIn, double duOut, double dvOut, double dz )
    \brief constructor with specified rotation matrix
    \param TBName name of the object
    \param center center position point object
    \param m rotation matrix wrt main reference system
    \param duIn internal half width  (x) in local reference system
    \param dvIn internal half height (y) in local reference system
    \param duOut outer half width  (x) in local reference system
    \param dvOut outer half height (y) in local reference system
    \param dz half depth  (z) in local reference system
  */
  Tunnel3D( std::string TBName, Point center, TMatrix rotM,
    double duIn, double dvIn, double duOut, double dvOut, double dz ): 
    Obj3D( TBName, tunnel, center, rotM ),
    duIn_(duIn), dvIn_(dvIn), duOut_(duOut), dvOut_(dvOut), dz_(dz)
    {}
    
  /*! \fn void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1)
    \brief create/draw nodes corresponding to object with given color and line width
    \param parent the parent node in wich new ones are to be created
    \param color color to be used for the object
    \param width line width
  */
  void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1 );
     
  double duIn_;    //!< inner half width DRS [mm]
  double dvIn_;    //!< inner half height DRS [mm]
  double duOut_;   //!< outerhalf width DRS [mm]
  double dvOut_;   //!< outer half height DRS [mm]
  double dz_;    //!< half depth [mm]

};

//__________________________________________________________________________________________
/*!
   \class   Circle3D
   \brief   rotated pave object (used for rectangle, dead zone/active area)
*/
class Circle3D: public Obj3D {
  public:
  
  /*! \fn Circle3D( std::string TBName, Point center, TMatrix m, double r, double dz )
    \brief constructor with specified rotation matrix
    \param TBName name of the object
    \param center center position point object
    \param m rotation matrix wrt main reference system
    \param r radius (in x,y plane)
    \param dz half depth  (z) in local reference system
  */
  Circle3D( std::string TBName, Point center, TMatrix rotM,
    double r, double dz ): 
    Obj3D( TBName, circle, center, rotM ),
    r_(r), dz_(dz)
    {}
  
  /*! \fn void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1)
    \brief create/draw nodes corresponding to object with given color and line width
    \param parent the parent node in wich new ones are to be created
    \param color color to be used for the object
    \param width line width
  */
  void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1 );
     
  double r_;     //!< radius perp to the beam [mm]
  double dz_;    //!< half depth [mm]

};

//__________________________________________________________________________________________
/*!
   \class   PolyLine3D
   \brief   polyline object, used for tracks
*/
class PolyLine3D: public Obj3D {
  public:
    
  /*! \fn PolyLine3D( std::string TBName )
    \brief empty constructor 
    \param TBName name of the object
  */
  PolyLine3D( std::string TBName ):Obj3D( TBName, polyline, Point::Orig ) {}
      
  /*! \fn PolyLine3D( std::string TBName,std::list< Point > points )
    \brief full constructor
    \param TBName name of the object
    \param points list of Point objects used in between which segments are drawn
  */
  PolyLine3D( std::string TBName, std::list< Point > points ):Obj3D( TBName, polyline, Point::Orig ), points_( points ) {}
  
  void AddPoint( Point p ) { points_.push_back( p ); } //!< add point to list of points
  void Sort( void )        { points_.sort(); }         //!< sort points according to z
  unsigned int Size()      { return points_.size(); }  //!< get the number of points
    
  /*! \fn void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1 )
    \brief create/draw nodes corresponding to object with given color and line width
    \param parent the parent node in wich new ones are to be created
    \param color color to be used for the object
    \param width line width
  */
  void MakeNodes( TNode* parent, unsigned int color = 0, unsigned int width = 1 );

  std::list< Point > points_; //!< list of polyline points

};

#endif

