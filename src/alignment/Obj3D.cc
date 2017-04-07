// $Id: Obj3D.cc,v 1.14 2010/05/07 09:14:11 suhl Exp $
/*!
   \file    Obj3D.cc
   \brief   base 3D object
   \author  Hugo Pereira
   \version $Revision: 1.14 $
   \date    $Date: 2010/05/07 09:14:11 $
*/

#include "Obj3D.h"
#include "Point.h"
#include <TMath.h>
#include <TRotMatrix.h>
#include <iostream>

#include <TNode.h>
#include <TShape.h>
#include <TTUBE.h>
#include <TBRIK.h>
#include <TTRD1.h>
#include <TGeometry.h>

using namespace std;

//______________________________________________________________________________
Obj3D::Obj3D( string TBName, Shape3D shape, Point center ):
  TBName_( TBName ), shape_(shape), center_( center ), drawn_( false )
{
  //=== get matrix elements
  double *rot = new double[9];
  rot[0] = 1; rot[1] = 0; rot[2] = 0;
  rot[3] = 0; rot[4] = 1; rot[5] = 0;
  rot[6] = 0; rot[7] = 0; rot[8] = 1;
  
  rotM_ = new TRotMatrix( "rotM", "rotM", rot );
  SafeDelete( rot );
}

//______________________________________________________________________________
Obj3D::Obj3D( string TBName, Shape3D shape, Point center, TMatrix m ):
  TBName_( TBName ), shape_(shape), center_( center ), drawn_( false )
{
  
  //=== get matrix elements
  double *rot = new double[9];
  if( m.GetNrows() == 3 && m.GetNcols() == 3 ) {
    TMatrix im( TMatrix::kInverted, m );
    rot[0] = im(0,0); rot[1] = im(0,1); rot[2] = im(0,2);
    rot[3] = im(1,0); rot[4] = im(1,1); rot[5] = im(1,2);
    rot[6] = im(2,0); rot[7] = im(2,1); rot[8] = im(2,2);
  } else {
    cout << "Obj3D::Obj3D - wrong size for rotation matrix.\n";
    rot[0] = 1; rot[1] = 0; rot[2] = 0;
    rot[3] = 0; rot[4] = 1; rot[5] = 0;
    rot[6] = 0; rot[7] = 0; rot[8] = 1;
  }
  
  rotM_ = new TRotMatrix( "rotM", "rotM", rot );
  SafeDelete( rot );
}

//______________________________________________________________________________
Obj3D::~Obj3D( void )
{ 
  DeleteTObjects(); 
  SafeDelete( rotM_ );
}

//______________________________________________________________________________
void Obj3D::DeleteTObjects( void )
{ 
  // delete TNodes
  for( unsigned int i=0; i<nodes_.size(); i++ ) 
  SafeDelete( nodes_[i] ); 
  nodes_.clear(); 

  // delete TShapes
  for( unsigned int i=0; i<shapes_.size(); i++ ) 
  SafeDelete( shapes_[i] ); 
  shapes_.clear(); 

  drawn_ = false;
}

//______________________________________________________
//! Used to dump informations needed for the alignment
ostream &operator << (ostream &out,const Obj3D &o) 
{
  const Pave3D*     po = (const Pave3D*)     &o;
  const Tunnel3D*   to = (const Tunnel3D*) &o;
  const Circle3D*   co = (const Circle3D*)   &o;
  const PolyLine3D* lo = (const PolyLine3D*) &o;
 
  out << o.center_ << " - ";
  switch( o.shape_ ) {
    case Obj3D::pave:     printf("%s pave     du=%13.4f, dv=%13.4f, dz=%13.4f", o.TBName_.c_str(), po->du_, po->dv_, po->dz_ ); break;
    case Obj3D::tunnel:   printf("%s tunnel   duIn=%13.4f, dvIn=%13.4f, duOut=%13.4f, dvOut=%13.4f, dz=%13.4f", 
      o.TBName_.c_str(), 
      to->duIn_,  to->dvIn_, 
      to->duOut_, to->dvOut_, 
      to->dz_ ); break;
    case Obj3D::circle:   printf("%s circle   r =%13.4f, dz=%13.4f",            o.TBName_.c_str(), co->r_,  co->dz_ ); break;
    case Obj3D::polyline: printf("%s polyline %zu", o.TBName_.c_str(), lo->points_.size() ); break;
    default: out << o.TBName_ << " unknown";
  }
  return out;
}

//____________________________________
void Pave3D::MakeNodes( TNode* parent, unsigned int color, unsigned int width )
{
  
  if( parent ) parent->cd();
  TBRIK *brik = new TBRIK( "brik", "brik", "void", du_,  dv_,  dz_ );
  shapes_.push_back( (TShape*) brik );
  
  TNode *node = new TNode( TBName_.c_str(), TBName_.c_str(), brik, center_.x, center_.y, center_.z );
  node->SetMatrix( rotM_ );
  if( color ) node->SetLineColor( color );
  node->SetLineWidth( width );
  nodes_.push_back( node );
  drawn_ = true;
  if( parent ) parent->cd();
  return;
}

//____________________________________
void Tunnel3D::MakeNodes( TNode* parent, unsigned int color, unsigned int width )
{
  if( parent ) parent->cd();
  TBRIK *brik = new TBRIK( "brik", "brik", "void", 0.1,  0.1,  0.1 );
  TTRD1* trap1 = new TTRD1( "lr","lr", "void", dvIn_, dvOut_, dz_, 0.5*(duOut_-duIn_) );
  TTRD1* trap2 = new TTRD1( "tb","tb", "void", duIn_, duOut_, dz_, 0.5*(dvOut_-dvIn_) );
  shapes_.push_back( (TShape*) trap1 );
  shapes_.push_back( (TShape*) trap2 );
  shapes_.push_back( (TShape*) brik );
    
  TNode* node = new TNode( TBName_.c_str(), TBName_.c_str(), brik, center_.x, center_.y, center_.z );
  nodes_.push_back( node );
  node->cd();
  
  TRotMatrix *rot;
  string name;
    
  name = TBName_+"_left";
  double rotl[] = { 0,-1,0,0,0, 1,-1,0,0 };
  rot =  new TRotMatrix( "rot", "rot", rotl );
  TNode *nodel = new TNode( name.c_str(), name.c_str(), trap1, -0.5*(duIn_+duOut_), 0, 0 );
  nodel->SetMatrix( rot );
  if( color ) nodel->SetLineColor( color );
  nodel->SetLineWidth( width );
  
  name = TBName_+"_right";
  double rotr[] = { 0,-1,0,0,0,-1,1,0,0 };
  rot =  new TRotMatrix( "rot", "rot", rotr );
  TNode *noder = new TNode( name.c_str(), name.c_str(), trap1,  0.5*(duIn_+duOut_), 0, 0 );
  noder->SetMatrix( rot );
  if( color ) noder->SetLineColor( color );
  noder->SetLineWidth( width );
  
  name = TBName_+"_top";
  double rott[] = { 1, 0, 0, 0, 0, -1, 0, 1, 0 };
  rot =  new TRotMatrix( "rot", "rot", rott );
  TNode *nodet = new TNode( name.c_str(), name.c_str(), trap2, 0, 0.5*(dvIn_+dvOut_), 0 );
  nodet->SetMatrix( rot );
  if( color ) nodet->SetLineColor( color );
  nodet->SetLineWidth( width );
  
  name = TBName_+"_bottom";
  double rotb[] = { 1, 0, 0, 0, 0, 1, 0, -1, 0 };
  rot =  new TRotMatrix( "rot", "rot", rotb );
  TNode *nodeb = new TNode( name.c_str(), name.c_str(), trap2, 0, -0.5*(dvIn_+dvOut_), 0 );
  nodeb->SetMatrix( rot );
  if( color ) nodeb->SetLineColor( color );
  nodeb->SetLineWidth( width );
  drawn_ = true;
  if( parent ) parent->cd();
  return;
}

//____________________________________
void Circle3D::MakeNodes( TNode* parent, unsigned int color, unsigned int width )
{
  if( parent ) parent->cd();
  TTUBE *tube = new TTUBE( "tube", "tube", "void", r_,  dz_ );
  TNode *node = new TNode( TBName_.c_str(), TBName_.c_str(), tube, center_.x, center_.y, center_.z );
  node->SetMatrix( rotM_ );
  if( color ) node->SetLineColor( color );
  node->SetLineWidth( width );
  nodes_.push_back( node );
  shapes_.push_back( (TShape*) tube );
  drawn_ = true;
  if( parent ) parent->cd();
  return;
}

//____________________________________
void PolyLine3D::MakeNodes( TNode* parent, unsigned int color, unsigned int width )
{
  if( parent ) parent->cd();
  if( points_.size() < 2 ) return;
  
  points_.sort();
  list< Point >::iterator p0 = points_.begin(); 
  list< Point >::iterator p1 = points_.begin(); p1++;
  
  while( p1 != points_.end() ) {
    if( p1->z == p0->z ) { 
      cout << "PolyLine3D::MakeNodes - ERROR: Two points with same z. point skipped.\n";
      p1++; 
      continue; }
    double l = sqrt( 
      pow(p1->x - p0->x, 2) + 
      pow(p1->y - p0->y, 2) + 
      pow(p1->z - p0->z, 2) );

    Point center( 0.5*(p1->x + p0->x), 0.5*(p1->y + p0->y),  0.5*(p1->z + p0->z) );
    Point vect( (p0->x - p1->x)/l, (p0->y - p1->y)/l,  (p1->z - p0->z)/l );
    
    //=== Build direct rotation matrix to change [0,0,1] into vect
    double x = vect.x,   x2 = pow( vect.x, 2 );
    double y = vect.y,   y2 = pow( vect.y, 2 );
    double z = vect.z+1, z2 = pow( vect.z+1, 2 );
    double n2 = x2 + y2+ z2;
    double rot[] = {
      1.0-2.0*x2/n2, -2.0*x*y/n2,   2.0*x*z/n2, 
      -2.0*x*y/n2,   1.0-2.0*y2/n2, 2.0*y*z/n2, 
      -2.0*x*z/n2,   -2.0*y*z/n2,    -1.0+2.0*z2/n2 };

    rotM_ = new TRotMatrix( "rotM", "rotM", rot );

    TTUBE *tube = new TTUBE( "tube", "tube", "void", 1,  0.5*l );
    TNode *node = new TNode( TBName_.c_str(), TBName_.c_str(), tube, center.x, center.y, center.z );
    node->SetMatrix( rotM_ );
    
    if( color ) node->SetLineColor( color );
    node->SetLineWidth( width );
    nodes_.push_back( node );
    shapes_.push_back( (TShape*) tube );
      
    p0 = p1;
    p1++;
  }
  drawn_ = true;
  if( parent ) parent->cd();
  return;
}
