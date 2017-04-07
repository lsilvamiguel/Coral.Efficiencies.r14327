// $Id: THlx.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef THlx_h
#define THlx_h


#include <cassert>
#include <cmath>
#include <iostream>

#include "TMtx.h"
#include "TOpt.h"

class THit;

/*!
  \class THlx
  \brief Track parameters
 
  Class for track state vector (Helix)
  with corresponding Cov matrix
  (lower triangle only)

*/

class CsHelix;

class THlx {
  
public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  THlx();
  THlx(const THlx&  H);
  ~THlx();
  
  // Overloaded operators
  THlx& operator = (const THlx&);  
  THlx  operator - (const THlx&);  
  double & operator () (const int i, const int j=0);       //!< "()" accessor
  const double & operator () (const int i, const int j=0) const {return const_cast<THlx*>(this)->operator()(i,j);} //!< "()" accessor
  THlx& operator *= (const double k);
  
  //Accessors / converters
  double dip(bool in_deg = false) const;
  double azi(bool in_deg = false) const;
  double DirCos(int) const;
  double Mom() const;
  void Set(double& x, TMtx& V, TMtx& M);
  void Get(double& x, TMtx& V, TMtx& M) const;
  CsHelix ExportHelix();
  void ImportHelix(const CsHelix& h);

  //Methods
  void Print(const char* str="") const;    //!< Print
  bool empty() const;                 //!< "true" it's empty helix (not yet defined)
  bool empty_cov() const;             //!< "true" if the helix cov. matrix is not defined
  bool with_mom() const;              //!< "true" if momentum had been defined

  void Rotate(double cosang, double sinang, THlx& Hout); //!< Rotate


  /*! 
    Extrapolate "this" helix and covariance matrix in magnetic field to the 
    position X = Hout(0) (has to be prepared before call to the function).
    By default, material map is used (if its switched ON by ReMode[10] option)
    If hidden parameter "mmap" is set to "false", maretial map is not used (to speed up)
    "path", "radLenFr" and "eLoss" data members are set by this method 
  */
  bool Extrapolate(THlx& Hout, bool mmap=true) const;       
  /*!
    Extrapolate helix w/ Multi-Scattering from material maps AND from detectors outside maps (Note that the latter is not taken into account in standard "THlx::Extrapolate" and may be signficant, e.g. in DRs equipped w/ lead).
  */
  bool Propagate(THlx &Hout) const;
  double Path()     const {return path;}       // Returns trajectory length of previous extrapolation
  double RadLenFr() const {return radLenFr;}   // Returns fraction of radiation length: it's either that corresponding to latest extrapolation or that explicitly assigned to the helix by "setRadLenFr".
  double ELoss() const {return eLoss;}         // Returns energy loss cumulated during latest extrapolation.
  //                                              Mutators
  void SetqP(double qP) { Hpar[5] = qP; }
  void SetPath(double l) { path = l; }
  void SetRadLenFr(double f) { radLenFr = f; }
  void SetELoss(double el) { eLoss = el; }

  // also "extrapolations" but not for general use.
  bool Extrap(double x, THlx& Hout) const;  //!< Extrapolate (new: single pass for cov matrix propagation)
  bool NpassExtrap(THlx& Hout) const;       //!< Extrapolate (old: 4 passes for cov matrix propagation)

  void AddNoise(float x, float RadL); //!< Add mult. scatt. "noise" to "this" helix 
  void AddNoise(int ipl);             //!< Add mult. scatt. "noise" on the plane ipl to "this" helix           

  double HitChi2(const THit& h);      //!< Estimate contribution of mesurement h to Chi2 
  double HlxChi2(const THlx& H);      //!< Chi2 of difference between "this" helix and helix H 

  bool Update(const THit& h, THlx& Hupd) const;                 //!< Update helix by mesurement h
  bool Update(const THlx& M, THlx& Hupd, double& dChi2) const;  //!< Update helix by another helix

  double Dist(THlx& H) const;    //!< Straight line distance from "this" helix (X,Y,Z) to helix H (X,Y,Z)

  /*! 
    Find "this" helix and helix H at point of closest distance of aproach (CDA).  
    Returns "false" in the CDA has not found within the range of search. 
    \param X0  0-th approximation of CDA X position
    \param Xmin, Xmax - range of search

  */
  bool FindCDA(THlx& H, float X0 = 0, float Xmin = -1000, float Xmax = +1000);

private:
  
  double Hpar[6]; //! Track parameters: x0(fixed), x1, x2, dx1/dx0, dx2/dx0, q/P
  double Hcov[15]; //! Cov matrix lower triangle (15 elements)

  mutable double path;   
  mutable double radLenFr;
  mutable double eLoss;    // Note (Y.B.): is the mutable attribute useful?

};

//! default constructor
inline THlx::THlx():
  path    (-1.),
  radLenFr(-1.),
  eLoss   (0)
{
  NobjCreated++;
  for (int i=0; i<6;  i++) Hpar[i]=0;

  Hcov[0]= 0;
  Hcov[1] =0; Hcov[2] =0;
  Hcov[3] =0; Hcov[4] =0; Hcov[5] =0;
  Hcov[6] =0; Hcov[7] =0; Hcov[8] =0; Hcov[9] =0;
  Hcov[10]=0; Hcov[11]=0; Hcov[12]=0; Hcov[13]=0; Hcov[14] =0;
};

//! copy constructor
inline THlx::THlx(const THlx& H)
{
  NobjCreated++;
  for (int i=0; i<6;  i++) Hpar[i]=H.Hpar[i];
  for (int i=0; i<15; i++) Hcov[i]=H.Hcov[i];
  
  path     = H.path;
  radLenFr = H.radLenFr;
  eLoss    = H.eLoss;
};

//! destructor
inline THlx::~THlx()
{
  NobjDestructed++;
};


//! "=" operator
inline THlx& THlx:: operator = (const THlx& H)
{
  for (int i=0; i<6;  i++) this->Hpar[i]=H.Hpar[i];
  for (int i=0; i<15; i++) this->Hcov[i]=H.Hcov[i];
  path     = H.path;
  radLenFr = H.radLenFr;
  eLoss    = H.eLoss;

  return(*this);
};

/*!
   Helix(i) - accessor to i-th element of thestate vector:

   i=0  is X0, i=1-5 are y0,z0,yp,zp,charge/p;  

   Helix(i,j) - accessor to i,j element of Cov matrix. I,j = 1...5

   (Helix(i,j)=Helix(j,i))
*/
inline double& THlx::operator () (const int i, const int j)
{
  if(j==0){
    if((i < 0) || (i > 5)) {
      std::cout<<"THlx::double& operator()  ==> Index out of range : "<<i<<std::endl;
      assert(false);
    }
  } else {
    if((i <= 0) || (i > 5)||(j <= 0) || (i > 5)) {
      std::cout<<"THlx::double& operator()  ==> Indecies out of range : "<<i<<" "<<j<<std::endl; 
      assert(false);
    }
  }
  if(j==0){
    return(Hpar[i]);
  }
  else{
    return ( j > i ? Hcov[(i-1)+((j-1)*j)/2] : Hcov[(j-1)+((i-1)*i)/2] );
  } 
};

/*!
 "*=" operator to scale cov matrix of the helix.
 Nondiagonal elements are set to 0
*/

inline THlx& THlx::operator *= (const double k)
{
  Hcov[0]*=k;
  Hcov[1] =0; Hcov[2]*=k;
  Hcov[3] =0; Hcov[4] =0; Hcov[5]*=k;
  Hcov[6] =0; Hcov[7] =0; Hcov[8] =0; Hcov[9]*=k;
  Hcov[10]=0; Hcov[11]=0; Hcov[12]=0; Hcov[13]=0; Hcov[14]*=k;
  return(*this);
};

//! - operator

inline THlx THlx::operator - (const THlx& H)
{
  THlx Res;
  for(int i=0; i < 6;  i++) Res.Hpar[i] = Hpar[i] - H.Hpar[i];
  for(int i=0; i < 15; i++) Res.Hcov[i] = Hcov[i] + H.Hcov[i];
  return(Res);
};

//! Functions return Dip angle
inline double THlx::dip(bool in_deg) const
{ 
  double     c = 1.;
  if(in_deg) c = 57.295779513;
  return(c * atan(Hpar[4]/sqrt(1.+Hpar[3]*Hpar[3])));
}
//! Functions return Azi angle
inline double THlx::azi(bool in_deg) const
{
  double     c = 1.;
  if(in_deg) c = 57.295779513;
  return(c * atan(Hpar[3]));
}

//! Functions return direction cosine i (i=1,2,3)
inline double THlx::DirCos(int i) const
{
  if(i < 1 || i > 3){
    std::cout<<"Wrong THlx::DirCos argument"<<std::endl;
    return(0);
  }
  double a = sqrt(1.+ Hpar[3]*Hpar[3] + Hpar[4]*Hpar[4]);
  if(i == 1){
    return(1./a);
  } else {
    return(Hpar[i+1]/a);
  }
};

//! Functions return momenum
inline double THlx::Mom() const
{
  if(Hpar[5] == 0.) return(10000.);
  else              return(fabs(1./Hpar[5]));
}

inline bool THlx::empty() const
{
  double sum=0;
  for (int i=0; i<6;  i++) sum+=Hpar[i];
  for (int i=0; i<15; i++) sum+=Hcov[i];
  if(sum == 0.) return (true);
  else          return(false);
}

inline bool THlx::empty_cov() const
{
  double sum=0.;
  for (int i=0; i<15; i++) sum+=Hcov[i];
  if(sum == 0.) return (true);
  else          return(false);
}

inline bool THlx::with_mom() const
{
  if(Hpar[5] != 0.) return (true);
  else              return(false);
}

inline double THlx::Dist(THlx& H) const
{
  return(sqrt((Hpar[0]-H.Hpar[0])*(Hpar[0]-H.Hpar[0])+
	      (Hpar[1]-H.Hpar[1])*(Hpar[1]-H.Hpar[1])+
	      (Hpar[2]-H.Hpar[2])*(Hpar[2]-H.Hpar[2])));
}

#endif















