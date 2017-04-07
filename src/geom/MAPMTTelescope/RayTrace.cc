// RayTrace.cc 
// 
// written by Andreas Mutter, 29.10.2005
// some changes (to avoid memory leaks) applied, 21.11.2005
// 
// added OptSyst::FindFirstLens, 15.12.2005 
//
// The function OptSyst::RayTrace to be called within CORAL 
// simulates the passage of photons through  
// the COMPASS RICH MAPMT "telescope".
// 
// 
//
// Please direct questions to andreas.mutter@cern.ch

#include <iostream>
#include <cmath>
#include "RayTrace.h"

using namespace std;

// global variables, should all be 0 if only Pin ID is required
int debug=0;
int printxy=0;
int verbose=0;

OptSyst::OptSyst()
{
  t[0]=new Sphere;
  t[1]=new OpticalPlane;
  t[2]=new CoordinateBreak;
  t[3]=new Asphere;
  t[4]=new Sphere;
  t[5]=new CoordinateBreak;
  t[6]=new OpticalPlane;
  t[7]=new OpticalPlane;

// first (spherical) surface of the field lens
  t[0]->SetLoc(1e-3,0,0);         // set a little off the origin to make sure that rays start before the surface
  t[0]->SetRadius(33.94);
  t[0]->SetSurf(1./54.93746,0,0);
  t[0]->SetSell(6.961663000E-001, 4.679148000E-003, 4.079426000E-001, 1.351206300E-002, 8.974794000E-001, 9.793400250E+001);

// back surface of the field lens (prismatic)
  t[1]->SetLoc(20.81269,0,0);
  t[1]->SetTilt(-5.,0,0);

// coordinate break
  t[2]->SetLoc(75.80508,0,-4);
  t[2]->SetTilt(5.4604757,0,0);

// first (aspherical) surface of the condensor lens 
  t[3]->SetLoc(0,0,0);
  t[3]->SetRadius(17.33948);
  t[3]->SetSurf(1./20.69576,0,-6.1388451e-005);
  t[3]->SetSell(6.961663000E-001, 4.679148000E-003, 4.079426000E-001, 1.351206300E-002, 8.974794000E-001, 9.793400250E+001);

// back surface of the condensor lens (spherical)
  t[4]->SetLoc(19.42022,0,0);
  t[4]->SetSurf(1./(-24.96372),0,0);
  t[4]->SetRadius(16.54558);

// coordinate break (PMT is decentered)
  t[5]->SetLoc(26.83498,0,1.5);
  

// the PMT window
  t[6]->SetLoc(0,0,0);
  t[6]->SetSell(6.961663000E-001, 4.679148000E-003, 4.079426000E-001, 1.351206300E-002, 8.974794000E-001, 9.793400250E+001);
  
  t[7]->SetLoc(1.5,0,0);
}


int OptSyst::RayTrace(LightRay &photon)
{
  int ipin=-1;
  double x, y, ax, ay, e;

  if(photon.getval(6)<=0){
    cerr << "Photon direction not valid!" << endl;
    return -1;
  }
  
  x=photon.getval(1)/10.; // convert to cm
  y=photon.getval(2)/10.;
  ax=photon.getval(4);
  ay=photon.getval(5);
  e=photon.getval(7);

  const double hc=1239.84186;    // conversion from wavelength to energy
  if(e>0.) e=hc/e;
  else
    {
      cerr << "Photon energy not valid!" << endl;
      return -1;
    }

  ipin=RayTrace(x, y, ax, ay, e);

  return ipin;
}


int OptSyst::RayTrace(double x, double y, double ax, double ay, double e)
{
  // x and y in cm!
  LightRay photon;
  int status=0;
  int ipin=-1, ipx, ipy;
  double xout, yout;

  const double pixelsize=4.5;    // size of PMT pixel in mm
                                 // no dead zone between pixels 
                                 // is implemented, yet.
  //  const int num_of_surf=8;       // total number of optical surfaces
  const double hc=1239.84186;    // conversion from energy to wavelength

  double start[3];
  double dir[3];
  double wave;



  if(e>0.) wave=hc/e;
  else
    {
      cerr << "Photon energy not valid!" << endl;
      return -1;
    }

// convert input variables x and y to mm
  start[0]=x*10.;
  start[1]=y*10.;
// find the corresponding z on the first surface
  start[2]=t[0]->SurfSag(start[0],start[1]);

  dir[0]=ax;
  dir[1]=ay;
// find z coordinate of the photon direction
  dir[2]=(1-dir[0]*dir[0]-dir[1]*dir[1]);
  if(dir[2]<0.) 
    {
      cerr << "Photon Direction not valid!" << endl;
      return -1;
    }
  dir[2]=sqrt(dir[2]);
  
  photon.SetWave(wave);
  photon.SetStart(start);
  photon.SetDir(dir);
  //  photon.norm();

  if(debug)
    {
      cout << endl << "===============================================================================================================" << endl;
      photon.print();
    }
  
// k should be <8...
  for (int k=0 ; k<7 ; k++)
    {
   t[k]->debug=debug;   
   if(debug)
	{      
          cout << "---------------------------------------------------------------------------------------------------------------" << endl;
          cout << "Surface no. " << k+1 << endl;
          cout << "---------------------------------------------------------------------------------------------------------------" << endl; 
	}
      status=t[k]->Refract(photon);
      if(status) 
	{
	  cerr << "Error in Refract at surface "<< k+1 << ", aborting function RayTrace! (status == " << status << ")"<< endl;
	  if(printxy) cout << "-999.  -999." << endl;
	  return ipin;
	}
    }
  
  if(debug) cout << "---------------------------------------------------------------------------------------------------------------" << endl;
  
// the coordinates on the PMT  
//  xout= photon.getval(1);
//  yout= photon.getval(2);
//
// Some deviations from the ZEMAX result (up to 0.2mm) have been found.
// Try to correct them by adding some linear terms in ax, ay and x, y:    
// This has been empirically tested...
  xout= photon.getval(1) - 1.25*ax - 0.004*x*10.;
  yout= photon.getval(2) -      ay - 0.006*y*10.;


  if(printxy) cout << xout << "\t" <<  yout << endl;
  
  if(debug)  cout << "===============================================================================================================" << endl;


  ipx=(int)floor((xout+2.*pixelsize)/pixelsize);
  ipy=(int)floor((yout+2.*pixelsize)/pixelsize);

//  if(printxy) cout << ipx << "\t" << ipy << endl;  

  if( ipx < 0 || ipx > 3 || ipy < 0 || ipy > 3 )
  {
    cerr << "Unexpected pins ipx = " << ipx <<" ipy = " << ipy << endl;
    cerr << " xout = " << xout << " yout = " << yout << endl;
    return ipin;
   }
  ipin = ipx + 4*ipy;
    
  return ipin;
}


// Function definitions for class OpticalSurface
// Constructor sets everything to 0
OpticalSurface::OpticalSurface()
{
  debug=0;
  dist_to_previous=0.;
  xshift=0.;
  yshift=0.;
  xtilt=0.;
  ytilt=0.;
  ztilt=0.;  
  radius=0.;  
  curv=0.;    
  asph1=0.;
  asph2=0;  
  for(int k=0 ; k<3; k++)
    {
      L[k]=0.;
      K[k]=0.;
    }              
}

// set the Sellmeier coefficients (cf. ZEMAX manual)
void OpticalSurface::SetSell(double a, double b, double c, double d, double e, double f)
{
  L[0]=b;
  L[1]=d;
  L[2]=f;
  K[0]=a;
  K[1]=c;
  K[2]=e;
}


// set different values cocerning geometry etc.
void OpticalSurface::SetLoc(double z, double x, double y)
{
  dist_to_previous=z;
  xshift=x;
  yshift=y;
}
void OpticalSurface::SetTilt(double x, double y, double z)
{
  xtilt=x;
  ytilt=y;
  ztilt=z;
}
void OpticalSurface::SetSurf(double c, double a1, double a2)
{
  curv=c;
  asph1=a1;
  asph2=a2;
}
void OpticalSurface::SetRadius(double r)
{
  radius=r;
}


// function for the shape of the surface
double OpticalSurface::SurfSag(double x, double y)
{
  double z, rsq;
  rsq=x*x+y*y;
  z=curv*rsq/(1+sqrt(1-curv*curv*rsq))+asph1*rsq+asph2*rsq*rsq;
  return z;
}


// calculate index of refraction from wavelength using the Sellmeier 1 formula
double OpticalSurface::RefInd(double wav)
{
  double refind=0.;
  double wl;
  if((wav>180.) && (wav<3000.))
    {
      wl=wav/1000.;  // convert nm to microns
      for(int k=0 ; k<3 ; k++) refind=refind+(K[k]*wl*wl)/(wl*wl-L[k]);
      refind=sqrt(refind+1);
    }
  return refind;
}


// calculates intersection of ray with the plane in which the lens edge lies
int OpticalSurface::CheckVig(LightRay photon)
{
  double s[3], d[3], x,y,z, l;

  if(debug) cout << "OpticalSurface::CheckVig" << endl;

  if(debug) cout << "     radius: " << radius << endl;

  if(radius>0.)
    {
      for(int k=0 ; k<3 ; k++) 
	{
	  s[k]=photon.getval(k+1);
	  d[k]=photon.getval(k+4);
	}

      if(d[2]<=0.) return 2;

      z=SurfSag(radius,0);

      l=(z-s[2])/d[2];
      x=l*d[0]+s[0];
      y=l*d[1]+s[1];

      if(debug) cout << "     coordinates of intersection: " << x << "\t" << y << "\t" << z << endl; 

      if((x*x+y*y)>radius*radius) return 1;

      return 0;
    }
  else
    return 0;

}


int OpticalSurface::DoRefraction(LightRay &photon, double normal[3])
{
  if(debug) cout << "OpticalSurface::DoRefraction" << endl;
  double n1, n2, cosa, rho, sigma, e[3], a[3], wav;

  for(int k=0; k<3 ; k++) e[k]=photon.getval(k+4);

  cosa=0;
  for(int k=0; k<3 ; k++) cosa=cosa+e[k]*normal[k];

  n1=photon.getval(8);
  wav=photon.getval(7);
  n2=RefInd(wav);

  if(n2>0)
    {
      rho=n1/n2;
      sigma=(1+rho*rho*(cosa*cosa -1));

      if(sigma<0) return 2;

      sigma=sqrt(sigma)-rho*cosa;
      
      for(int k=0; k<3 ; k++) a[k]=rho*e[k]+sigma*normal[k];

      photon.SetDir(a);
      photon.SetInd(n2);

      if(debug) photon.print();
      
      return 0;
    }
  else return 1;
}


int OpticalSurface::FindNormal(double x, double y, double n[3])
{
  if(debug) cout <<  "OpticalSurface::FindNormal" << endl;
  double rsq;
  double length;
  double w;
  rsq=x*x+y*y;
  w=(1-curv*curv*rsq);
  if(w<0.) return 2;
  w=sqrt(w);
  n[0]= -curv*x - 2*asph1*x*w -4*asph2*x*rsq*w;
  n[1]= -curv*y - 2*asph1*y*w -4*asph2*y*rsq*w;
  n[2]= w;
  length=0;
  for(int k=0; k<3 ; k++) length=length+n[k]*n[k];
  length=sqrt(length);
  if(length>0.)
    {
      for(int k=0; k<3 ; k++) n[k]=n[k]/length;
      return 0;
    }
  else return 1;
}



int OpticalSurface::TransformRay(LightRay &photon)
// The "back" transformation function only supports rotations around x...
// The definition of the rotation differs from the ZEMAX convention:
// 1. Rotation around the "old" x-axis
// 2. Rotation around the "old" y-axis!
{
  if(debug) cout << "OpticalSurface::TransformRay" << endl;
  double z[3];
  long double x[3], y[3], ang, angy;
  const long double pi=3.14159265358979323846;

  ang=-(xtilt*pi)/180.;
  angy=-(ytilt*pi)/180.;

  x[0]=photon.getval(1)-xshift;
  x[1]=photon.getval(2)-yshift;
  x[2]=photon.getval(3)-dist_to_previous;
  
//   y[0]=(double)(x[0]);
//   y[1]=(double)(cos(ang)*x[1]-sin(ang)*x[2]);
//   y[2]=(double)(sin(ang)*x[1]+cos(ang)*x[2]);

  y[0]=(x[0]);
  y[1]=(cos(ang)*x[1]-sin(ang)*x[2]);
  y[2]=(sin(ang)*x[1]+cos(ang)*x[2]);
  
  z[0]=(double) y[0]*cos(angy) - sin(angy)*cos(ang)*y[2] + sin(angy)*sin(ang)*y[1];
  z[1]=(double) y[1]*cos(angy) - sin(angy)*sin(ang)*y[0] + cos(ang)*( 1-cos(angy) )*( cos(ang)*y[1] + sin(ang)*y[2] );
  z[2]=(double) y[2]*cos(angy) + sin(angy)*cos(ang)*y[0] + sin(ang)*( 1-cos(angy) )*( cos(ang)*y[1] + sin(ang)*y[2] );
  
  photon.SetStart(z);
  
  x[0]=photon.getval(4);
  x[1]=photon.getval(5);
  x[2]=photon.getval(6);
  
  //  y[0]=(double)x[0];
  //  y[1]=(double)(cos(ang)*x[1]-sin(ang)*x[2]);
  //  y[2]=(double)(sin(ang)*x[1]+cos(ang)*x[2]);

  y[0]=(x[0]);
  y[1]=(cos(ang)*x[1]-sin(ang)*x[2]);
  y[2]=(sin(ang)*x[1]+cos(ang)*x[2]);
  
  z[0]=(double) y[0]*cos(angy) - sin(angy)*cos(ang)*y[2] + sin(angy)*sin(ang)*y[1];
  z[1]=(double) y[1]*cos(angy) - sin(angy)*sin(ang)*y[0] + cos(ang)*( 1-cos(angy) )*( cos(ang)*y[1] + sin(ang)*y[2] );
  z[2]=(double) y[2]*cos(angy) + sin(angy)*cos(ang)*y[0] + sin(ang)*( 1-cos(angy) )*( cos(ang)*y[1] + sin(ang)*y[2] );
  
  
  
  photon.SetDir(z);
  
  if(debug) photon.print();
  
  return 0;
}

int OpticalSurface::TransformRayBack(LightRay &photon)
{
  if(debug) cout << "OpticalSurface::TransformRayBack" << endl;

  double y[3];
  long double x[3], ang;

  const long double pi=3.14159265358979323846;

  ang=(xtilt*pi)/180.;

  x[0]=photon.getval(1);
  x[1]=photon.getval(2);
  x[2]=photon.getval(3);
  
  y[0]=(double)(x[0]+xshift);
  y[1]=(double)(cos(ang)*x[1]-sin(ang)*x[2]+yshift);
  y[2]=(double)(sin(ang)*x[1]+cos(ang)*x[2]);

  photon.SetStart(y);

  x[0]=photon.getval(4);
  x[1]=photon.getval(5);
  x[2]=photon.getval(6);
  
  y[0]=(double)x[0];
  y[1]=(double)(cos(ang)*x[1]-sin(ang)*x[2]);
  y[2]=(double)(sin(ang)*x[1]+cos(ang)*x[2]);
 
  photon.SetDir(y);

  if(debug) photon.print();

  return 0;
}


// the following are different for the different types of surfaces
// the coordinate break only transforms the ray coordinate and does nothing else
int CoordinateBreak::Refract(LightRay &photon)
{
  int status;
  if(debug) cout << "CoordinateBreak::Refract" << endl;
  status=TransformRay(photon);
  if(status) 
    {
      cerr << "Error in TransformRay, aborting function Refract! (status == " << status << ")"<< endl;
      return status;
    }
  return status;
}


// Refract is the same for OpticalPlane, Sphere and Asphere, only FindIntersection will differ. 
// refraction consist of the following steps
// 1. transform the ray to the local coordinate system of the current surface 
//    (if no tilts and decenters are involved, this only changes the z coordinate
// 2. check if it will pass through the surface at all
// 3. find the point of intersection 
// 4. find the normal vector of the surface at the point of intersection
// 5. actually refract it (i.e. change its direction)
// 6. transform the ray coordinates to the previous coordinate system (except for the z-coordinate) 

int OpticalPlane::Refract(LightRay &photon)
{
  double normal[3]={0, 0, 0},x,y;
  int status=0;
  if(debug) cout << "OpticalPlane::Refract" << endl;

  status=TransformRay(photon);
  if(status) 
    {
      cerr << "Error in TransformRay, aborting function Refract! (status == " << status << ")"<< endl;
      return status;
    }

  status=CheckVig(photon);
  if(status) 
    {
      cerr << "Ray will be vignetted, aborting function Refract! (status == " << status << ")"<< endl;
      return status;
    }

  status=FindIntersection(photon);
  if(status) 
    {
      cerr << "Error in FindIntersection, aborting function Refract! (status == " << status << ")"<< endl;
      return status;
    }

  x=photon.getval(1);
  y=photon.getval(2);

  status=FindNormal(x,y,normal);
  if(status) 
    {
      cerr << "Error in FindNormal, aborting function Refract! (status == " << status << ")"<< endl;
      return status;
    }

  if(debug) cout << "     normal at " << x << "\t" << y << ":\t";
  for(int k=0; k<3; k++) if(debug) cout << normal[k] << "\t"; 
  if(debug) cout << endl;

  status=DoRefraction(photon, normal);
  if(status) 
    {
      cerr << "Error in DoRefraction, aborting function Refract! (status == " << status << ")"<< endl;
      return status;
    }

  status=TransformRayBack(photon);
  if(status) 
    {
      cerr << "Error in TransformRayBack, aborting function Refract! (status == " << status << ")"<< endl;
      return status;
    }
 
  return status;
}

// the following calculate the intersection point for the three different surfaces (plane, sphere, asphere)
// the direction of the ray is not changed, only its starting point
int OpticalPlane::FindIntersection(LightRay &photon)
{
  if(debug) cout << "OpticalPlane::FindIntersection" << endl;
  
  double s[3], d[3], x[3], l;

  for(int k=0 ; k<3 ; k++) 
    {
      s[k]=photon.getval(k+1);
      d[k]=photon.getval(k+4);
    }
  
  if(d[2]<=0.) return 2;
  
   l=-s[2]/d[2];
  x[0]=l*d[0]+s[0];
  x[1]=l*d[1]+s[1];
  x[2]=0.;
  
  if(debug) cout << "     coordinates of intersection: " << x[0] << "\t" << x[1] << endl; 
   
  photon.SetStart(x);

  if(debug) photon.print();

  return 0;
  
  }


int Sphere::FindIntersection(LightRay &photon)
{
  if(debug) cout << "Sphere::FindIntersection" << endl;

  double s[3], d[3], x[3], l, l1, l2, rs, sq;

  if(curv==0.) return 4;

  for(int k=0 ; k<3 ; k++) 
    {
      s[k]=photon.getval(k+1);
      d[k]=photon.getval(k+4);
    }
  
  if(d[2]<=0.) return 3;

  rs=0;
  for(int k=0 ; k<3 ; k++) rs=rs+s[k]*d[k];

  sq=0;
  for(int k=0 ; k<3 ; k++) sq=sq+s[k]*s[k];

  l=(-rs+d[2]/curv)*(-rs+d[2]/curv)-sq+(2*s[2])/curv;

  if(l<0) return 2;

// the quadratic equation in general has two solutions  
  l1=sqrt(l)-rs+d[2]/curv;
  l2=-sqrt(l)-rs+d[2]/curv;

  if(debug) cout << "     l1 and l2: " << l1 << "\t" << l2 << endl;

// we are only interested in the forward direction and the first hit
  if((l1<=l2))
    {
      if(l1>=0.) l=l1;
      else if(l2>=0.) l=l2; 
      else return 1;
    }
  if((l2<l1))
    {
      if(l2>=0.) l=l2;
      else if(l1>=0.) l=l1; 
      else return 1;
    }   

  if(debug)  cout << "     l: " << l <<  endl;

  for(int k=0 ; k<3 ; k++)  x[k]=l*d[k]+s[k];

  photon.SetStart(x);

  if(debug) photon.print();  

  return 0;
}


// this one uses a iterative method to find the intersection
int Asphere::FindIntersection(LightRay &photon)
{
  if(debug) cout << "Asphere::FindIntersection" << endl;


  double s[3], d[3], x[3], l, z, dist;
  int iter=0;

// finding the intersection within 1e-6 mm = 1nm should be precise enough...
  const double max_dist=1e-6;
// end the iterative process if it doesn't coverge
  const int max_iter=1000;

  // Asphere must have finite size!
  if(radius==0.) return 4;
  z=SurfSag(radius,0);

  for(int k=0 ; k<3 ; k++) 
    {
      s[k]=photon.getval(k+1);
      d[k]=photon.getval(k+4);
    }

// first step: find the intersection with the "edge plane"  
  if(d[2]<=0.) return 3;
  l=(z-s[2])/d[2];  
  for(int k=0 ; k<3 ; k++) x[k]=l*d[k]+s[k];

  
  do
    {
      iter++;
      if(iter > max_iter) return 2;    // emergency exit
      dist=x[2]-SurfSag(x[0],x[1]);    // distance in z from the preliminary "intersection point" to the surface
      if(debug) cout << "     iteration " << iter << ", l=" << l << ", dist=" << dist << endl;  
 // change the length of the vector from the photon starting point to the preliminary "intersection point" by
// the ration of the according z distances
      if(abs(x[2]-s[2])>0.) l=l*(1-dist/(x[2]-s[2]));  
      else return 1;
      for(int k=0 ; k<3 ; k++) x[k]=l*d[k]+s[k];
    } 
  while((abs(dist) > max_dist));
    
  photon.SetStart(x);

  if(debug) photon.print();


  return 0;
}



// the default ray propagates in vacuum along the z axis and has a wavelenght of 300nm
LightRay::LightRay() 
{ 
  for(int k=0 ; k<3 ; k++) sp[k]=0.;
  for(int k=0 ; k<2 ; k++) dir[k]=0.;
  dir[2]=1.;
  wl=300.;
  n=1.;
}


void LightRay::SetStart(double s[3])
{
  for(int k=0 ; k<3 ; k++) sp[k]=s[k];
} 

void LightRay::SetDir(double d[3])
{
  for(int k=0 ; k<3 ; k++) dir[k]=d[k];
} 

void LightRay::SetWave(double wav)
{
  wl=wav;
} 

void LightRay::SetInd(double ind)
{
  n=ind;
} 


void LightRay::print()
{
      for(int k=0 ; k<3 ; k++ ) cout << sp[k] << "\t";
      for(int k=0 ; k<3 ; k++ ) cout << dir[k] << "\t";
      cout << wl << "\t" << n << endl;
}

// normalize direction vector of light ray (usually this function is not needed)
void LightRay::norm()
{
  double length=0.;

  for(int k=0 ; k<3; k++) length=length+dir[k]*dir[k];
  
  length=sqrt(length);
  
  if(length>0)
    for(int k=0 ; k<3; k++) dir[k]=dir[k]/length;
  else
    cerr << "LightRay::norm: Direction vector has zero length!" << endl;
}

// return individual values 
double LightRay::getval(int k)
{
  if((k>0) && (k<4)) return sp[k-1];
  if((k>3) && (k<7)) return dir[k-4];
  if(k==7) return wl;
  if(k==8) return n;
  else return -99999.9;
}


int OptSyst::FindFirstLens(double xr, double yr, double zr, double xm, double ym, double zm, double xd, double yd, double e, LightRay &Photon)
// all lengths are in cm!
{

  LightRay testray;

  int lens=-1;
  int status=0;
  int ilx, ily; // naive lens numbers from CsI coordinates

  int top=0, jura=0, bottom=0, saleve=0;

  const double topend=192.; // borders of the CsI cathodes in MRS and cm
  const double sideend=60.5;
  const double distic=5.;  // z-distance in mm of lens center to CsI cathode

  double lx, ly;
  double d[3];
  double s[3];

  if(zm > 0 && zm <  topend) top=1;
  if(zm < 0 && zm > -topend) bottom=1;
  if(ym > 0 && ym <  sideend) jura=1;
  if(ym < 0 && ym > -sideend) saleve=1;

  if(!(top+bottom) || !(jura+saleve)){
    cerr << "OptSyst::FindFirstLens: Photon not on one of the four central cathodes!" << endl;
    if(verbose){      
      cout << "---------------------------------------------------------------------------------------" << endl;
      cout << "Photon not on one of the central cathodes" << endl;
      cout << "---------------------------------------------------------------------------------------" << endl;
    }
    return lens;
  }

  if(verbose){
    cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Cathode:                     \t";
    if(top) cout << "top ";
    if(bottom) cout << "bottom ";
    if(jura) cout << "jura" << endl;
    if(saleve) cout << "saleve" << endl;
  }




  const double pi=3.14159265358979323846;
  const double ls=4.8;         // size of the "pseudo-lens"
  const double hc=1239.84186;
  const double auxdist=50.;    // distance between auxiliary plane and CsI cathode in mm  

  const double lsx=47.7;
  const double lsy=44.8;
    
  ilx=(int)floor((xd+6.*ls)/ls);
  ily=(int)floor((yd+6.*ls)/ls);
  lx=xd+ls*(5.5-ilx);
  ly=yd+ls*(5.5-ily);
  if(verbose)  cout << "Lens and coord. from xd,yd:  \t" << (ilx+1)+ily*12 << "; \t" <<	\
    -lx*top+lx*bottom << ", " << -ly*bottom+ly*top << endl;

  static CoordinateBreak FirstLensLoc[12][12];
  for(int k=-6; k<6; k++){
    for(int i=-6 ; i<6 ; i++){
// distic in this case is assumed to be perpendicular to the CsI
// which it is probably not. Probaby we have to set it below with
// FirstLens.SetLoc(distic,0,0);
      FirstLensLoc[k+6][i+6].SetLoc(auxdist-distic,25.09+k*48.2,25.56+i*46.7);
        }
  }
  

  static CoordinateBreak FirstLensOri[2][2][2][2]; // top, bottom, jura, saleve
// the original values (ytilt wrong!!!)
//  const double xtilt=15.5;   // in deg
//  const double ytilt=5.;  // in deg
// to play around...
  const double xtilt=15.5;   // in deg
  const double ytilt=-5.;  // in deg


  // these are the orientations after the 180 degree rotations!
  FirstLensOri[1][0][1][0].SetTilt(-xtilt,ytilt,0); // top jura
  FirstLensOri[1][0][0][1].SetTilt(-xtilt,-ytilt,0); // top saleve 
  FirstLensOri[0][1][1][0].SetTilt(-xtilt,-ytilt,0); // bottom jura
  FirstLensOri[0][1][0][1].SetTilt(-xtilt,ytilt,0); // bottom saleve

  static Sphere FirstLens;
//  FirstLens.SetLoc(distic,0,0); // distance to CsI plane
  FirstLens.SetRadius(33.94);
  FirstLens.SetSurf(1./54.93746,0,0);

  const double csitilt=0.15*180./pi;
  static CoordinateBreak CsIupper;
  static CoordinateBreak CsIlower;
  CsIupper.SetTilt(csitilt,0,0);
  CsIlower.SetTilt(-csitilt,0,0);


  static OpticalPlane AuxiliaryOpticalPlane;
  AuxiliaryOpticalPlane.SetLoc(auxdist,0,0); //  auxdist mm above current CsI

  if(e>0.) Photon.SetWave(hc/e);
  else 
    {
      cerr << "OptSyst::FindFirstLens: Energy must be > 0!" << endl;
      return(-1);
    }
    

  if(verbose){
    d[0] = xm-xr;  
    d[1] = ym-yr;
    d[2] = zm-zr;
    s[0] = xm;
    s[1] = ym;
    s[2] = zm;
    Photon.SetDir(d);
    Photon.norm();
    Photon.SetStart(s);
    cout << "Photon parameters (in MRS):  \t";
    Photon.print();
  }
  


  // "backward" photon direction in MRS (variables renamed to correspond to DRS)
  d[2] = -xm+xr;  
  d[0] = -ym+yr;
  d[1] = -zm+zr;

  Photon.SetDir(d);
  if(debug) {
    cout << endl << "Unnormalised direction of photon starting on CsI plane:" << endl;
    Photon.print();
  }
  Photon.norm();
  if(debug) {
    cout << endl << "Normalised direction of photon starting on CsI plane:" << endl;
    Photon.print();
  }



  // Here, the backward photon direction will be transformed
  // to the DRS (= CsI frame) 
  if(top){
    CsIupper.Refract(Photon);
    if(debug) cout << endl << "Transforming to top cathode" << endl;
  }
  
  if(bottom){
    CsIlower.Refract(Photon);
    if(debug) cout << endl << "Transforming to bottom cathode" << endl;
  }
  
  s[0] = xd*10.;  // convert to mm
  s[1] = yd*10.;
  s[2] = 0;
  
  Photon.SetStart(s);

  if(debug) {
    cout << "Direction and origin of photon starting on CsI plane:" << endl;
    Photon.print();
  }


  // trace the photon to the auxil. plane and invert its direction
  // AuxiliaryOpticalPlane.debug=debug;
  AuxiliaryOpticalPlane.Refract(Photon);
  if(debug) {
    cout << endl << "Direction and point of incidence on auxiliary plane:" << endl;
    Photon.print();
  }
  for(int k=0 ; k<3 ; k++) d[k]= -Photon.getval(k+4);
  Photon.SetDir(d);

  if(debug) {
    cout << endl << "Direction and origin of photon starting on auxiliary plane:" << endl;
    Photon.print();
  }  

  int itemplens[9][3]; // stores lenses x and y number and if it is hit or not
  double tempdist[9];      // stores the distance the photon has travelled from the aux plane to the lens
  LightRay tempray[9];   // stores the photon coordinates on the aux plane and on the first lens
  LightRay saveray[9];
  int count=-1;

  for(int k=ilx-1 ; k<ilx+2 ; k++){
    for(int i=ily-1 ; i<ily+2 ; i++){
      
      count++;
      itemplens[count][1]=k;
      itemplens[count][2]=i;
      itemplens[count][3]=0;
      tempdist[count]=9999.;
      
      // check only lenses within valid range
      if(k >-1 && k < 12 && i >-1 && i < 12){
	
	if(debug) cout << endl << "Checking lens " << k << ", " << i << endl;
	
	testray=Photon;
	
	if(debug) cout << "-- photon starting at:               ";
	if(debug) testray.print();
	
	//	FirstLensLoc[k][i].debug=debug;
	FirstLensLoc[k][i].Refract(testray);
	
	if(debug) cout << "-- shifted to lens center:           ";
	if(debug) testray.print();		
	
	
	
	for(int j=0 ; j<3 ; j++) 
	  {
	    s[j]=testray.getval(j+1);
	    d[j]=testray.getval(j+4);
	  }
	s[0]=s[0]*bottom-s[0]*top;
	s[1]=s[1]*top-s[1]*bottom;
	//	s[2]=-s[2];
	testray.SetStart(s);
	d[0]=d[0]*bottom-d[0]*top;
	d[1]=d[1]*top-d[1]*bottom;
	d[2]=-d[2];
	testray.SetDir(d);
	if(debug) cout << "-- renamed coordinates (180deg rot): ";
	if(debug) testray.print();
	
	//	FirstLensOri[top][bottom][jura][saleve].debug=debug;
	FirstLensOri[top][bottom][jura][saleve].Refract(testray);
	// if(debug) Photon.print();

	if(debug) cout << "-- rotated to local lens system:     ";
	if(debug) testray.print();

	saveray[count]=testray;

	status=FirstLens.FindIntersection(testray);

	double xtemp, ytemp;
	xtemp=testray.getval(1);
	ytemp=testray.getval(2);

	// cout << xtemp << "\t" << ytemp << endl;

	if(status || abs(xtemp) > lsx/2. || abs(ytemp) > lsy/2.){
	  if(debug) cout << "-- no intersection!" << endl; 
	}
	else{
	  if(debug) cout << "-- intersection with first lens:     ";
	  if(debug) testray.print();
	  itemplens[count][3]=1;
	  tempray[count]=testray;
	  for(int j=0 ; j<3 ; j++) 
	    {
	      tempdist[count]=tempdist[count] + \
		(tempray[count].getval(j+1)-saveray[count].getval(j+1))*\
		(tempray[count].getval(j+1)-saveray[count].getval(j+1));
	    } 
	  tempdist[count]=sqrt(tempdist[count]);
	  if(debug) cout << "Distance travelled: " << tempdist[count] << endl;
	}
      }
    }
  }
  
  // Now, sort for minimum travelling distance:

  ilx=-1;
  ily=-1;
  double dist=9999.;

  for(int k=0 ; k< 3; k++){
    s[k]=0;
    d[k]=0;
  }
  Photon.SetDir(d);
  Photon.SetStart(s);
  
  
  for(count=0 ; count < 9 ; count++){
    if(debug){
      cout << count << "\t";
      cout << "Distance travelled to lens " << itemplens[count][1]	\
	   << ", " << itemplens[count][2] << ": " << tempdist[count];
      cout << endl;
    }
    
    if(itemplens[count][3]){
      if(tempdist[count]<dist){
	ilx=itemplens[count][1];
	ily=itemplens[count][2];
	dist=tempdist[count];
	Photon=tempray[count];
      }
    }
    
    
    
  }    

  if( ilx>=0 && ily >= 0) lens=(ily)*12+(ilx+1);
  
  if(verbose){
    cout << "Photon parameters (lens " << lens << "): \t" ;
    Photon.print();
    cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;
  }

  return lens;
}
