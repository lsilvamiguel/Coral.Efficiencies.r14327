#include "THlx.h"

using namespace std;

/*
*----------------------------------------------------------------------*
* Sort the elements of the array into itself in ascending order        *
*                                                                      *
* if IMOD = 1                                                          *
* The procedure sorts the elements of the array A(1:N)                 *
* into ascending order. Algorithm 271, Comm., ACM(1965) p.669          *
*                                                                      *
* if IMOD = 2                                                          *
* The procedure sorts the element almost sorted array A(1:N)           *
*                                                                      *
* Input parameters  : A    - one-dimension array of elements           *
*                            to be sorted                              *
*                     IN   - one-dimension array of elements           *
*                            to be sorted together with A              *
*                     N    - number of words to be sorted              *
* Author : I.Gavrilenko                                                *
*----------------------------------------------------------------------*
*/
void thesort(double a[], int in[], int n, int imod)
{
  //cout<<"INSIDE thesort...."<<endl<<flush;  
  /* Local variables */
  static double b, c;
  static int i, j, k, m;
  static double t, x;
  static int k1, k2, ib, ip, iq, it, lt[200], mt[200], ix, iq1;
  
  /* Parameter adjustments */
  --in;
  --a;
  
  /* Function Body */
  if (n <= 1) {
    return;
  }
  if (imod == 2) {
    goto L100;
  }
  
  /*   First algorithm */
  
  j = n;
  i = 1;
  m = 1;
 L6:
  if (j - i <= 1) {
    goto L7;
  }
  ip = (j + i) / 2;
  t = a[ip];
  it = in[ip];
  a[ip] = a[i];
  in[ip] = in[i];
  iq = j;
  
  for (k = i + 1; k <= iq; ++k) {
    if (a[k] <= t) {
      goto L71;
    }
    for (iq1 = iq; iq1 >= k; --iq1) {
      if (a[iq1] >= t) {
	goto L2;
      }
      x = a[k];
      ix = in[k];
      a[k] = a[iq1];
      in[k] = in[iq1];
      a[iq1] = x;
      in[iq1] = ix;
      iq = iq1 - 1;
      goto L71;
    L2:
      ;
    }
    iq = k - 1;
    goto L3;
  L71:
    ;
  }
 L3:
  a[i] = a[iq];
  in[i] = in[iq];
  a[iq] = t;
  in[iq] = it;
  if (iq + iq <= i + j) {
    goto L4;
  }
  lt[m - 1] = i;
  mt[m - 1] = iq - 1;
  i = iq + 1;
  goto L5;
 L4:
  lt[m - 1] = iq + 1;
  mt[m - 1] = j;
  j = iq - 1;
 L5:
  ++m;
  if(m > 200) cout<<" thesort.cc ==> mt, lt index is out of range "<< m <<endl;
  goto L6;
  
 L7:
  if (i >= j || a[i] <= a[j]) {
    goto L8;
  }
  x = a[i];
  ix = in[i];
  a[i] = a[j];
  in[i] = in[j];
  a[j] = x;
  in[j] = ix;
 L8:
  --m;
  if (m <= 0) {
    return;
  }
  i = lt[m - 1];
  j = mt[m - 1];
  goto L6;
  
  /*   Second algorithm */
  
 L100:
  c = a[1];
  for (k = 2; k <= n; ++k) {
    if (a[k] >= c) {
      goto L11;
    }
    b = a[k];
    ib = in[k];
    k1 = k - 1;
    a[k] = c;
    in[k] = in[k1];
  L9:
    if (k1 == 1) {
      goto L10;
    }
    k2 = k1 - 1;
    if (a[k2] <= b) {
      goto L10;
    }
    a[k1] = a[k2];
    in[k1] = in[k2];
    k1 = k2;
    goto L9;
  L10:
    a[k1] = b;
    in[k1] = ib;
  L11:
    c = a[k];
  }
  //cout<<"End of thesort...."<<endl<<flush;
  return;
}





void lindist(THlx *HTr1, THlx *HTr2, double xmin[], double *dist)
{
  double x1, x2, y1, y2, z1, z2;
  double dyPrdX, dzPrdX, dyPidX, dzPidX;
  double aa, bb, cc;
  double l1,n1,m1,l2,n2,m2,t1,t2,dd,lam1,lam2,xv1,xv2,yv1,yv2,zv1,zv2;
 
    x1 = (*HTr1)(0);
    y1 = (*HTr1)(1);
    z1 = (*HTr1)(2);
    dyPrdX = (*HTr1)(3);
    dzPrdX = (*HTr1)(4);
    t1 = sqrt(dyPrdX*dyPrdX+dzPrdX*dzPrdX+1.);
    if(t1!=0.){
      l1 = 1./t1;
      n1 = dyPrdX/t1;
      m1 = dzPrdX/t1;
    }
    else {
      l1 = 1.;
      n1 = 0.;
      m1 = 0.;
    }
 
    x2 = (*HTr2)(0);
    y2 = (*HTr2)(1);
    z2 = (*HTr2)(2);
    dyPidX = (*HTr2)(3);
    dzPidX = (*HTr2)(4);
    t2 = sqrt(dyPidX*dyPidX+dzPidX*dzPidX+1.);
    if(t2!=0.){
      l2 = 1./t2;
      n2 = dyPidX/t2;
      m2 = dzPidX/t2;
    }
    else{
      l2 = 1.;
      n2 = 0.;
      m2 = 0.;
    }
 
    aa = l2*l1 + n2*n1 + m2*m1;
    bb = l1*(x1-x2) + n1*(y1-y2) + m1*(z1-z2);
    cc = l2*(x1-x2) + n2*(y1-y2) + m2*(z1-z2);
 
    dd = (aa*aa -1.);
 
    lam1 = (bb - aa*cc)/dd;
    lam2 = (bb*aa - cc)/dd;
 
    xv1 = x1 + l1*lam1;
    yv1 = y1 + n1*lam1;
    zv1 = z1 + m1*lam1;
    xv2 = x2 + l2*lam2;
    yv2 = y2 + n2*lam2;
    zv2 = z2 + m2*lam2;
 
    xmin[0] = (xv1+xv2)/2.;
    xmin[1] = (yv1+yv2)/2.;
    xmin[2] = (zv1+zv2)/2.;
 
    *dist=sqrt((xv1-xv2)*(xv1-xv2)+(yv1-yv2)*(yv1-yv2)+(zv1-zv2)*(zv1-zv2));
    
    return;
}



inline void dist2(double xc, THlx *HTr1, THlx *HTr2, 
			   double& fc, double& gc)
{
  double yTr1, yTr2, zTr1, zTr2;
  double dydx1, dydx2, dzdx1, dzdx2;
  double dy12, dz12;
 
  const double two = 2.;
 
  THlx HTrExtr1, HTrExtr2;
 
  HTrExtr1(0) = xc;
  HTrExtr2(0) = xc;
 
  (*HTr1).Extrapolate(HTrExtr1,0);
  (*HTr2).Extrapolate(HTrExtr2,0);
 
  yTr1 = HTrExtr1(1);
  zTr1 = HTrExtr1(2);
  dydx1 = HTrExtr1(3);
  dzdx1 = HTrExtr1(4);
 
  yTr2 = HTrExtr2(1);
  zTr2 = HTrExtr2(2);
  dydx2 = HTrExtr2(3);
  dzdx2 = HTrExtr2(4);
 
  dy12 = yTr1 - yTr2 ;
  dz12 = zTr1 - zTr2 ;
 
  fc = dy12*dy12+dz12*dz12;
  gc = two*( dy12*(dydx1-dydx2) + dz12*(dzdx1-dzdx2) );
 
  //   cout<<" func xc="<<xc << " fc="<<*fc<<" gc="<<*gc<<endl;
 
}


 
void helixdist(THlx *H1, THlx *H2, double a,double b, 
			       double x0, double ver[], double *dist, 
			       THlx *H1v, THlx *H2v)
{
  double f0, g0;
  double xcurr,fcurr,gcurr,xnext,fnext,gnext,xprev,fprev,gprev;
  double x1,x2,x3,y1,y2,y3,aPar,bPar,XvrtParabolic;
  double XVert1,YVert1,ZVert1,XVert2,YVert2,ZVert2;
  double dYdXVert1, dYdXVert2, dZdXVert1, dZdXVert2;
  int NIter, EndFlag;
 
 
  dist2(x0, H1, H2, f0, g0);
 
  //   cout<<"x0="<<x0<<" f0="<<f0<<" g0="<<g0<<endl;
 
  xcurr = x0;
  xnext = x0;
  xprev = x0;
  fcurr = f0;
  fnext = f0;
  fprev = f0;
  gcurr = g0;
  gnext = g0;
  gprev = g0;
  EndFlag = 0;
 
  for (NIter=0 ; (EndFlag==0 && NIter<20) ; NIter++){
    if(gcurr!=0.){
 
 
      if(gcurr*gprev >= 0.) {
        xnext = xcurr - fcurr/gcurr;
        if(xnext>b){xnext = x0 + double(NIter)*5.;}
        if(xnext<a){xnext = x0 - double(NIter)*5.;}
      }
      else{
        if(fabs(xprev-xcurr)<25.) { EndFlag = 1;}
        xnext = (xprev + xcurr)/2.;
      }
 
      dist2(xnext, H1, H2, fnext, gnext);
 
      //      cout<<" "<<endl;
      //      cout<<NIter<<" xprev="<<xprev<<" fprev="<<fprev<<" gprev="<<gprev<<endl;
      //      cout<<NIter<<" xnext="<<xnext<<" fnext="<<fnext<<" gnext="<<gnext<<endl;
      //      cout<<NIter<<" xcurr="<<xcurr<<" fcurr="<<fcurr<<" gcurr="<<gcurr<<endl;
 
      if(EndFlag == 0){
        xprev = xcurr;
        fprev = fcurr;
        gprev = gcurr;
        xcurr = xnext;
        fcurr = fnext;
        gcurr = gnext;
      }
 
    }
    else { EndFlag = 2;}
  }
 
  //  cout<<"Xvrt="<<xnext<<" fnext="<<fnext<<" NIter="<<NIter<<" Eflg="<<EndFlag<<endl;
 
  if(EndFlag==1){
    x1 = xprev; x2 = xnext; x3 = xcurr;
    y1 = fprev; y2 = fnext; y3 = fcurr;
 
    //Det= x1*x1*(x2-x3)-x2*x2*(x1-x3)+x3*x3*(x1-x2);
    //if(Det==0.) {cout<<"Det=0!!!!"<<endl;}
 
    aPar= y1*(x3-x2)+y2*(x1-x3)-y3*(x1-x2);
    bPar= x1*x1*(y3-y2)-x2*x2*(y3-y1)+x3*x3*(y2-y1);
 
    XvrtParabolic = -bPar/(2.*aPar);
  }
  else{
    if (EndFlag==2) {
      XvrtParabolic = xcurr;
    }
    else{
      XvrtParabolic = x0;
    }
  }
 
  //cout<<"Xvrt_parabolic="<<XvrtParabolic<<endl;
 
  (*H1v)(0) = XvrtParabolic ;
  (*H2v)(0) = XvrtParabolic ;
 
  (*H1).Extrapolate(*H1v,0);
  XVert1 = (*H1v)(0);
  YVert1 = (*H1v)(1);
  ZVert1 = (*H1v)(2);
  dYdXVert1 = (*H1v)(3);
  dZdXVert1 = (*H1v)(4);
 
  (*H2).Extrapolate(*H2v,0);
  XVert2 = (*H2v)(0);
  YVert2 = (*H2v)(1);
  ZVert2 = (*H2v)(2);
  dYdXVert2 = (*H2v)(3);
  dZdXVert2 = (*H2v)(4);
 
  ver[0] = (XVert1 + XVert2)/2.;
  ver[1] = (YVert1 + YVert2)/2.;
  ver[2] = (ZVert1 + ZVert2)/2.;
 
  *dist = sqrt((XVert1-XVert2)*(XVert1-XVert2)+
	       (YVert1-YVert2)*(YVert1-YVert2)+
	       (ZVert1-ZVert2)*(ZVert1-ZVert2));
 
  //  cout<<"ParabolicVer x="<<ver[0]<<" y="<<ver[1]<<" z="<<ver[2]<<" dist ="<<*dist<<endl;
}


 

double Q2( double pmu0, double pmu, double Cos )
{
  const double mmu = 0.105658357;
  double pp = pmu0*pmu*Cos;
  double E0 = sqrt( pmu0*pmu0 + mmu*mmu );
  double E  = sqrt( pmu *pmu  + mmu*mmu );
  double q2 = 2 * ( - mmu * mmu - pp + E0 * E );
  return q2;
}

double xbj( double pmu0, double pmu, double Cos )
{
  const double m_P = 0.93827231;
  const double mmu = 0.105658357;
  double q2 = Q2( pmu0, pmu, Cos );
  double E0 = sqrt( pmu0*pmu0 + mmu*mmu );
  double E  = sqrt( pmu *pmu  + mmu*mmu );
  double x = q2 / ( 2 * m_P * ( E0 - E ) );
  return x;
}
