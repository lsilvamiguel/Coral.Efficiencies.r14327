
#include <iostream>

int setN1N2(int N, int&N1, int&N2) {
  int sr=int(sqrt(N));
  if      ( sr   * sr   ==N){N1=sr  ;N2=sr  ;}
  else if ((sr-1)*(sr+1)>=N){N1=sr-1;N2=sr+1;}
  else if ( sr   *(sr+1)>=N){N1=sr;  N2=sr+1;}
  else                      {N1=sr+1;N2=sr+1;}
  return N1*N2;
}

double eS(double x, double *p) {
  double xx=x-p[1];
  return 0.5*( (p[2]+p[3]) *       xx  +
	       (p[2]-p[3]) * (sqrt(xx*xx+p[4]*p[4])-fabs(p[4])) ) + p[5];
}

double eSprime(double x, double *p) {
  double xx=x-p[1];
  return 0.5*( (p[2]+p[3]) -
	       (p[2]-p[3]) * xx / sqrt(xx*xx+p[4]*p[4]) );
}

double eSi(double s, double *p) {
  double apc=p[2]+p[3];
  double amc=p[2]-p[3];
  double f=s-p[5]+0.5*p[4]*amc;
  double ac=p[2]*p[3];
  double x=0.5/ac*(apc*f-amc*sqrt(f*f+p[4]*p[4]*ac));
  return x+p[1];
}

double APV_r(double x, double *p) {
  return p[0]*exp(-exp(-eS(x, p)));
}

double APV_r_(double *x, double *p) {
  return APV_r(x[0], p);
}

double APV_t(double x, double *p) {
  return eSi(-log(-log(x/p[0])), p);
}

double APV_t_(double *x, double *p) {
  return APV_t(x[0], p);
}

double APV_SI_ampl(double t, int ratio,
		       double *pp,
		       double delta, int N) {
  t+=delta;
  double ampl = 1;
  double tau2 = delta/log(pp[0]);
  for (int i=0; i<N; i++) {
    ampl *= (ratio==0) ? APV_r(t+double(i)*delta, pp)
                       : APV_r(t+double(i)*delta, pp);
  }
  return ampl * exp(-(t+double(N)*delta)/tau2);
}

double APV_SI_ampl_f(double *t, double *pp) {
  return pp[7]+pp[6]*APV_SI_ampl(t[0], 0, pp, 1.0, 8);
}




