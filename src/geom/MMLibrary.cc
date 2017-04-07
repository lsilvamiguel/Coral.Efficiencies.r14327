#include "MMLibrary.h"
#include <map>
#include <math.h>
#include <stdio.h>

using namespace std;
using namespace MM;


// Constants -------------------------------------------------------
// const float Constant::sTotnorm = 100; 
const int   Constant::sRetainMaxHits = 7; 
const int   Constant::sRetainMaxHitsTime = 2; 

// Digit -----------------------------------------------------------

int Digit::avstotmode = STANDARD;

Digit::Digit(int ch, double t, double tot, bool largePitch) 
  : ch(ch), t(t), tot(tot) {

  switch(Digit::avstotmode) {
  case REBOURGEARD:
    // a is in 1000*e-, tot in ns. 
    a = 6.241*exp( (tot/134)*(tot/134) );
    break;
  case REBOURGEARD2:
    // a is in 1000*e-, tot in ns. 
    a = 6.241*exp( tot/134. );
    break;
  case STANDARD:
    if (largePitch) a = exp(tot/100);
    else            a = exp(tot/100);
    //a = exp(tot/100.);
    break;
  case NEW:
    a = 6.241*exp(tot/150.);
    break;
  default:
    cerr<<"MM::Digit::Digit(ch,t,tot), Aargh ! this should never happen !"<<endl;
    a = 0;
    break;
  }
}

void Digit::SetAvsTot(unsigned mode) {
  switch(mode) {
  case REBOURGEARD:
  case REBOURGEARD2:
  case STANDARD:
  case NEW:
    Digit::avstotmode = mode;
    break;
  default:
    cerr<<"MM::Digit::SetAvsTot, undefined mode, using STANDARD."<<endl;
    Digit::avstotmode = STANDARD;
  }  
}


ostream& MM::operator<<(ostream& out, Digit& d) {
  if(!out) return out;
  out<<"digit\t"<<d.ch<<"\t"<<d.t<<"\t"<<d.tot<<"\t"<<d.a;
  return out;
}


// Cluster --------------------------------------------------------

Cluster::Cluster(double wire,double t,double tot,int s) 
  : wire(wire),t(t), tot(tot), a(0), s(s) {}


Cluster::Cluster(const Cluster& c)  {
  wire = c.wire;
  t = c.t;
  tot = c.tot;
  s = c.s;
  a = c.a;

  fDigits = c.fDigits;
}


Cluster::Cluster(const vector<Digit>& v)  {
  fDigits = v;
  Update();
}


Cluster::~Cluster() {
  fDigits.clear();
}


void Cluster::Dump() {
  cout.setf(ios::fixed);
  cout<<"cluster\t"<<"\t"<<wire<<"\t"<<t<<"\t"<<tot<<"\t"<<a<<"\t"<<s<<endl;
  for(unsigned i=0; i<fDigits.size(); i++) {
    cout<<"\t";
    MM::operator<<(cout,fDigits[i]);
    cout<<endl;
  }
  cout.unsetf(ios::fixed);
}
  
 
vector<Cluster*> Cluster::Split() {


  vector<Cluster*> products;
  const int splitminsize=5;
  if( s < splitminsize) {
    products.push_back(this);
    return products;
  }

  vector<Digit> v;
  v.push_back(fDigits[0]);
  for(unsigned i=1; i<fDigits.size()-1; i++) {
    double d1 = fDigits[i].a - fDigits[i-1].a;
    double d2 = fDigits[i+1].a - fDigits[i].a;
    
    v.push_back(fDigits[i]);
    
    if(d1<0 && d2>0) {
      v.back().a /= 2;
      Digit digit = v.back();
      products.push_back( new Cluster(v) );
      v.clear();
      v.push_back(digit);
    }
  }
  v.push_back(fDigits[fDigits.size()-1]);
  this->SetDigits(v);
  this->Update();
  products.push_back( this );
  
  return products;
}


void Cluster::Update() {

  // cluster size is by def the number of digits in fDigits 
  s = fDigits.size();
       
  map<double, Digit*, greater<double> > digits; 
  for(unsigned i=0; i<fDigits.size(); i++) {
    digits[ fDigits[i].a ] = & fDigits[i];
  }

  double numwire = 0, denomwire = 0;
  int wiremin = 0xffff, wiremax = -1;
  double numtime = 0, denomtime = 0;

  int n=0;
  a = 0;
  for(map<double, Digit*, greater<double> >::iterator id=digits.begin();
      id != digits.end(); id++) {
    
    double weight = id->second->a;
    a += weight; // recalculating amplitude 

    if(n < Constant::sRetainMaxHits) {
      wiremin = min(wiremin,id->second->ch);
      wiremax = max(wiremax,id->second->ch);
      
      numwire += weight*id->second->ch;
      denomwire += weight;
    }
    
    if(n < Constant::sRetainMaxHitsTime) {
      numtime += weight*id->second->t;  
      denomtime += weight;	  
    }
    n++;
  }

  wire = numwire / denomwire;
  t = numtime/denomtime;
  tot = digits.begin()->second->tot;
}


