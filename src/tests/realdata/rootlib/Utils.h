// $Id: Utils.h,v 1.13 2002/12/09 08:51:06 hpereira Exp $
#ifndef Utils_h
#define Utils_h

#include <string>
#include <TROOT.h>
#include <TObject.h>

class TMatrix;
class TCanvas;
class TH1;
class TH1S;
class TH1D;
class TH2;
class TH2S;
class TH2D;
class ostream;

class Utils: public TObject {
public:
  static char*  GetTimeStamp( void );
  static string Convert(string sIn, string c1, string c2); 
  static void SetVisibility( TCanvas *tcv );
  static double HDiv(TH2* h1, TH2* h2, TH2D* h3);
  static double HDiv(TH1* h1, TH1* h2, TH1D* h3);
  static double HDiv(TH1* h1, TH1* h2, TH1D* h3, unsigned int i1, unsigned int i2); 
  static void HIntegrate( TH1S* h1, TH1S* h2 );
  static void HDerive( TH1S* h1, TH1S* h2 );
  static void DumpMatrix( TMatrix m, ostream &out = cout );
  
  enum WStatus { NO, YES, EXIT, ALWAYS_W, NEVER_W };
  static WStatus Writable( string fileName );
  inline static void SetWritableStatus( const WStatus wStat ) { wStat_ = wStat; }
 
  private:
  static WStatus wStat_;
  ClassDef(Utils,0)
};

#endif
