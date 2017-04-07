
#ifndef CompassSoft_PixelMMDecoding__include
#define CompassSoft_PixelMMDecoding__include

#include<vector>
#include<map>
#include<complex>
#include "DaqError.h"

using namespace std;

// template <typename T> class compsort : public complex<T> {
//   bool operator< (const complex<T>& rhs) const;
// };

// template <typename T> class cmp_complex {
//  public:
//   bool operator() (const complex<T>& lhs, const complex<T>& rhs) {
//     if (lhs.real()==rhs.real()) return lhs.imag()>rhs.imag();
//     return lhs.real()>rhs.real();
// //     if (lhs.imag()==rhs.imag()) return lhs.real()<rhs.real();
// //     return lhs.imag()<rhs.imag();
//   }
// };


// template<typename T> inline bool operator< (const complex<T>& lhs, const complex<T>& rhs) {
//   if (lhs.real()==rhs.real()) return lhs.imag()<rhs.imag();
//   return lhs.real()<rhs.real();
// }
// 
// template<typename T> inline bool operator() (const complex<T>& lhs, const complex<T>& rhs) {
//   return complex<T>(lhs.real()<rhs.real();
// }

// template<typename T> inline bool compsort<T>::operator< (const complex<T>& rhs) const {
//   if (this->real()==rhs.real()) return this->imag()<rhs.imag();
//   return this->real()<rhs.real();
// }


// Private class for PixelMM detector channel to X,Y mapping, square version
class pixmm {

 private:
  static const int nbChannel=128;
  int nbConn;
  int NVSConn[20], ConnVSN[20];
  double pitchU, pitchV, firstU, firstV;
  int nWirU, nWirV;
//   int pix[nbConn*nbChannel];
//   float X[nbConn*nbChannel], Y[nbConn*nbChannel];
  vector<int> pix;
  vector<float> X, Y;
//   map< complex<int>, int, cmp_complex<int> > pix_address;
  map< pair<int,int>, int > pix_address;
  int version;

  void initconndata_v1(void);
  void initconndata_v2(void);
  void initconndata_v3(void);

 public:
  bool dataloaded;

  pixmm(int _version = 1);
  void initconndata(void);

  int GetNbPix()  {return nbConn*nbChannel;}
  int GetNbConn()  {return nbConn;}
  int GetIdxConn(int indice)  { if (indice<0 || indice>=20) return -1; return ConnVSN[indice]; }
  int GetPixNb(int indice)  { if (indice<0 || indice>=GetNbPix()) return -1; return pix[indice]; }
  int GetPixNb(int conn, int chan)  { 
                             if (conn<0 || conn>=20 || chan<0 || chan>=nbChannel) return -1;
                             if (ConnVSN[conn]<0) return -1;
                             return pix[ConnVSN[conn]*nbChannel+chan];
                           }
  // X, Y positions in mm, center of detector at 0.,0.
  float GetXPix(int pixnb) { if (pixnb>=0 && pixnb<(nbChannel*nbConn)) return X[pixnb];
                                    throw CS::Exception("pixmm::GetXPix(): bad pixel number %d", pixnb); }
  float GetYPix(int pixnb) { if (pixnb>=0 && pixnb<(nbChannel*nbConn)) return Y[pixnb];
                                    throw CS::Exception("pixmm::GetYPix(): bad pixel number %d", pixnb); } 
  bool GetLargePix(int pixnb) { if (version==2) return pixnb>=384; return pixnb>=640; }
  int GetNConn(int conn)  { if (conn>=0 && conn<20) return ConnVSN[conn]; return 0; }
//   map< complex<int>, int, cmp_complex<int> > PixelAddresses() { return pix_address; }
  map< pair<int,int>, int > PixelAddresses() { return pix_address; }
  int PixelAddr(int i, int j) const  { return pix_address.at(make_pair(i,j)); }
  int GetNbPixU() const   { return nWirU; }
  int GetNbPixV() const   { return nWirV; }
  double GetPitchU() const   { return pitchU; }
  double GetPitchV() const   { return pitchV; }
  double GetFirstWireU() const   { return firstU; }
  double GetFirstWireV() const   { return firstV; }

 private:
  static pixmm* pixmmv1;
  static pixmm* pixmmv2;
  static pixmm* pixmmv3;

 public:
  static pixmm* Instance(int version);


}; // class pixmm


#endif // CompassSoft_PixelMMDecoding__include
