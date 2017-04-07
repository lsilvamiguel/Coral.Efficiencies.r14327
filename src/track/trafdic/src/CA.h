#if !defined CA_h
#define CA_h

#include <list>

typedef float ftype;

ftype const PI=3.1415925636;

extern "C"{
  void ranlux_(float *,int *);
}

extern void myexit(const char *c);

inline ftype sqr(ftype x){ return x*x; }

extern void fill_zero(ftype G[6][6]);
extern void mirror(ftype G[6][6]);
extern void copy_matrix(ftype G[6][6], ftype G1[6][6]);
extern void copy_vector(ftype R[6], ftype R1[6]);
extern void count_rot(const ftype G[6][6], const ftype R[6], ftype G1[6][6], ftype R1[6], ftype C, ftype S);
extern void count_shift(const ftype G[6][6], const ftype R[6], ftype G1[6][6], ftype R1[6], ftype dx);
extern void count_sp(ftype x1, ftype u1, ftype C1, ftype S1,
		     ftype x2, ftype u2, ftype C2, ftype S2,
		     ftype R[6], ftype x0, int hallo);

extern void count_segment(ftype x1, ftype u1, ftype C1, ftype S1, ftype sigma1, // == sigma1^2
			  ftype x2, ftype u2, ftype C2, ftype S2, ftype sigma2, // == sigma2^2
			  ftype x3, ftype u3, ftype C3, ftype S3, ftype sigma3, // == sigma3^2
			  ftype x4, ftype u4, ftype C4, ftype S4, ftype sigma4, // == sigma4^2
			  ftype G[6][6], ftype R[6], ftype x0, 
			  ftype G1[6][6], ftype R1[6], ftype x01, 
			  ftype sigma_ms);

void IPLm(int, float*, float*); //substitution for HIGZ IPL (higz_mods.cc)
void IPMm(float, float);        //substitution for HIGZ IPM (higz_mods.cc)

static int 
  White      =   0,
  Black      =   1, 
  Red        =   2, 
  Green      =   3,  
  Blue       =   4, 
  Yellow     =   5,
  Magenta    =   6, 
  LightBlue  =   7,
  Brown      = 111; 

// ==  define classes ==
  
class TLayer;

class TCAHit{
 public:
  int origin_index; // index of corresponding THit object in event.vHit[]
  TLayer *Layer;
  ftype x, u, v;
  ftype sigma_u, sigma_v;
  int used;	    // is hit used by reconstructed track
  int num4;       // number of long segments which are used this hit
  TCAHit(){
    Layer = NULL;
  }
  std::vector<int> mirror;
};

class TLayer{
 public:
  int index;
  int global_plane_index;
  int projection;         // Projection number in TRAFFIC program
  int type;               // detector type in TRAFFIC
  ftype resol;            // ?
  ftype pitch;            // pitch between wires
  ftype range_u;          // length of active area in measurement direction   cm
  ftype range_v;          // length of active area perpendicular to measurement direction   cm
  ftype holl;             // radius of non-active area cm
  ftype holl2;            // holl^2
  ftype size_x;           // size in beam direction 
  ftype Sa, Ca, X0;       // destination in space
  TCAHit Hits[200];
  int NumberOfHits;
  ftype SEGMENT_CUT_TARGET; // cut angle for segments with target center direction tg(angle)
  ftype SEGMENT_CUT_HALLO;  // cut angle for segments with beam direction tg(angle)

  TLayer(){
    NumberOfHits = 0;
  }

  int active_zone(ftype y, ftype z, ftype cut=0){
    ftype u =  y*Ca+z*Sa;
    ftype v = -y*Sa+z*Ca;
    return (fabs(u)-cut<=range_u && fabs(v)-cut<=range_v && sqr(u)+sqr(v)+cut>=sqr(holl));
  }
  int active_zone(ftype x, ftype R[6]){
    ftype dx = X0-x;
    return active_zone(R[1] + dx*R[3], R[2] + dx*R[4]);
  }
  TCAHit *find_nearest_hit(ftype u);
};

class TSegment;

class TSuperLayer{
 public:
  int NumberOfLayers;
  TLayer *Layers[10];
  int CrossedLayers[20];
  int NumberOfCrossedLayers;
  int RefSegmentNumberOfHits;
  int ExtraSegmentNumberOfHits;
  int f_seg, l_seg;
  TSuperLayer(){
    NumberOfLayers = NumberOfCrossedLayers = 0;
  }
};

class TCATrack{
 public:
  std::list<TCAHit*> Hits;
  std::list<ftype> crosspoints;
  ftype vx, vR[6], vG[6][6];
  ftype Vx, VR[6], VG[6][6];
  ftype chi2, pinv;
  TCATrack(){}

  int NumberOfHits(){
    return Hits.size();
  }
  void add_hit(TCAHit *hit);
  void CountPQ();
  void Fit();
  void Fit(ftype x1, ftype x2);
  void Propagate(ftype x, ftype R[6], ftype G[6][6]);
  void FitQK();
  void FitQK(ftype x1, ftype x2);
  void PropagateQK(ftype x, ftype R[6], ftype G[6][6]);
  void Copy(const TCATrack& t);
};

extern void Fill_Layers(TLayer Layers[], int &NumberOfLayers);
extern void Fill_Super_Layers(TSuperLayer SuperLayers[], int &NumberOfSuperLayers, 
			      TLayer Layers[], int NumberOfLayers);
#endif
