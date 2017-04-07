#ifndef __CS_DST_04_H__
#define __CS_DST_04_H__

#include "CsDstObject.h"

class CsDstCluster_04 : public CsDstCluster
{
  float32 u;
  float32 v;
  float32 w;
  float32 cov[6];
  float32 time;
  float32 analog;
  uint8   hasTime;
  uint8   hasAnalog;
  int16   detectorPosition;
 public:
           CsDstCluster_04();
  virtual ~CsDstCluster_04() {};

  virtual void        Init(const CsCluster& cl);
  virtual CsCluster*  GetCluster() const;
  virtual bool        Load(std::istream& );
  virtual bool        Load(FlatFile& );
  virtual bool        IsEquil(const CsDstObject*) const;
  virtual void        Print(std::ostream& ) const;
  virtual void        GetBuffer(std::vector<uint8>& buffer) const;
};

class CsDstHelix_04 : public CsDstHelix
{
  float32 x;
  float32 y;
  float32 z;
  float32 dxdz;
  float32 dydz;
  float32 cop;
  float32 cov[15];
 public:
           CsDstHelix_04();
  virtual ~CsDstHelix_04() {};

  virtual void     Init(const CsHelix& hl);
  virtual CsHelix* GetHelix() const;
  virtual bool     Load(std::istream& );
  virtual bool     Load(FlatFile& ); 
  virtual bool     IsEquil(const CsDstObject*) const;
  virtual void     Print(std::ostream& ) const;
  virtual void     GetBuffer(std::vector<uint8>& buffer) const;
};

class CsDstTrack_04 : public CsDstTrack
{
  float32 chi2;
  float32 time;
  float32 timeError;
  float32 XX0;
  float32 rich1Probs[CSTRACK_RICHDATASIZE];
  uint16  flags;
  uint32  expDets[CSTRACK_MAPSIZE];
  uint32  firDets[CSTRACK_MAPSIZE];
  
 public:
           CsDstTrack_04();
  virtual ~CsDstTrack_04() {};

  virtual void     Init(const CsTrack& tr);
  virtual CsTrack* GetTrack() const; // Returns CsTrack w/o helix, clusters etc associations
  virtual bool     Load(std::istream& );
  virtual bool     Load(FlatFile& );
  virtual bool     IsEquil(const CsDstObject*) const;
  virtual void     Print(std::ostream& ) const;
  virtual void     GetBuffer(std::vector<uint8>& buffer) const;
  virtual void     SetAsBeam() { flags |= 1; }
  virtual int      GetNClusters() const { return 0; }
};

class CsDstVertex_04 : public CsDstVertex
{
  uint8   isPrimary;
  float32 x;
  float32 y;
  float32 z;
  float32 chi2;
 public:
           CsDstVertex_04();
  virtual ~CsDstVertex_04() {};

  virtual void      Init(const CsVertex& vt);
  virtual CsVertex* GetVertex() const; // Returns CsVertex w/o tracks etc associations
  virtual bool      Load(std::istream& );
  virtual bool      Load(FlatFile& );
  virtual bool      IsEquil(const CsDstObject*) const;
  virtual void      Print(std::ostream& ) const;
  virtual void      GetBuffer(std::vector<uint8>& buffer) const;
};

class CsDstCalobject_04 : public CsDstCalobject
{
  float32 e;
  float32 x;
  float32 y;
  float32 z;
  float32 eErr;
  float32 xErr;
  float32 yErr;
  float32 zErr;
 public:
           CsDstCalobject_04();
  virtual ~CsDstCalobject_04() {};

  virtual void                       Init(const Reco::CalorimeterParticle& vt);
  virtual Reco::CalorimeterParticle* GetCalobject() const; 
  virtual bool                       Load(std::istream& );
  virtual bool                       Load(FlatFile& );
  virtual bool                       IsEquil(const CsDstObject*) const;
  virtual void                       Print(std::ostream& ) const;
  virtual void                       GetBuffer(std::vector<uint8>& buffer) const;
};

class CsDstParticle_04 : public CsDstParticle
{
  int16 partType;
  int32 PDGid;
 public:
           CsDstParticle_04();
  virtual ~CsDstParticle_04() {};

  virtual void        Init(const CsParticle& pt);
  virtual CsParticle* GetParticle() const;
  virtual bool        Load(std::istream& );
  virtual bool        Load(FlatFile& );
  virtual bool        IsEquil(const CsDstObject*) const;
  virtual void        Print(std::ostream& ) const;
  virtual void        GetBuffer(std::vector<uint8>& buffer) const;
  virtual int16       GetPartType() const { return partType; }
  virtual std::string GetPartName() const;
};

class CsDstEventHeader_04 : public CsDstEventHeader
{
  float32 extVars[3];

 public:
           CsDstEventHeader_04();
  virtual ~CsDstEventHeader_04() {};

  virtual void   Init(const CsRecoEvent&);
  virtual void   GetHeader(CsRecoEvent&)  const;
  virtual bool   Load(std::istream& );
  virtual bool   Load(FlatFile& );
  virtual bool   IsEquil(const CsDstObject*) const;
  virtual void   Print(std::ostream& ) const;
  virtual void   GetBuffer(std::vector<uint8>& buffer) const;

  virtual float32 GetExternalVariable(int n) const { return extVars[n]; }
  virtual void    SetExternalVariable(int n, float32 val) { extVars[n] = val; }
};

struct VertexToTrackAssociations
{
  int         vertex;
  std::vector<int> tracks;
  std::vector<int> trackParams;
  std::vector<int> cov;
};

struct TrackToClusterAssociation
{
  int         track;
  std::vector<int> clusters;
};

struct ParticleToTrackCaloAssociation
{
  int         particle;
  int         track;
  std::vector<int> calos;
};

class CsDstEvent_04 : public CsDstEvent
{
  std::vector<VertexToTrackAssociations>      vta;
  std::vector<TrackToClusterAssociation>      tca;
  std::vector<ParticleToTrackCaloAssociation> ptca;

  int32  nVtxToTrk;
  int32  nTrkToClu;
  int32  nPartTkCo;
  int32  nClusters;
  int32  nHelices;
  int32  nTracks;
  int32  nVertices;
  int32  nCalobjects;
  int32  nParticles;
  int32  nTrkParAtVtx;
  int32  nVtxErrMatrices;
  int32  nReserved;

  int32   scalers[8];
  int32   nHodoData;
  float32 triggerTime;

  uint8   fatness;

  std::vector<int16>         vtxToTrk;
  std::vector<int16>         trkToClu;
  std::vector<int16>         partTkCo;
  std::vector<uint32>        hodoData;
  std::vector<float32>       reserved;

  std::vector<CsDstCluster*>   clusters;
  std::vector<CsDstHelix*>     helices;
  std::vector<CsDstTrack*>     tracks;
  std::vector<CsDstVertex*>    vertices;
  std::vector<CsDstCalobject*> calobjects;
  std::vector<CsDstParticle*>  particles;
  std::vector<CsDst3Vector*>   trkParAtVtx;
  std::vector<CsDst3Matrix*>   vtxErrMatrices;

  std::vector<CsDstObject*>    objects;
 public:
           CsDstEvent_04(int fat);
  virtual ~CsDstEvent_04();

  virtual void    Init(const CsRecoEvent& );
  virtual void    GetEvent(CsRecoEvent&) const;
  virtual bool    Load(std::istream& );
  virtual bool    Load(FlatFile& );
  virtual bool    IsEquil(const CsDstObject*) const;
  virtual void    Print(std::ostream& ) const;
  virtual void    GetBuffer(std::vector<uint8>& buffer) const;
  virtual float32 GetExternalVariable(int n) const { return ((n<3)?extVars[n]:0); }
  virtual void    SetExternalVariable(int n, float32 val) { (n<3)?(extVars[n]=val):0;};
 protected:
  virtual bool    SetAssociations();
};

class CsDstChunk_04 : public CsDstChunk
{
 public:
  CsDstChunk_04();
  virtual ~CsDstChunk_04() {};

  virtual void        Init(const char* chunkName );
  virtual CsDstEvent* GetNextEvent(std::istream& ) const;
  virtual CsDstEvent* GetNextEvent(FlatFile& ) const;
  virtual CsDstEvent* GetNextEvent(std::vector<uint8>&) const;
  virtual bool        Load(std::istream& );
  virtual bool        Load(FlatFile& );
  virtual bool        IsEquil(const CsDstObject*) const;
  virtual void        Print(std::ostream& ) const;
  virtual void        GetBuffer(std::vector<uint8>& buffer) const;
  virtual void        SetExtVarName(int n, const char* name);
  virtual const char* GetExtVarName(int n) const 
    { 
      return ((n>=0) && (n<3))?extVarNames[n].c_str():NULL; 
    }
};


#endif
