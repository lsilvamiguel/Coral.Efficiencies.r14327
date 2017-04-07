/*!
   \file    CsRolandPattern.h
   \brief   Compass Vertex Pattern Class.
   \author  Alexandre Korzenev
   \version $Revision: 1.11 $
   \date    $Date: 2010/04/20 00:13:48 $ 

*/

#ifndef CsRolandPattern_h
#define CsRolandPattern_h



/*! \class CsRolandPattern 
    \brief Coral interface to Vertex Reconstruction Pattern.
*/

#include "CsVrtPattern.h"
#include "CsHistograms.h"
#include "CsBeam.h"

class CsRolandPattern : public CsVrtPattern {

 private:
  bool hist_;
  int NSpec_;
  std::list<CsTrack*> beams_;
  std::list<CsTrack*> tracks_;
  std::list<CsVertex*> vrts_;
  double c_[4];  // CUTS

 private:
  /*! \fn int filtervtx(int nin, int & nout,int* list,float toldis,float **vtk,int &jflag); 
    \brief Function filter used in vtxfit 
  */
  int filtervtx(int nin, int & nout,int* list,float toldis,float **vtk,int &jflag);

  /*! \fn int ntrvtx(int nrtr,float **vtk,float zin,float* vtxout,float &disav); 
    \brief Function ntrvtx calc. the closest  
    app. used in vtxfit 
  */
  int ntrvtx(int nrtr,float **vtk,float zin,float* vtxout,float &disav);

  /*! \fn int trackinfield(float* parin,float zend,float* parout); 
    \brief Function trackinfield used in vtxfit 
  */
  int trackinfield(float* parin,float zend,float* parout);

  /*! \fn int trackbeam(float* parin,float zend,float* parout); 
    \brief Function trackinfield used for the beam in vtxfit. 
  */
  int trackbeam(float* parin,float zend,float* parout);

  /*! \fn int trackinfield2(float* parin1,float*parin2,float& zend,float* parout1,float *parout2,float & dist); 
    \brief Function trackinfield used in vtxfit 
  */
  int trackinfield2(float* parin1,float*parin2,float& zend,float* parout1,float *parout2,float & dist);

  /*! \fn int quickvtx(float* vec1,float* vec2,float &zapp,float &xapp,float &yapp,float &disapp); 
    \brief Function quickvtx calc closest. 
    app of 2 lines, used in vtxfit 
  */
  int quickvtx(float* vec1,float* vec2,float &zapp,float &xapp,float &yapp,float &disapp);

  /*! \fn int vtxfit(int & mlassc,float **par,int & jmu,float* prbeam,float* tol,float* vtxbest,float* vtx2,int & jret); 
    \brief Function vtxfit :calls filter,trackin..,ntrvtx,quick 
  */
  int vtxfit(int & mlassc,float **par,int & jmu,float* prbeam,float* tol,float* vtxbest,float* vtx2,int & jret);

  /*! \fn CsVertex* getPrimaryVertex(float* vtx2,int *nflag,int& vertflag,int & bad, 
                                     int & less,int & muon,int & muon2); 
    \brief <b><i>Rolands</i></b> vertex finding routine 
  */
  CsVertex* getPrimaryVertex(float* vtx2,int *nflag,int& vertflag,int & bad,
                             int & less,int & muon,int & muon2);
  
 public:

  CsRolandPattern();             //!< constructor
  virtual ~CsRolandPattern() {}; //!< destructor

  /*! \fn bool doPattern( vector<CsParticle*> &particles )
    \brief This method performs the patter recognition on a given
    collection of particles. Returns \c true if the operation ended correctly;
    returns \c false otherwise.
    \param particles Vector of pointers to the particles to be used in the 
    pattern recognition procedure.
  */
  bool doPattern( std::vector<CsParticle*> &particles, double *reTrackT0);

  /*! \fn bool getPatterns( list<CsVertex*>& vrts, map<CsTrack*,bool>& specials )
    \brief This method returns the list of found track candidates after
    the prepattern procedure. These tracks could be a not completely 
    fitted ones.
    \param The list of track where to add the found patterns
  */
  bool getPatterns( std::list<CsVertex*>& vrts, std::map<CsTrack*,bool>& specials );
  
  /*! \fn const CsVertex *getT0SettingVertex() const
    \brief Dummy method.
   */
  const CsVertex *getT0SettingVertex() const { return 0; }
};

#endif //CsRolandPattern_h
