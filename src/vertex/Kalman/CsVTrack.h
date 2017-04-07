/*!
   \file    CsVTrack.h
   \brief   Reconstructed track at vertex class.
   \version $Revision: 1.6 $
   \author  C.Ulvegren, A.Korzenev
   \date    $Date: 2007/11/11 00:21:26 $
*/

#ifndef CsVTrack_h
#define CsVTrack_h

#include <math.h>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "coral_config.h"
#include "CsTrack.h"
#include "TMtx.h"

class CsTrack;

/*! \class CsVTrack
 *  \brief Compass Track at Vertex Class.
 *  
 *  Matrices for Kalman filter notation as in R.Fruhwirth et al CPC 96(1996) pp 190-192
 */
class CsVTrack {
  
  private :
    
    CsTrack* TrkRef_;  //!< reference to associated reconstructed track
  /*! \fn int status_;
    \brief Information about track status in vertex finding process.

    0-in prim. vert; 
    1-rejected by dist. cut; 
    2-rejected by IKF.; 
    3-rejected as track which goes not from prim. vrt;
    4-rejected as track which has not MC associated trk;
    5-apriori in fit (for track beam and for scattered muon).
  */
    int status_;      
    double Chi2Tr_;    //!< Chi2 of track in inverse Kalman filter
    double RadLength_; //!< The thickness of the scattering medium in radiation length.
    double ELoss_;     //!< The loss of energy in target.
    double Res_[5];    //!< Reduced smoothed residuals ("pulls"). It contains pulls values for \a x, \a y, \a dx/dz, \a dy/dz, \a q/p.
    
  public :

    // Vk - covariance matrices of track parameters Pk
  
    TMtx Pk;  //!< Vector (5) of track parameters at reference point
    TMtx Gk;  //!< Weight matrix (5x5) at reference point. \b G <sub>k</sub>= \b V <sub>k</sub><sup>-1</sup>.

    TMtx Ak;  //!< Matrix (5x3) of derivatives at expansion point <b>A</b><sub>k</sub>=[d<b>h</b><sub>k</sub>/d<b>x</b>]<sub>e</sub>.
    TMtx AkT; //!< Matrix(3x5). \b A <sub>k</sub> transposed.
    TMtx Bk;  //!< Matrix (5x3) of derivatives at expansion point <b>B</b><sub>k</sub>=[d<b>h</b><sub>k</sub>/d<b>q</b>]<sub>e</sub>.
    TMtx BkT; //!< Matrix(3x5). <b>B</b><sub>k</sub> transposed.
    TMtx Cke; //!< Vector (5)
    TMtx GkB; //!< Matrix (5x5)
    TMtx Wk;  //!< Matrix (3x3)

    TMtx Qnk; //!< Vector (3) of smoothed track parameters at vertex
    TMtx Pnk; //!< Vector (5) of smoothed track parameters at ref. plane 

    TMtx Dnk; //!< Matrix (3x3) of covariance of smoothed track parameters at vertex
    TMtx Enk; //!< Matrix (3x3) of covariance of vertex position and smoothed track parameters at vertex
    
    CsVTrack();                 //!< Constructor
    CsVTrack(const CsVTrack& ); //!< Constructor
    
    // Destructor
    ~CsVTrack(){};

    // Operators
    CsVTrack& operator = (const CsVTrack& );
 
    // Methods

    /*! \fn bool primary()const
      \brief Returns true if track is from primary vertex.
    */
    bool primary()const {if(status_ == 0 || status_ == 5 )return true;else return false;};

    double getELoss(){ return ELoss_; };

    double getRadLength(){ return RadLength_; };

    /*! \fn int getStatus()const
      \brief Returns track status variable \a status_ .
    */
    int getStatus()const {return status_;}

    /*! \fn void setChi2(double Chi2Tr);
      \brief Sets chi-square of track.
    */
    void setChi2(double Chi2Tr){ Chi2Tr_=Chi2Tr; };

    void setELoss(double el){ ELoss_=el; };

    void setRadLength(double rl){ RadLength_=rl; };

    /*! \fn void setRes(TMtx pulls)
      \brief Sets reduced smoothed residuals ("pulls"). 
      \param pulls Vector(5). It contains pulls values for \a x, \a y, 
      \a dx/dz, \a dy/dz, \a q/p.
    */
    void setRes(TMtx pulls) {Res_[0]=pulls(1);Res_[1]=pulls(2);
                             Res_[2]=pulls(3);Res_[3]=pulls(4);
                             Res_[4]=pulls(5);}

    /*! \fn void setTrkRef(CsTrack* TrkRef)
      \brief Sets reference to associated reconstructed track.
    */
    void setTrkRef(CsTrack* TrkRef){ TrkRef_=TrkRef; };
    
    /*! \fn void SetStatus(int st);
      \brief Sets status of track.
    */
    void setStatus(int st){status_=st;};
    
    /*! \fn CsTrack* getAssociatedTrk()const 
      \brief Return reference to associated reconstructed track
    */
    CsTrack* getAssociatedTrk() const {return TrkRef_;}
    
    /*! \fn double getChi2() const
      \brief Return Chi-square which track has in inverse Kalman filter.
    */
    double getChi2() const {return Chi2Tr_;}; 

     
    /*! \fn double getDXDZ()
      \brief Return value for dx/dz
    */
    double getDXDZ() {return Qnk(1);}
    
    /*! \fn double getDYDZ() 
      \brief Return value for dy/dz
    */
    double getDYDZ() {return Qnk(2);}
    
    /*! \fn double getCop() 
      \brief Return value for Cop
    */
    double getCop() {return Qnk(3);}


    /*! \fn double getMom() 
      \brief Return track's |momentum|.
    */
    double getMom(); 
    
    /*! \fn Hep3Vector ec()
      \brief Return Vector @ Vertex.
    */
    CLHEP::Hep3Vector Vec() const;

    /*! \fn HepLorentzVector LzVec(const double mass)
      \brief Return Lorentz Vector @ Vertex.
    */
    CLHEP::HepLorentzVector LzVec (const double mass) const;

    /*! \fn double getResX()const
      \brief Return pull value for x
    */
    double getResX()const {return Res_[0];}
    
    /*! \fn double getResY()const
      \brief Return pull value for y
    */
    double getResY()const {return Res_[1];}
    
    /*! \fn double getResDXDZ()const
      \brief Return pull value for dx/dz
    */
    double getResDXDZ()const {return Res_[2];}
    
    /*! \fn double getResDYDZ()const 
      \brief Return pull value for dy/dz
    */
    double getResDYDZ()const {return Res_[3];}
    
    /*! \fn double getResCop()const 
      \brief Return pull value for Cop
    */
    double getResCop()const {return Res_[4];}
    
};

#endif














