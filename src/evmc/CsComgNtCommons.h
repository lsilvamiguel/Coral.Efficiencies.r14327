// $Id: CsComgNtCommons.h,v 1.3 1999/08/25 08:52:47 benigno Exp $
 
/*!
   \file    CsComgNtCommons.h
   \brief   Comgeant Ntuples Commons.
   \author  Benigno Gobbo
   \version $Revision: 1.3 $
   \date    $Date: 1999/08/25 08:52:47 $
*/

#ifndef CsComgNtCommons_h
#define CsComgNtCommons_h

//! \struct Qbea
//! \brief the beam ntuple
typedef struct {      
  int ibtyp;         //!< type of beam particle
  int ibfla;         //!< 0: beam starting point, 1: interaction point
  float bpara[6];    //!< Vertex(X,Y,Z), Slope(X,Y), Momentum
	      } QbeaType;

//! \struct Qhea 
//! \brief the ntuple head
typedef struct { 
  int ieve;  //!< event number
  int irun;  //!< run number
  int iend;  //!< IEOTRI (see Comgeant)
	      } QheaType;

//! \struct Qhit
//! \brief the hit and digit ntuple
const int mxhitnt  = 10000; //!< max number of hits
const int mxdignt  = 12000; //!< max number of digitising
const int mxhdignt = 16000; //!< max number of digitising->hit links
const int mxddig   = 2;
typedef struct { 
  int lmxhitnt;             //!< actual max number of hits
  int lmxdignt;             //!< actual max number of digitising
  int lmxhdignt;            //!< actual max number of digitising->hit links
  int nhit;                 //!< number of hits in detectors
  int nhitall;              //!< full number of hits
  int ip1hit[mxhitnt];      //!< packed ithitnt, jthitnt (2 B each)
  int ip2hit[mxhitnt];      //!< packed idhitnt, itmhitnt+32768 (2 B each)
  float hit[mxhitnt][2];    //!< transversal hit coordinates in MRS
  float xhit[mxhitnt][3];   //!< coordinates in MRS
  int ithitnt[mxhitnt];     //!< track number
  int jthitnt[mxhitnt];     //!< 0: original trk, N: secondary track
  int idhitnt[mxhitnt];     //!< detector number
  int itmhitnt[mxhitnt];    //!< the time of the hit ( ns * 10 )
  int ndig;                 //!< number of digitisings
  int ndigall;              //!< full number of digits
  int ip1dig[mxdignt];      //!< packed iddignt, jpdig ( 2 B each )
  int ip2dig[mxdignt];      //!< packed jddignt[0,], jddignt[1,] ( 2 B each ) 
  int iddignt[mxdignt];     //!< detector number
  int jddignt[mxddig][mxdignt]; //!< digitisation valur ( mxddig=2 at present)
  int nhdig;                //!< arrays for hits
  int ntdig;                //!< arrays for tracks
  int npdig;                
  int ihdig[mxhdignt];
  int itdig[mxhdignt];
  int ipdig[mxhdignt];
  int jhdig[mxhdignt];
  int jtdig[mxhdignt];
  int jpdig[mxhdignt];
  int nhdiglnt;
  int nhitant;
  int ndigant;
	      } QhitType;

//! \struct Qkin; 
//! \brief the kinematics ntuple
const int mxtraknt = 250;   //!< max number of tracks
const int mxverknt = 120;   //!< max number of vertices
typedef struct {   
  int nver;                 //!< number of vertices         
  int igev[mxverknt];       //!< Geant vertex number
  float vert[mxverknt][3];  //!< Vertex Position
  int ltimv[mxverknt];      //!< Vertex time ( ns * 10 )
  int imov[mxverknt];       //!< Mother track
  int ntdv[mxverknt];       //!< Number of daughter tracks
  int itdv[mxverknt];       //!< Number of first daughter track
  int ntra;                 //!< Number of tracks
  int iget[mxtraknt];       //!< Geant track number
  float ptra[mxtraknt][3];  //!< Momentum
  int itra[mxtraknt];       //!< Particle type
  int itvb[mxtraknt];       //!< Origin Vertex
  int itve[mxtraknt];       //!< End Vertex (if any)
  int nkinc;                //!< Number of daughters
  float xkinc[2];           //!< Values of X Feynman
  float pkinc[2][3];        //!< CM momentum
	      } QkinType;

#endif // CsComgNtCommons_h
