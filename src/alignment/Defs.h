// $Id: Defs.h,v 1.11 2008/06/05 13:13:12 rgazda Exp $

#ifndef Defs_h
#define Defs_h
 
/*!
   \file    Defs.h
   \brief   Some global parameters used everywhere
   \author  Hugo Pereira
   \version $Revision: 1.11 $
   \date    $Date: 2008/06/05 13:13:12 $
*/

#define BUFFER_SIZE 32000  /*!< root buffer max size before flushing to file */
#define NPLAN 370         /*!< max number of detectors. */
#define NCLUST 20         /*!< max number of clusters per plane. may bias the tree*/
#define NDIGIT 30         /*!< max number of digits per cluster*/
#define DURANGE 20        /*!< max residual (in number of pitches). may bias the tree*/

/*! number of parameters/track
  0: x           ( perp to the beam horizontal) 
  1: tan(teta_x) ( in Ox,Oz plane)
  2: y           ( perp to the beam vertical
  3: tan(teta_y) (in Ox,Oz plane)
*/ 
#define NPARTRCK 4	      

/*! number of alignement parameters per plane:
  0: U offset perp to the wire: U+=\alpha_U
  1: Z offset along the beam: Z+=\alpha_Z
  2: T rotational offset, perp to the beam: T+=\alpha_T
  3: P pitch: U*=(1+\alpha_P)
  4: R R0 offset for the distance to the wire, proportional to T0 for drift like detectors: U+=sign(R)*\alpha_R
  5: L lorentz angle RT scaling for drift like detectors: R*=(1+\alpha_L) i.e. U+=\alpha_L*R
*/
#define NPARPLAN 6	

#define NGLB NPLAN*NPARPLAN /*!< maximum number of global parameters */ 
#define PI 3.14159          /*!< pi */

#endif




