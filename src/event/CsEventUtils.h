// $Id: CsEventUtils.h,v 1.9 2010/02/15 10:26:49 tnagel Exp $

/*!
   \file    CsEventUtils.h
   \brief   Compass Event Utilities Class.
   \author  Benigno Gobbo
   \version $Revision: 1.9 $
   \date    $Date: 2010/02/15 10:26:49 $
*/


#ifndef CsEventUtils_h
#define CsEventUtils_h

#include "CsSTD.h"
#include "CsDetector.h"
#include "CsCluster.h"

/*! \class CsEventUtils 
    \brief Event Utilities Class.

    Contains some method to easy associate event quantities.
*/

class CsEventUtils{
 public:

  /* \fn static void associateCluster()
     \brief Perform cluster association between staggered planes (for DC, DR,
     DW, MB, ST) and compute Left/Right probabilities, based upon an option
     specified assumption on the origin (infinity or target pointing) of the
     particles.
     \warning Called only if doClusterAssociation_ is true.
  */
  static void associateClusters(void);

  /* \fn static void useLRFilter(list<CsCluster*> &clusters)
     \brief Remove elements from a list of pointers to clusters which does not
     pass the left/right probability cut. LRProbCut is taken from the
     CsDetector class.
     \param clusters The cluster list to be filtered 
     \warning This variable is changed!
  */
  static void useLRFilter(std::list<CsCluster*> &clusters);
	
  /* \fn static bool isALRGoodCluster(const CsCluster cluster)
     \brief Returns false if the Left/Right robability of a drift-like cluster is below LRProbCut (LRProbCut is taken from the CsDetector class). Returns true otherwise.
     \param clus The cluster to be tested
  */
  static bool isALRGoodCluster(const CsCluster cluster);

  /*! \fn static void LRmonitor(list<CsCluster*> clusters);
    \brief Build and fill histograms for Left/Right ambiguities.
    \param clusters The cluster list used for the monitoring
  */		
  static void LRmonitor(std::list<CsCluster*> clusters);

  /*! \fn static void SetLRToMatchAngle(CsCluster *c0, CsCluster *c1, const double a);
    \brief Performs the LR raising for argument clusters, based upon argument
    angle for the particle incidence: i.e. set the good LRProb values and put
    the pointers to the associated clusters in the proper order.
    \param c0, c1 Pointers to the cluster to be associated;
    \param a Tangent of the angle used for LR ambiguity raising;
  */
  static void setLRToMatchAngle(CsCluster *c0, CsCluster *c1, const double a);

 private:
   static double TPQ(CsCluster *c0, CsCluster *c1, const int LRMode);  //!< Target Pointing Quality

};

#endif // CsEventUtils_h
