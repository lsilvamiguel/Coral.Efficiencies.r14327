
// $Id: Defs.h,v 1.3 2003/03/05 12:52:26 hpereira Exp $

#ifndef Defs_h
#define Defs_h
 
/*!
   \file    Defs.h
   \brief   Some global parameters used everywhere
   \author  Hugo Pereira
   \version $Revision: 1.3 $
   \date    $Date: 2003/03/05 12:52:26 $
*/

#define NPLAN 450         /*!< max number of detectors. */
#define BUFFER_SIZE 32000  /*!< root buffer max size before flushing to file */
#define PI 3.14159        /*!< pi */

#define ROOT_PATH "/CsEfficiency"  /*!< Root TDirectory where calibration objects are to be saved */
#define TREE_NAME "T_eff_"         /*!< Name of the root TTree to be saved */
#define LIST_NAME "TBName_list_"   /*!< Name of the root TObjString used to store the probed detectors TBNames */

#endif
