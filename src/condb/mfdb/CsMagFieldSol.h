/*!
   \file CsMagFieldSol.h
   \brief COMPASS Magnetic Filed Database Class
   \author  Takeaki TOEDA
   \version $Revision: 1.1 $
   \date    $Date: 2000/05/23 14:24:43 $
*/

#ifndef CsMagFieldSol_h
#define CsMagFieldSol_h
#include "CsSTD.h"
#include "CsFieldGridSol.h"

/*! \class CsMagFieldSol
 *  \brief Magnetic Field Class for Solenoid.
 */

class CsMagFieldSol{
public:
  CsMagFieldSol();                                  //!< Default Constructor

  /*! \fn  list<CsFieldGridSol> getList();
    \brief Get grid list
  */
  list<CsFieldGridSol> getList(){return _list;}

  /*! \fn  void addGrid(CsFieldGridSol grid);
    \brief Add grid to list
  */
  void addGrid(CsFieldGridSol grid);


  /*! \fn  int getNoGrid();
    \brief Get number of grid
  */
  int getNoGrid();
private:
  list<CsFieldGridSol> _list;    //!< list of grid
};
#endif
