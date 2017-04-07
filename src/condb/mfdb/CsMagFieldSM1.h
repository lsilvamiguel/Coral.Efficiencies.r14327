/*!
   \file CsMagFieldSM1.h
   \brief COMPASS Magnetic Filed Database Class
   \author  Takeaki TOEDA
   \version $Revision: 1.2 $
   \date    $Date: 2000/05/25 16:33:32 $
*/

#ifndef CsMagFieldSM1_h
#define CsMagFieldSM1_h
#include "CsSTD.h"
#include "CsFieldGridSM1.h"

/*! \class CsMagFieldSM1
 *  \brief Magnetic Field Class for SM1.
 */


class CsMagFieldSM1{
public:
  CsMagFieldSM1();                                  //!< Default Constructor

  /*! \fn  list<CsFieldGridSol> getList();
    \brief Get grid list
  */
  list<CsFieldGridSM1> getList(){return _list;}

  /*! \fn  void addGrid(CsFieldGridSM1 grid);
    \brief Add grid to list
  */
  void addGrid(CsFieldGridSM1 grid);

  /*! \fn  int getNoGrid();
    \brief Get number of grid
  */
  int getNoGrid();
private:
  list<CsFieldGridSM1> _list;                       //!< list of grid
};
#endif
