/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/GUICellsMatrix.h,v $
   $Date: 2011/02/01 16:20:24 $
   $Revision: 1.5 $
   -------------------------------------------------------------------------

   This file is part of cellular calorimeter reconstruction program.

   Authors:
     Vladimir  Kolosov   ( Kolosov@mx.ihep.su )
     Alexander Zvyagin   ( Alexander.Zviagine@cern.ch, Zvyagin@mx.ihep.su )

   Copyright(C): 1999-2000  V.Kolosov,A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef COLORCELLS_H
#define COLORCELLS_H

#include <vector>

#include "Reco_config.h"

#if USE_Qt

#include <QtGui/QGraphicsView>

class QMouseEvent;

namespace Reco {

class Calorimeter;

////////////////////////////////////////////////////////////////////////////////

class GUICellsMatrix : public QGraphicsView
{
  public:

    /// Base constructor
                        GUICellsMatrix          (QWidget* parent=0);

    /*! \brief Set calorimeter that will be shown
         The second argument is how many screen points will be in one unit of length.
         If point size of your display monitor is 0.23 mm, you use cm to represent lengths,
         than a rectangular cell with side size of 10 (cm) and screen_points_per_length_unit=10
         will be shown in 23x23 points box of your display.
    */
    virtual void        SetCalorimeter          ( Reco::Calorimeter &c );

    /*! \brief Set calorimeter SubSet that will be shown */
    virtual void        SetCalorimeterSubSets   ( Reco::Calorimeter &c,size_t subset_cells_number );

    void                Update                  (int what, int when, float min_value, float max_value);

  protected:

    void                mousePressEvent(QMouseEvent *e);
    void                resizeEvent(QResizeEvent* e);
    void                wheelEvent(QWheelEvent* e);

    Reco::Calorimeter *calorimeter;

    QGraphicsScene     *scene;

    std::vector < std::pair<size_t,QGraphicsItem *> >  pointers_to_items;

  private:

    QColor  color_to_restore;

};

////////////////////////////////////////////////////////////////////////////////

} // namespace Reco

#endif // USE_Qt

#endif // COLORCELLS_H
