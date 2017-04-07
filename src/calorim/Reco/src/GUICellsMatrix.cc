#include <cstdio>
#include <cmath>

#include "Calorimeter.h"
#include "CellDataRaw.h"
#include "StatInfo.h"

#include "GUICellsMatrix.h"
#include "GUICalorimeter.h"

#if USE_Qt

#include <QtGui/QGraphicsItem>
#include <QtGui/QMouseEvent>

using namespace std;

Reco::GUICellsMatrix::GUICellsMatrix(QWidget* parent) :
  QGraphicsView   (parent),
  calorimeter   (NULL),
  scene         (NULL)
{
  scene = new QGraphicsScene();
  scene->setBackgroundBrush(QBrush(Qt::white, Qt::SolidPattern));
  setScene(scene);

  setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

  pointers_to_items.clear();
  color_to_restore=Qt::white;
}

////////////////////////////////////////////////////////////////////////////////

void Reco::GUICellsMatrix::SetCalorimeter( Reco::Calorimeter &c )
{
  if( calorimeter!=NULL )
    throw "GUICellsMatrix::SetCalorimeter(): calorimeter has been set already";
  calorimeter = &c;

  for( size_t icell =0; icell < calorimeter->NCells(); icell++)
  {
    const Reco::Cell &cell = calorimeter->GetCells()[icell];

    const double x1 =   cell.GetX() - cell.GetCellType().GetSizeX()/2.;
    const double y1 = -(cell.GetY() + cell.GetCellType().GetSizeY()/2.);

    const double xw = cell.GetCellType().GetSizeX();
    const double yw = cell.GetCellType().GetSizeY();

    QGraphicsItem* item = scene->addRect(QRectF(x1,y1,xw,yw));
    pointers_to_items.push_back(pair<size_t,QGraphicsItem *>(icell,item));
    item->show();
  }
}

////////////////////////////////////////////////////////////////////////////////

void Reco::GUICellsMatrix::SetCalorimeterSubSets( Reco::Calorimeter &c,size_t subset_number )
{
  if( calorimeter!=NULL )
    throw "GUICellsMatrix::SetCalorimeter(): calorimeter has been set already";
  calorimeter = &c;

  for( size_t i =0;  i< calorimeter->GetSubSets()[subset_number].Size();  i++ )
  {
    size_t icell = calorimeter->GetSubSets()[subset_number].GetCells()[i];
    const Reco::Cell &cell = calorimeter->GetCells()[icell];

    const double x1 =   cell.GetX() - cell.GetCellType().GetSizeX()/2.;
    const double y1 = -(cell.GetY() + cell.GetCellType().GetSizeY()/2.);

    const double xw = cell.GetCellType().GetSizeX();
    const double yw = cell.GetCellType().GetSizeY();

    QGraphicsItem *item = scene->addRect(QRectF(x1,y1,xw,yw));
    pointers_to_items.push_back(pair<size_t,QGraphicsItem *>(icell,item));
    item->show();
  }
}

////////////////////////////////////////////////////////////////////////////////

void Reco::GUICellsMatrix::mousePressEvent(QMouseEvent *e)
{
  QList<QGraphicsItem*> l=scene->items(mapToScene(e->pos()));

  if (l.isEmpty ()) return;
  string name;
  for( QList<QGraphicsItem*>::Iterator it=l.begin(); it!=l.end(); ++it )
  {
    char s[111];
    sprintf(s,"cell number %d\n",scene->items().indexOf(*it));
    name += s;
    QGraphicsRectItem* r = dynamic_cast<QGraphicsRectItem*>(*it);
    if( r==NULL )
      throw "GUICellsMatrix::contentsMousePressEvent():  internal error";
    QColor it_was_color = r->brush().color();
    QColor color_to_be_restored=Qt::white;
    if(it_was_color != Qt::gray )
    {
      color_to_be_restored=color_to_restore;
      color_to_restore=it_was_color;
    }
//    r->setBrush( QBrush(r->brush().color()==Qt::blue?Qt::white:Qt::blue,Qt::SolidPattern) );
    r->setBrush( QBrush(Qt::gray, Qt::SolidPattern) );
    for( vector<pair<size_t,QGraphicsItem*> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
    {
      QGraphicsRectItem* rr = dynamic_cast<QGraphicsRectItem*>(item->second);
      if (r == rr )
      {
        size_t i = item->first;
        if (e->button () == Qt::LeftButton)
        {
//           cout << " This is cell # " << i << " Of calorimeter " << calorimeter->GetName() << endl;
          if(calorimeter->GetGUI().selection_2->currentIndex() == 5)
	    calorimeter->DrawMonitorCellHisto(i);
          else if(calorimeter->GetGUI().selection_2->currentIndex() == 8)
	    calorimeter->DrawMonitorCellHistoInSpill(i);
          else
	    calorimeter->DrawCellHisto(i);
          calorimeter->GetGUI ().cell_text->setText (tr (calorimeter->GetCellName (i).c_str ()));
        }
        else
        if (e->button () == Qt::RightButton)
        {
          calorimeter->GetGUI ().cell_text->setText (tr (calorimeter->GetCellName (i).c_str ()));
        }
      }
      else
      {
        if(rr->brush().color()==Qt::gray) rr->setBrush( QBrush(color_to_be_restored,Qt::SolidPattern) );
      }
    }
    break; // comment this line to allow selection of several cells
  }
  scene->update();
    //cout << name.c_str() << "\n";

//     QMessageBox("Selected cells",name.c_str(),QMessageBox::Information,
//                  QMessageBox::NoButton,QMessageBox::NoButton,QMessageBox::NoButton)
//       .show();
}

////////////////////////////////////////////////////////////////////////////////

void Reco::GUICellsMatrix::resizeEvent(QResizeEvent* e) {
    fitInView(scene->sceneRect(), Qt::KeepAspectRatio);
}

////////////////////////////////////////////////////////////////////////////////

void Reco::GUICellsMatrix::wheelEvent(QWheelEvent* e) {
    if ((e->modifiers()&Qt::ShiftModifier) == 0) {
        QGraphicsView::wheelEvent(e);
        return;
    }

    const QPointF before = mapToScene(e->pos());

    const double reScale = std::pow(1.1, e->delta()/120.);
    scale(reScale, reScale);

    centerOn(before);
}

////////////////////////////////////////////////////////////////////////////////

void Reco::GUICellsMatrix::Update (int ct, int vt, float min_value, float max_value)
  {
    float value;
    Reco::Calorimeter::CellInfoType what;
    Reco::Calorimeter::CalibTimeType when = Reco::Calorimeter::NEW;

    switch (ct) {
      case 0: what = Reco::Calorimeter::RAW;  break;
      case 1: what = Reco::Calorimeter::CALIB;  break;
      case 2: what = Reco::Calorimeter::LED;    break;
      case 3: what = Reco::Calorimeter::PED;    break;
      case 4: what = Reco::Calorimeter::TIME;   break;
      case 5: what = Reco::Calorimeter::NOISE;  break;
      default: what = Reco::Calorimeter::RAW;  break;
//       default: return;
    }
    switch (vt) {
      case 0:                                          //  Draw statistic
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          size_t i = item->first;
          value = (float)calorimeter->GetCellInfo (what, when, i).GetEntries();
          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
          if (value <= min_value) r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
          else if (value <= max_value) r->setBrush (QBrush (Qt::yellow, Qt::SolidPattern));
          else r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
        }
        break;
      case 1:                                          //  Draw mean
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          size_t i = item->first;
          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
          double  entries = calorimeter->GetCellInfo (what, when, i).GetEntries();
          if (entries <= 0 )
          {
            r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
            continue;
          }
          value = (float)calorimeter->GetCellInfo (what, when, i).GetMean();
          if (value <= min_value ) r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
          else if (value <= max_value) r->setBrush (QBrush (Qt::yellow, Qt::SolidPattern));
          else r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
        }
        break;
      case 2:                                          //  Draw sigma
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          size_t i = item->first;
          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
          double  entries = calorimeter->GetCellInfo (what, when, i).GetEntries();
          if (entries <= 0 )
          {
            r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
            continue;
          }
          value = (float)calorimeter->GetCellInfo (what, when, i).GetSigma();
          if (value <= min_value) r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
          else if (value <= max_value) r->setBrush (QBrush (Qt::yellow, Qt::SolidPattern));
          else r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
        }
        break;
      case 3:                                          //  Draw event
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          int i = item->first;
          value = 0;
          if( what == Reco::Calorimeter::LED || what == Reco::Calorimeter::PED )
          {
            for(vector<Reco::CellDataRaw>::const_iterator it=calorimeter->GetSignals().begin();
                                                  it!=calorimeter->GetSignals().end(); it++)
            {
               int address = it->GetCellIdx();
               if( address == i ) value = it->GetAmplitude();
            }
          }
          else
          {
            for(vector<Reco::CellDataRaw>::const_iterator it=calorimeter->GetSignals().begin();
                                                  it!=calorimeter->GetSignals().end(); it++)
            {
               int address = it->GetCellIdx();
               if( address == i ) value = 100.*it->GetEnergy();
            }
          }
          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
          if (value <= min_value) r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
          else if (value <= max_value) r->setBrush (QBrush (Qt::yellow, Qt::SolidPattern));
          else r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
        }
        break;
      case 4:                                          //  Draw fit validity
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          size_t i = item->first;
          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
          bool  fit_ok = calorimeter->GetFitInfo(what,i).fit_ok;
          if (fit_ok )
          {
            r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
            continue;
          }
          else r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
        }
        break;
      case 5:                                          //  Draw monitoring statistic
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          size_t i = item->first;
	  Reco::StatInfoStore * store = calorimeter->GetStatInfoStore(what, when, i);
          if( store != NULL) value = (float)store->CalculateAverage().GetEntries();
	  else value =0;

          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
          if (value <= min_value) r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
          else if (value <= max_value) r->setBrush (QBrush (Qt::yellow, Qt::SolidPattern));
          else r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
        }
        break;
      case 6:                                          //  AutoFit
        break;
      case 7:                                          //  Draw monitoring statistic
        when = Reco::Calorimeter::OLD;
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          size_t i = item->first;
          double  value = calorimeter->GetCellInfo (what, when, i).GetMean();
          if( ct == 0 )  value = calorimeter->GetEnergyCutInCell( i); // No RAW calibrations, return ecut
          if( ct == 5 )  value = calorimeter->GetMinEnergyCutInCell( i); // No OLD NOISE calibrations, return min ecut

          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
          if (value <= min_value) r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
          else if (value <= max_value) r->setBrush (QBrush (Qt::yellow, Qt::SolidPattern));
          else r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
        }
        break;
      case 8:                                          //  Draw Mean value NEW/PRIM
//        when = Reco::Calorimeter::NEW;
//        when = Reco::Calorimeter::OLD;
        when = Reco::Calorimeter::MONITOR;
//	cout <<" Development /Ref value  ??????????? " << endl;
//         for (size_t i = 0; i < calorimeter->NCells (); i++)
        for( vector<pair<size_t,QGraphicsItem *> >::const_iterator item = pointers_to_items.begin();
                                                             item!= pointers_to_items.end(); item++ )
        {
          size_t i = item->first;
          double  value = calorimeter->GetCellInfo (what, when, i).GetMean();
          double  value_ref = calorimeter->GetCellInfo (what, Reco::Calorimeter::PRIM, i).GetMean();
//	  cout << i << " value " << value << " value_ref " << value_ref << endl;
          if( ct == 0 )  value = calorimeter->GetEnergyCutInCell( i); // No RAW calibrations, return ecut
          if( ct == 5 )  value = calorimeter->GetMinEnergyCutInCell( i); // No OLD NOISE calibrations, return min ecut
 	  bool refok = true;
	  double ratio = value;
          if( value_ref > 20. )
	  {
	    ratio = 100.*(1.- value/value_ref);
	    refok = true;
	  }
	  else
	  {
	    ratio = 0.;
	    refok = false;
	  }

	  ratio = fabs(ratio);

          QGraphicsRectItem *r = dynamic_cast<QGraphicsRectItem*>(item->second);
	  if( refok )
	  {
            if ( value < 20. ) r->setBrush (QBrush (Qt::black, Qt::SolidPattern));
            else if (ratio <= min_value) r->setBrush (QBrush (Qt::green, Qt::SolidPattern));
            else if (ratio <= max_value) r->setBrush (QBrush (Qt::yellow, Qt::SolidPattern));
            else r->setBrush (QBrush (Qt::red, Qt::SolidPattern));
	  }
	  else
	  {
            if ( value > 20. )
	      r->setBrush (QBrush (Qt::white, Qt::SolidPattern));
	    else
              r->setBrush (QBrush (Qt::blue, Qt::SolidPattern));
	  }
        }
        break;
    }
    scene->update ();
  }

////////////////////////////////////////////////////////////////////////////////

#endif
