#include "Calorimeter.h"
#include "GUICalorimeter.h"
#include "GUICellsMatrix.h"

using namespace std;

namespace Reco {

#if USE_Qt

////////////////////////////////////////////////////////////////////////////////

GUICalorimeter::GUICalorimeter(Calorimeter &c, QWidget* parent, const char* name, Qt::WFlags fl,bool xy_select) :
  GUIBaseCalorimeter(parent,name,fl,xy_select),
  cal(c),
  xy_selection(xy_select),
  qt_timer(this)
{
  setWindowTitle(name);

  gui_cells_table -> setSortingEnabled(false); // Sorting is too slow!

  //gui_cells_table -> setLeftMargin(QFontMetrics(font()).width("88888"));

  gui_cells_table -> setSelectionMode(QAbstractItemView::NoSelection);

  gui_cells_table->verticalHeader()->setResizeMode(QHeaderView::Fixed);
  gui_cells_table->verticalHeader()->setClickable(false);

  gui_cells_table->horizontalHeader()->setResizeMode(QHeaderView::Interactive);
  gui_cells_table->horizontalHeader()->setClickable(false);
//    h->setSortIndicator(0,true);

  gui_cells_table -> setColumnCount(5);

  QTableWidgetItem* hdrCol0 = new QTableWidgetItem("Type");
  QTableWidgetItem* hdrCol1 = new QTableWidgetItem("Position");
  QTableWidgetItem* hdrCol2 = new QTableWidgetItem("HV");
  QTableWidgetItem* hdrCol3 = new QTableWidgetItem("PED");
  QTableWidgetItem* hdrCol4 = new QTableWidgetItem("Calib");

  gui_cells_table->setHorizontalHeaderItem(0, hdrCol0);
  gui_cells_table->setHorizontalHeaderItem(1, hdrCol1);
  gui_cells_table->setHorizontalHeaderItem(2, hdrCol2);
  gui_cells_table->setHorizontalHeaderItem(3, hdrCol3);
  gui_cells_table->setHorizontalHeaderItem(4, hdrCol4);

  if(c.GetOptions().show_all_cells)
    cells_matrix->SetCalorimeter(cal);

//   cells_matrix_0->SetCalorimeter(cal,screen_points_per_length_unit_x,
//                                     screen_points_per_length_unit_y );
//
  cut1_value_from->setMaximum(cal.NCells());
  cut2_value_from->setMaximum(cal.NCells());
  cut1_value_to->setMaximum(cal.NCells());
  cut2_value_to->setMaximum(cal.NCells());


//  connect( &qt_timer, SIGNAL(timeout()), this, SLOT(ShowEvent()) );

//  connect( &qt_timer, SIGNAL(timeout()), this, SLOT(Update()) );

  connect ((QObject*)PushButton3, SIGNAL (clicked ()), this, SLOT (ShowEvent ()));

  connect ((QObject*)PushButton4, SIGNAL (clicked ()), this, SLOT (Update ()));

  connect ((QObject*)PushButton5_2, SIGNAL (clicked ()), this, SLOT (ShowCells ()));

  connect ((QObject*)PushButton7, SIGNAL (clicked ()), this, SLOT (Draw ()));

  connect ((QObject*)PushButton8, SIGNAL (clicked ()), this, SLOT (Fit ()));

  connect ((QObject*)PushButton6, SIGNAL (clicked ()), this, SLOT (PrintPS ()));

  connect ((QObject*)button1, SIGNAL (clicked ()), this, SLOT (ResetHisto ()));

  connect ((QObject*)button3, SIGNAL (clicked ()), this, SLOT (FitOK ()));

  connect ((QObject*)button2, SIGNAL (clicked ()), this, SLOT (SaveFit ()));

  connect (selection, SIGNAL (activated (int)), this, SLOT (ChangeRange ()));
  connect (selection_2, SIGNAL (activated (int)), this, SLOT (ChangeRange ()));
//   connect (selection_2, SIGNAL (highlighted (int)), this, SLOT (ChangeRange ()));
//   connect (selection, SIGNAL (highlighted (int)), this, SLOT (ChangeRange ()));

  SetAutoUpdateTime(2); // Update every two seconds

  Update();
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::Update(void)
{
  if( gui_cells_table->rowCount()!=(int)cal.NCells() )
    gui_cells_table -> setRowCount(cal.NCells());

  for( size_t i=0; i<cal.NCells(); i++ )
  {
    const Cell &c = cal.GetCells()[i];

    QTableWidgetItem* itemCol0 = new QTableWidgetItem(c.GetCellType().GetName().c_str());
    QTableWidgetItem* itemCol1 = new QTableWidgetItem(QString().sprintf("(%g;%g;%g)",c.GetX(),c.GetY(),c.GetZ()));
    QTableWidgetItem* itemCol2 = new QTableWidgetItem(QString().sprintf("?"));
    QTableWidgetItem* itemCol3 = new QTableWidgetItem(QString().sprintf("N=%.1e m=%.2e s=%.2e",
                                                                        (float)cal.GetCellInfo(Calorimeter::PED,   Calorimeter::NEW, i).GetEntries(),
                                                                        (float)cal.GetCellInfo(Calorimeter::PED,   Calorimeter::NEW, i).GetMean(),
                                                                        (float)cal.GetCellInfo(Calorimeter::PED,   Calorimeter::NEW, i).GetSigma()));
    QTableWidgetItem* itemCol4 = new QTableWidgetItem(QString().sprintf("N=%.1e m=%.2e s=%.2e",
                                                                        (float)cal.GetCellInfo(Calorimeter::CALIB, Calorimeter::NEW, i).GetEntries(),
                                                                        (float)cal.GetCellInfo(Calorimeter::CALIB, Calorimeter::NEW, i).GetMean(),
                                                                        (float)cal.GetCellInfo(Calorimeter::CALIB, Calorimeter::NEW, i).GetSigma()));

    gui_cells_table->setItem(i, 0, itemCol0);
    gui_cells_table->setItem(i, 1, itemCol1);
    gui_cells_table->setItem(i, 2, itemCol2);
    gui_cells_table->setItem(i, 3, itemCol3);
    gui_cells_table->setItem(i, 4, itemCol4);
  }

  gui_cells_table->resizeColumnsToContents();
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::UpdateForPreviousSpill(void)
{
  const EventID &evid = cal.GetEventID();
  if( evid.Initialized() )
  {
//    cout << " EventID Run " << evid.GetRunNumber() <<" Spill " << evid.GetBurstNumber() << endl;
//    char s[100];
//    sprintf(s," Run %d Spill %d ",evid.GetRunNumber(),evid.GetBurstNumber() );
//    update_time_last->setText( tr( s  ) );
    update_time_last->setText( QString().sprintf(" Run %d Spill %d ",evid.GetRunNumber(),evid.GetBurstNumber() ) );
//    update_time_last->setText(i,1,QString().sprintf("(%g;%g;%g)",c.GetX(),c.GetY(),c.GetZ()));
  }
  else
  {
    update_time_last->setText( QString().sprintf(" Run ?? Spill ?? " ) );
//    cout << " EventID not initialized " << endl;
  }

  ShowCells();
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::SetAutoUpdateTime(unsigned t)
{
  qt_timer.start(t*1000);
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::ShowCells (void)
{
  int ct = selection->currentIndex();
  int vt = selection_2->currentIndex();
  float min_value = (float)cut1_value_min->value ();
  float max_value = (float)cut1_value_max->value ();
  int units = selection_units->currentIndex();
  float units_scale=1.;
  if(units == 1 ) units_scale=0.1;
  else if(units == 2 ) units_scale=0.01;
  else if(units == 3 ) units_scale=0.001;
  else if(units == 4 ) units_scale=10.;
  else if(units == 5 ) units_scale=100.;
  else if(units == 6 ) units_scale=1000.;
  float fmin = min_value*units_scale;
  float fmax = max_value*units_scale;
  cells_matrix->Update (ct, vt, fmin, fmax);
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::Draw (void)
{
  bool debug=true;
  if( debug ) std::cout << " GUICsCalorimeter::Draw() " << cal.GetName() << endl;
#if USE_Qt

  int min_value = cut1_value_min->value ();
  int max_value = cut1_value_max->value ();
  int units = selection_units->currentIndex();
  double units_scale=1.;
  if(units == 1 ) units_scale=0.1;
  else if(units == 2 ) units_scale=0.01;
  else if(units == 3 ) units_scale=0.001;
  else if(units == 4 ) units_scale=10.;
  else if(units == 5 ) units_scale=100.;
  else if(units == 6 ) units_scale=1000.;
  double fmin = min_value*units_scale;
  double fmax = max_value*units_scale;

  int selection_type= selection->currentIndex();
//   cout << " GUICsCalorimeter::Draw selection " << selection_type  << endl;
  if(selection_type <= 0) return;
  Reco::Calorimeter::CellInfoType info_type=Reco::Calorimeter::CALIB;
  if(selection_type == 1) info_type = Reco::Calorimeter::CALIB;
  else if(selection_type == 2) info_type = Reco::Calorimeter::LED;
  else if(selection_type == 3) info_type = Reco::Calorimeter::PED;
  else if(selection_type == 4) info_type = Reco::Calorimeter::TIME;
  else if(selection_type == 5) info_type = Reco::Calorimeter::NOISE;

  if( xy_selection )
  {
    int xmin = cut1_value_from->value ();
    int xhow_much = cut1_value_to->value ();
    int xstep = cut1_value_step->value ();

    int ymin = y1_value_from->value ();
    int yhow_much = y1_value_to->value ();
    int ystep = y1_value_step->value ();

    cal.DrawFitCellsHistoXY(info_type,"DRAW",xmin,xhow_much,xstep,ymin,yhow_much,ystep,fmin,fmax);
  }
  else
  {
    int min = cut1_value_from->value ();
    int how_much  = cut1_value_to->value ();
    int step = cut1_value_step->value ();
    if(step<=0) step=1;
    if( selection_2->currentIndex() == 5 )
    {
      cal.DrawMonitorHisto(info_type,"DRAW",min,how_much,step,fmin,fmax);
    }
    else
    {
      cal.DrawFitCellsHisto(info_type,"DRAW",min,how_much,step,fmin,fmax);
    }
  }

#endif
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::Fit (void)
{
#if USE_Qt

  int min_value = cut1_value_min->value ();
  int max_value = cut1_value_max->value ();
  int units = selection_units->currentIndex();
  double units_scale=1.;
  if(units == 1 ) units_scale=0.1;
  else if(units == 2 ) units_scale=0.01;
  else if(units == 3 ) units_scale=0.001;
  else if(units == 4 ) units_scale=10.;
  else if(units == 5 ) units_scale=100.;
  else if(units == 6 ) units_scale=1000.;
  double fmin = min_value*units_scale;
  double fmax = max_value*units_scale;

  int selection_type= selection->currentIndex();
//   cout << " GUICsCalorimeter::Draw selection " << selection_type  << endl;
  if(selection_type <= 0) return;
  Reco::Calorimeter::CellInfoType info_type=Reco::Calorimeter::CALIB;
  if(selection_type == 1) info_type = Reco::Calorimeter::CALIB;
  else if(selection_type == 2) info_type = Reco::Calorimeter::LED;
  else if(selection_type == 3) info_type = Reco::Calorimeter::PED;
  else if(selection_type == 4) info_type = Reco::Calorimeter::TIME;
  else if(selection_type == 5) info_type = Reco::Calorimeter::NOISE;

  if( selection_2->currentIndex() == 6 )   //  AutoFit
  {
    cal.FitAllCellsHisto(info_type,"AUTO");
    return;
  }

  {

    if( xy_selection )
    {
      int xmin = cut1_value_from->value ();
      int xhow_much = cut1_value_to->value ();
      int xstep = cut1_value_step->value ();

      int ymin = y1_value_from->value ();
      int yhow_much = y1_value_to->value ();
      int ystep = y1_value_step->value ();

      cal.DrawFitCellsHistoXY(info_type,"FIT",xmin,xhow_much,xstep,ymin,yhow_much,ystep,fmin,fmax);
    }
    else
    {
      int min = cut1_value_from->value ();
      int how_much  = cut1_value_to->value ();
      int step = cut1_value_step->value ();
      if(step<=0) step=1;
      if( selection_2->currentIndex() == 5 )
      {
        cal.DrawMonitorHisto(info_type,"FIT",min,how_much,step,fmin,fmax);
      }
      else
      {
        cal.DrawFitCellsHisto(info_type,"FIT",min,how_much,step,fmin,fmax);
      }
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::FitOK (void)
{
#if USE_Qt
  int selection_type= selection->currentIndex();
//   cout << " GUICsCalorimeter::Draw selection " << selection_type  << endl;
  if(selection_type <= 0) return;
  Reco::Calorimeter::CellInfoType info_type=Reco::Calorimeter::CALIB;
  if(selection_type == 1) info_type = Reco::Calorimeter::CALIB;
  else if(selection_type == 2) info_type = Reco::Calorimeter::LED;
  else if(selection_type == 3) info_type = Reco::Calorimeter::PED;
  else if(selection_type == 4) info_type = Reco::Calorimeter::TIME;
  else if(selection_type == 5) info_type = Reco::Calorimeter::NOISE;

  if( xy_selection )
  {
    int xmin = cut1_value_from->value ();
    int xhow_much = cut1_value_to->value ();
    int xstep = cut1_value_step->value ();

    int ymin = y1_value_from->value ();
    int yhow_much = y1_value_to->value ();
    int ystep = y1_value_step->value ();

    cal.DrawFitCellsHistoXY(info_type,"FIT_ON_OFF",xmin,xhow_much,xstep,ymin,yhow_much,ystep);
  }
  else
  {
    int min = cut1_value_from->value ();
    int how_much  = cut1_value_to->value ();
    int step = cut1_value_step->value ();
    if(step<=0) step=1;
    cal.DrawFitCellsHisto(info_type,"FIT_ON_OFF",min,how_much,step);
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::SaveFit (void)
{
  int selection_type= selection->currentIndex();
  if(selection_type <= 0) return;
  Reco::Calorimeter::CellInfoType info_type=Reco::Calorimeter::CALIB;
  if(selection_type == 1) info_type = Reco::Calorimeter::CALIB;
  else if(selection_type == 2) info_type = Reco::Calorimeter::LED;
  else if(selection_type == 3) info_type = Reco::Calorimeter::PED;
  else if(selection_type == 4) info_type = Reco::Calorimeter::TIME;
  cal.SaveFit2CellInfo(info_type);
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::ChangeRange (void)
{
 int max;
  switch ( selection->currentIndex() )
  {
    case 1:
        switch ( selection_2->currentIndex() )
        {
          case 0:
            max = 1000;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
          case 1:
            max = 1000;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
          case 2:
            max = 100;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
          case 3:
            max = 1000;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
        }
        break;
    case 2:
        switch ( selection_2->currentIndex() )
        {
          case 0:
            max = 1000;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
          case 1:
            max = 4096;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);

            break;
          case 2:
            max = 100;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
          case 3:
            max = 4096;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
        }
        break;
    case 3:
        switch ( selection_2->currentIndex() )
        {
          case 0:
            max = 1000;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
          case 1:
            max = 1000;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);

            break;
          case 2:
            max = 100;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
          case 3:
            max = 4096;
            cut1_value_min->setMaximum(max);
            cut2_value_min->setMaximum(max);
            cut1_value_max->setMaximum(max);
            cut2_value_max->setMaximum(max);
            break;
        }
        break;
    default: return;
  }
  cout << " Change Range " << endl;
}

////////////////////////////////////////////////////////////////////////////////

#else // USE_Qt=0

////////////////////////////////////////////////////////////////////////////////

// GUICalorimeter::GUICalorimeter     (Calorimeter &c,double screen_points_per_length_unit_x,
//               double screen_points_per_length_unit_y,QWidget* parent, const char* name, WFlags fl,bool xy_select)
//   : cal(c),
//     xy_selection(xy_select)

GUICalorimeter::GUICalorimeter     (Calorimeter &c, QWidget* parent, const char* name,
                                                                  Qt::WFlags fl,bool xy_select)
  : cal(c)
  {}

void GUICalorimeter::Update (void) {}

void GUICalorimeter::ShowCells (void) {}

void GUICalorimeter::Draw (void) {}

void GUICalorimeter::Fit (void) {}

void GUICalorimeter::FitOK (void) {}

void GUICalorimeter::SaveFit (void) {}

void GUICalorimeter::ChangeRange (void) {}

// void GUICalorimeter::PrintPS (void) {}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::SetAutoUpdateTime  (unsigned t) {}

void GUICalorimeter::UpdateForPreviousSpill  (void) {}

#endif  // USE_Qt

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::ShowEvent (void)
{
  cal.ShowEvent();
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::Dummy (void)
{
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::ResetHisto (void)
{
  cal.ResetHisto();
  cal.ResetStatistic();
}

////////////////////////////////////////////////////////////////////////////////

void GUICalorimeter::PrintPS (void)
{
  cal.PrintPS();
}

////////////////////////////////////////////////////////////////////////////////


}
