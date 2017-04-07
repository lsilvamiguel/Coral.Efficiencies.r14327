#ifndef GUICalorimeter___include
#define GUICalorimeter___include

#include "Reco_config.h"
#if USE_Qt
#include "GUIBaseCalorimeter.h"

#include <Qt/qobject.h>
#include <Qt/qtimer.h>
#include <QtGui/qspinbox.h>

#else
typedef int QWidget;
namespace Qt {
  typedef int WFlags;
}
#endif

namespace Reco {

class Calorimeter;

class GUICalorimeter
#if USE_Qt
: public GUIBaseCalorimeter
#endif
{
#if USE_Qt
  Q_OBJECT
#endif

  // =========================================================================
  // Constructors and destructor
  // =========================================================================

 public:

  GUICalorimeter          ( Calorimeter &c, QWidget* parent = 0, const char* name = 0, Qt::WFlags fl = 0,
			    bool xy_select = false );

  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~GUICalorimeter()  {}

 private:
  GUICalorimeter          (const GUICalorimeter &g);

  // =========================================================================
  // Operators
  // =========================================================================

  GUICalorimeter     &operator =              (const GUICalorimeter &g);

  // =========================================================================
  // Methods
  // =========================================================================

 public:
#if USE_Qt
    public slots:
#endif
  virtual void        UpdateForPreviousSpill  (void);
 protected:
#if USE_Qt
  protected slots:
#endif

  /// Set update time in seconds.
  virtual void        SetAutoUpdateTime       (unsigned t);

  virtual void        Update                  (void);

  virtual void        ShowEvent               (void);

  virtual void        ShowCells               (void);

  virtual void        Draw                    (void);

  virtual void        Fit                     (void);

  virtual void        FitOK                   (void);

  virtual void        SaveFit                 (void);

  virtual void        PrintPS                 (void);

  virtual void        ResetHisto              (void);

  virtual void        ChangeRange             (void);

  void                Dummy                   (void);
  // =========================================================================
  // Attributes
  // =========================================================================

 protected:

  /// Reference to calorimeter
  Calorimeter        &cal;

  ///  Define cells selection(for fit and drawing) policy
  bool               xy_selection;

 private:

#if USE_Qt
  QTimer              qt_timer;
#endif
};

}

#endif // GUICalorimeter___include
