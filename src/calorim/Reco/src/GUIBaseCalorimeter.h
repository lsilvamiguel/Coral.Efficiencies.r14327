#ifndef RECO_GUIBASECALORIMETER_H
#define RECO_GUIBASECALORIMETER_H

#include "Reco_config.h"

#include <QtGui/QWidget>

#include "GUIBaseCalorimeter_uic.h"

namespace Reco {

class GUIBaseCalorimeter : public QWidget, public Ui::GUIBaseCalorimeter {
    Q_OBJECT

    public:
                        GUIBaseCalorimeter(QWidget* parent = 0, const char* name = 0, Qt::WFlags fl = 0,bool xy_selection = false);
        virtual        ~GUIBaseCalorimeter();

    private slots:
        virtual void    selMinChanged();
        virtual void    selMaxChanged();
        virtual void    fitXMinChanged();
        virtual void    fitXMaxChanged();
        virtual void    fitYMinChanged();
        virtual void    fitYMaxChanged();
};

} // namespace Reco

#endif // RECO_GUIBASECALORIMETER_H
