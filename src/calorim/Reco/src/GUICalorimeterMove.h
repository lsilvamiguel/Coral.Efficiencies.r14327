#ifndef RECO_GUICALORIMETERMOVE_H
#define RECO_GUICALORIMETERMOVE_H

#include "Reco_config.h"

#include <QtGui/QWidget>

#include "GUICalorimeterMove_uic.h"

namespace Reco {

class GUICalorimeterMove : public QWidget, public Ui::GUICalorimeterMove {
    public:
        GUICalorimeterMove(QWidget* parent = 0);
};

} // namespace Reco

#endif // RECO_GUICALORIMETERMOVE_H
