#include <qapplication.h>
#include "DtpImp.h"

int main( int argc, char* argv[] )
{
  QApplication app( argc, argv );

  DtpImp dtp;
  app.setMainWidget( &dtp );
  dtp.show();

  int ret = app.exec();
  return ret;
}
