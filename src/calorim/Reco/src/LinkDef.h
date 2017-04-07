/***************************************************************************
    $Header: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/LinkDef.h,v 1.1 2000/07/25 15:33:05 zvyagin Exp $
    This file is part of software for data analysing from GAMS experiment.
    ------------------------------------------------------------------------
    last modified   $Date: 2000/07/25 15:33:05 $
    copyright            : (C) 2000 by Zvyagin Alexander
    email                : Zvyagin@mx.ihep.su, Alexander.Zviagine@cern.ch
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//#pragma link C++ class Reco;
//#pragma link C++ class Run;
#pragma link C++ class Event-!;
#pragma link C++ class EventPhysics;
#pragma link C++ class EventCalibration;
#pragma link C++ class EventNoises;
#pragma link C++ class EventPedestals;
#pragma link C++ class RunID;
#pragma link C++ class EventID;

#endif
