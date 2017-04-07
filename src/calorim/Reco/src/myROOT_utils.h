/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/myROOT_utils.h,v $
   $Date: 2010/06/18 10:44:20 $
   $Revision: 1.3 $
   -------------------------------------------------------------------------

   Some useful utilites to work with ROOT

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch,Kolosov@mx.ihep.su )

   Copyright(C): 1999-2001  V.Kolosov, A.Zvyagin, D.Murashev

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

#ifndef myROOT_utils___include
#define myROOT_utils___include

#include "Reco_config.h"
#include <cassert>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TTree.h"
#include "TObject.h"
#include "TKey.h"
#include "TList.h"

class myROOT_utils
{

public:

static TDirectory * TDirectory_checked(const char* name);

static TH1D * TH1D_checked(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup);

static TH2D * TH2D_checked(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup, Int_t nbinsy, Axis_t ylow, Axis_t yup);

static TProfile * TProfile_checked(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup);

static void ResetRecursive(TDirectory* node);

};

#endif  // myROOT_utils___include
