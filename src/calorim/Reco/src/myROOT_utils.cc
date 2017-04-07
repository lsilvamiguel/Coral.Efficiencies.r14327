/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/myROOT_utils.cc,v $
   $Date: 2010/06/18 10:44:20 $
   $Revision: 1.4 $
   -------------------------------------------------------------------------

   Some useful utilites to work with ROOT

   Authors:
     Vladimir  Kolosov   ( Vladimir.Kolosov@cern.ch,Kolosov@mx.ihep.su )

   Copyright(C): 2002  V.Kolosov

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

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cstdlib>   // for exit()

#include "Reco_config.h"

#include "myROOT_utils.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TF1.h"

////////////////////////////////////////////////////////////////////////////////

TDirectory * myROOT_utils::TDirectory_checked(const char* name)
{
  TDirectory *dir_ok=NULL;
  TDirectory *dir_check;
  dir_check = (TDirectory*) gDirectory->Get(name);
  if(dir_check == NULL)
  {
    dir_ok=gDirectory->mkdir(name);
    if ( NULL==dir_ok )
    {
      assert(false);
    }
  }
  else
  {
    dir_ok = dir_check;
  }
  return dir_ok;
}

////////////////////////////////////////////////////////////////////////////////

TH1D * myROOT_utils::TH1D_checked(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup)
{
  TH1D *h1_ok=NULL;
  TH1D *h1_check;
  h1_check = (TH1D*) gDirectory->Get(name);
  if(h1_check != NULL)
  {
    h1_ok = h1_check;
  }
  else
  {
    h1_ok  = new TH1D(name,title, nbinsx,xlow,xup);
    if( h1_ok == NULL )
    {
      std::cerr << " TH1D_checked memory problems ?? for " << name << " Title " << title << std::endl;
      exit(1);
    }
  }
  return h1_ok;
}

////////////////////////////////////////////////////////////////////////////////

TH2D * myROOT_utils::TH2D_checked(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup, Int_t nbinsy, Axis_t ylow, Axis_t yup)
{
  TH2D *h2_ok=NULL;
  TH2D *h2_check;
  h2_check = (TH2D*) gDirectory->Get(name);
  if(h2_check != NULL)
  {
    h2_ok = h2_check;
  }
  else
  {
    h2_ok  = new TH2D( name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
    if( h2_ok == NULL )
    {
      std::cerr << " TH2D_checked memory problems ?? for " << name << " Title " << title << std::endl;
      exit(1);
    }
  }
  return h2_ok;
}

////////////////////////////////////////////////////////////////////////////////

TProfile * myROOT_utils::TProfile_checked(const char* name, const char* title, Int_t nbinsx, Axis_t xlow, Axis_t xup)
{
  TProfile *p1_ok=NULL;
  TProfile *p1_check;
  p1_check = (TProfile*) gDirectory->Get(name);
  if(p1_check != NULL)
  {
    p1_ok = p1_check;
  }
  else
  {
    p1_ok  = new TProfile(name,title, nbinsx,xlow,xup);
    if( p1_ok == NULL )
    {
      std::cerr << " TProfile_checked memory problems ?? for " << name << " Title " << title << std::endl;
      exit(1);
    }
  }
  return p1_ok;
}

////////////////////////////////////////////////////////////////////////////////

// Macro to reset histograms starting from current directory recrusively

TObject *obj;
TKey    *key;

void myROOT_utils::ResetRecursive(TDirectory* node)
{
  printf(" Enter ResetRecursive node name=%s \n",node->GetName());
//  TDirectory *dirsav;
  //We create an iterator to loop on all objects(keys) of first file
//   node->cd();
  TList* list=node->GetListOfKeys();
  if(list->IsEmpty()) printf(" List is empty \n");
  else printf(" List is NOT empty \n");
  TIter nextkey(list);
  int nobject=0;
  while ( (key = (TKey*)nextkey()) )
  {
    node->cd();
    printf("Object in the loop # %d  \n",nobject);
    obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom("TTree"))
    { //case of a TTree or TNtuple
//VK      TTree* t1 = (TTree*)obj;
      // this part still to be implemented
      // use TChain::???? instead
    }
    else if(obj->IsA()->InheritsFrom("TH1"))
    { //case of TH1 or TProfile
      printf("Found TH1 name=%s title=%s \n",
          obj->GetName(),obj->GetTitle());
      TH1* h1 = (TH1*)obj;
      h1->Reset();
    }
    else if(obj->IsA()->InheritsFrom("TDirectory"))
    { //case of TDirectory
      // recursion
      printf("Found TDirectory name=%s title=%s \n",
          obj->GetName(),obj->GetTitle());
      TObject *objsave = obj;
      TKey    *keysave = key;
      ResetRecursive((TDirectory*)obj);
      obj = objsave;
      key = keysave;
    }
    else
    { //another object
      printf("another obj name=%s, title=%sn",obj->GetName(),obj->GetTitle());
    }
    nobject++;
  }
}

////////////////////////////////////////////////////////////////////////////////
