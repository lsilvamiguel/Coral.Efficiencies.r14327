// $Id: TEvDumpTrackList.cc 13148 2011-12-28 16:55:25Z kbicker $

#include <iostream>
#include <cstdio>
#include <CsSTD.h>
#include "Coral.h"
#include "TEv.h"

  /*!
    \file    TEvDumpTrackList.cc
    \brief   Dump for debug purpose
    \author  Y.B
  */

using namespace std;

void  TEv::DumpTrackList()
{
  //********** RECONSTRUCTED TRACKS

  list<TTrack>::iterator it;
  for (it = listTrack.begin(); it!=listTrack.end(); it++) {
    TTrack& t = (*it);
    t.DumpHits(1);
  }
}













