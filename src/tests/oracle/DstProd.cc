// DST Production User's Coral routines (ORACLE DB)

#include <coral_config.h>
#if !(USE_ORACLE)
#warning("This program works only with Oracle")
#define USE_ORACLE 1
#endif
#undef  __STRICT_ANSI__
#include "CsEvent.h"
#include "CsOraStore.h"

#include "CoralUser.h"

const double M_Pi   = 0.139567;
const double M_P    = 0.93827231;
const double M_N    = 0.93956563;
const double M_K    = 0.493677;
const double M_mu   = 0.105658357;
const double M_D0_  = 1.8645;
const double M_K0_  = 0.497672;
const double M_nucl = (M_P+M_N)/2.;

using namespace std;

bool GetMyExternalVars(double& Q2, double& Y, int& nTracksWithMomentum)
{

  const list<CsTrack*> &Trks = CsEvent::Instance()->getTracks();
  list<CsTrack*>::const_iterator it;
  for( it=Trks.begin(); it!=Trks.end(); it++ ) 
    {
      CsTrack *trk = (*it);
      
      const vector<CsHelix> v = trk->getHelices();
      
      if( v.size()!=0 && v[0].getCop()!=0 ) nTracksWithMomentum++;   // number of tracks with momentum
      
    }

  const vector<CsParticle* > &parts = CsEvent::Instance()->getParticles();
  
  double P_i = 0;
  double P_s = 0;
  double dxdz_i = 0;
  double dxdz_s = 0;
  double dydz_i = 0;
  double dydz_s = 0;
  bool isFound = false;
  for(unsigned i = 0; i < parts.size(); i++)
    {
      if(parts[i]->getName() != "mu") // not scattered muon
	continue;
      const CsTrack* tr = parts[i]->getTrack();
      if(tr)
	{
	  const list<CsVertex*> &v = tr->getVertices();
	  list<CsVertex*>::const_iterator it;
	  for(it = v.begin(); it != v.end(); it++)
	    {
	      if((*it)->isPrimary()) // it is primary vertex
		{
		  isFound = true;
		  const vector<CsHelix>& hlx = tr->getHelices();
		  for(unsigned j = 0; j < hlx.size(); j++)
		    {
		      if(hlx[j].getCop() != 0)
			{
			  P_s = 1./hlx[j].getCop();
			  dxdz_s = hlx[j].getDXDZ();
			  dydz_s = hlx[j].getDYDZ();
			  break;
			}
		    }
		  if(P_s <= 0.) return false;
		  // looking for beam track:
		  const list<CsTrack* > &trks = (*it)->getTracks();
		  list<CsTrack*>::const_iterator itr;
		  for( itr = trks.begin(); itr != trks.end(); itr++ )
		    {
		      if((*itr)->IsBeamTrack())
			{
			  const vector<CsHelix>& hlx = (*itr)->getHelices();
			  for(unsigned j = 0; j < hlx.size(); j++)
			    {
			      if(hlx[j].getCop() != 0)
				{
				  P_i = 1./hlx[j].getCop();
				  dxdz_i = hlx[j].getDXDZ();
				  dydz_i = hlx[j].getDYDZ();
				  break;
				}
			    }
			  if(P_i <= 0.)
			    return false;
			}
		    }
		  break;
		}
	    }
	}
    }

  if(!isFound) return false;

  if(P_i <= P_s)
    return false;

  double E_i = sqrt(P_i*P_i+M_mu*M_mu);
  double E_s = sqrt(P_s*P_s+M_mu*M_mu);

  double pz_i = P_i*sqrt(1+dxdz_i*dxdz_i+dydz_i*dydz_i);
  double px_i = dxdz_i*pz_i;
  double py_i = dydz_i*pz_i;
  double pz_s = P_s*sqrt(1+dxdz_s*dxdz_s+dydz_s*dydz_s);
  double px_s = dxdz_s*pz_s;
  double py_s = dydz_s*pz_s;

  double dE  = E_i-E_s;
  double dPx = px_i-px_s;
  double dPy = py_i-py_s;
  double dPz = pz_i-pz_s;

  Q2 = -(dE*dE-dPx*dPx-dPy*dPy-dPz*dPz);

  if(Q2 < 0) 
    return false;

  Y = dE/E_i;

   return true;
}

bool CoralUserSetup(int argc, char* argv[])
{
  CsOraStore::setExtVarName(0,string("Q_2"));      // Q square
  CsOraStore::setExtVarName(1,string("Y_bjork"));  // Y Bjorken
  CsOraStore::setExtVarName(2,string("Ntr_with_P"));  // Number of tracks with momentum
  return true;
}

void DST_PreuploadFunc()
{

  double Q2 = 0.;
  double Y  = 0.;
  int nTrWithMom  = 0;

  if(!GetMyExternalVars(Q2,Y,nTrWithMom))
    Q2 = Y = -1.;

  CsOraStore::setExtVariable(0, Q2);
  CsOraStore::setExtVariable(1, Y);
  CsOraStore::setExtVariable(2, nTrWithMom);

}

bool CoralUserInit() 
{

  CsOraStore* store = CsOraStore::Instance();
  if(store)
    store->setPreuploadFunc(DST_PreuploadFunc); // the function for filling DST4 ext.vars
  else
    return false;
  return true;
}

bool CoralUserEvent()
{
  return true;
}

bool CoralUserEnd()
{
  return true;
}
