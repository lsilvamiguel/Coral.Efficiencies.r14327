#include "GroupPGEM.h"
#include "TThread.h"
#include "GroupPanel.h"
#include "PlanePGEM.h"
#include "TriggerTime.h"
#include "DaqEvent.h"
#include "GeomPlanePGEM.h"

#include <cmath>

ClassImp(GroupPGEM);

void GroupPGEM::Init() {
  std::string name0 = fName + "_U_vs_V";
  fHistList.push_back(new TH2F(name0.c_str(), name0.c_str(), 110, -5.5, 5.5, 110, -5.5, 5.5));
  fHistList[0]->GetXaxis()->SetTitle("V_{strip} / cm");
  fHistList[0]->GetYaxis()->SetTitle("U / cm");
  fHistList[0]->SetOption("COLZ");
  fHistList[0]->SetStats(false);

  std::string name1 = fName + "_Y_vs_X";
  fHistList.push_back(new TH2F(name1.c_str(), name1.c_str(), 110, -5.5, 5.5, 110, -5.5, 5.5));
  fHistList[1]->GetXaxis()->SetTitle("X / cm");
  fHistList[1]->GetYaxis()->SetTitle("Y / cm");
  fHistList[1]->SetOption("COLZ");
  fHistList[1]->SetStats(false);

  std::string name2 = fName + "_U_V_AmplCorr";
  fHistList.push_back(new TH2F(name2.c_str(), name2.c_str(), 50, 0, 1000, 50, 0, 1000));
  fHistList[2]->SetOption("COLZ");

  std::string name3 = fName + "_Y_X_AmplCorr";
  fHistList.push_back(new TH2F(name3.c_str(), name3.c_str(), 50, 0, 1000, 50, 0, 1000));
  fHistList[3]->SetOption("COLZ");

  std::string name4 = fName + "_A2U/A2V";
  fHistList.push_back(new TH1F(name4.c_str(), name4.c_str(), 50, -2.45, 2.55));

  std::string name5 = fName + "_A2Y/A2X";
  fHistList.push_back(new TH1F(name5.c_str(), name5.c_str(), 50, -2.45, 2.55));
}

void GroupPGEM::EndEvent(const CS::DaqEvent &event) {
  if(fPlanes.size()<3) return;

  if (thr_flag) TThread::Lock();

  for (int k=0; k<2; k++) { // k: 0-VU 1-XY
    const GeomPlane *pixgeom = NULL;
    const GeomPlane *xgeom   = NULL;
    const GeomPlane *ygeom   = NULL;

    for (register unsigned int pl=0; pl<fPlanes.size(); pl++) {
      char planeLetter = fPlanes[pl]->GetName()[4];
      if (planeLetter=='P') { // plane is a pixel plane, find the correct one
                              // VU->GP??P1, XY->GP??P2
        char planeNumber = fPlanes[pl]->GetName()[5];
        if (planeNumber=='1' && k==0) pixgeom = fPlanes[pl]->GetGeometry();
        if (planeNumber=='2' && k==1) pixgeom = fPlanes[pl]->GetGeometry();
      }
      if (planeLetter=='V' && k==0) xgeom = fPlanes[pl]->GetGeometry();
      if (planeLetter=='U' && k==0) ygeom = fPlanes[pl]->GetGeometry();
      if (planeLetter=='X' && k==1) xgeom = fPlanes[pl]->GetGeometry();
      if (planeLetter=='Y' && k==1) ygeom = fPlanes[pl]->GetGeometry();
    }

    if (pixgeom && xgeom && ygeom) {
      const std::set<CCluster1>& xclu = xgeom->GetClusters();
      const std::set<CCluster1>& yclu = ygeom->GetClusters();

      int xsize = xclu.size();
      int ysize = yclu.size();
      float xamp2[xsize];
      float yamp2[ysize];
      float xpos[xsize];
      float ypos[ysize];
      bool xused[xsize], yused[ysize];

      int nx=0;
      for (std::set<CCluster1>::iterator xcl=xclu.begin(); xcl!=xclu.end(); xcl++)
        if (xcl->dt.size()>=3) {
          xamp2[nx] = xcl->dt[2];
	      xpos[nx] = xcl->pos;
	      xused[nx] = false;
	      nx++;
        }

      int ny=0;
      for (std::set<CCluster1>::iterator ycl=yclu.begin(); ycl!=yclu.end(); ycl++)
        if (ycl->dt.size()>=3) {
	      yamp2[ny] = ycl->dt[2];
	      ypos[ny] = ycl->pos;
	      yused[ny] = false;
	      ny++;
	    }
      
      int msize = nx*ny;
      int xindex[msize], yindex[msize];
      float distance[msize];
      for (int i=0; i<msize; i++) distance[i]=10000;

      for (int ix=0; ix<nx; ix++)
    	for (int iy=0; iy<ny; iy++) {

	      // Distance to yamp2 = xamp2
    	  float this_dist = fabs(yamp2[iy]-xamp2[ix])/sqrt(2.);
	  
	      // Maximum deviation (this_dist < slope*average+offset)
    	  this_dist = this_dist - 0.05*(xamp2[ix]+yamp2[iy])/2.;
	    
	      // Sort according to this_dist
    	  for (int i=0; i<msize; i++) {
	        if (this_dist<distance[i]) {
	          for (int ii=msize-2; ii>=i; ii--) {
                distance[ii+1] = distance[ii];
                xindex[ii+1] = xindex[ii];
                yindex[ii+1] = yindex[ii];
	          }
	          distance[i]=this_dist;
	          xindex[i]=ix;
	          yindex[i]=iy;
	          break;
	        }
	      }
	    }
      
      // Associate clusters accordingly
      for (int i=0; i<msize; i++) {
        // Check if within road and intersection point not in pixels
        if ( distance[i]<80. && !(fabs(xpos[xindex[i]])<pixgeom->GetDY()/2. && fabs(ypos[yindex[i]])<pixgeom->GetDZ()/2.) ) {
	  
    	  // Clusters already associated
	      if (!xused[xindex[i]] && !yused[yindex[i]]) {
	        fHistList[k  ]->Fill(xpos[xindex[i]], ypos[yindex[i]]);
    	    fHistList[k+2]->Fill(xamp2[xindex[i]], yamp2[yindex[i]]);
	        fHistList[k+4]->Fill(yamp2[yindex[i]]/xamp2[xindex[i]]-1.);
	        xused[xindex[i]]=true;
    	    yused[yindex[i]]=true;
	      }
	    } else 
          break;
      }

      // fill the position of the pixel clusters
      std::set<CCluster1>::iterator pcluit;
      for (pcluit=pixgeom->GetClusters().begin(); pcluit!=pixgeom->GetClusters().end(); pcluit++) {
        if (pcluit->dt.size()>=5) {
          double posX=pcluit->er[3];
          double posY=pcluit->er[4];
          if (pixgeom->GetName()[5] == '1') { // UV detector => (strips:v->x, u->y) pixels: v->-x, u->y
            double temp=posX;
            posX = -posY;
            posY = temp;
          }
          fHistList[k]->Fill(posX, posY);
        }
      }
    }
  }
  
  if (thr_flag) TThread::UnLock();
}

void GroupPGEM::ControlPanel(const TGWindow *p, const TGWindow *main) {
  if (!fControlPanel) fControlPanel = new GroupPanel(p, main, 100, 100, this);
}





