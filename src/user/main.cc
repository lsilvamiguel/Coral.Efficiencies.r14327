// $Id: main.cc,v 1.53 2010/09/07 18:27:02 tnagel Exp $

#include "DaqDataDecoding/Exception.h"
#include "CsInit.h"
#include "CsEvent.h"
#include "CsGeom.h"
#include "CsStopwatch.h"
#include "CsRegistrySing.h"
#include "CsErrLog.h"
#include "CoralUser.h"
#include "Coral.h"
#if USE_ObjectsCounter
#  include "ObjectsCounter.h"
#endif

#ifdef TIME_CHECK  
bool TiTrack, TiClust, TiDecod;

double find_proj[5]; // TAlgo::FindProj
int n_hit[5];

//double t_algo1, t_algo2;  //TAlgo::FindSpace

double clu_times[8]; //CsEvent::_clusterize clustering time for given detector
                      //0-BMS;1-MicroMega;2-DriftChamber;3-Strow;4-Muon Wall
		      //5-GEMS; 6-SiFi; 7-MWPC
int    clu_nb[8];

double dec_times, dec1, dec2, dec3, dec4, dec5, dec6, dec7;//CsEvent::_decode
int n_dig;

#endif

// Analysis Coral main program
int main( int argc, char *argv[] ) {

  try
  {  
#ifdef TIME_CHECK 
TiTrack=TiClust=TiDecod=false;
if(char* TiTrack_tmp= getenv("TITRACK")) TiTrack=true;
if(char* TiClust_tmp= getenv("TICLUST")) TiClust=true;
if(char* TiDecod_tmp= getenv("TIDECOD")) TiDecod=true;

//    t_algo1=t_algo2=0;   

    double endofini, endofep, endofcoral;    
    CsStopwatch stopwatch;
    int chrono = stopwatch.start(); // start a chronometer
#endif  
                 
    // Package Initialization 
    ///\todo "Coral init is not needed. But some packages instantiate it. They should not."
    Coral::init( argc, argv ); 
    // this should be preferred...
    //CsInit::Instance( argc, argv );
    CsEvent* event = CsEvent::Instance();
    // User Initialization 
    CoralUserInit();
    
#ifdef TIME_CHECK    
//booking histograms
    string pathname =  "/TIME_DEBUG";
    CsHistograms::SetCurrentPath(pathname);

    CsHist2D* td_1 = new CsHist2D("time1", "FindProj_time vs zone", 5, 0., 5., 250, 0.,0.5);
    CsHist2D* td_2 = new CsHist2D("time2", "FindProj_time vs nb_of_hits", 100, 0., 1000., 250, 0.,0.5);
    CsHist2D* td_3 = new CsHist2D("time3", "nb_of_hits vs zone", 5., 0., 5., 100, 0.,1000.); 
   
    CsHist2D* td_4 = new CsHist2D("time4", "clusterization_time vs detector_type", 8., 0., 8., 100, 0.,0.2);   
    CsHist2D* td_5 = new CsHist2D("time5", "clusterization_time vs nb_of_clusters", 100., 0., 1000., 100, 0.,0.2);        
    CsHist2D* td_6 = new CsHist2D("time6", "nb_of_clusters vs detector_type", 8., 0., 8., 100, 0.,1000);   	
 
    CsHist2D* td_7 = new CsHist2D("time7", "decoding_time vs nb_of_digits", 500., 1000., 6000., 40, 0.,0.4);

    CsHist2D* td_8 = new CsHist2D("time8", "clusterization_time vs nb_of_clusters for BMS", 10., 0., 10., 10, 0.,0.1);
    CsHist2D* td_9 = new CsHist2D("time9", "clusterization_time vs nb_of_clusters for MicroMegas",100., 0., 1000., 20, 0.,0.2);     
    CsHist2D* td_10 = new CsHist2D("time10", "clusterization_time vs nb_of_clusters for DC", 100., 0., 1000., 20, 0.,0.2);   
    CsHist2D* td_11 = new CsHist2D("time11", "clusterization_time vs nb_of_clusters for Strows ", 100., 0., 1000., 20, 0.,0.2);    	
    CsHist2D* td_12 = new CsHist2D("time12", "clusterization_time vs nb_of_clusters for MW ", 100., 0., 1000., 20, 0.,0.2);
    CsHist2D* td_13 = new CsHist2D("time13", "clusterization_time vs nb_of_clusters for GEM-s ", 100., 0., 1000., 20, 0.,0.2); 
    CsHist2D* td_14 = new CsHist2D("time14", "clusterization_time vs nb_of_clusters for SiFi ", 100., 0., 1000., 20, 0.,0.2);   
    CsHist2D* td_15 = new CsHist2D("time15", "clusterization_time vs nb_of_clusters for MWPC ", 100., 0., 1000., 20, 0.,0.2);    
    
#endif 
    int nevt=0;
    
#ifdef TIME_CHECK    
    endofini=stopwatch.stop(chrono);   
    chrono = stopwatch.start(); // start a chronometer   
#endif 

    // Loop on events    
    while( event->getNextEvent() ) { 
       
#ifdef TIME_CHECK
      if(TiTrack){
//+++ Find Projection(Zone*) stuff+++
	for(int i=0; i<4; i++){
	  td_1->Fill((double)i,find_proj[i]); 
	  td_2->Fill((double)n_hit[i],find_proj[i]);
	  td_3->Fill((double)i,(double)n_hit[i]);     
	  //     std::cout<<i<<" "<<find_proj[i]<<std::endl; 
	  find_proj[i]=0;
	  n_hit[i]=0;
	}
      }

      if(TiClust){
	//+++ Clusterization (detector) stuff+++
	for(int i=0; i<8; i++){
	  td_4->Fill((double)i,clu_times[i]);      
	  td_5->Fill((double)clu_nb[i],clu_times[i]);
	  td_6->Fill((double)i,(double)clu_nb[i]);
      
	  switch(i){
	  case 0: {td_8-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  case 1: {td_9-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  case 2: {td_10-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  case 3: {td_11-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  case 4: {td_12-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  case 5: {td_13-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  case 6: {td_14-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  case 7: {td_15-> Fill((double)clu_nb[i],clu_times[i]); break;}
	  default:
	    break;
	  }
	  
	  clu_times[i]=0;
	  clu_nb[i]=0;
	}    
      } 

      if(TiDecod){
	td_7->Fill((double)n_dig,dec_times);
	dec_times=0;
	n_dig=0;   
      }
           
#endif
      nevt++;
      if((nevt)%5 == 0 )
	std::cout << "Event: " << nevt << std::endl;
      CoralUserEvent();
      if( !(CsRegistrySing::Instance()->callEoeMethods()) ) {
	break;
      }
    }
#ifdef TIME_CHECK
    endofep=stopwatch.stop(chrono);
    endofep/=nevt;
    chrono = stopwatch.start(); // start a chronometer    
#endif
    
    // End session
    CoralUserEnd();
    CsRegistrySing::Instance()->callEndMethods();
    CsErrLog::Instance()->dump( elDebugging );
    
#ifdef TIME_CHECK   
    endofcoral=stopwatch.stop(chrono);    

    std::cout<<" Run nb. "<<event->getRunNumber()<<std::endl; 
    std::cout<<" nb. of events "<<nevt<<std::endl;
    std::cout << "Timing ( " << nevt << " events ): \n" 
         << "   Coral Initialization Time:  " << endofini << std::endl
         << "   Event Reading and Building: " << endofep << std::endl
         << "   Coral end  session:         " << endofcoral << std::endl
         << "   Total time:   "<<endofini+endofep*nevt+endofcoral<<std::endl; 


//    std::cout<<"Times for Calling IActive in TAlgo::FindSpace: "<<std::endl
//        <<" in the 5-th loop "<<t_algo1<<" "<<std::endl
//	<<" in the loop over combinations "<<t_algo2<<std::endl; 


if(TiDecod){
    std::cout<<"Times spent in decoding CsEvent: "<<std::endl
        <<" getDaqMaps " <<dec1/nevt<<std::endl
	<<" getDaqOptions "<<dec2/nevt<<std::endl
	<<" ReadChips "<<dec3/nevt<<std::endl
	<<" (*chip)->Decode "<<dec4/nevt<<std::endl
	<<" _decodeTriggerTime "<<dec5/nevt<<std::endl
	<<" decodeChipDigit measured loop over digits "<<dec6/nevt<<std::endl
	<<" Rich decoding "<<dec7/nevt<<std::endl;
}
#endif 	 
	               
  }  
  catch(std::exception &e ) { std::cerr << "Exception:\n" << e.what() << std::endl; }
  catch(std::string &e) { std::cerr << "Exception:\n" << e << "\n"; }
  catch(const char *e) { std::cerr << "Exception:\n" << e << "\n"; }
  catch( ... ) { std::cerr << "Unknown exception!!\n"; }

  //  CS::Exception::PrintStatistics();

#if USE_ObjectsCounter
  ObjectsCounterMaster::SetStream(&std::clog);
  ObjectsCounterMaster::Print();
#endif

  return 0;
}
