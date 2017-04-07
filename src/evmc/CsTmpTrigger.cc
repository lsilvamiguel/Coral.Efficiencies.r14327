#include "coral_config.h"
#include <cmath>
#include "CsGeant3.h"
#include "CsMCParticle.h"
#include "CsMCTrack.h"
#include "CsTmpTrigger.h"
#include "CsErrLog.h"
#include "CsEvent.h"
#include "CsOpt.h"

using namespace std;

bool CsTmpTrigger::ReadTriggerMap(){
  noTrig=1; //number of notrigered events
  conf=1;
  string triggerFile;
  string tag = ""; 
  string key = "temporary trigger file";
  string str;

  if(CsOpt::Instance()->getOpt( tag, key, str ) ) {
      triggerFile = str; 
 cout<<" It's a good time for a coffe, reading tmp trigger map... "<<endl  
     <<triggerFile.c_str();

   const int lineSize = 256;  
   char   line[lineSize];

   ifstream f(triggerFile.c_str(), ios::in ); //open file temporary here!
   if(!f){
	  string str = "Trigger map file "; 
       	  str.append(triggerFile);
	  str.append(" not opened." );
	  CsErrLog::Instance()->mes( elFatal, str );
          exit(0);
   }

//chossing setup
   string tag1="";
  string key1="temporary trigger configuration";

  if(!CsOpt::Instance()->getOpt( tag1, key1, conf ))
   CsErrLog::mes( elInfo, " --Trigger configuration not specified: all sys assumed" );
  //    cout<<" --  configuration "<<conf<<endl; 
  TrigVar trygger;      
  do {
    f.getline( line, lineSize, '\n' );
    if( f.eof() ) continue;
    istringstream s(line);

    s>>trygger.qq2;
    s>>trygger.y;
    s>>trygger.phi;
    s>>trygger.pp;
    s>>trygger.l;
    s>>trygger.p;
    s>>trygger.unp;    
  
         tryg.push_back(trygger); 
     
  }while( !f.eof() );

  return true;
    } else{
    return false;
    }
}

bool CsTmpTrigger::CheckTmpTrigger(){
 list<CsMCTrack*> tracks=CsEvent::Instance()->getMCTracks();
 list<CsMCTrack*>::iterator track;
   
 int k=0;
 double Px_in=0;
 double Py_in=0; 
 double Pz_in=0;
 double E_in=0;
 double Px_sc=0;
 double Py_sc=0; 
 double Pz_sc=0;
 double E_sc=0;

  for(track=tracks.begin(); track!=tracks.end(); track++){
    if(k==0){    
        Px_in=-(*track)->getPX();
        Py_in=-(*track)->getPY();
        Pz_in=-(*track)->getPZ();  //beam convention
        E_in=(*track)->getE();         
    } else  if(k==1){
        Px_sc=(*track)->getPX();
        Py_sc=(*track)->getPY();
        Pz_sc=(*track)->getPZ();
        E_sc=(*track)->getE();
    } else break; //endif k
    k++;  
  } // MC tracks loop 

  double q2=(E_in-E_sc)*(E_in-E_sc)-(Px_in-Px_sc)*(Px_in-Px_sc)
                                   -(Py_in-Py_sc)*(Py_in-Py_sc) 
                                   -(Pz_in-Pz_sc)*(Pz_in-Pz_sc);
  double Q2=-q2;

  double y=(E_in-E_sc)/E_in;

  if(Px_sc == 0.) Px_sc = 0.00001; //just in case, to prevent infinites
  double phi=atan(fabs(Py_sc/Px_sc));
  if(Py_sc>=0 && Px_sc>0)  phi=phi;
  if(Py_sc<=0 && Px_sc>0)  phi=phi+3*M_PI_2;
  if(Py_sc<=0 && Px_sc<0)  phi=phi+M_PI; 
  if(Py_sc>=0 && Px_sc<0)  phi=phi+M_PI_2;
  
  //  cout<<" Q2 y phi "<<Q2<<" "<<y<<" "<<phi<<endl;

  if(!tryg.size()) {return true;} // if one do not want to use trigger...

  float q2Tol=0.01;
  if(Q2>=0 && Q2<0.2) q2Tol=0.01;
  if(Q2>=0.2 && Q2<0.5) q2Tol=0.03;
  if(Q2>=0.5 && Q2<1.)q2Tol=0.05;
  if(Q2>=1. && Q2<2.) q2Tol=0.1;
  if(Q2>=2. && Q2<5.) q2Tol=0.3;
  if(Q2>=5. && Q2<10.) q2Tol=0.5;
  if(Q2>=10 && Q2<20.) q2Tol=1.;
  if(Q2>=20 && Q2<50.) q2Tol=3.;
  if(Q2>=50.) q2Tol=5.;

  bool trigOk=true;
        list<TrigVar>::iterator it;
      // trigger species
      switch (conf) {
      case 2001:    
       for(it=tryg.begin(); it!=tryg.end(); it++){
        if((it->pp) || (it->p)  || (it->l)){            
          if(fabs((it->qq2)-Q2)<q2Tol && fabs((it->y)-y)<0.01 && 
                                        fabs((it->phi)-phi)<0.063 ){
	    //	   cout<<" trigger map 2001 == event: Q2, y, phi "<<endl;
	    //           cout<<it->qq2<<" == "<<Q2<<" "<<it->y<<" == "<<y<<" "
	    //               <<it->phi<<" == "<<phi<<endl;
           trigOk=true;
           break;  // pp or l or p
          }  
	}
	 trigOk=false;
       }//end of loop over tryg
      break;
      default:  

       for(it=tryg.begin(); it!=tryg.end(); it++){
        if((it->pp) || (it->p)  || (it->l) || (it->unp)){            
          if(fabs((it->qq2)-Q2)<q2Tol && fabs((it->y)-y)<0.01 && 
                                       fabs((it->phi)-phi)<0.063 ){
	    //	    cout<<" trigger map  == event: Q2, y, phi "<<endl;
	    //            cout<<it->qq2<<" == "<<Q2<<" "<<it->y<<" == "<<y<<" "
	    //		<<it->phi<<" == "<<phi<<" tolq2 "<<q2Tol<<endl;
           trigOk=true;
           break;  // any kind of trigger (pp or l or p or unp)
         }  
	}
	 trigOk=false;
       } //end of loop over tryg
	break;
      } //switch
   
      if(!trigOk) cout<<noTrig++<<" notrig event ";
       return trigOk;
}






