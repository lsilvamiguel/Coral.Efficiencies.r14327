#include "CsEvent.h"
#include "CsPPI_EC02time.h"
#include "CsOpt.h"
#include "CsInit.h"
#include "DaqDataDecoding/DaqEvent.h"
#include "DaqDataDecoding/ChipSADC.h"
#include <sstream>
#include <cstdio>
#include <sstream>
#include <cassert>
#include <stdlib.h>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <fstream>
#include <sys/mman.h>

using namespace std;
using CS::DetID;

/****************************
 ** Error code definitions **
 ****************************/

#define INVALID_CALL          -10
#define NOT_MAPPED            -11

#define UNKNOWN_EXCEPT        -100


//                            -3xx for db errors
#define DB_CONN_FAILED        -300
#define NO_DB_HOST            -301
#define NO_DB_DB              -302
#define NO_DB_USER            -303
#define NO_DB_PASS            -304
#define NO_DB_TABLE           -305
#define DB_QUERY_FAILED       -306
#define NO_DB_OUTPUT_TABLE    -307
//                            -4xx for cfd errors
#define NO_CFD_THR            -400
#define NO_CFD_DELAY          -401
#define NO_CFD_AMPL           -402
#define NO_HIST_WIDTH         -403
#define NO_CFD_HITS           -404
#define MULTI_CFD_HITS        -405
#define TO_BIG_TIMEDIFF       -406
#define NO_BASELINE           -407
#define WRONG_XY_COORDINATE   -408
#define TO_HIGH_PORT_NMB      -409
#define TO_HIGH_CHANNEL_NMB   -410
#define NO_MIN_SRCID          -411


const float CsPPI_EC02time::cycle_length = 12.86;


/****************************
 ** CsPPI_EC02time section **
 ****************************/


  CsPPI_EC02time::CsPPI_EC02time(const DetID& id, const std::string& TBname) 
:  CsPPI(id,TBname), fTime(NULL), fCTime(NULL), fConn(NULL), fOutConn(NULL),
  run_nmb(-9999), spill_nmb(-9999), 
  nentries(200*6*8*8),filesize((nentries+41)*sizeof(ec02_tcorr_data_t))
{

  testcount=0; //used to abort after N events, only for debugging

    const string cfd_tag = "CsPPI_EC02time_TIMEpar";

      if( !CsOpt::Instance()->getOpt( cfd_tag, "output_folder",fdir)) {
       MSG( WARNING, NO_CFD_THR ,
          "No Output folder  specified, switching to database mode", 
          __FILE__, __LINE__);

       output_mode=database; 
    }
    else
      output_mode=file;

    // get cfd_threshold
    if ( !CsOpt::Instance()->getOpt( cfd_tag , "cfd_threshold" , settings.cfd_threshold ) ) {
      MSG( WARNING, NO_CFD_THR ,
          "No CFD threshold  specified, setting to 20", 
          __FILE__, __LINE__);
      settings.cfd_threshold=20;
    }
    // get cfd_delay
    if ( !CsOpt::Instance()->getOpt( cfd_tag , "cfd_delay" , settings.cfd_delay ) ) {
      MSG( WARNING, NO_CFD_DELAY ,
          "No CFD delay  specified, setting to 2", 
          __FILE__, __LINE__);
      settings.cfd_delay=2;
    }
    // get cfd_amplification
    if ( !CsOpt::Instance()->getOpt( cfd_tag , "cfd_ampl" , settings.cfd_ampl ) ) {
      MSG( WARNING, NO_CFD_DELAY ,
          "No CFD amplification  specified, setting to 2", 
          __FILE__, __LINE__);
      settings.cfd_ampl=2;
    }
    // get baseline
    if ( !CsOpt::Instance()->getOpt( cfd_tag , "baseline" , settings.baseline ) ) {
      MSG( WARNING, NO_BASELINE ,
          "No baseline  specified, setting to 50", 
          __FILE__, __LINE__);
      settings.baseline=50;
    }

    // get width of the array containing the timing values
    if ( !CsOpt::Instance()->getOpt( cfd_tag , "time_hist_width" , settings.time_hist_width ) ) {
      MSG( FATAL, NO_HIST_WIDTH ,
          "No width of the timing histogramm specified", 
          __FILE__, __LINE__);
      ClearMSG(0);
      exit(NO_HIST_WIDTH);
    }


    if(output_mode == database) {
    if ( !CsOpt::Instance()->getOpt( "CsPPI_EC02time_DB" , 
          "maxcombined_queries" , 
          maxcombined ) ) {
      MSG( WARNING, 0 ,
          "No number for max combined queries given: using 1", 
          __FILE__, __LINE__);
      maxcombined = 1;
    }



    // Start the DBwriter thread
    fDBThreadState = THREAD_RUN;
    pthread_create(&fDBThread, NULL, 
        start_db_thread, (void*)this);
    }


    vector<int> empty1;
    vector<vector<int> > empty2;

    for(int i=0;i<settings.time_hist_width;i++)
      empty1.push_back(0);
    for(int i=0;i<8;i++)
      empty2.push_back(empty1);
    for(int i=0;i<8;i++)
      empty_time_hist.push_back(empty2);




  if ( int status=OpenConn() ) {
    ClearMSG(0);
    exit(status);
  }
  // open db connection if it is not open by now
  //   if ( fOutConn==NULL ) {
  //     if ( int status=OpenOutConn() ) {
  //       ClearMSG(0);
  //     exit(status);
  //     }
  //   }
  // 

};

/****************
 ** Destructor **
 ****************/

CsPPI_EC02time::~CsPPI_EC02time() 
{

  //write out the time of the last processed spill
  FillGlobalHistogram();

  if( output_mode == database) {
  write_out_time2db();
  // join thread this may take until all db queries are executed
  cout<<"Writing ECAL timing corrections to database, please be patient.\n";
  pthread_mutex_lock(&fDBmutex);
  fDBThreadState = THREAD_TERM;
  pthread_mutex_unlock(&fDBmutex);
  void *status;
  pthread_join(fDBThread, &status);

  ClearMSG(1);
  cout<<"Writing ECAL timing corrections to database, finished.\n";
  }
  else {
    write_out_time2file();
  }
}

int
CsPPI_EC02time::end()  {
  this->~CsPPI_EC02time();
  return 0;

}

/*********************************************************************************
 ** OpenOutConn() opens the connection to the mysql server and return 0 on success **  
 *********************************************************************************/
int
CsPPI_EC02time::OpenOutConn() {

  const string tag = "CsPPI_EC02time_DB";
  string host, db, user, pass;

  // Only one connection allowed
  if ( fOutConn!=NULL ) return 0;

  // get host
  if ( !CsOpt::Instance()->getOpt( tag , "output_host" , host ) ) {
    MSG( FATAL, NO_DB_HOST ,
        "No out db host specified", 
        __FILE__, __LINE__);
    ClearMSG(0);
    exit( NO_DB_HOST);
  }
  // get db
  if ( !CsOpt::Instance()->getOpt( tag , "output_db" , db ) ) {
    MSG( FATAL, NO_DB_DB ,
        "No out db specified", 
        __FILE__, __LINE__);
    ClearMSG(0);
    exit(NO_DB_DB);
  }  
  // get user
  if ( !CsOpt::Instance()->getOpt( tag , "output_user" , user ) ) {
    MSG( FATAL, NO_DB_USER ,
        "No out db user specified", 
        __FILE__, __LINE__);
    ClearMSG(0);
    exit( NO_DB_USER);
  }
  // get password
  if ( !CsOpt::Instance()->getOpt( tag , "output_pass" , pass ) ) {
    pass = "";
    MSG( ERROR, NO_DB_PASS ,
        "No out db password specified using \"\" ", 
        __FILE__, __LINE__);
  }

  // open the connection
  try {
    fOutConn = new MySQLInterface(host.c_str(), user.c_str(), pass.c_str(), db.c_str());
  } catch ( const char* e ) {
    char *msg;
    asprintf( &msg, 
        "Catch mysqlpp::Exception while opening mysqlpp::Connection"
        "to '%s'@'%s' using '%s' with '%s': %s",
        db.c_str(), host.c_str(), user.c_str(), pass.c_str(), e);
    MSG( FATAL, DB_CONN_FAILED, msg, __FILE__,__LINE__);
    free(msg);
    return DB_CONN_FAILED;
  } catch  ( ... ) {
    char *msg;
    asprintf( &msg, 
        "Catch unknown exception while opening mysqlpp::Connection"
        "to '%s'@'%s' using '%s' with '%s'!",
        db.c_str(), host.c_str(), user.c_str(), pass.c_str());
    MSG( FATAL, UNKNOWN_EXCEPT, msg, __FILE__,__LINE__);
    free(msg);
    return  UNKNOWN_EXCEPT;
  }
  ostringstream ost;
  ost << "USE " << db; 
  fOutConn->query(ost.str());

  return 0;

};



/*********************************************************************************
 ** OpenConn() opens the connection to the mysql server and return 0 on success **  
 *********************************************************************************/
int
CsPPI_EC02time::OpenConn() {

  const string tag = "CsPPI_EC02time_DB";
  string host, db, user, pass;

  // Only one connection allowed
  if ( fConn!=NULL ) return 0;

  // get host
  if ( !CsOpt::Instance()->getOpt( tag , "host" , host ) ) {
    MSG( FATAL, NO_DB_HOST ,
        "No db host specified", 
        __FILE__, __LINE__);
    ClearMSG(0);
    exit( NO_DB_HOST);
  }
  // get db
  if ( !CsOpt::Instance()->getOpt( tag , "db" , db ) ) {
    MSG( FATAL, NO_DB_DB ,
        "No db specified", 
        __FILE__, __LINE__);
    ClearMSG(0);
    exit(NO_DB_DB);
  }  
  // get user
  if ( !CsOpt::Instance()->getOpt( tag , "user" , user ) ) {
    MSG( FATAL, NO_DB_USER ,
        "No db user specified", 
        __FILE__, __LINE__);
    ClearMSG(0);
    exit( NO_DB_USER);
  }
  // get password
  if ( !CsOpt::Instance()->getOpt( tag , "pass" , pass ) ) {
    pass = "";
    MSG( ERROR, NO_DB_PASS ,
        "No db password specified using \"\" ", 
        __FILE__, __LINE__);
  }
  if(output_mode==database) {
    // get time correction output table name
    if ( !CsOpt::Instance()->getOpt( tag , "corr_tbl" , corr_tbl ) ) {
      MSG( FATAL, NO_DB_OUTPUT_TABLE ,
          "No table to write out the timing correction parameters specified ", 
          __FILE__, __LINE__);
      ClearMSG(0);
      exit(NO_DB_OUTPUT_TABLE);
    }
  }


  // open the connection
  try {
    fConn = new MySQLInterface(host.c_str(), user.c_str(), pass.c_str(), db.c_str());
  } catch ( const char* e ) {
    char *msg;
    asprintf( &msg, 
        "Catch mysqlpp::Exception while opening mysqlpp::Connection"
        "to '%s'@'%s' using '%s' with '%s': %s",
        db.c_str(), host.c_str(), user.c_str(), pass.c_str(), e);
    MSG( FATAL, DB_CONN_FAILED, msg, __FILE__,__LINE__);
    free(msg);
    return DB_CONN_FAILED;
  } catch  ( ... ) {
    char *msg;
    asprintf( &msg, 
        "Catch unknown exception while opening mysqlpp::Connection"
        "to '%s'@'%s' using '%s' with '%s'!",
        db.c_str(), host.c_str(), user.c_str(), pass.c_str());
    MSG( FATAL, UNKNOWN_EXCEPT, msg, __FILE__,__LINE__);
    free(msg);
    return  UNKNOWN_EXCEPT;
  }
  ostringstream ost;
  ost << "USE " << db; 

  fConn->query(ost.str());

  return 0;

};


/***************************************************
 ** Retrive time calibration for srcid from DB    **
 ** Return error code (negative) or number of     **
 ** time calibrations                             **
 ***************************************************/

int
CsPPI_EC02time::DBGetTCalib(const unsigned int srcid) {

  // open db connection if it is not open by now
  if ( fConn==NULL ) {
    if ( int status=OpenConn() ) {
      ClearMSG(0);
      exit(status);
    }
  }

  if ( !fCTime )
    convertTime();

  string table;

  // get table
  if ( !CsOpt::Instance()->getOpt( "CsPPI_EC02time_DB" , 
        "table_tcalib" , 
        table ) ) {
    MSG( FATAL, NO_DB_TABLE ,
        "No table for time calibration specified!", 
        __FILE__, __LINE__);
    ClearMSG(0);
    exit(NO_DB_TABLE);
  }

  ostringstream ost;
  ost << "SELECT port,ch,mean FROM " << table 
    << " WHERE srcid=" << srcid
    << " AND valid_from<'" << fCTime 
    << "' AND (valid_to>'" << fCTime 
    << "' OR ISNULL(valid_to) )";


  if (fConn->query(ost.str())) {
    float time_corr;

    if ( !CsOpt::Instance()->getOpt( "CsPPI_EC02time" , 
          "time_corr" , 
          time_corr ) )
      time_corr = 0;

    const unsigned int port = atoi(fConn->getCol(0));
    const unsigned int ch   = atoi(fConn->getCol(1));
    const float mean = atof(fConn->getCol(2)) + time_corr;

    SetTCalib(srcid, port, ch, mean);

    while (fConn->getNextRow()) {
      const unsigned int port = atoi(fConn->getCol(0));
      const unsigned int ch   = atoi(fConn->getCol(1));
      const float mean = atof(fConn->getCol(2)) + time_corr;

      SetTCalib(srcid, port, ch, mean);
    }
  }

  fConn->endQuery();

  return 0;


}


/**********************************************
 ** process function for the interface class **
 **********************************************/

int 
CsPPI_EC02time::process() {
  CsEvent* event = CsEvent::Instance();
  // on change of spill writeout calculated times
  if( (unsigned int)spill_nmb!=event->getBurstNumber() ||(unsigned int)run_nmb!=event->getRunNumber() ){
    
    //first time process is called writeout is skipped
    if(spill_nmb!=-9999&& run_nmb!=-9999){
      FillGlobalHistogram();

    }
    spill_nmb=event->getBurstNumber();
    run_nmb=event->getRunNumber();
    time_hist.clear();
  }
  //is it a physics event?
  if(event->getDaqEvent().GetType()!=7)
    return -1;

  testcount++;
  const CS::ChipSADC::Digits& Digs = event->getChipDigits();
  CS::ChipSADC::Digits::const_iterator it;
  //   float cfd_time;
  vector<float> cfd_time_v;

  for(it = Digs.begin(); it != Digs.end(); ++it) { 
//     if( it->first.GetName() != "EC02P1__" )
//      if( it->first.GetName() == "EC02P1__" )
     if( it->first != GetID() )
      continue;
    CS::ChipSADC::Digit* chd=(CS::ChipSADC::Digit*)(it->second);


    assert(chd != NULL);
    calc_time((vector<float>)(chd->GetNtupleData()),&cfd_time_v);
//     for(unsigned int i=0;i<chd->GetNtupleData().size();i++)


    GetTCalib(((CS::ChipSADC::DataID&)chd->GetDataID()).GetSourceID(),
        ((CS::ChipSADC::DataID&)chd->GetDataID()).GetPort(),
        ((CS::ChipSADC::DataID&)chd->GetDataID()).GetChannel()+(16*(((CS::ChipSADC::DataID&)chd->GetDataID()).GetChip())),
        &time_reference);

    for(unsigned int i=0; i<cfd_time_v.size();i++){

      if(cfd_time_v.at(i)>0){
        int timediff  = (int)(
            0.5+time_reference/cycle_length
            + settings.time_hist_width
            -(cfd_time_v.at(i) +  ((event->getTCSPhaseTime()) / cycle_length)  )
            );

        if( timediff >= 0 && timediff < settings.time_hist_width ) {
          increment_time_bin((CS::ChipSADC::DataID&)chd->GetDataID(),(unsigned int)timediff);
        }
        //         else {
        // //           MSG(INFO,TO_BIG_TIMEDIFF,
        //               "Called CsPPI_EC02time::process(), negativ time value.",
        //               __FILE__,__LINE__);
        //         }
      }
    }
    cfd_time_v.clear();

  }
  return 0;
}




/**************************************************
 ** calculate the time of the sample using a CFD **
 *************************************************/

float 
CsPPI_EC02time::calc_time(const vector<float>& data,vector<float>* cfd_time_v) {
  double data_delayed=0;
  double data_orig=0;
  double difference[2];
  int cfd_hit=false;
  double ltime=0;
  difference[0]=0;
  difference[1]=0;
  double diff=0;
  double diff_delayed=0;

  /*********************************************************************/
float  basl1=data.at(2);
float  basl2=data.at(3);
  for( unsigned int i=2+settings.cfd_delay; i<data.size();i++) {
    if(i%2){
    data_delayed=((data.at(i-settings.cfd_delay))-basl2)*settings.cfd_ampl;
    data_orig=(data.at(i))-basl2;
    }
    else {

    data_delayed=((data.at(i-settings.cfd_delay))-basl1)*settings.cfd_ampl;
    data_orig=(data.at(i))-basl1;
    }

    diff_delayed=diff;
    diff= data_orig-data_delayed;
    if(i>3) {
      if(diff_delayed -diff >settings.cfd_threshold &&diff<0 && diff_delayed>=0) {
        cfd_time_v->push_back((i-5) + ((float)(diff_delayed)/(diff_delayed-diff)));
        cfd_hit=true;
      }

    }
  }

  /*********************************************************************/


  if(cfd_hit==true) {
    return ltime;
  }
  else {
    //     MSG(INFO,NO_CFD_HITS,
    //         "Called CsPPI_EC02time::calc_time(),  sample returned no timing value.",
    //         __FILE__,__LINE__);
    return NO_CFD_HITS;
  }
};

/***************************************************
 ** creates the query to write the time to a file **
 **************************************************/
int
CsPPI_EC02time::write_out_time2file(void) {
  unsigned int                             Port;
  float                                    tmp_sum;
  int                                      maxbin;
  float                                    qualty_ratio;
  unsigned int                             checksum=0;
  unsigned int                             fname_counter=0;
  vector <vector <vector<int> > >          *port_vec;

   map< unsigned int,
    map< unsigned int,
    map< unsigned int, vector <vector <vector<int> > > > > >::iterator Run;
  map< unsigned int,
    map< unsigned int, vector <vector <vector<int> > > > >::iterator Spill;
  map< unsigned int, vector <vector <vector<int> > > >::iterator SrcID;
 char*                                     pre_fname_tmp;



  for(Run=fGlobalTimeHist.begin();Run!=fGlobalTimeHist.end();Run++) {
    fname_counter=0;
    run_nmb=Run->first;
//    cout<< *CsInit::Instance()->getDateFilesList().front()<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
   size_t found;
    found=CsInit::Instance()->getDateFilesList().front()->find("cdr");
//    cout<< CsInit::Instance()->getDateFilesList().front()->substr(found,8)<<"lkjkljlkjkjljljjlk\n";





    asprintf(&pre_fname_tmp,"%s/EC02_timecorr_run%u_%s",fdir.c_str(),Run->first,CsInit::Instance()->getDateFilesList().front()->substr(found,8).c_str());
    //check if file exists
    //----------------------------------------
    {
      bool flag = false;
      do {
        fstream fin;
        asprintf(&pre_fname,"%s_%u",pre_fname_tmp,fname_counter);
        fin.open(pre_fname,ios::in);
        if( fin.is_open() )
        {
          fname_counter++;
          flag=true;
          free(pre_fname);
        }
        else {
          flag=false;
        }
        fin.close();
      } while(flag==true);
      free(pre_fname_tmp);
    }
    //----------------------------------------
    fname=pre_fname;
    OpenOutputFile(fname);
    free(pre_fname);
    for(Spill=(Run->second).begin();Spill!=(Run->second).end();Spill++) {
      checksum=0;
      //----------------------------------------
      //loop  over srcids
      for(SrcID=(Spill->second).begin();SrcID!=(Spill->second).end();SrcID++)
      {
        port_vec=&(SrcID->second);
        //----------------------------------------
        //loop  over ports
        for(Port=0;Port < port_vec->size();Port++)
        {
          //----------------------------------------
          //loop over channel-groups 
          for( unsigned int chgrp=0;chgrp<(port_vec->at(Port).size());chgrp++)
          {

            //find bin with the most entries and determine the ratio between maxbin and all bins

            unsigned int tempbin[5];
            for(int ii=0;ii<5;ii++) {
              tempbin[ii]=port_vec->at(Port).at(chgrp).at(ii);
              checksum+=tempbin[ii];
            }
            Write2OutputFile(((Spill->first)-1)*384+64*(SrcID->first-616)+8*Port+chgrp,tempbin);
          }
          // end loop over channel groups
          //----------------------------------------
        }
        // end loop over ports
        //----------------------------------------

      }
      // end loop over srcids
      //----------------------------------------
      WriteChecksum(Spill->first-1,checksum);
    }
    // end loop over spills
    //----------------------------------------
    CloseOutputFile();
  }
  // end loop over runs
  //----------------------------------------
  return 0;
}


/***************************************************
 ** creates the query to write the time to the db **
 **************************************************/
int
CsPPI_EC02time::write_out_time2db(void) {
  char * query = NULL;
  char * query_vals = NULL;
  bool first_chan=true;
  unsigned int Port;
  float tmp_sum;
  int  maxbin;
  float qualty_ratio;
  unsigned int checksum=0;
  map< unsigned int,
    map< unsigned int,
    map< unsigned int, vector <vector <vector<int> > > > > >::iterator Run;
  map< unsigned int,
    map< unsigned int, vector <vector <vector<int> > > > >::iterator Spill;
  map< unsigned int, vector <vector <vector<int> > > >::iterator SrcID;
  vector <vector <vector<int> > >  *port_vec;
  for(Run=fGlobalTimeHist.begin();Run!=fGlobalTimeHist.end();Run++) {
    for(Spill=(Run->second).begin();Spill!=(Run->second).end();Spill++) {
      //----------------------------------------
      //loop  over srcids
      for(SrcID=(Spill->second).begin();SrcID!=(Spill->second).end();SrcID++)
      {
        port_vec=&(SrcID->second);
        //----------------------------------------
        //loop  over ports
        first_chan=true;
        for(Port=0;Port < port_vec->size();Port++)
        {
          //----------------------------------------
          //loop over channel-groups 
          for( unsigned int chgrp=0;chgrp<(port_vec->at(Port).size());chgrp++)
          {

            //find bin with the most entries and determine the ratio between maxbin and all bins
            maxbin=0;
            qualty_ratio=0;
            tmp_sum=0;
            for(int bin=0;bin<settings.time_hist_width;bin++) { 
              if( port_vec->at(Port).at(chgrp).at(maxbin) < port_vec->at(Port).at(chgrp).at(bin)  ) {
                maxbin=bin;
              }
              tmp_sum+=port_vec->at(Port).at(chgrp).at(bin);
            }

            if(tmp_sum>0)
              qualty_ratio=((  (float)(port_vec->at(Port).at(chgrp).at(maxbin))/tmp_sum ));
            else {
              if(port_vec->at(Port).at(chgrp).at(maxbin) == 0)
                qualty_ratio=0;
              else if( port_vec->at(Port).at(chgrp).at(maxbin) >100)
                qualty_ratio=100;
              else
                qualty_ratio=50;
            }

            maxbin-=(settings.time_hist_width>>1);
            unsigned int tempbin[5];
            for(int ii=0;ii<5;ii++) {
              tempbin[ii]=port_vec->at(Port).at(chgrp).at(ii);
              checksum+=tempbin[ii];
            }

            if(first_chan==true) {
              //                     rn,sn,sr,po,ch,b0,b1,b2,b3,b4,bm
              asprintf(&query_vals,"(%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%i,%f)",
                  Run->first,Spill->first,
                  SrcID->first,Port,chgrp,
                  port_vec->at(Port).at(chgrp).at(0),port_vec->at(Port).at(chgrp).at(1),
                  port_vec->at(Port).at(chgrp).at(2),port_vec->at(Port).at(chgrp).at(3),
                  port_vec->at(Port).at(chgrp).at(4),(int)maxbin,qualty_ratio);
              first_chan=false;
              continue;
            }
            //                        rn,sn,sr,po,ch,b0,b1,b2,b3,b4,bm
            asprintf(&query_vals,"%s,(%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,%i,%f)",query_vals,
                Run->first,Spill->first,
                SrcID->first,Port,chgrp,
                port_vec->at(Port).at(chgrp).at(0),port_vec->at(Port).at(chgrp).at(1),
                port_vec->at(Port).at(chgrp).at(2),port_vec->at(Port).at(chgrp).at(3),
                port_vec->at(Port).at(chgrp).at(4),(int)maxbin,qualty_ratio);
                      }
          // end loop over channel groups
          //----------------------------------------
        }
        // end loop over ports
        //----------------------------------------
        asprintf(&query,"INSERT INTO %s (run_nmb,spill,srcid,port,ch_grp,bin0,bin1,bin2,bin3,bin4,shift,qual) VALUES %s",corr_tbl.c_str(),query_vals);
        Write2DB(query);
        free(query_vals);

      }
      // end loop over srcidss
      //----------------------------------------
    }
    // end loop over spills
    //----------------------------------------

  }
  // end loop over runs
  //----------------------------------------
  return 0;
}


/*****************************************************************
 ** GetTCalib try to find the time calibration for the channel, **
 ** which is identified using (srcid,port,ch). If it returns    **
 ** 0 the value of *tcalib is set to the calibration. Otherwise **
 ** an error code is returned.                                  **
 *****************************************************************/

int                     
CsPPI_EC02time::GetTCalib(const unsigned int srcid, 
    const unsigned int port, 
    const unsigned int ch, 
    float* tcalib) {

  // No NULL allowed for tcalib
  if ( !tcalib ) { 
    MSG( FATAL, INVALID_CALL, 
        "While retieving time calibration: tcalib = NULL. Fix your code!", 
        __FILE__, __LINE__);
    exit( INVALID_CALL);
  }

  // Get calibrations for srcid if not initialized yet
  if ( !fInitSrcId.count(srcid) ) {
    InitSrcId(srcid);
  }

  // Find the port map (match srcid)
  map< unsigned int, 
    map< unsigned int, 
    map< unsigned int, float> > >::iterator SrcID 
      = fTCalib.find(srcid);
  if ( SrcID == fTCalib.end() ) {
    char *msg;
    asprintf( &msg, 
        "While retieving time calibration: srcid %u not mapped!", 
        srcid );
    MSG( ERROR, NOT_MAPPED, msg, __FILE__, __LINE__);
    free(msg);
    return NOT_MAPPED; 
  }

  // Find channel map (match port)
  map< unsigned int, 
    map< unsigned int, float> >::iterator Port = 
      (SrcID->second).find(port);
  if (  Port == (SrcID->second).end() ) {
    char *msg;
    asprintf( &msg, 
        "While retieving time calibration: "
        "for srcid %u port %u not mapped!", 
        srcid, port );
    MSG( ERROR, NOT_MAPPED, msg, __FILE__, __LINE__);
    free(msg);
    return NOT_MAPPED;
  }

  // Time calib (match ch)
  map< unsigned int, float>::iterator Ch 
    = (Port->second).find(ch);
  if (  Ch == (Port->second).end() ) {
    char *msg;
    asprintf( &msg, "While retieving time calibration: "
        "for srcid %u, port %u ch %u not mapped!", 
        srcid, port, ch );
    MSG( ERROR, NOT_MAPPED, msg, __FILE__, __LINE__);
    free(msg);
    return NOT_MAPPED;
  }

  (*tcalib) = Ch->second;

  return 0;

};

/****************************************************
 ** Make the initialization of srcid (if not done) **
 ****************************************************/

int
CsPPI_EC02time::InitSrcId(const unsigned int srcid) {

  if ( fInitSrcId.count(srcid) ) return 0;

  fInitSrcId.insert(srcid);

  if ( int status = DBGetTCalib(srcid) < 0 ) {
    return status;
  }

  return 0;

};

/************************************************************
 ** Set the current event time for retrieving calibrations **
 ************************************************************/

void
CsPPI_EC02time::getTime() {


  if ( fTime ) return;

  unsigned int epoch = CsEvent::Instance()->getEventTime().secFrEpoch();

  fTime = new time_t( *(time_t*)(&epoch) );

  return;

};

/****************************************************
 ** Convert the fTime to a char string for queries **
 ****************************************************/
void
CsPPI_EC02time::convertTime() {

  if ( fCTime ) return;

  if ( !fTime ) getTime();

  fCTime = new char[20];

  strftime(fCTime, 20, "%Y-%m-%d %H:%M:%S", gmtime(fTime)  );

  return;


};

/***************************************************************
 ** SetTCalib sets the time calibration value of the channel, **
 ** which is identified by (srcid,port,ch) to tcalib.         ** 
 ***************************************************************/
void                    
CsPPI_EC02time::SetTCalib(unsigned int srcid, unsigned int port, 
    unsigned int ch, float tcalib) {

  // Find the port map (match srcid)
  map< unsigned int, 
  map< unsigned int, 
  map< unsigned int, float> > >::iterator SrcID = fTCalib.find(srcid);
  if ( SrcID == fTCalib.end() ) {
    map< unsigned int, map< unsigned int, float> > Port;
    map< unsigned int, float> Ch;
    fTCalib[srcid] = Port;
    fTCalib[srcid][port] = Ch;
    fTCalib[srcid][port][ch] = tcalib;
    return;
  }

  // Find the ch map (match srcid)
  map< unsigned int, map< unsigned int, float> >::iterator Port = 
    (SrcID->second).find(port);
  if (  Port == (SrcID->second).end() ) {
    map< unsigned int, float> Ch;
    fTCalib[srcid][port] = Ch;
    fTCalib[srcid][port][ch] = tcalib;
    return;
  }

  // Does not matter if channel exist
  fTCalib[srcid][port][ch] = tcalib;

  return;

};

void
CsPPI_EC02time::increment_time_bin(CS::ChipSADC::DataID& id ,unsigned int bin) {
  unsigned int srcid=id.GetSourceID();
  unsigned int port=id.GetPort();
  unsigned int ch=id.GetChannel()+(16*id.GetChip());
  //   unsigned int ch=(id.GetChannel()*(id.GetChip()+1));


  // Find the port map (match srcid)
  map< unsigned int, vector<vector<vector<int> > > >::iterator SrcID = time_hist.find(srcid);     
  if ( SrcID == time_hist.end() ) {
    time_hist[srcid]= empty_time_hist;
  }

  //increment coresponding bin
  if(port>7) {
    char *msg;
    asprintf( &msg, 
        "CsPPI_EC02time::increment_time_bin, port number to high: "
        "%i > 7",port );
    MSG( WARNING, TO_HIGH_PORT_NMB, msg, __FILE__,__LINE__);
    free(msg);

  }
  else if(ch>63) {
    char *msg;
    asprintf( &msg, 
        "CsPPI_EC02time::increment_time_bin, channel number to high: "
        "%i > 63",ch );
    MSG( WARNING,TO_HIGH_CHANNEL_NMB, msg, __FILE__,__LINE__);
    free(msg);
  }
  else {

    time_hist[srcid].at(port).at(ch/8).at(bin)++;
  }

  return;

}

/********************************************************
 ** Write2DB will enque a query to be submitted to DB. **
 ********************************************************/

void
CsPPI_EC02time::Write2DB(char * query) {

  if ( fOutConn==NULL ) {
    if ( int status=OpenOutConn() ) {
      ClearMSG(0);
      exit(status);
    }
  }


  pthread_mutex_lock(&fDBmutex);
  fQuery.push_back(query);
  pthread_mutex_unlock(&fDBmutex);
  return;

}

/***************************
 ** Thread to write to db **
 ***************************/
void* start_db_thread(void* arg) {

  pthread_exit(((CsPPI_EC02time*)arg)->DB_Thread(NULL));

};

void*
CsPPI_EC02time::DB_Thread(void*) {

  pthread_mutex_unlock(&fDBmutex);

  // Enter a eternal loop until the end of the run
  while ( fDBThreadState&THREAD_RUN  || fQuery.size()  ) {
    //       cout<<"still in thread\t"<<  fQuery.size()<<'\t'<<fDBThreadState<<endl ;

    pthread_mutex_unlock(&fDBmutex);
    if ( maxcombined < 1 ) maxcombined = 1;

    // lock to get size
    pthread_mutex_lock(&fDBmutex);

    while ( !fDBThreadState&THREAD_TERM ||
        fQuery.size() ) {

      // unlock from size
      pthread_mutex_unlock(&fDBmutex);

      int count = 0;
      stringstream Query;

      pthread_mutex_lock(&fDBmutex);

      // Make sure that all queued queries are executed
      while ( ( !fDBThreadState&THREAD_TERM || 
            fQuery.size() ) &&
          count < maxcombined) {

        Query << *fQuery.begin() << "; ";

        ++count;

        free(*fQuery.begin());
        fQuery.erase(fQuery.begin());


      }

      // unlock for execution
      pthread_mutex_unlock(&fDBmutex);

      if (!fOutConn->query(Query.str())) {
        char *msg;
        asprintf( &msg, 
            "Error writing to CsPPIoutput database: %s",
            Query.str().c_str());
        MSG( FATAL, UNKNOWN_EXCEPT, msg, __FILE__,__LINE__);
        free(msg);
      }

      fOutConn->endQuery();

      // lock to get size
      pthread_mutex_lock(&fDBmutex);      

    }

    // lock to get size
    pthread_mutex_unlock(&fDBmutex);

    sleep(10);
    pthread_mutex_lock(&fDBmutex);      

  }

  pthread_mutex_unlock(&fDBmutex);
  //   cout<<"leaving thread\n";
  return NULL;
};

void 
CsPPI_EC02time::FillGlobalHistogram() {
  map< unsigned int, 
    map< unsigned int, 
    map< unsigned int, vector <vector <vector<int> > > > > >::iterator Run = fGlobalTimeHist.find(run_nmb);

  if( Run == fGlobalTimeHist.end() ) { 
    map< unsigned int, 
      map< unsigned int, vector <vector <vector<int> > > > > Spill;
    fGlobalTimeHist[run_nmb] = Spill;
    fGlobalTimeHist[run_nmb][spill_nmb] = time_hist;
    return;
  }

  map< unsigned int, 
    map< unsigned int, vector <vector <vector<int> > > > >::iterator Spill =(Run->second).find(spill_nmb);
  if( Spill == (Run->second).end() ) {
    fGlobalTimeHist[run_nmb][spill_nmb] = time_hist;
    return;
  }
  else {
    map< unsigned int, vector <vector <vector<int>  > > >::iterator SrcID;

    for(SrcID = time_hist.begin(); SrcID != time_hist.end(); SrcID++)
    {
      if((Spill->second).find(SrcID->first) != (Spill->second).end()){
        fGlobalTimeHist[run_nmb][spill_nmb][SrcID->first]= SrcID->second;
        return;
      }
      else {
        for(unsigned int port=0;port<SrcID->second.size();port++ ) {
          for(unsigned int ch_grp=0; ch_grp<SrcID->second.at(port).size();ch_grp++) {
            for(unsigned int bin=0; bin<SrcID->second.at(port).at(ch_grp).size();bin++) {
              fGlobalTimeHist[run_nmb][spill_nmb][SrcID->first].at(port).at(ch_grp).at(bin) +=
                SrcID->second.at(port).at(ch_grp).at(bin);
            }
          }
        }
      }
    }

    return;
  }

};


void
CsPPI_EC02time::OpenOutputFile(string filename) {
  int result;
  ec02_tcorr_data_t* ec02_tcorr;

  output_file = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0644);
  if (output_file == -1) {
  MSG( WARNING, 0 ,
          "Error opening output file", 
          __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }
  result = lseek(output_file, filesize-1, SEEK_SET);
  if (result == -1) {
    close(output_file);
    MSG( WARNING, 0 ,
        "Stretching output file to requested size", 
        __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  result = write(output_file, "", 1);
  if (result != 1) {
    close(output_file);
    MSG( WARNING, 0 ,
        "Error writing to output file", 
        __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  time_t seconds;
  seconds = time (NULL);
  ec02_tcorr_header_t header;
  header.version=CsPPI_EC02time_VERSION;
  header.time=seconds;
  header.run=run_nmb;
  header.nu_1=0;
  header.nu_2=0;
    /* mmapp.
     *      */
  ec02_tcorr = (ec02_tcorr_data_t*)mmap(0, filesize, PROT_READ | PROT_WRITE, MAP_SHARED, output_file, 0);
  if (ec02_tcorr == MAP_FAILED) {
    close(output_file);
    MSG( WARNING, 0 ,
        "Error mmapping the output file", 
        __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  ec02_tcorr[0]=*((ec02_tcorr_data_t*)&header);

  output_target =ec02_tcorr;
  return;

}


void
CsPPI_EC02time::WriteChecksum(const unsigned int spill,const unsigned int checksum) {
  ((unsigned int*)output_target)[5+spill]=checksum; 
}

void 
CsPPI_EC02time::Write2OutputFile(const unsigned int index,const unsigned int bin[5]) {
  output_target[index+41]=*(ec02_tcorr_data_t*)bin;
  return;
}
  
void                    
CsPPI_EC02time::CloseOutputFile() {
  if (munmap(output_target, filesize) == -1) {
    MSG( WARNING, 0 ,
        "Error un-mmapping the output file", 
        __FILE__, __LINE__);
  }
  close(output_file);
  return;

}
