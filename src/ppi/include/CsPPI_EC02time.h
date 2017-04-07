/*!
  \file      CsPPI_EC02time.h
  \brief     Compass PreProcessor Class to determine ECAL2 timinig corrections
  \author    Stefan Huber
  \version   $Revision: 1.5 $
  \date      $Date: 2010/11/12 18:25:52 $

  \par       History:
  20091201   Design start of this class
*/

#ifndef __CsPPI_EC02time_H__
#define __CsPPI_EC02time_H___

#include <string>
#include <time.h>
#include <vector>
#include "DaqDataDecoding/DetID.h"
#include "CsPPI.h"
#include "CsEvent.h"
#include "DaqDataDecoding/ChipSADC.h"
#include "MySQLInterface.h"
#include <pthread.h>
#include "CsTimingInfoFormat.h"

#define THREAD_RUN 1
#define THREAD_TERM 2
#define PRE_THREAD_TERM 3
#define CsPPI_EC02time_VERSION 1

void* start_db_thread(void*);
enum output_mode_t {
    database = 0,
    file     =1
};

/*************************************************************
 ** The CsPP_EC02time, derived from the class CsPP,         **
 ** provides the functionality to calculate the mean timing **
 ** per channel and per spill                               **
 *************************************************************/

class CsPPI_EC02time : public CsPPI {

 protected:

                                  
   time_t                              *fTime; /**< time of the run (to select calibrations) */
   char                                *fCTime; /**< a string representation of fTime to use with mysql */
   std::set< unsigned int >             fInitSrcId; /**< Allready initzialized srcids */

   std::map< unsigned int, 
     std::map< unsigned int, 
     std::map< unsigned int, float> > >      fTCalib;  /**< Storage for time calibration */

   std::map< unsigned int,
     std::map< unsigned int,
     std::map< unsigned int, std::vector <std::vector <std::vector<int> > > > > > fGlobalTimeHist;


   /* Variables */
  int                     maxcombined;
  MySQLInterface*         fConn;    /**< Connection to the MYSQL server for input*/
  MySQLInterface*         fOutConn;    /**< Connection to the MYSQL server for output*/
  pthread_mutex_t         fDBmutex;  /**< Mutex for the query queue */
  std::vector<char*>      fQuery;   /**< Query queue */
  pthread_t               fDBThread; /**< Thread to do the insert queries */
  unsigned int            fDBThreadState; /**< runstate of the thread */

  void*                   DB_Thread(void*); /**< Worker funczion of the db thread */
  friend  void*           start_db_thread(void*);

  void                    Write2DB(char *query); /** Function to add a query to fQuery, query will be freed*/


  void                    CloseOutputFile();
  void                    OpenOutputFile(std::string filename) ;
  void                    Write2OutputFile(const unsigned int index,const unsigned int bin[5]);
  void                    WriteChecksum( const unsigned int spill,const unsigned int checksum) ;

  int                     GetTCalib(unsigned int srcid, 
                                    unsigned int port, 
                                    unsigned int ch, 
                                    float *tcalib); /** Retrieve time calibration */
  void                    SetTCalib(unsigned int srcid, 
                                    unsigned int port, 
                                    unsigned int ch, 
                                    float tcalib); /** Store time calibration */

  int                     OpenConn(); /** opens the connection to the MYSQL Server */
  int                     OpenOutConn(); /** opens the connection to the output MYSQL Server */
  int                     DBGetTCalib(unsigned int srcid); /** Returns the number of srcids which were mapped uppon call (call during runtime) */
 public:
                          CsPPI_EC02time(const CS::DetID &id,  
                                         const std::string &TBname);   
//   virtual           
                          ~CsPPI_EC02time();   


  void                    getTime(); /**< set fTime */
  void                    convertTime(); /**< set fCTime */
  int                     InitSrcId(unsigned int srcid); /**< Prepare for SrcID */
 
  // function, which is called during preprocessing 
  int                     process() ; 
  int                     end();
  void                    FillGlobalHistogram();
  int                     write_out_time2file(void);
  int                     write_out_time2db(void);
  void                    increment_time_bin(CS::ChipSADC::DataID& , unsigned int bin);
  std::vector<int>        GetData(void);
 private:

  /***************
   * CFD settings
   **************/
  struct settings_struct {
    int cfd_ampl;
    int cfd_delay;
    int cfd_threshold;
    int time_hist_width;
    int baseline;
  } settings;
  //name of the table to store the time corrections in
  std::string corr_tbl;

  int                     testcount;
  int                     run_nmb;
  int                     spill_nmb;
  output_mode_t           output_mode;
  std::map< unsigned int, std::vector <std::vector <std::vector<int> > > >       time_hist;

  std::vector <std::vector <std::vector<int> > >        empty_time_hist;

  float                   time_reference;

//   float                   tcs_phase;

  static const float      cycle_length;
  int                     min_srcid;
  int                     output_file;
  ec02_tcorr_data_t*      output_target;
  std::string             fdir;
  char*                   pre_fname;
  std::string             fname;
  CS::DetID               *fakeID;



  float                   calc_time(const std::vector<float>& data,std::vector<float>* cfd_time_v);
  int                     nentries;
  int                     filesize;

};

#endif
