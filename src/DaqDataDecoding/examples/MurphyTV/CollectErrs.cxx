#include "CollectErrs.h"

bool operator <(const CollectErrs::VPortID &id1,const CollectErrs::VPortID &id2) {
  return (id1.id<id2.id);
}

CollectErrs::CollectErrs(const char *mapfile) {
  uint32 a;
  pthread_mutex_init(&sync,0);
  RunCount=0;
  Rnr=true;
  Rne=true;
  mapsread=false;
  for (a=0;a<MAX_SOURCE;a++) Catches[a]=0;
  Reset();
  if (mapfile) ReadMaps(mapfile);
}

CollectErrs::~CollectErrs(){
  for (uint32 a=0;a<MAX_SOURCE;a++) {
    delete Catches[a];
    Catches[a]=0;
  }
}

void CollectErrs::ReadMaps(const char *mapfile) {

  static char mapfn[1000]="";
  if (mapfile) strncpy(mapfn,mapfile,1000);
  if (!strlen(mapfn)) strncpy(mapfn,".",1000);
  if (RunNr) {
    vector <string> dn;
    maps.clear();
    mapsread=false; 
    try {Chip::ReadMaps(RunNr,mapfn,maps,opts,dn);}
    catch(...) {cout << "MurphyTV: Error reading mapping file: " << mapfn << endl; return;}
    dn.clear();
    mapsread=true;
  }
}

void CollectErrs::Reset(){
  int a,b;
  for (a=0;a<HISTBUF_SIZE;a++) {
    RecEventsHist[a]=0;
    EventSizeHist[a]=0;
    for (b=0;b<R_MAX_TYPE;b++) ErrHist[a][b]=0;
  }
  for (a=0;a<R_MAX_TYPE;a++) {
    Err[a]=0;
    NrEvntsWithErr[a]=0;
  }
  for (a=0;a<MAX_SOURCE;a++) {
    delete Catches[a];
    Catches[a]=0;
  }
  
  HistBufPos=0;
  NrRecEvents=0;
  ParsedSize=0;
  TriggerNr=0;
  EventTime=0;
  EventSpillNb=0;
  RunNr=0;

}

void CollectErrs::Check4NewRun(DaqEvent *event) {
  uint32 a=event->GetRunNumber();
  bool nr=(a!=RunNr);
  if (nr) RunCount++;
  if (Rne || (nr && Rnr)) Reset();
  RunNr=a;
  if ((Rne || nr) && !mapsread) ReadMaps(); 
  Rne=false;
}

void CollectErrs::DecodeEvent(DaqEvent *evnt) {
  try {
    Check4NewRun(evnt);
    evnt->GetTT().InitAndClear(opts.GetTTConfig());	
    evnt->ReadChips(opts);
    if (!mapsread) {
    	cout << "MurphyTV: Mapping file not available." << endl;
        cout << "MurphyTV: Exiting." << endl;
	exit(1);
    }
    Chip::Digits digits;
    for( vector<Chip*>::iterator it=evnt->GetChips().begin();
	 it!=evnt->GetChips().end(); it++ ) (*it)->Decode(maps,digits,opts);
//    evnt->CheckSrcIDs(opts);

    evnt->GetTT().Decode(digits,opts,evnt->GetDaqErrors());
    digits.Clear();

  }
  catch ( DaqError err) {
    err.Print(cout,"MurphyTV: DaqError: ");
    return;            
  }
  catch( Exception &e ) {
    cerr << "MurphyTV: Chip Decoding Exception: " << e.what() << "\n";
    return;
  }
  catch (...) {
    cerr << "MurphyTV: unknown exception in DecodeEvent.\n Exiting. \n";
    exit(1);
  }
}

void CollectErrs::HandleNewError(uint32 trnr,const DaqError &errit){

  try {
    DaqError::Type typ=errit.GetType();
    if(  (typ>=DaqError::MAX_TYPE) ) {
      //unknown error
      cout << "MurphyTV: Error in event: " << trnr << endl;
      errit.Print(cout,"MurphyTV: Unknown Daq Error: ");
      return;
    }
        
    int sourceid=0;
    switch( errit.GetArgs().count(DaqError::SOURCE_ID) )
    {
        case 0: sourceid=0; break; // defaul value is used
        case 1: sourceid=(int)errit.GetArg(DaqError::SOURCE_ID); break;
	default:
	{
          cout << "MurphyTV: Error in event: " << trnr << endl;
          errit.Print(cout,"MurphyTV: Error without source ID:  ");
          return;
	}
    }
    
    if (sourceid>=MAX_SOURCE||sourceid<0) {
      cout << "MurphyTV: Error in event: " << trnr << endl;
      cout << "MurphyTV: Got SourceID " <<sourceid<<" >1024 from DaqError"<<endl;
      return;
    }

    // TIGER monitoring sources (not handled in MurphyTV)
    if (sourceid>=880 && sourceid<=881) return;
    // GANDALF
    if (sourceid>=800 && sourceid<=881) {
        // conflict in format type for first event in spill
        if (typ>=DaqError::ChipHotGeSiCA_WRONG_EVENT1 && typ<=DaqError::ChipHotGeSiCA_WRONG_EVENT) return;
    }

  
    ErrHist[HistBufPos][typ]++;
    Err[typ]++;
    ErrEx[typ]=true; 
 
    bool hasport=(errit.GetArgs().count(DaqError::PORT)>0);
    bool hasgeoid=(errit.GetArgs().count(DaqError::GEO_ID)>0);
    bool hasvalue=(errit.GetArgs().count(DaqError::VALUE)>0);
    int port=INVALID;
    int geoid=INVALID;
    int value=INVALID;
    if (hasport) {
      port=(int)errit.GetArg(DaqError::PORT);
      if (port>=MAX_PORT) {
        cout << "MurphyTV: Error in event: " << trnr << endl;
	cout << "MurphyTV: Got Port# " <<port<<" >="<<MAX_PORT<<" from DaqError"<<endl;
	return;
      }
    }
    if (hasgeoid) {
      geoid=(int)errit.GetArg(DaqError::GEO_ID);
      if (geoid>=MAX_GEOID) {
        cout << "MurphyTV: Error in event: " << trnr << endl;
	cout << "MurphyTV: Got GeoID " <<geoid<<" >="<<MAX_GEOID<<" from DaqError"<<endl;
	return;
      }
    }
    // use GeoID field for parameter VALUE, but only for errors 
    //   without GeoID or port
    if (!hasgeoid && !hasport && hasvalue) {
       value=(int)errit.GetArg(DaqError::VALUE);
       if (value!=INVALID && value>0 && value<MAX_GEOID) geoid=value;
    }
    CatchStat *ctch=GetCatch(sourceid);
    ctch->AddError(trnr,EventSpillNb,typ,port,geoid);    
  }
  catch( Exception &e ) { 
    cerr << "MurphyTV: Decoding Library Exception: " << e.what() << "\n";
    return;
  }
  catch(...) {
     cerr << "MurphyTV: Unknown exception in decoding library. \n Exiting. \n";
     exit(1);
  }
}

void CollectErrs::HandleNewEvent(DaqEvent *event){
  uint32 a,b,c,d;
  const uint32 *p=0,*p2=0;
  
  try {
  Check4NewRun(event);
  
  EventTime=event->GetTime().first;
  EventSpillNb=event->GetBurstNumber();
  ParsedSize+=event->GetLength();
  uint32 trnr=event->GetEventNumberInBurst();
  c=GetNrRecEvents()/HISTBUF_RES;
  
  if (c < HISTBUF_SIZE) HistBufPos=c;  
  
  for (a=0;a<R_MAX_TYPE;a++) ErrEx[a]=false; 
  for (a=0; a<MAX_SOURCE; a++) Catch_in_event[a]=false;
  for (a=0;a<MAX_SOURCE;a++) if (CatchExists(a)) GetCatch(a)->BeginEvent();
  
  vector <Chip *>::iterator ci=event->GetChips().begin();
  uint32 tstenr=0xffffffff;
  for (;ci!=event->GetChips().end();ci++){
    if (tstenr==0xffffffff)
      tstenr=(*ci)->GetSLink().GetEventNumber();

    a=(*ci)->GetSourceID();

   
    
    //Checks the status bit of the TCS (important in case of event skips!)
    //    0x18: normal physics event without problems
    
    // Gesicas have a different TCS bit setting!
      if ((*ci)->GetSLink().GetStatus() !=0x18 && ((a<510) && (a!=110) && (a!=111) && (a!=380) && (a!=381) && (a!=382))) {
     		DaqError::Type typ;
		if ((*ci)->GetSLink().GetStatus() & 0x80) typ=DaqError::TCS_SYNC;
		else if (((*ci)->GetSLink().GetStatus() & 0x03) == 0x03) typ=DaqError::TCS_FIFO;
		else if (((*ci)->GetSLink().GetStatus() & 0x03) == 0x02) typ=DaqError::TDC_ERR12;
		else typ=DaqError::TCS_UNDEF;
     		int port=INVALID;
     		int geoid=INVALID;
     		ErrHist[HistBufPos][typ]++;
     		Err[typ]++;
     		ErrEx[typ]=true;   
     		CatchStat *ctch=GetCatch(a);
     		ctch->AddError(trnr,EventSpillNb,typ,port,geoid);   
  	}

    if (a >= MAX_SOURCE||a<0) {
      cout<< "MurphyTV : got SourceId " << a << ">1024 from DaqDataDecoding"<<endl;
      continue;
    }
    
    Catch_in_event[a]=true;
    

    
    //the rest of this loop counts header and data words of the different modes,
    //maybe this should also be moved to decoding library some day.

    if ((*ci)->GetSLink().IsFirstEventInRun()) continue; //todo
    CatchStat *ctch=GetCatch(a);
    ctch->FillSizeHist((*ci)->GetLength());
    if ((*ci)->GetSLink().GetFormat()&1) {
      if ((*ci)->GetSLink().GetFormat()&2) {
        p=(*ci)->GetDataStart();
        p2=(*ci)->GetDataEnd();
        if (p>=p2) continue; 
        a=*p >> 22;
        c=0;
        for (;p<p2;p++){
          if ((*p >>20) !=0xFFF) {
            b= *p >> 22;
            if (b!=a) {
              ctch->AddHD(INVALID,a,0,c); 
              a=b;
              c=0;
            } 
            c++;    
          } else p++;
        }
	ctch->AddHD(INVALID,a,0,c);  
      }
      else {
        p2=(*ci)->GetDataEnd();
        for (p=(*ci)->GetDataStart();p<p2;p++){
          if (((*p >>20) & 0xFFF )==0xFFF) continue;
          a=(*p >> 4) & 0xf;
          if (*p & (1<<31)) ctch->AddHD(a,INVALID,0,1);          
          else ctch->AddHD(a,INVALID,1,0);              
        }
      }
    }
    else {
      if ((((*ci)->GetSLink().GetFormat() >> 4) & 7 ) ==0) {
        ctch->AddData((*ci)->GetDataEnd()-(*ci)->GetDataStart());
        continue;
      }
	

      p2=(*ci)->GetDataEnd();
      a=INVALID;
      b=c=0;

      // GandalfADC
      if ((((*ci)->GetSLink().GetFormat() >> 4) & 7) == 1) {
    	  int channel=0,total=0;
    	  bool inChannel=false;

    	  // debug mode
    	  if ((((*ci)->GetSLink().GetFormat() >> 2) & 1) == 1) {
			  for (p=(*ci)->GetDataStart();p<(*ci)->GetDataEnd();p++){
				  // header/trailer
				  if ( ((*p) >> 31) == 0) {
					  c++;
					  if (!inChannel) {
						inChannel=true;
					  	channel = ((*p)>>20)&0x0f;
					  }
					  else {
						  // TODO: if (channel != ((*p)>>20)&0x0f)
						  ctch->AddHD(channel,a,2,b);
						  total+=b;
						  b=0;
						  inChannel=false;
					  }
				  }
				  // count data
				  else b++;
			  }
			  ctch->AddHD(INVALID,a,c,total);
    	  } else {
    		  ctch->AddHD(INVALID,a,0,(*ci)->GetDataEnd()-(*ci)->GetDataStart());
    	  }
    	  continue;
      }

      
      uint16 hpos=0;
      for (p=(*ci)->GetDataStart();p<p2;p++){
        if (((*p >>20) & 0xFFF )==0x7FF) {
          p++;
          continue;
        }
        if (*p & (1<<31)) c++;          
        else {
          if (*p & (1<<30)) {
            d=(*p >> 20) & 0x3ff;
            if (d!=a) ctch->AddHD(INVALID,d,1,0);  
            else b++;
            if (b>0 || c>0) ctch->AddHD(INVALID,a,b,c);  
            c=b=0;
            a=INVALID;
          } 
          else {
            if (b>0 || c>0) {
              ctch->AddHD(INVALID,a,b,c);   
              c=0;
            }
            a=(*p >> 20) & 0x3ff;
            b=1;
            
            hpos++;
          }
        }
      }
      if (b>0 || c>0) ctch->AddHD(INVALID,a,b,c);      
    } 
  }
 
  const vector<DaqError> &errs = event->GetDaqErrors().errors;  // this is just a short name.
  for( vector<DaqError>::const_iterator e=errs.begin(); e!=errs.end(); e++ )
      HandleNewError(trnr,*e);

         
  for (a=0;a<MAX_SOURCE;a++) if (CatchExists(a)) GetCatch(a)->EndEvent();
 
  for (a=0;a<R_MAX_TYPE;a++) if (ErrEx[a]) NrEvntsWithErr[a]++;
  NrRecEvents++;
  RecEventsHist[HistBufPos]++; 
  EventSizeHist[HistBufPos]+=event->GetLength();

} 
 catch ( DaqError err) {
    err.Print(cout,"MurphyTV : DaqError: ");
    return;
    }  
  catch( Exception &e ) { 
    cerr << "MurphyTV: Decoding Library Exception: " << e.what() << "\n";
    return;
    }
  catch (...) {
    cerr << "MurphyTV: Unknown exception in event handling.\n Exiting \n.";
    exit(1);
    }
}


uint32 CollectErrs::GetErrSum() {
  uint32 a,b;
  for (a=b=0;a<R_MAX_TYPE;a++) b+=Err[a];
  return b;
}
void CollectErrs::GetDetsAtCatch(uint16 source,set <string> &detnames){
	Chip::DataID did0;
	Chip::DataID did1;
	// TODO: this gandalf specific DataID should be joined with the default one
	if (source >=800 && source <=840) {
		did0 = ChipGandalf::DataID(source,0,0);
		did1 = ChipGandalf::DataID(source+1,0,0)-1;
	} else {
		did0 = Chip::CreateDataID5(source,0,0,0,0);
		did1 = Chip::CreateDataID5(source+1,0,0,0,0)-1;
	}
  	for (Chip::Maps::iterator it=maps.lower_bound(did0);it!=maps.upper_bound(did1);it++) 
    		detnames.insert(it->second->GetDetID().GetName());
}

CollectErrs::CatchStat *CollectErrs::GetCatchNC(uint16 source) {
  if (source>=MAX_SOURCE) return &defctch; 
  else if (!Catches[source]) return &defctch;
  else return Catches[source];
}





string CollectErrs::GetErrorName(DaqErrorType err) {
  string en;
  if (err==DaqError::EXCEPTION) en="DaqDataDecoding EXCEPTION";
  else if (err<DaqError::MAX_TYPE) en=CS::DaqErrorType::GetDaqErrorType(err).GetName();
  else en="Unknown error type";
  return en;
}

void CollectErrs::PrintErr(ostream &out,const string &suffix,int nr,int bad_events){
  int n=GetNrRecEvents();
  if (n>0) {
    ostringstream tmp;
        if ((float(bad_events)/float(n)*100.0)<10) 
    	  tmp << setw(7) << nr << " (" << setprecision(2) <<
	  (float(bad_events)/float(n)*100.0) << "% bad events) "; 
        else
	  tmp << setw(7) << nr << " (" << setprecision(3) <<
	  (float(bad_events)/float(n)*100.0) << "% bad events) "; 
    out.setf(ios::left, ios::adjustfield);
    out  << setw(30) << tmp.str() << suffix <<endl;
  }
}


void CollectErrs::PrintReport(ostream &out, set<int> srcids) {

 multimap <uint32,int> errl;

 out << endl;
 out << "MurphyTV error report for COMPASS data:" << endl;
 out << "=======================================" << endl << endl;
 out << "Run number: " << GetRunNr() << endl;
 out << "Number of Events analyzed: " << GetNrRecEvents()  << endl;
 out << endl << endl;
 for (int a=DaqError::EXCEPTION;a<R_MAX_TYPE;a++) {
   uint32 nr=GetErr(DaqError::Type(a));
   if (nr>0) errl.insert(make_pair(nr,a));  
 }    
 if (errl.size()==0) {
   out << "No errors found!" << endl;
 } else {
  out << "Errors found from the following SourceIDs:" << endl;
  out << "==========================================" << endl << endl;
  
  multimap <uint32,int> catchl;
  for (int a=0;a<MAX_SOURCE;a++) if (CatchExists(a)) {
    if (srcids.size()>0 && srcids.count(a)==0) continue;
    uint32 nr=GetCatchErrSum(a);
    if (nr>0) catchl.insert(make_pair(nr,a));
  }
  for (multimap <uint32,int>::reverse_iterator cit=catchl.rbegin();cit!=catchl.rend();cit++) {
    int a=cit->second;
    out << "SourceID " << a <<":";
    set <string> atdet;
    GetDetsAtCatch(a,atdet);
    if (atdet.size()) {
      out << " (";
      for (set<string>::iterator dni=atdet.begin();dni!=atdet.end();dni++) {
	out << *dni << " ";
      }
      out << ")";
    }
    atdet.clear();
    out << endl;
    out << "-------------" << endl;
    PrintErr(out,"Sum of errors on this SourceID",GetCatchErrSum(a),GetCatchNrBadEvnts(a));
    multimap <uint32,DaqError::Type> errs;
    for (int err=DaqError::EXCEPTION;err<R_MAX_TYPE;err++) {
      uint32 nerr=GetErrOnCatch(a,DaqError::Type(err));
      if (nerr>0) errs.insert(make_pair(nerr,(DaqError::Type)err));
    }
    for (multimap<uint32,DaqError::Type>::reverse_iterator eit=errs.rbegin();eit!=errs.rend();eit++) {
      int err=eit->second;
      PrintErr(out,GetErrorName(DaqError::Type(err)),
                   GetErrOnCatch(a,DaqError::Type(err)),
	           GetEvntsWithErrOnCatch(a,DaqError::Type(err)));
      multimap <uint32,int> portl;
      multimap <uint32,int> geoidl;
      
      map <uint16,uint32> pl0;
      map <uint16,uint32> gl0;
      const VPortStatMap &pst=GetVPortStats(a);
      for (VPortStatMap::const_iterator it=pst.begin();it!=pst.end();it++) {
	uint32 nr=it->second.GetNumErr((DaqError::Type)err);
	if (nr==0) continue;
	if (it->first.Port()!=INVALID) pl0[it->first.Port()]+=nr;
	if (it->first.GeoID()!=INVALID) gl0[it->first.GeoID()]+=nr;
      }
      for (multimap<uint16,uint32>::iterator it=pl0.begin();it!=pl0.end();it++) 
	portl.insert(make_pair(it->second,it->first));
      pl0.clear();
      for (multimap<uint16,uint32>::iterator it=gl0.begin();it!=gl0.end();it++) 
	geoidl.insert(make_pair(it->second,it->first));
      gl0.clear();
      if (portl.size() || geoidl.size()) {
        out << "                            ";
        if (portl.size()) {
	  out << "  Port ";
	  for (multimap <uint32,int>::reverse_iterator pit=portl.rbegin();pit!=portl.rend();pit++) 
	    out << pit->second << " ";
        }
        if (geoidl.size()) {
          if (portl.size()) {
            out << "  GeoID "; 
            } else {
	        out << "  Value ";
            }
	  for (multimap <uint32,int>::reverse_iterator pit=geoidl.rbegin();pit!=geoidl.rend();pit++) 
	    out << pit->second << " ";
        }
        out << endl;
      }
      portl.clear();
      geoidl.clear();
    }
    errs.clear();
    out << endl;
  }
  catchl.clear(); 
 
  out << endl; 
  out << "Summary of all errors from all SourceIDs in this run:" << endl;
  out << "=====================================================" << endl;
  for (multimap <uint32,int>::reverse_iterator eit=errl.rbegin();eit!=errl.rend();eit++){
    DaqError::Type err=(DaqError::Type)eit->second;
    PrintErr(out,GetErrorName(err),eit->first,GetNrEvntsWithErr(err));
    multimap <uint32,int> catchl;
    for (int a=0;a<MAX_SOURCE;a++) if (CatchExists(a)) {
      if (srcids.size()>0 && srcids.count(a)==0) continue;
      uint32 nr=GetErrOnCatch(a,err);
      if (nr>0) catchl.insert(make_pair(nr,a));
    }
      out <<"                              SourceIds ";
      for (multimap <uint32,int>::reverse_iterator cit=catchl.rbegin();cit!=catchl.rend();cit++) {
	int a=cit->second;
	out  << a << " ";
    }
    catchl.clear();
    out << endl;
  }
  errl.clear();
 }
  
 out << endl << endl; 
 out << "List of all SourceIDs in this run:" << endl;
 out << "==================================" << endl; 
 out << right << setw(10) << "SourceID" << setw(20) << "BadEvents" << setw(20) << "#Errors" << setw(20);
 out << "#Headers" << setw(20) << "#Data" << endl;
 for (int a=0;a<MAX_SOURCE;a++) if (CatchExists(a)&&(GetNrRecEvents()>0)) {
    if (srcids.size()>0 && srcids.count(a)==0) continue;
    float nr=(float)GetCatchErrSum(a)/GetNrRecEvents();
    float hdr=(float)GetNrHeaderOnCatch(a)/GetNrRecEvents();
    float data=(float)GetNrDataOnCatch(a)/GetNrRecEvents();
    float nrbe=(float)100.0*GetCatchNrBadEvnts(a)/GetNrRecEvents();
    out << right <<setw(10) << a   <<fixed<<showpoint <<setprecision(1);
    out <<setw(19) << nrbe << "%" << setprecision(4) <<setw(20) << nr;
    out <<setprecision(2) <<setw(20) << hdr << setw(20) << data << endl;
 }
 out << endl;
}


CollectErrs::CatchStat::CatchStat(){
  uint32 a;

  BadEvnts=Data=Header=ErrSum=0;
  NbSpillWithErrGen = LastSpillWithErrGen = 0;
 
  for ( a=0 ; a<R_MAX_TYPE ; a++ ) {
    ErrEx[a]=false;
    EvntsWithErr[a]=0;
    Err[a]=FirstErrPos[a]=0;
    NbSpillWithErr[a] = 0;
    LastSpillWithErr[a] = 0;
  }

  for (a=0;a<SIZEHIST_SIZE;a++) SizeHist[a]=0;
}

CollectErrs::VPortStatMap::iterator  CollectErrs::CatchStat::FindVPort(uint16 port,uint16 geoid){
  VPortStatMap::iterator vpi=VPortStats.find(VPortID(port,geoid));
  if (vpi!=VPortStats.end()) return vpi;
  if (VPortStats.size() < 50 ) {    
      return VPortStats.insert(VPortStats.begin(),make_pair(VPortID(port,geoid),VPortStat()));    
  }
  return vpi;
}

void CollectErrs::CatchStat::AddError(uint32 eventnr, uint32 spillnb, DaqError::Type typ,uint16 port,uint16 geoid){
  ErrSum++;    
  Err[typ]++;
  ErrEx[typ]=true;
  if (0==FirstErrPos[typ]) FirstErrPos[typ]=eventnr;

  if (spillnb > LastSpillWithErrGen) {
    LastSpillWithErrGen = spillnb;
    NbSpillWithErrGen++;
  }
  if (spillnb > LastSpillWithErr[typ]) {
    LastSpillWithErr[typ] = spillnb;
    NbSpillWithErr[typ]++;
  }

  if (port==INVALID && geoid==INVALID) return;
  map <VPortID,VPortStat>::iterator vpi=FindVPort(port,geoid);
  if (vpi==VPortStats.end()) return;
  vpi->second.NumAnyErr++;
  vpi->second.NumErr[typ]++;
}

void CollectErrs::CatchStat::AddHD(uint16 port,uint16 geoid,uint32 header,uint32 data){
  Data+=data;
  Header+=header;
  VPortStatMap::iterator vpi=FindVPort(port,geoid);
  if (vpi==VPortStats.end()) return;
  vpi->second.NumData+=data;
  vpi->second.NumHeader+=header;
}

void CollectErrs::CatchStat::BeginEvent() {
  for (uint32 a=0;a<R_MAX_TYPE;a++) ErrEx[a]=false;
}


void CollectErrs::CatchStat::EndEvent() {
  bool bad=false;
  for (uint32 a=0;a<R_MAX_TYPE;a++) if (ErrEx[a]) {
    EvntsWithErr[a]++;
    bad=true;
  }
  if (bad) BadEvnts++;
}


void CollectErrs::CatchStat::FillSizeHist(uint32 size) {
  uint32 bin=size/SIZEHIST_RES;
  if (bin<SIZEHIST_SIZE) SizeHist[bin]++;  
}

CollectErrs::VPortStat& CollectErrs::VPortStat::operator +=(const VPortStat &ps) {
  NumHeader+=ps.NumHeader;
  NumData+=ps.NumData;
  NumAnyErr+=ps.NumAnyErr;
  for (map <DaqError::Type,uint32>::const_iterator it=ps.NumErr.begin();it!=ps.NumErr.end();it++) 
    NumErr[it->first]+=it->second;
  return *this;
}

uint32 CollectErrs::VPortStat::GetNumErr(DaqError::Type typ) const {
  map <DaqError::Type,uint32>::const_iterator it=NumErr.find(typ);
  if (it!=NumErr.end()) return it->second;
  return 0;
}
