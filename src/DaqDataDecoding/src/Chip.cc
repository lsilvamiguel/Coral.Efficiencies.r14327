#include <cstdio>
#include <fstream>
#include <sstream>
#include <dirent.h>
#include <sys/stat.h> 
#include <unistd.h>

#include "utils.h"
#include "Chip.h"
#include "ChipSinica.h"
#include "ChipF1.h"
#include "ChipADC.h"
#include "ChipSADC.h"
#include "ChipAPV.h"
#include "ChipAPVRICH.h"
#include "ChipGassiplex.h"
#include "ChipHotGeSiCA.h"
#include "ChipGandalf.h"
#include "Scaler.h"
#include "DaqOption.h"
#include "DaqEvent.h"
#include "Event1Run.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

bool                    Chip::fill_stat_flag            = true;
map<int,int>            Chip::map__slink_source_id;
map<int,int>            Chip::map__slink_type;
map<int,int>            Chip::map__slink_format;
map<pair<int,int>,int>  Chip::map__src_fmt;

////////////////////////////////////////////////////////////////////////////////

Chip::Chip(void const * const buffer,bool copy_buf,DaqOption &opt,DaqEvent &ev) :
  flag_buffer_delete (copy_buf),
  flag_can_modify    (false),
  buf                ((uint32*)buffer),
  is_scaned          (false),
  event              (&ev)
{
  if( buf==NULL )
    throw Exception("Chip::Chip():  buffer=NULL");

  if( flag_buffer_delete )
  {
    const size_t len = GetLength();
    buf = (uint32*) realloc(NULL,len);
    if( buf==NULL )
      throw Exception("Chip::Chip(): Can not allocate %d bytes of memory.",len);
    memcpy(buf,buffer,len);
  }
  
  if( fill_stat_flag )
  {
    map__slink_source_id [GetSLink().GetSourceID()     ] ++; 
    map__slink_type      [GetSLink().GetEventType()    ] ++;
    map__slink_format    [GetSLink().GetFormat()       ] ++;
    map__src_fmt[pair<int,int>(GetSLink().GetSourceID(),GetSLink().GetFormat())]++;
  }
  
  GetSLink().AnalyseErrors(*this,opt);
}

////////////////////////////////////////////////////////////////////////////////

Chip *Chip::Create(void const * const buffer,bool copy_buf,DaqOption &opt,DaqEvent &ev)
{
  Chip *eq=NULL;

  const SLink *slink = (const SLink*)(buffer);

  if( slink->IsFirstEventInRun() )
  {
    if( ev.GetEvent1Run()!=NULL )
        ev.GetEvent1Run()->Add(ev.GetRunNumber(),slink);
    return NULL;
  } 
  
  int32 format = slink->GetFormat()&127;

  for( vector<Chip::CreationRule>::const_iterator it=opt.GetChipCreationRules().begin();
       it!=opt.GetChipCreationRules().end(); it++ )
    if( it->Accept(slink->GetEventType(),slink->GetSourceID(),format) )
    {
      if     ( it->chip_name=="ChipF1" )
        eq = new ChipF1 (buffer,copy_buf,opt,ev);
      else if( it->chip_name=="ChipSinica" )
	eq = new ChipSinica(buffer,copy_buf,opt,ev);
      else if( it->chip_name=="ChipADC" )
        eq = new ChipADC(buffer,copy_buf,opt,ev);
      else if( it->chip_name=="ChipSADC" )
        eq = new ChipSADC(buffer,copy_buf,opt,ev);
      else if( it->chip_name=="ChipAPV" )
        eq = new ChipAPV(buffer,copy_buf,opt,ev);
      else if( it->chip_name=="ChipAPVRICH" )
        eq = new ChipAPVRICH(buffer,copy_buf,opt,ev);
      else if( it->chip_name=="ChipGassiplex" )
        eq = new ChipGassiplex(buffer,copy_buf,opt,ev);
      else if( it->chip_name=="Scaler" )
        eq = new Scaler(buffer,copy_buf,opt,ev);
      else if( it->chip_name=="ChipGandalf" )
        eq = new ChipGandalf(buffer,copy_buf,opt,ev);
      else
      {
        ev.GetDaqErrors().Add(Exception("Chip::Create(): unknown chip \"%s\"",it->chip_name.c_str()));
        return NULL;
      }
      break;
    }

  if( eq==NULL && opt.Use_SLinkFormat() )
  {
    // We failed to guess a chip type from the CreationRules.
    // Let's try to use information from the SLink format field.

    if( format&1 )
    {
      // This is TDC-F1 readout
      eq = new ChipF1(buffer,copy_buf,opt,ev);
    }
    else if( format&2 ){
      // This is Academia Scinica Format
      eq = new ChipSinica(buffer,copy_buf,opt,ev);
    }
    else
    {
      // This is ADC readout.

      const uint8 who=(format>>4);
      //uint8 mode=(format>>1)&7;
      //Exception("Chip::Create():II: ADC: who=%d mode=%d",int(who),int(mode));

      switch( who )
      {
        case 0:  // GeSiCA (APV)
          eq = new ChipAPV(buffer,copy_buf,opt,ev);
          break;

        case 1:
          eq = new ChipHotGeSiCA(buffer,copy_buf,opt,ev);
          break;

        case 2:  // Pattern-CMC
          throw Exception("Chip::Create(): no code for Pattern-CMC data");
          break;

        case 3:  // Scaler-CMC
          eq = new Scaler(buffer,copy_buf,opt,ev);
          break;

        case 4:  // FI-ADC HOTLink
          eq = new ChipADC(buffer,copy_buf,opt,ev);
          break;

        case 5:  // RICH HOTFibre
          eq = new ChipGassiplex(buffer,copy_buf,opt,ev);
          break;

        case 7:  // Scaler HOTLink
          throw Exception("Chip::Create(): no code for Scalar HOTLink data");
          break;

        default:
          Exception("Chip::Create(): format %d is not specified!!!",format);
      }
    }
  }

  if( eq==NULL ) // Chip was not identifed... Create a dummy version!
  {
    ev.GetDaqErrors().Add(Exception("Chip::Create():  Dummy chip. DATE:[Type=%d EquipmentID=%d] SLink:[SourceID=%d EventType=%d Format=%d]",
              -1, -1,
              (int)(slink==NULL?-1:(int)slink->GetSourceID()),
              (int)(slink==NULL?-1:(int)slink->GetEventType()),
              (int)(slink==NULL?-1:format) ));

    eq = new Chip(buffer,copy_buf,opt,ev);
  }

  return eq;
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Scan(DaqOption &opt)
{}

////////////////////////////////////////////////////////////////////////////////

void Chip::AddError(const DaqError &e,DaqOption &o)
{
    /// First we register the error
    event->GetDaqErrors().Add(e);

    /// And second we take the action
    switch( e.GetDaqErrorType().GetAction() )
    {
        case DaqError::Action::NOTHING:
            break;

        case DaqError::Action::DISCARD_EVENT_ON_GEOID:
        {
            int srcID=-1, geoID=-1, port=-1;

            try
            {
                srcID = (int) e.GetArg(DaqError::SOURCE_ID);
            }
            catch(...)
            {}

            try
            {
                geoID = (int) e.GetArg(DaqError::GEO_ID);
            }
            catch(...)
            {}

            try
            {
                port = (int) e.GetArg(DaqError::PORT);
            }
            catch(...)
            {}
            
            if( srcID==-1 || (geoID==-1 && port==-1) )
            {
                printf("Chip::AddError(): I need both srcID and geoID/port for the error:\n");
                e.Print(cout,"(missing srcID or geoID):   ");
            }
            else
            {
                // Add geoID,port or both as a bad data source.

                if( geoID!=-1)
                    event->AddBadDataSource(srcID,geoID);
                if( port!=-1)
                    event->AddBadDataSource(srcID,port);
            }
            
            break;
        }

        case DaqError::Action::DISCARD_EVENT_ON_PORT:
        {
            int srcID=-1, port=-1;

            try
            {
                srcID = (int) e.GetArg(DaqError::SOURCE_ID);
            }
            catch(...)
            {}

             try
            {
                port = (int) e.GetArg(DaqError::PORT);
            }
            catch(...)
            {}
            
            if( srcID==-1 || port==-1 )
            {
                printf("Chip::AddError(): I need both srcID and port for the error:\n");
                e.Print(cout,"(missing srcID or port):   ");
            }
            else
            {
	      // Add port or both as a bad data source.
	      
	      if( port!=-1)
                    event->AddBadDataSource(srcID,port);
            }
            
            break;
        }


        case DaqError::Action::DISCARD_EVENT_ON_SRCID:
        {
            int srcID=-1;

            try
            {
                srcID = (int) e.GetArg(DaqError::SOURCE_ID);
            }
            catch(...)
            {}
            
            if( srcID==-1 )
            {
                printf("Chip::AddError(): I need to know srcID for the error:\n");
                e.Print(cout,"(missing srcID):   ");
            }
            else
                event->AddBadDataSource(srcID);
            
            break;
        }

        case DaqError::Action::DISCARD_EVENT:
            // The DaqEvent already knows about it...
            break;

        default:
            throw Exception("Chip::AddError(): internal error");
    }
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Add(const void *buffer,size_t n)
{
  // n is size in bytes

  if( !CanModify() )
    throw Exception("Chip::Add(): You do not have rights to modify the chip data.");
  
  buf = reinterpret_cast<uint32*>(realloc(buf,GetLength()+n*4));
  if( buf==NULL )
    throw Exception("Chip::Add(): Can not allocate %d bytes of memory.",n*4);

  memcpy(reinterpret_cast<uint8*>(buf)+GetLength(),buffer,n*4);

  GetSLink().SetEventSize(GetSLink().GetEventSize()+n);
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Print(ostream &o,const string &prefix) const
{
  GetSLink().Print(o,prefix +" S-Link ");
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Dump(ostream &o,const string &prefix) const
{
  o << prefix << "Chip::Dump():  This is the S-Link header:\n";
  GetSLink().Print(o,prefix+" S-Link header: ");

  for( unsigned i=0; i<sizeof(SLink)/4; i++ )
  {
    o << prefix;
    bits_print(o,reinterpret_cast<const uint32*>(&GetSLink())[i],0,32,"%8x\n");
  }

  o << prefix << "Chip::Dump():  These are " << int(GetDataEnd()-GetDataStart())
    << " raw data words:\n" ;

  for( const uint32 *p=GetDataStart(); p<GetDataEnd(); p++ )
  {
    o << prefix;
    char s[111];
    sprintf(s,"%4d ",int(p-GetDataStart()));
    o<<s;
    bits_print(o,*p,0,32,"%8x\n");
  }
  
  o << prefix << "Chip::Dump():  This is the end.\n";
}

////////////////////////////////////////////////////////////////////////////////

static void my_print(const map<int,int> &m,const string &s,ostream &o,const string &prefix)
{
  for( map<int,int>::const_iterator it=m.begin(); it!=m.end(); it++ )
  {
    o << prefix;
    char ss[1111];
    sprintf(ss,"%s %5d was found %8d time%c\n",s.c_str(), it->first, it->second,  it->second>1?'s':' ' );
    o << ss;
  }
  o << "\n";
}

void Chip::PrintStatistics(ostream &o,const string &prefix)
{
  o << "\n\n"
    << prefix << "Chip::PrintStatistics:\n\n";

  my_print(map__slink_source_id,"S-Link source ID",     o,prefix);
  my_print(map__slink_type,     "S-Link type",          o,prefix);
  my_print(map__slink_format,   "S-link format",        o,prefix);

  for( map<pair<int,int>,int>::const_iterator it=map__src_fmt.begin(); it!=map__src_fmt.end(); it++ )
  {
    o << prefix;
    char s[1111];
    sprintf(s,"SrcID,format=(%3d,%3d)    was found %8d time%c\n",
              it->first.first, it->first.second, it->second,  it->second>1?'s':' ' );
    o << s;
  }

  o << "\n\n";
}

////////////////////////////////////////////////////////////////////////////////

bool Chip::ReadMaps(uint32 run,const string &path,Maps &maps,DaqOption &daq_option,vector<string> &dets)
{
  bool ok=true;
  //cout << "Chip::ReadMaps(): run=" << run << "   path is \"" << path << "\"\n";

  // First of all we try to open
  struct stat info;
  if( 0!=stat(path.c_str(),&info) )
    throw Exception("Chip::ReadMaps(): bad path \"%s\"",path.c_str());

  if( S_ISLNK(info.st_mode) )
    Exception("Chip::ReadMaps():WW: LINK!").Print();

  if( S_ISDIR(info.st_mode) )
  {
    // The path is a directory.
    // Scan that directory for files with the extension ".xml"

    cout << "Chip::ReadMaps(): reading files from the directory \"" << path << "\"\n";

    struct dirent **namelist;
    int n = scandir(path.c_str(),&namelist,0,0);

    while( n-- > 0 )
    {
      const char *name=namelist[n]->d_name;  // file name
      if( strlen(name)>4 &&
          0==strcmp(name+strlen(name)-4,".xml") )
      {
        string s=path+'/'+name;
        try
        {
          if( !Chip::ReadMaps(run,s,maps,daq_option,dets) )
            ok = false;
        }
        catch( const std::exception &e )
        {
          cerr << "Chip::ReadMaps: std::exception: "<< e.what() << " while reading map file "<<name<< "\n";
          ok = false;
        }
        catch( const DaqError &e )
        {
          cerr << "Chip::ReadMaps: DaqError while reading map file "<<name<<"\n";
          e.Print();
          ok = false;
        }
        catch( const char *s )
        {
          cerr << "Chip::ReadMaps: exception : "<< s << " while reading map file "<<name<<"\n";
          ok = false;
        }
        catch( ... )
        {
          cerr << "Chip::ReadMaps: Unknown problem in map file \"" << name << "\" reading.\n";
          ok = false;
        }
      }
      free(namelist[n]);
    }
    free(namelist);
  }
  else
  if( S_ISREG(info.st_mode) )
  {
    //Regular file
    ObjectXML objs("all");
    ObjectXML::Parse(path.c_str(),objs);
    
    if( objs.GetChildren().size()!=1 )
        throw Exception("Chip::ReadMaps(): Bad maps?");

    for( list<ObjectXML*>::iterator it=objs.GetChildren().front()->GetChildren().begin(); it!=objs.GetChildren().front()->GetChildren().end(); it++ )
    {
      // Check that run number is the correct runs range
      string runs;
      if( NULL!=(*it)->GetAttribute("runs",runs) )
      {
        DaqMap m(**it);
        m.SetParent(NULL);
        if( !m.TestRunsRange(run) )
          continue;
      }

      if( (*it)->GetName()=="Map" )
        continue;

      if( (*it)->GetName()=="ChipCreation" || (*it)->GetName()=="ChipsCreation" )
      {
        // must be considered by DaqOption.
        if( (*it)->GetParent()==NULL ||
            ((*it)->GetParent()->GetName()!="DaqOption" &&
             (*it)->GetParent()->GetName()!="DaqOptions") )
            throw Exception("Chip::ReadMaps(): section \"ChipsCreation\" must be inside \"DaqOptions\"");
        continue;
      }

      if( (*it)->GetName()=="DaqOption" || (*it)->GetName()=="DaqOptions" )
      {
        daq_option.ReadXML(**it);
        continue;
      }

      if( (*it)->GetName()=="ChipSADC_HT" )
      {
        daq_option.ReadChipSADC_HT(**it);
        continue;
      }
      
      if( (*it)->GetBody().empty() )
      {
        Exception("Chip::ReadMaps(): file \"%s\" contains an empty mapping section.\n",path.c_str()).Print();
        continue;
      }

      // Proceed all body lines
      ObjectXML o(**it);
      o.SetParent(NULL);
      for( list<string>::const_iterator s=(*it)->GetBody().begin(); s!=(*it)->GetBody().end(); s++ )
      {
        // Create an ObjectXML with only ONE string in the body.
        o.GetBody().clear();
        o.GetBody().push_back(*s);

        // Create a new map
        Map *m = NULL;
        
        try
        {
            m = Map::Create(o);
        }
        catch( const std::exception &e )
        {
          cerr << "Chip::ReadMaps: std::exception "<< e.what() << " in Map::Create while treating line: " << *s <<"\n";
        }
        catch( const char *st )
        {
          cerr << "Chip::ReadMaps: exception "<< st << " in Map::Create while treating line: " << *s <<"\n";
        }
        catch( ... )
        {
          cerr << "Chip::ReadMaps: Unknown exception in Map::Create while treating line: " << *s <<"\n";
          ok = false;
        }
     
        if( m==NULL )
        {
          Exception("Chip::ReadMaps(): Failed to load the map \"%s\" from the file %s",(*it)->GetName().c_str(),path.c_str()).Print();
          ok = false;
          continue;
        }

        RegisterMap(m,maps,dets,daq_option);
        delete m;
      }
    }
  }
  else
    throw Exception("Chip::ReadMaps(): not a file or directory: \"%s\"",path.c_str());

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

void Chip::RegisterMap(Map *m,Maps &maps,vector<string> &dets,DaqOption &options)
{
  // Fill detector's list
  size_t det_n;
  for( det_n=0; det_n<dets.size(); det_n++ )
    if( dets[det_n]==m->GetDetID().GetName() )
      break;

  if( det_n==dets.size() )  // New detector was found.
    dets.push_back(m->GetDetID().GetName());
  m->GetDetID() = DetID(m->GetDetID().GetName(),det_n);

  // Register the sourceID
  maps.GetDetSrcIDes()[m->GetDetID()].insert(m->GetSourceID());
  options.GetEventSrcIDs().insert(m->GetSourceID());

  // Add newly created map to the list
  m->AddToMaps(maps,options);
}

////////////////////////////////////////////////////////////////////////////////

Chip::Map::Map(const ObjectXML &o)
:   DaqMap        (o),
    id            (DetID("Unknown",uint32(-1))),
    chanF         (1),
    chanL         (0),
    chanS         (1),
    wireF         (1),
    wireL         (0),
    wireS         (1)
{
    if( o.GetBody().size()!=1 )
        throw Exception("Chip::Map::Map(): internal error");

    map_line = o.GetBody().front();
    
    vector<string> parts;
    
    // Convert all hexadecimal numbers (which starts from "0x") internally to decimal ones
    istringstream s(map_line.c_str());

    string ss;
    while( s>>ss)
    {
        if( ss.substr(0,2)=="0x" )
        {
            int i;
            if( 1!=sscanf(ss.c_str(),"%x",&i) )
                throw Exception("Chip::Map::Map(): bad hexadecimal number in line \"%s\"",map_line.c_str());
            char ss[33];
            sprintf(ss,"%d",i);
            dec_line += ss;
        }
        else
            dec_line += ss;

        dec_line += ' ';
    }
    map_line=dec_line;

    // Set the default version to '1'
    if( version==0 )
        version=1;
}

////////////////////////////////////////////////////////////////////////////////

Chip::Map* Chip::Map::Create(const ObjectXML &o)
{
  Chip::Map *m=NULL;

  if( o.GetName()=="ChipF1" )
    m = new ChipF1::Map(o);
  else
    if( o.GetName()=="ChipSinica" )
      m = new ChipSinica::Map(o);
  else
  if( o.GetName()=="ChipAPV" )
    m = new ChipAPV::Map(o);
  else
  if( o.GetName()=="ChipAPVRICH" )
    m = new ChipAPVRICH::Map(o);
  else
  if( o.GetName()=="ChipADC" )
    m = new ChipADC::Map(o);
  else
  if( o.GetName()=="ChipSADC" )
    m = new ChipSADC::Map(o);
  else
  if( o.GetName()=="ChipGassiplex" )
    m = new ChipGassiplex::Map(o);
  else
  if( o.GetName()=="Scaler" )
    m = new Scaler::Map(o);
  else
  if( o.GetName()=="ChipHotGeSiCA" )
    m = new ChipHotGeSiCA::Map(o);
  else
  if( o.GetName()=="ChipGandalf" )
    m = new ChipGandalf::Map(o);

  if( m!=NULL )
    m->SetParent(NULL);

  return m;
}

////////////////////////////////////////////////////////////////////////////////

int32 Chip::Map::CalcWire(int32 chan) const
{
  if( (chan-GetChanF())*GetChanS()<0 ||
      (chan-GetChanL())*GetChanS()>0 ||
      (chan-GetChanF())%GetChanS()!=0 )
    throw int(0);

  // transfer channel to wire
  return GetWireF()+(chan-GetChanF())/GetChanS()*GetWireS();
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Map::Check(void)
{
  if( chanS==0 )
    throw Exception("Chip::Map::Check(): %s: channel_step=0  %s",GetName().c_str(),map_line.c_str());

  int chanN=1+(chanL-chanF)/chanS;

  if( (chanL-chanF)%chanS!=0 || chanN<=0 )
    throw Exception("Chip::Map::Check(): %s: bad channels range or bad step: %s",
                    GetName().c_str(),map_line.c_str());

  // -----

  if( wireS==0 )
    throw Exception("Chip::Map::Check(): %s: wire_step=0  %s",GetName().c_str(),map_line.c_str());

  if( wireS<0 && wireF<wireL )
  {
    int w=wireF;
    wireF=wireL;
    wireL=w;
  }

  int wireN=1+(wireL-wireF)/wireS;

  if( (wireL-wireF)%wireS!=0 || chanN<=0 )
    throw Exception("Chip::Map::Check(): %s: bad wires range or bad step: %s",
                    GetName().c_str(),map_line.c_str());

  // -----

  if( chanN!=wireN )
    throw Exception("Chip::Map::Check(): %s: the number of channels and wires are not equal: chanN=%d wireN=%d : %s",
                    GetName().c_str(),chanN,wireN,map_line.c_str());
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Map::Print(ostream &o,const string &prefix) const
{
  DaqMap::Print(o,prefix);
  o<<prefix;
  char s[1111];
  sprintf(s,"  det=\"%s\" n=%d (chanF=%d chanL=%d chanS=%d chanN=%d)   (wireF=%d wireL=%d wireS=%d wireN=%d)\n",
               id.GetName().c_str(),id.GetNumber(),
               chanF,chanL,chanS,(chanL-chanF)/chanS+1,
               wireF,wireL,wireS,(wireL-wireF)/wireS+1);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Maps::Clear(void)
{
  for( Maps::const_iterator it=begin(); it!=end(); it++ )
    delete it->second;
  multimap<DataID,Digit*>::clear();
  det_srcIDes.clear();
  channels.clear();
}

////////////////////////////////////////////////////////////////////////////////

Chip::Maps::~Maps(void)
{
    Clear();
}

////////////////////////////////////////////////////////////////////////////////

uint32 Chip::Maps::GetWires(const DetID &id) const
{
    map<DetID,uint16>::const_iterator it=channels.find(id);
    return it==channels.end() ? 0 : it->second;
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Maps::GetSrcIDs(const string &pattern,set<uint16> &srcIDs) const
{
    for( map<DetID,set<uint16> >::const_iterator it=det_srcIDes.begin();
         it!=det_srcIDes.end(); it++ )
        if( string_match(it->first.GetName(),pattern) )
            GetSrcIDs(it->first,srcIDs);
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Maps::GetSrcIDs(const DetID& id,set<uint16> &srcIDs) const
{
    map<DetID,set<uint16> >::const_iterator it=det_srcIDes.find(id);
    if( it==det_srcIDes.end() )
        return;
    
    srcIDs.insert(it->second.begin(),it->second.end());
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Maps::PrintSrcIDs(const DetID& id) const
{
    set<uint16> srcIDs;
    GetSrcIDs(id,srcIDs);
    cout << "Source IDes for \"" << id.GetName() << "\":";
    for( set<uint16>::const_iterator it=srcIDs.begin(); it!=srcIDs.end(); it++ )
        cout << "  " << *it;
    cout << "\n";
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Maps::Print(ostream &o,const string &prefix) const
{
  for( Maps::const_iterator it=begin(); it!=end(); it++ )
  {
    o << hex << it->first << "   ";
    it->second->Print(o,prefix);
  }
}

////////////////////////////////////////////////////////////////////////////////

Chip::CreationRule::CreationRule(const string &ss)
{
  istringstream s(ss.c_str());

  s >> Dtype >> Did >> Stype >> Sid >> format >> chip_name >> srcID;

  if( s.fail() )
  {
    Exception e("CreationRule::CreationRule(): bad format in line: %s",ss.c_str());
    cerr<<e.what();
    throw e;
  }
}

////////////////////////////////////////////////////////////////////////////////

bool Chip::CreationRule::Accept(int32 _Dtype,int32 _Did,int32 _Stype,int32 _Sid,int32 _format) const
{
  return
    (Dtype ==-1 || Dtype ==_Dtype ) &&
    (Did   ==-1 || Did   ==_Did   ) &&
    (Stype ==-1 || Stype ==_Stype ) &&
    (Sid   ==-1 || Sid   ==_Sid   ) &&
    (format==-1 || format==_format);
}

////////////////////////////////////////////////////////////////////////////////

bool Chip::CreationRule::Accept(int32 _Stype,int32 _Sid,int32 _format) const
{
  return
    (Stype ==-1 || Stype ==_Stype ) &&
    (Sid   ==-1 || Sid   ==_Sid   ) &&
    (format==-1 || format==_format);
}

////////////////////////////////////////////////////////////////////////////////

void Chip::CreationRule::Print(ostream &o,const string &prefix) const
{
  o<<prefix;
  char s[222];
  sprintf(s,"chip=%14s Dtype=%4d Did=%4d Stype=%4d Sid=%4d",
             chip_name.c_str(), Dtype, Did, Stype, Sid);
  o<<s;

  o<<"  format=";
  if( format==-1 )
    o<<"any";
  else
    bits_print(o,format,0,8," %3d");

  o<<"\n";
}

////////////////////////////////////////////////////////////////////////////////

void Chip::Digit::Print(ostream &o,const string &prefix) const
{
  o << prefix << "Digit from " << GetDetID() << "\n";
}

////////////////////////////////////////////////////////////////////////////////

Chip::ErrorWord::ErrorWord(uint32 word,uint32 error_marker,uint32 source,Chip &chip,DaqOption &o)
{
  data.all = word;
  is_error = GetMarker()==error_marker;

  if( IsError() )
  {
    if( 0!=GetErrorCode() )
      chip.AddError(DaqError(DaqError::Type(DaqError::TDC_ERR1+GetErrorCode()-1),
                             DaqError::SOURCE_ID,source,
                             DaqError::PORT,GetPort(),
                             DaqError::GEO_ID,GetGeographicID()),o);

    if(  0==GetSpecialErrorCode() )
    {
      if (0==GetErrorCode())
        chip.AddError(DaqError(DaqError::TDC_WRGEOID,DaqError::SOURCE_ID,source,DaqError::PORT,GetPort(),DaqError::GEO_ID,GetGeographicID()),o);
    }
    else if (2==GetSpecialErrorCode())
      chip.AddError(DaqError(DaqError::TDC_SERR2,DaqError::SOURCE_ID,source,DaqError::PORT,GetPort(),DaqError::GEO_ID,GetGeographicID()),o);
    else if (3==GetSpecialErrorCode())
      chip.AddError(DaqError(DaqError::TDC_SERR3,DaqError::SOURCE_ID,source,DaqError::PORT,GetPort(),DaqError::GEO_ID,GetGeographicID()),o);
    else if (1==GetSpecialErrorCode() && (GetErrorCode()==0 || GetErrorCode()==3 || GetErrorCode()==7 )) 
      chip.AddError(DaqError(DaqError::TDC_SERR1,DaqError::SOURCE_ID,source,DaqError::PORT,GetPort(),DaqError::GEO_ID,GetGeographicID()),o);  
  }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
