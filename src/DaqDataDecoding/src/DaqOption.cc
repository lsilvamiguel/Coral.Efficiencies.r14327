#include <cstdio>
#include <cstdlib>
#include "DaqOption.h"
#include "ChipSADC_huffman.h"

using namespace std;

namespace CS {

namespace {

bool get_bool(const string &s)
{
    char opt[5+1];

    if( 1!=sscanf(s.c_str(),"%5s",opt) )
        throw Exception("DaqOption::Set(): bad option %s",s.c_str());

    if( 0==strcasecmp(opt,"yes") || 0==strcasecmp(opt,"true") || 0==strcasecmp(opt,"1") )
        return true;
    else
        if( 0==strcasecmp(opt,"no") || 0==strcasecmp(opt,"false") || 0==strcasecmp(opt,"0") )
            return false;
        else
            throw Exception("DaqOption::Set(): unknown value %s",s.c_str());
}

void get_set(const string &s,set<uint16> &myset)
{
    istringstream str(s.c_str());

    string run;  // Examples:   "1"   "2-5"
    while( str>>run )
    {
        int runF,runL;
        size_t p = run.find('-');
        if( p!=string::npos )
        {
            runF=atoi(run.substr(0  ,p           ).c_str());
            runL=atoi(run.substr(p+1,string::npos).c_str());
        }
        else
            runF=runL=atoi(run.c_str());

        if( runF>runL || runF<0 )
            throw Exception("get_set():  NFirst>NLast or negative numbers:  NFirst=%d  NLast=%d",
                             runF,runL);

        for( int i=runF; i<=runL; i++ )
            myset.insert(i);
    }
}

}

////////////////////////////////////////////////////////////////////////////////

DaqOption::DaqOption(void)
:   DaqMap(ObjectXML("")),
    use_ERROR_ECC_checksum(true),
    use_TDC_WRPLL(true),
    use_SLinkFormat(true),
    fixSLinkMultiplexerSize(true),
    online_filter_srcID(-1),
    chip_SADC_ht(NULL)
{
    ObjectXML o("DaqOption");
    o.GetAttributes().push_back(ObjectXML::Attribute("runs","1-999999"));
    *this=DaqOption(o);
}

////////////////////////////////////////////////////////////////////////////////

DaqOption::DaqOption(const ObjectXML &o)
:   DaqMap(o),
    use_ERROR_ECC_checksum(true),
    use_TDC_WRPLL(true),
    use_SLinkFormat(true),
    fixSLinkMultiplexerSize(true),
    online_filter_srcID(-1),
    chip_SADC_ht(NULL)
{
    Clear();
    
    ReadXML(o);
}

////////////////////////////////////////////////////////////////////////////////

void DaqOption::ReadXML(const ObjectXML &o)
{
    bool master_time_found=false;

    if( o.GetName()!="DaqOption" && o.GetName()!="DaqOptions" )
        throw Exception("DaqOption::ReadXML(): bad ObjectXML: \"%s\"",o.GetName().c_str());

    for( list<ObjectXML::Attribute>::iterator it=GetAttributes().begin(); it!=GetAttributes().end(); it++ )
        if( it->first=="use_ERROR_ECC_checksum" )
            use_ERROR_ECC_checksum = get_bool(it->second);
        else
        if( it->first=="use_TDC_WRPLL" )
            use_TDC_WRPLL = get_bool(it->second);
        else
        if( it->first=="use_SLinkFormat" )
            use_SLinkFormat = get_bool(it->second);
        else
        if( it->first=="FixSLinkMultiplexerSize" )
            fixSLinkMultiplexerSize = get_bool(it->second);
        else
        if( it->first=="srcID_scan" )
        {
            try
            {
                get_set(it->second,GetSrcIDScan());
            }
            catch( const std::exception &e )
            {
                cerr << e.what() << "\n";
            }
        }

    for( list<ObjectXML*>::const_iterator it=o.GetChildren().begin(); it!=o.GetChildren().end(); it++ )
        if( (*it)->GetName()=="ChipCreation" || (*it)->GetName()=="ChipsCreation" )
            for( list<string>::const_iterator l=(*it)->GetBody().begin(); l!=(*it)->GetBody().end(); l++ )
                chip_creation_rules.push_back(*l);
        else
        if( (*it)->GetName()=="Error" )
        {
            try
            {
                string s;
                if( NULL==(*it)->GetAttribute("name",s) )
                {
                    Exception e("DaqOption::ReadXML(): missing attribute \"name\" in the definition of <Error>");
                    e.Print(cerr);
                    continue;
                }

                DaqErrorType &e=DaqErrorType::GetDaqErrorType(s);

                int level;
                if( NULL!=(*it)->GetAttribute("level",level) )
                    e.SetSeverityLevel(DaqError::SeverityLevel(level));

                int action;
                if( NULL!=(*it)->GetAttribute("action",action) )
                    e.SetAction(DaqError::Action::Type(action));
            }
            catch( const std::exception &e )
            {
                cerr << e.what() << "\n";
            }
        }
        else
        if( (*it)->GetName()=="TriggersIgnore" )
        {
            string s;
            if( NULL!=(*it)->GetAttribute("bit",s) )
            {
                istringstream str(s.c_str());
                
                uint16 bit;
                while(str>>bit)
                    trigger_bits_ignore.insert(bit);
            }

            if( NULL!=(*it)->GetAttribute("name",s) )
            {
                istringstream str(s.c_str());
                
                string name;
                while(str>>name)
                {
                    try
                    {
                        trigger_bits_ignore.insert(GetTriggerBit(name));
                    }
                    catch(...)
                    {
                        printf("TriggersIgnore: unknown trigger %s\n",name.c_str());
                    }
                }
            }
        }
        else
        if( (*it)->GetName()=="Trigger" )
        {
            const ObjectXML &xml = **it;    // short name
            for( list<ObjectXML::Attribute>::const_iterator it=xml.GetAttributes().begin();
                 it!=xml.GetAttributes().end(); it++ )
            {
                // Get the bit number.
                char *end=NULL;
                int n = strtol(it->second.c_str(), &end, 10);
                if( *end!=0 || n<0 )
                    throw Exception("DaqOption:: bad trigger bit setting: \"%s\"=\"%s\"",
                                     it->first.c_str(),it->second.c_str());

                Trigger t;
                t.SetBit(n);
                t.SetName(it->first);

                // Be sure that we do not have two triggers with the same bit settings.
                if( triggers.count(t.GetBit())>0 )
                    throw Exception("DaqOption::ReadXML(): second trigger with bit %d",t.GetBit());
                triggers[t.GetBit()] = t;
            }
        }
        else
        if( (*it)->GetName()=="TriggerSettings" )
        {
            Trigger t(**it);

            if( triggers.count(t.GetBit())>0 )
                if( triggers[t.GetBit()].GetName()!=t.GetName() )
                    throw Exception("Trigger name conflict for the trigger bit %d:  \"%s\"!=\"%s\"",
                                    t.GetBit(),triggers[t.GetBit()].GetName().c_str(),t.GetName().c_str());

            triggers[t.GetBit()] = t;
        }
        else
        if( (*it)->GetName()=="TriggerTimesCheck" )
        {
            unsigned buffer_size,check_entries_min,check_on_fill;
            float cut_ratio,allowed_shift;

            if( NULL!=(*it)->GetAttribute("buffer_size",buffer_size) )
                tt_stats.SetBufferSize(buffer_size);

            if( NULL!=(*it)->GetAttribute("buffer_cut_ratio",cut_ratio) )
                tt_stats.SetBufferCutRaitio(cut_ratio);

            if( NULL!=(*it)->GetAttribute("check_entries_min",check_entries_min) )
                tt_stats.SetCheckOnFill(check_entries_min);
                
            if( NULL!=(*it)->GetAttribute("check_on_fill",check_on_fill) )
                tt_stats.SetCheckEntriesMin(check_on_fill);

            if( NULL!=(*it)->GetAttribute("allowed_shift",allowed_shift) )
                tt_stats.SetAllowedShift(allowed_shift);
        }
        else
        if( (*it)->GetName()=="OnlineFilter" )
        {
            if( NULL==(*it)->GetAttribute("srcID",online_filter_srcID) )
                throw Exception("DaqOption::ReadXML(): bad online filter settings");
            GetEventSrcIDs().insert(online_filter_srcID);
        }
        else
        if( (*it)->GetName()=="TriggerTimeCorrection" )
        {
            for( list<string>::const_iterator s=(*it)->GetBody().begin(); s!=(*it)->GetBody().end(); s++ )
            {
                try
                {
                    char ss[1111], *end=NULL, stu[1111]="0.12892312,0.06446156";

                    float t;
                    
                    if( 2>sscanf(s->c_str(),"%g %s %s",&t,ss,stu) )
                        throw Exception("DaqOption::ReadXML(): bad TriggerTimeCorrection: %s",s->c_str());

                    // Get the bit number.
                    int m = strtol(ss, &end, 10);
                    if(*end==0)
                        throw "Bad format for TriggerTimeCorrection!";
                    m = GetTriggerBit(ss);

                    triggers[m].SetCorrection(t);
                    triggers[m].SetBit(m);
                }
                catch(const std::exception &e)
                {
                    cerr << e.what() << "\n";
                }
            }
        }
        else
        if( (*it)->GetName()=="LastTriggerTicks" )
        {
            string s;
            if( NULL!=(*it)->GetAttribute("IgnoreSrcID",s) )
            {
                istringstream str(s.c_str());
                
                int srcID;
                while(str>>srcID)
                    last_trigger_ticks_srcID_ignore.insert(srcID);
            }
            if( NULL!=(*it)->GetAttribute("IgnoreLocalSrcID",s) )
            {
                istringstream str(s.c_str());
                
                int srcID;
                while(str>>srcID)
                    local_trigger_ticks_srcID_ignore.insert(srcID);
            }
        }
        else
        if( (*it)->GetName()=="MasterTime" )
        {
            master_time_found=true;
            string s;
            if( NULL!=(*it)->GetAttribute("srcID",s) )
            {
                istringstream ss(s.c_str());
                uint16 n;
                while( ss>>n )
                    master_time_srcID.insert(n);
            }
            if( NULL!=(*it)->GetAttribute("no_check",s) )
            {
                istringstream ss(s.c_str());
                uint16 n;
                while( ss>>n )
                    master_time_srcID_dont_check.insert(n);
            }
        }
        else
        {
            Exception("DaqOption::ReadXML():WW: Unrecognized object \"%s\"",(*it)->GetName().c_str()).Print();
        }
        
    // automatically add catch number 2 to the MasterTime, if no option is found!
    if( !master_time_found )
        master_time_srcID.insert(2);

    if(0)
        for( map<uint16,Trigger>::const_iterator tr=triggers.begin(); tr!=triggers.end(); tr++ )
            tr->second.Print();
}

////////////////////////////////////////////////////////////////////////////////

DaqOption::~DaqOption(void)
{
    if( chip_SADC_ht!=NULL )
        free_tree(chip_SADC_ht);
}

////////////////////////////////////////////////////////////////////////////////

DaqOption & DaqOption::operator = (const DaqOption &o)
{
    if( chip_SADC_ht!=NULL )
        free_tree(chip_SADC_ht);

    if( o.GetChipSADC_HT()!=NULL )
    {
        chip_SADC_ht = copy_tree(o.GetChipSADC_HT());
        if( chip_SADC_ht==NULL )
            throw "DaqOption::operator = (): failed to copy the Huffman tree!";
    }
    else
        chip_SADC_ht = NULL;
    
    use_ERROR_ECC_checksum              = o.use_ERROR_ECC_checksum;
    use_TDC_WRPLL                       = o.use_TDC_WRPLL;
    use_SLinkFormat                     = o.use_SLinkFormat;
    fixSLinkMultiplexerSize             = o.fixSLinkMultiplexerSize;
    chip_creation_rules                 = o.chip_creation_rules;
    srcID_scan                          = o.srcID_scan;
    triggers                            = o.triggers;
    trigger_bits_ignore                 = o.trigger_bits_ignore;
    online_filter_srcID                 = o.online_filter_srcID;
    last_trigger_ticks_srcID_ignore     = o.last_trigger_ticks_srcID_ignore;
    calibrations                        = o.calibrations;
    master_time_srcID                   = o.master_time_srcID;
    master_time_srcID_dont_check        = o.master_time_srcID_dont_check;
    event_srcIDs                        = o.event_srcIDs;
    f1_cmc_data                         = o.f1_cmc_data;
    cmc_apv_always                      = o.cmc_apv_always;
    sadc_data_ver2                      = o.sadc_data_ver2;
    sadc_data_ver3                      = o.sadc_data_ver3;
    tt_config                           = o.tt_config;

    return *this;
}

////////////////////////////////////////////////////////////////////////////////

void DaqOption::SetCalibration(const Chip::DataID &data_id,const Chip::Calibration &calib)
{
    const Chip::Calibration *&c=calibrations[data_id];

    if( c!=NULL )
        throw Exception("DaqOption::SetCalibration(): calibration constants with given DataID already exists.");
    
    c = &calib;
}

////////////////////////////////////////////////////////////////////////////////

const Chip::Calibration* DaqOption::FindCalibration(const Chip::DataID &data_id) const
{
    map<Chip::DataID,const Chip::Calibration*>::const_iterator it=calibrations.find(data_id);

    if( it==calibrations.end() )
        return NULL;
    else
        return it->second;
}

////////////////////////////////////////////////////////////////////////////////

void DaqOption::ReadChipSADC_HT(const ObjectXML &o)
{
    std::string hufftree_str;
  
    for( list<string>::const_iterator s=o.GetBody().begin(); s!=o.GetBody().end(); s++ )
    {
        hufftree_str += *s + "\n";
    }
   
    if( chip_SADC_ht!=NULL )
         free_tree(chip_SADC_ht);
    chip_SADC_ht = hufftree_read_from_string(hufftree_str.c_str());

    if( !chip_SADC_ht)
        throw "DaqOption::ReadChipSADC_HT(): failed reading huffman tree!";
}

////////////////////////////////////////////////////////////////////////////////

bool DaqOption::NeedScanSrcID(uint16 srcID) const
{
    return GetSrcIDScan().size()==0 || GetSrcIDScan().count(srcID)>0;
}

////////////////////////////////////////////////////////////////////////////////

Trigger & DaqOption::GetTrigger(uint16 bit)
{
    map<uint16,Trigger>::iterator m=triggers.find(bit);
    if( m==triggers.end() )
        throw Exception("DaqOption::GetTrigger(): there is no trigger with bit %d",bit);
    
    return m->second;
}

////////////////////////////////////////////////////////////////////////////////

unsigned DaqOption::GetTriggerBit(const string &name) const
{
    for( map<uint16,Trigger>::const_iterator m=triggers.begin(); m!=triggers.end(); m++ )
        if( m->second.GetName()==name )
            return m->second.GetBit();

    throw Exception("DaqOption::GetTriggerBit(): there is no trigger with the name \"%s\"",name.c_str());
}


////////////////////////////////////////////////////////////////////////////////

void DaqOption::Clear(void)
{
    use_ERROR_ECC_checksum  = true;
    use_TDC_WRPLL           = true;
    use_SLinkFormat         = true;
    fixSLinkMultiplexerSize = false;
    srcID_present_ports.clear();
}

////////////////////////////////////////////////////////////////////////////////

void DaqOption::Print(ostream &o,const string &prefix) const
{
    DaqMap::Print(o,prefix);
    o << prefix << "use_ERROR_ECC_checksum ... " << (Use_ERROR_ECC_checksum()?"yes":"no") << "\n";
    o << prefix << "use_TDC_WRPLL ............ " << (Use_TDC_WRPLL()         ?"yes":"no") << "\n";
    o << prefix << "use_SLinkFormat .......... " << (Use_SLinkFormat()       ?"yes":"no") << "\n";
    for( vector<Chip::CreationRule>::const_iterator it=chip_creation_rules.begin();
         it!=chip_creation_rules.end(); it++ )
      it->Print(o,prefix+"  chips creation rule:  ");

    o << prefix << "srcID to scan:";
    for( set<uint16>::const_iterator it=srcID_scan.begin(); it!=srcID_scan.end(); it++ )
        o << " " << (*it);
    o << "\n";
}

////////////////////////////////////////////////////////////////////////////////

void DaqOption::AddSrcIDPresentPorts(uint16 srcID, uint32 port) {
    srcID_present_ports[srcID].insert(port);
}

////////////////////////////////////////////////////////////////////////////////

} // namesapce CS
