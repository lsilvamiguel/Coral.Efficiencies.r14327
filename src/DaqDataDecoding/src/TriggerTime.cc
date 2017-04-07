#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <sstream>

#include "TriggerTime.h"
#include "ChipF1.h"
#include "ChipSinica.h"
#include "ChipGandalf.h"
#include "DaqError.h"
#include "DaqOption.h"
#include "Scaler.h"

/*! \file TriggerTime.cc
  \author Alexander Zvyagin
*/

namespace CS {

  using namespace std;

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTime::InitAndClear(const TriggerTimeConfig &ttc)
  {
    config = &ttc;

    ltt                    = -1;
    event_trigger          = 0;
    trigger_of_master_time = NULL;
    // Users wanting to avoid using bogus TCS phase must consult DaqErrors.
    TCS_phase              = -1e6;
    // Users wanting to avoid using bogus TimeInSpill must consult DaqErrors.
    TimeInSpill            = -1.;

    trigger_mask_data.clear();
    TCS_phase_data.clear();
    shifts.clear();

    data.clear();
    for( std::map<int,TriggerTimeConfig::TTC>::const_iterator it=ttc.tt.begin(); it!=ttc.tt.end(); it++ )
      {
        assert( it->first == it->second.index );
        data[it->second.index] = Data(&it->second);
      }
  }

  ////////////////////////////////////////////////////////////////////////////////

  bool TriggerTime::SetLastTriggerTicks(int t)
  {
    if( ltt==-2 )
      return true;    // at least one wrong measurement:
    // don't trust LastTriggerTicks at all (for the event)!

    if( ltt==-1 )
      {
        // New event, first measurement
        ltt=t;
        return true;
      }
    
    if( ltt!=t )
      {
        // That is bad: wrong measurement!
        //printf("TriggerTime::SetLastTriggerTicks(): %d!=%d\n",ltt,t);
        ltt=-2;
        return false;
      }

    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////

#define TT_DEBUG 0

  void TriggerTime::Decode(const DaqOption &options,DaqErrors &errors)
  { 
    if( config==NULL )
      throw Exception("TriggerTime::Decode(): config==NULL");

    assert(event_trigger==0);
    vector<uint16> trigger_bits;

    // Check the trigger mask
    if( config->trigger_mask_sources_max>0 )
      {
        if( trigger_mask_data.size()<config->trigger_mask_sources_min
            ||
            trigger_mask_data.size()>config->trigger_mask_sources_max )
            
	  errors.Add(DaqError(DaqError::WRONG_HITS_IN_TRIGGER_MASK,
			      DaqError::SOURCE_ID,config->trigger_mask_srcID,
			      DaqError::PORT,config->trigger_mask_port,
			      DaqError::COUNTER,trigger_mask_data.size()));
        
        for( vector<const ChipF1::Digit*>::const_iterator it=trigger_mask_data.begin();
             it!=trigger_mask_data.end(); it++ )
	  {
            uint16 bit = (*it)->GetChannel();
            event_trigger |= 1<<bit;
            trigger_bits.push_back(bit);
	  }
      }
    else
      {
        // Print a warning message once...

        static bool init=true;
        if( init )
	  {
            init=false;
            Exception("TriggerTime::Decode(): Trigger mask is not configured!").Print();
	  }
      }

    for( std::map<int,Data>::iterator it=data.begin(); it!=data.end(); it++ )
      {
        Data &d=it->second;
        if (d.GetConfig().srcID!=800)
	  d.decoded=false; //this is set to false for unknown reason

        // Do not decode if:
        //   - there are no data
        //   - it was already decoded
        //   - ID is empty
        //	 - it is a GandalfDigit that was already decoded
        if( d.data.empty() || d.decoded || d.GetConfig().id.GetName().length()==0 )
	  continue;

        if( d.data.size()!=(unsigned int)d.GetConfig().channels )
	  {
            errors.Add(DaqError(DaqError::TT_WRONG_HITS_IN_MASTER_TIME,
                                DaqError::SOURCE_ID,d.GetConfig().srcID,
                                DaqError::PORT,d.GetConfig().port,
                                DaqError::COUNTER,d.data.size(),
                                DaqError::COUNTER,d.GetConfig().channels,
                                DaqError::VALUE,d.GetConfig().time_unit));
            continue;
	  }

#if TT_DEBUG
        printf("%s: Data index=%d   ",d.config->id.GetName().c_str(),d.config->index);
#endif
        uint32 ref=d.data.begin()->second;  // This is the refernce time
        double sum=0, sum_sqr=0;

        for( map<Chip::DataID,uint16>::const_iterator it=d.data.begin();
             it!=d.data.end(); it++ )
	  {
            double diff=ChipF1::TimeDifference(it->second,ref,d.GetConfig().overolling,0);
#if TT_DEBUG
            cout << it->second << "(diff=" << diff << ") ";
#endif
            sum     += diff;
            sum_sqr += diff*diff;
	  }

        d.time = sum/d.data.size();
        bool sigma_ok=true;
        double sigma=0;
        
        if( d.data.size()>1 )
	  {
            // fabs() is used to overcome the precision problem (negative argument for sqrt())
            sigma = sqrt( fabs(sum_sqr/(d.data.size()-1) - double(d.time)*d.time) );
            if( sigma>d.GetConfig().sigma )
	      {
                errors.Add(DaqError(DaqError::TT_BIG_SIGMA,
                                    DaqError::SOURCE_ID,d.GetConfig().srcID,
                                    DaqError::PORT,d.GetConfig().port,
                                    DaqError::VALUE,sigma,DaqError::VALUE,d.GetConfig().time_unit));
                sigma_ok=false;
	      }
	  }
        
        d.time += ref;  // add the reference
        d.decoded=sigma_ok;
      }
    
    CorrectMasterTime(options,errors);
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTime::CorrectMasterTime(const DaqOption &options,DaqErrors &errors)
  {
    const bool debug = 0;

    if( config->trigger_mask_sources_max==0 )
      {
        static Exception e("TriggerTime::CorrectMasterTime(): can not correct master time: missing configuration!");
        return;
      }

    // We have to fill this structure:
    multimap<float,const Trigger*> triggers;

    // Work only with the 'base' resolution.
    const Data *ttd = Find(config->trigger_mask_TT_index);
    if( ttd==NULL )
      throw Exception("TriggerTime::CorrectMasterTime(): TT index %d is not found",config->trigger_mask_TT_index);

    if( !ttd->decoded )
      throw Exception("TriggerTime::CorrectMasterTime(): TT index %d was not decoded",config->trigger_mask_TT_index);

    const float &mt = ttd->time; // Master time

    for( vector<const ChipF1::Digit*>::const_iterator it=trigger_mask_data.begin();
         it!=trigger_mask_data.end(); it++ )
      {
        uint16 bit=(*it)->GetChannel();
        uint32 time=(*it)->GetTime();
        float tt_mt_diff = ChipF1::TimeDifference(time,mt,ttd->GetConfig().overolling,0);

        shifts.push_back(TriggerTimeShift(bit,event_trigger,mt,tt_mt_diff));

        const Trigger &trig = options.GetTrigger(bit);
        const float time_diff = (tt_mt_diff-trig.GetOffset())*ttd->GetConfig().time_unit;
        triggers.insert( pair<float,const Trigger*>(time_diff,&trig) );
      }

    // OK, everything is ready....
    
    if( trigger_mask_data.size()>=1 && trigger_mask_data.size()<=2 )
      {
        TriggerTimeStat &tt_stats = const_cast<TriggerTimeStat&>(options.GetStatTT());

        // Fill some statistics

        uint16 bit1  = trigger_mask_data[0]->GetChannel();
        uint32 time1 = trigger_mask_data[0]->GetTime();

        if( trigger_mask_data.size()==1 )
	  {
            float tt_mt_diff = ChipF1::TimeDifference(time1,mt,ttd->GetConfig().overolling,0);
	    tt_stats.Add(bit1,tt_mt_diff);
	  }
        
        if( trigger_mask_data.size()==2 )
	  {
            uint16 bit2  = trigger_mask_data[1]->GetChannel();
            uint32 time2 = trigger_mask_data[1]->GetTime();
            float diff = ChipF1::TimeDifference(time1-time2,0,ttd->GetConfig().overolling,0);
            tt_stats.Add(bit1,bit2,diff);
	  }
        
        // Check trigger time shifts.
        CheckShiftTT(tt_stats,options,errors);
      }
    
    float cor=0; // Master time correction in [ns]
    if( triggers.size()>0 )
      {
        map<float,const Trigger*>::const_iterator fired_trigger=triggers.end();
        GetFiredTrigger(triggers,fired_trigger,cor);
        // Now 'cor' is a time difference (in [ns]) between the fired trigger
        // and the trigger which we want to use as a master time.
        // This correction 'cor' must be non negative!
        assert(cor>=0);
        assert(fired_trigger!=triggers.end());
        
        trigger_of_master_time = fired_trigger->second;
        
        if( debug )
	  if( cor!=0 )
	    printf("Correction from triggers:   %g\n",cor);

        // Apply a trigger time correction from the fired trigger.
        cor += fired_trigger->second -> GetCorrection();
      }
    else
      {
        Exception("TriggerTime::Decode(): no triggers?!");
      }

    // Now the correction to be applied for all master times
    if( cor!=0 )
      for( std::map<int,Data>::iterator it=data.begin(); it!=data.end(); it++ )
        {
	  Data &d = it->second;
	  d.time += cor/d.GetConfig().time_unit;  // compensation for the trigger type
        }
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTime::GetFiredTrigger(const multimap<float,const Trigger*> &triggers,
				    multimap<float,const Trigger*>::const_iterator &trig_fired,
				    float &mt_correction)
  {
    const bool debug = 0;

    if( triggers.size()==0 )
      throw Exception("TriggerTime::GetFiredTrigger(): no triggers!");
    
    if( triggers.size()==1 )
      {
        trig_fired = triggers.begin();
        return;
      }

    // ============================================================================

    if(debug)
      {
        // This is a debug printout.
        for( multimap<float,const Trigger*>::const_iterator it=triggers.begin(); it!=triggers.end(); it++ )
	  {
            printf("(%2d,%+6.1f) ",it->second->GetBit(),it->first);
	  }
        printf("\n");
      }

    trig_fired = triggers.begin();
    multimap<float,const Trigger*>::const_iterator current_trig=triggers.begin();

    for( current_trig++; current_trig!=triggers.end(); current_trig++)
      {
        const Trigger
	  &trig1 = *trig_fired->second,
	  &trig2 = *current_trig->second;

        const bool is_precise_trigger = trig1.GetPrecision()<config->time_precise;
        if( is_precise_trigger )
	  {
            // Precise trigger is found. There is no need to check other triggers.
            return;
	  }

        // Situation: the trigger is not precise!

        float
	  tdiff           = current_trig->first - trig_fired->first,   // in [ns]!
	  sigma1          = trig1.GetPrecision(),
	  sigma2          = trig2.GetPrecision(),
	  sigma           = sqrt(sigma1*sigma1+sigma2*sigma2),
	  diff_in_sigma   = tdiff/sigma;

        assert(tdiff>=0);

        if( diff_in_sigma > config->triggers_diff_sigma )
	  {
            if(debug)
	      printf("Big difference: %s-%s = %g sigmas\n",trig1.GetName().c_str(),trig2.GetName().c_str(),diff_in_sigma);

            return;
	  }

        // This is the situation: we have an unprecise trigger and time difference to a next trigger is small!
        if(debug)
	  printf("diff = %g ns (%g sigmas) with:  %s %s\n",tdiff,diff_in_sigma,trig1.GetName().c_str(),trig2.GetName().c_str());

        // Choose a most precise trigger
        if( current_trig->second->GetPrecision() < trig_fired->second->GetPrecision() )
	  {
            trig_fired = current_trig;

            if( tdiff>config->time_jitter )
	      mt_correction += tdiff;

            // We need to change the triggers order!!
            if(debug)
	      printf("Triggers order changed!  old:%s  new:%s   cor=%g\n",
		     trig_fired->second->GetName().c_str(),
		     current_trig->second->GetName().c_str(),
		     tdiff);
        
	  }
        else
	  if(debug)
	    printf("That must be: two unprecize triggers!\n");
      }
  }

  ////////////////////////////////////////////////////////////////////////////////

  bool TriggerTime::CheckShiftTT(const TriggerTimeStat &tt_stats,const DaqOption &options,DaqErrors &errors) const
  {
    bool check_ok = true;

    if( (tt_stats.GetFillsCounter()%tt_stats.GetCheckOnFill()) == 0 )
      {
        for( map<uint32,Stat>::const_iterator it=tt_stats.GetTriggersStat().begin(); it!=tt_stats.GetTriggersStat().end(); it++ )
	  {
            const Stat &stat = it->second;

            // Check do we have enough statistics for analysis or not.
            if( stat.Size()<tt_stats.GetCheckEntriesMin() )
	      continue;

            // Get the trigger
            const int mask = it->first;
            vector<int> bits=bits_from_mask(mask);
            if( bits.size()!=1 )
	      continue;  // it is not a pure trigger
            const unsigned int bit = bits[0];
            const Trigger &trig=options.GetTrigger(bit);

            // Do not check triggers for which precision is not set.
            if( trig.GetPrecision()<=0 )
	      continue;

            // Get mean value and sigma with 90% of events
            const float
	      mean  = stat.Mean(stat.GetCutRatio()),
	      sigma = stat.Sigma(stat.GetCutRatio());

            // Now check it!
            if( fabs(mean-trig.GetOffset())>(trig.GetPrecision()+tt_stats.GetAllowedShift()) )
	      {
                check_ok = false;
                errors.Add(DaqError(DaqError::TT_MT_SHIFT_CHANGED,
                                    DaqError::SOURCE_ID,config->trigger_mask_srcID,
                                    DaqError::VALUE,bit));
                if( trig.IsVerbose()>0 )
		  {
                    printf("----------- Trigger bit problem with MT-TT shift. -----------\n");
                    trig.Print();
                    printf("stat_size=%u stat_mean=%g stat_sigma=%g trig_offset=%g trig_prec=%g diff=%g diff_sigma=%g\n",
			   stat.Size(),mean,sigma,trig.GetOffset(),trig.GetPrecision(),mean-trig.GetOffset(),
			   fabs(mean-trig.GetOffset())/trig.GetPrecision());
                    if( trig.IsVerbose()>100 )
		      {
                        printf("buffer(size %d):",stat.Size());
                        for( std::set<float>::const_iterator v=stat.GetBuffer().begin(); v!=stat.GetBuffer().end(); v++ )
			  printf(" %g",*v);
                        printf("\n");
		      }
                    printf("----------- End of the verbose trigger problem printout -----------\n");

		  }
	      }
	  }
      }
    
    return check_ok;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTime::Decode(const Chip::Digits &digits,const DaqOption &options,DaqErrors &errors)
  {
    typedef Chip::Digits::const_iterator m_it; // iterator type

    if( config->TCS_phase_DetID && config->TCS_phase_DetID!=DetID("") )
      {
        pair<m_it,m_it> m_range = digits.equal_range(config->TCS_phase_DetID);
        for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
	  {
            const ChipF1::Digit *f1 = dynamic_cast<const ChipF1::Digit*>(d_it->second);
	    if( f1==NULL )
	      throw Exception("TriggerTime::Decode(): digit is not originated from ChipF1!");
	   
            TCS_phase_data.push_back(f1);
	  }
      }
    
    if( TCS_phase_data.size()==0 )
      errors.Add(DaqError(DaqError::TCS_PHASE_NO_DIGITS,
			  DaqError::SOURCE_ID,config->TCS_phase_srcID,
			  DaqError::PORT,config->TCS_phase_port,
			  DaqError::CHANNEL,config->TCS_phase_channel));


    if( TCS_phase_data.size()>1 )
      errors.Add(DaqError(DaqError::TCS_PHASE_MANY_DIGITS,
			  DaqError::SOURCE_ID,config->TCS_phase_srcID,
			  DaqError::PORT,config->TCS_phase_port,
			  DaqError::CHANNEL,config->TCS_phase_channel));

    if( config->trigger_mask_DetID  && config->trigger_mask_DetID!=DetID("") )
      {
        pair<m_it,m_it> m_range = digits.equal_range(config->trigger_mask_DetID);
        for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
	  {
            const ChipF1::Digit *f1 = dynamic_cast<const ChipF1::Digit*>(d_it->second);
	    if( f1==NULL )
	      throw Exception("TriggerTime::Decode(): digit is not originated from ChipF1!");
	    

            // Check do we need to ignore some triggers?
            if( options.GetTriggerBitsIgnore().count(f1->GetChannel())>0 )
	      continue;
       
            trigger_mask_data.push_back(f1);
	  }
      }

    for( std::map<int,Data>::iterator it=data.begin(); it!=data.end(); it++ )
      {
        if( it->second.decoded )
	  throw Exception("TriggerTime::Decode(): attempt to add data to already decoded TT");

        const DetID &id = it->second.GetConfig().id;
        // get all digits for detector my_det, fast call
        pair<m_it,m_it> m_range = digits.equal_range(id);

        // loop on all digits for the given DetID
        for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
	  {
            const ChipF1::Digit *f1 = dynamic_cast<const ChipF1::Digit*>(d_it->second);
	    const ChipGandalf::DigitGADC *gandalf = dynamic_cast<const ChipGandalf::DigitGADC*>(d_it->second);
	    if( f1==NULL && gandalf==NULL )
	      throw Exception("TriggerTime::Decode(): digit is not originated from ChipF1 or ChipGandalf!");

	    // we only accept the first time measurement!
	    if (f1)
	      {
		if( it->second.data.find(f1->GetDataID()) == it->second.data.end() )
		  it->second.data[f1->GetDataID()] = f1->GetAmplitude();
	      }
	    // TODO: what does that mean ?? " we only accept the first time measurement! "
	    if (gandalf)
	      {  
		if( it->second.data.find(gandalf->GetDataID()) == it->second.data.end() ){
		  // only channel 0 has the triggertime; also, ignore frame digits
		  if (gandalf->getChannel()==0 && gandalf->getNumSamples()==0) {
		    it->second.time = gandalf->getTime();
		    it->second.decoded = true;
		  }
		}
	      }
	  }
      }

    // read TimeInSpill from scalers
    {
      typedef multimap<CS::DetID,CS::Chip::Digit*>::const_iterator mIt;
      const set<string> &tbnames = options.GetTTConfig().tis_tbnames;
      
      if ( ! tbnames.size() )
	throw Exception("Error: time_in_spill tbnames have not been defined! "
			"Please make sure to add a time_in_spill attribute which "
			"lists the tbnames of scalers to be used for time in spill "
			"measurement to the SCALER.xml mapping file!");

      long long sum = 0; unsigned count = 0;
      for ( set<string>::const_iterator tbname = tbnames.begin(); tbname != tbnames.end(); tbname++ ) {
	pair<mIt,mIt> mRange = digits.equal_range(*tbname);
	for (mIt dIt = mRange.first; dIt!=mRange.second; dIt++) {
	  CS::Scaler::Digit *sd = dynamic_cast<CS::Scaler::Digit*>(dIt->second);
	  assert ( sd != NULL );
	  if ( sd->GetChannel() != 33 )
	    continue;
	  // (Pseudo)-channel #33 records the time in spill on all modules (in
	  // fact it records the time of arrival of the TCS clock => hence
	  // small differences among the various modules. TimeInSpill is
	  // computed as average of all given tbnames for maximal stability
	  // and resilience against missing data.
	  sum += sd->GetValue();
	  count++;
	}
      }
      if ( count > 0 ) {
	TimeInSpill = (double) sum / count / CS::ChipF1::GetTCSFrequency();
      } else {
	errors.Add( DaqError(DaqError::TIS_NO_DATA) );
      }
    }

    Decode(options,errors);
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTime::DecodeAndSubtract(Chip::Digits &digits,const DaqOption &options,DaqErrors &errors)
  {
    // Create 
    list<ChipF1::Digit*> d;
    list<ChipSinica::Digit*> ds;
    list<ChipGandalf::DigitGADC*> dg;
    for( Chip::Digits::const_iterator it=digits.begin();
         it!=digits.end(); it++ )
      {
        ChipF1::Digit *f1 = dynamic_cast<ChipF1::Digit*>(it->second);
        if( f1!=NULL )
	  d.push_back(f1);
        ChipSinica::Digit *si = dynamic_cast<ChipSinica::Digit*>(it->second);
        if( si!=NULL )
	  ds.push_back(si);
        ChipGandalf::DigitGADC *gandalf = dynamic_cast<ChipGandalf::DigitGADC*>(it->second);
        if( gandalf!=NULL )
	  dg.push_back(gandalf);
      }

    Decode(digits,options,errors);

    // And now subtract the trigger time.
    for( list<ChipF1::Digit*>::iterator it=d.begin(); it!=d.end(); it++ )
      {
        const double tu=(*it)->GetTimeUnit();
        const Data *d = Find(tu);
        if( d==NULL )
	  {
            Exception("TriggerTime::DecodeAndSubtract(): bad time_unit=%g for det=\"%s\"",
                      (*it)->GetTimeUnit(),(*it)->GetDetID().GetName().c_str()).Print();
            continue;
	  }

        ChipF1::DataID id = (*it)->GetDataID();
        double MT=d->time;
     
        if( id.u.s.src_id!=d->config->srcID )
	  MT  += d->config->MT_shift/tu;

	TriggerTimeStat &tt_stats = const_cast<TriggerTimeStat&>(options.GetStatTT());

        double time = ChipF1::TimeDifference((*it)->GetAmplitude(),MT,d->GetConfig().overolling,(*it)->GetTimeReference());

	//            if( id.u.s.src_id!=d->config->srcID && d->config->MT_shift!=0 )
	//                 printf("srcID: %d    %s   cor=%g    t1=%g   t2=%g\n",
	//                             int(id.u.s.src_id),(*it)->GetDetID().GetName().c_str(),d->config->MT_shift,time,time2);

	(*it)->SetTimeDecoded(tu*time);
      }

    // now that the F1 digits have been corrected for the trigger time, the
    // TCS phase can be calculated
    if( TCS_phase_data.size() > 0 )
      {
        TCS_phase = TCS_phase_data.front()->GetTimeDecoded();
        if ( TCS_phase < -200. || TCS_phase > 200. )
          {
            errors.Add( DaqError(DaqError::TCS_PHASE_OUT_OF_RANGE,
                                 DaqError::VALUE, TCS_phase) );
          }
      }

    // and now subtract the trigger time for digits of type 'ChipSinica'
    for( list<ChipSinica::Digit*>::iterator it=ds.begin(); it!=ds.end(); it++ )
      {
        const double tu=(*it)->GetTimeUnit();

        // digits, in DC05 Chip with 1.07.... tme unit
        // on other chips triggertime may be double because it is averaged/adjusted before
        // but on sinica it is just a digit, like hit time, on front end's clock
        const double time = ChipSinica::TimeDifference((*it)->GetHitTime(), (*it)->GetTriggerTime(), (*it)->GetOverolling(), (*it)->GetTimeReference());

        (*it)->SetTimeDecoded(tu*time + TCS_phase);
      }


    for( list<ChipGandalf::DigitGADC*>::iterator it=dg.begin(); it!=dg.end(); it++ )
      {
        const double tu=(*it)->GetTimeUnit();
        // TODO: if Gandalf has no time unit --> what do we do ?
        if (tu==-1) {
	  Exception("TriggerTime::DecodeAndSubtract(): time_unit not found for det=\"%s\"",
		    (*it)->GetDetID().GetName().c_str()).Print();
	  continue;
        }
        const Data *d = Find(tu);
        if( d==NULL )
	  {
            Exception("TriggerTime::DecodeAndSubtract(): bad time_unit=%g for det=\"%s\"",
                      (*it)->GetTimeUnit(),(*it)->GetDetID().GetName().c_str()).Print();
            continue;
	  }

        ChipGandalf::DataID id = (*it)->GetDataID();
        double MT=d->time;

	//        TODO: do we need this for Gandalf ?
	//        if( id.u.s.src_id!=d->config->srcID )
	//            MT  += d->config->MT_shift/tu;

        double time = (*it)->getTime() - MT;
        //cout.precision(15);
        //std::cout << (*it)->GetDetID() << "  " << MT << " " <<  (*it)->getTime() << " " << time << " " << tu << "  " << time*tu << std::endl;
	//        cin.get();

	//            if( id.u.s.src_id!=d->config->srcID && d->config->MT_shift!=0 )
	//                 printf("srcID: %d    %s   cor=%g    t1=%g   t2=%g\n",
	//                             int(id.u.s.src_id),(*it)->GetDetID().GetName().c_str(),d->config->MT_shift,time,time2);

        (*it)->SetTimeDecoded(tu*time);

      }
  }

  ////////////////////////////////////////////////////////////////////////////////

  float TriggerTime::GetPhaseTCS(void) const
  {
    return TCS_phase;
  }

  ////////////////////////////////////////////////////////////////////////////////

  double TriggerTime::GetTimeInSpill(void) const
  {
    return TimeInSpill;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeConfig::SetTriggerMaskDataLimit (int min,int max)
  {
    if( min<=0 || min>max )
      throw Exception("TriggerTimeConfig::SetTriggerMaskDataLimit(): bad min/max:  min=%d max=%d",min,max);

    trigger_mask_sources_min = min;
    trigger_mask_sources_max = max;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeConfig::Add(const TriggerTimeConfig::TTC &ttc)
  {
    map<int,TriggerTimeConfig::TTC>::const_iterator m=tt.find(ttc.index);
    if( m!=tt.end() )
      {
        if( m->second!=ttc )
	  throw Exception("TriggerTime::Add(): bad trigger time settings for the index=%d",ttc.index);
      }
    else
      tt[ttc.index] = ttc;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeConfig::Print(const string &prefix) const
  {
    const char *p = prefix.c_str();
    printf("%sTrigger mask id=\"%s\"  hits: [%d,%d]\n",
	   p,trigger_mask_DetID.GetName().c_str(),trigger_mask_sources_min,trigger_mask_sources_max);
  }

  ////////////////////////////////////////////////////////////////////////////////

  const TriggerTime::Data * TriggerTime::Find(double time_unit) const
  {
    for( std::map<int,Data>::const_iterator it=data.begin(); it!=data.end(); it++ )
      if( it->second.GetConfig().time_unit==time_unit )
	return &it->second;
    return NULL;
  }

  ////////////////////////////////////////////////////////////////////////////////

  const TriggerTime::Data * TriggerTime::Find(int index) const
  {
    map<int,Data>::const_iterator m = data.find(index);
    return m==data.end() ? NULL : &m->second;
  }

  ////////////////////////////////////////////////////////////////////////////////

  double TriggerTime::GetTime(int index) const
  {
    const Data *d = Find(index);
    if( d==NULL )
      throw Exception("TriggerTime::GetTime(): TT with index %d was not found",index);

    if( d->decoded )
      return d->time;

    // Ooops! The time was not decoded! Do we have a recovery procedure?
    if( d->GetConfig().index_recover>=0 )
      return GetTime(d->GetConfig().index_recover);

    throw Exception("TriggerTime::GetTime(): trigger time was not decoded! index=%d",index);
  }

  ////////////////////////////////////////////////////////////////////////////////

  TriggerTimeConfig::TTC *TriggerTimeConfig::Find(double time_unit)
  {
    for( std::map<int,TriggerTimeConfig::TTC>::iterator it=tt.begin(); it!=tt.end(); it++ )
      if( it->second.time_unit==time_unit )
	return &it->second;
    return NULL;
  }

  ////////////////////////////////////////////////////////////////////////////////

  bool TriggerTimeConfig::TTC::operator ==(const TriggerTimeConfig::TTC &d) const
  {
    return
      d.id            == id               &&
      d.index         == index            &&
      d.index_recover == index_recover    &&
      d.channels      == channels         &&
      d.time_unit     == time_unit        &&
      d.overolling    == overolling       &&
      d.sigma         == sigma;
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeConfig::TTC::Print(const char *prefix) const
  {
    printf("%sid=%s index=%d index_recover=%d channels=%d unit=%g overolling=%d",prefix,id.GetName().c_str(),
	   index,index_recover,channels,time_unit,overolling);
    //printf(" overolling=%d sigma=%g(%gns) decoded=%d time=%g data_size=%d\n",
    //        overolling,sigma,sigma*time_unit,decoded,time,data.size());
  }

  ////////////////////////////////////////////////////////////////////////////////

  TriggerTimeStat::TriggerTimeStat(void) :
    buffer_size(1000),
    buffer_cut_ratio(0.5),
    check_entries_min(50),
    check_on_fill(100),
    allowed_shift(1),
    fills_counter(0)
  {
    SetBufferSize(buffer_size);
    SetBufferCutRaitio(buffer_cut_ratio);
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeStat::SetBufferSize(unsigned s)
  {
    buffer_size = s;
    for( map<uint32,Stat>::iterator it=triggers_stat.begin(); it!=triggers_stat.end(); it++ )
      it->second.SetBufferMaxSize(s);
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeStat::SetBufferCutRaitio(float r)
  {
    buffer_cut_ratio=r;
    for( map<uint32,Stat>::iterator it=triggers_stat.begin(); it!=triggers_stat.end(); it++ )
      it->second.SetCutRatio(r);
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeStat::Add(unsigned bit,float tt_mt_diff)
  {
    fills_counter++;
    triggers_stat[1<<bit].Add(tt_mt_diff);
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeStat::Add(unsigned bit1,unsigned bit2,float diff)
  {
    if( bit1==bit2 )
      return;

    fills_counter++;

    unsigned long int mask = (1UL<<bit1) | (1UL<<bit2);
    if( bit2<bit1 )
      diff = -diff;

    triggers_stat[mask].Add(diff);
  }

  ////////////////////////////////////////////////////////////////////////////////

  void TriggerTimeStat::Print(void) const
  {
    for( map<uint32,Stat>::const_iterator it=triggers_stat.begin(); it!=triggers_stat.end(); it++ )
      {
        const int mask = it->first;
        const Stat &stat = it->second;
        if( stat.Size()<2 )
	  continue;

        ostringstream o;
        o << " ";
        vector<int> bits=bits_from_mask(mask);
        for( vector<int>::iterator it=bits.begin(); it!=bits.end(); it++ )
	  o << *it << " ";

        printf("Trigger mask=%x bits:(%s) %s   entries=%d   mean(%g)=%g   sigma(%g)=%g\n",
               mask,o.str().c_str(),bits_string(mask).c_str(),stat.Size(),
               stat.GetCutRatio(),stat.Mean(stat.GetCutRatio()),stat.GetCutRatio(),stat.Sigma(stat.GetCutRatio()));
      }
  }

  ////////////////////////////////////////////////////////////////////////////////

}
