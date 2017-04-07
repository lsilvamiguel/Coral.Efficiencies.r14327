#include <cassert>

#include "ChipSADC.h"
#include "DaqEvent.h"
#include "ChipSADC_huffman.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

ChipSADC::ChipSADC(void const * const buf,bool copy_buf,DaqOption &opt,DaqEvent &ev)
:   Chip(buf,copy_buf,opt,ev)
{
    Clear();
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::ScanADCBlockFormat2(const HeaderADC & adc,const uint32 *start,DaqOption &opt)
{
    Digit digit(DataID(0,0,0,0),DetID(""));

    int chan_prev=-1;
    bool expect_channel_info=true, expect_data=false, data_end = false;
    const ChannelInfo *channel_info=NULL;

    for( int i=0; i<int(adc.size)-1; i++ )
    {
        Data &d = *(Data*)(start+i);
        bool word_data=false, word_channel_info=false;

        if( expect_data && (d.id==2 || (d.id==1 && !expect_channel_info)) )
        {
            word_data = true;
            //printf("This is a data word!\n");

            // We may expect another channel info word after a data word!
            expect_channel_info = true;

            digit.GetSamples().push_back(1023-d.sample0);
            digit.GetSamples().push_back(1023-d.sample1);
            digit.GetSamples().push_back(1023-d.sample2);
            
            data_end = i==int(adc.size)-2;
        }

        if( (!word_data && d.id==1 && expect_channel_info) || data_end )
        {
            //printf("This is a channel info!\n");
            word_channel_info=true;
            
            // Do we need to finalize a previous digit?
            // For this I check for a GetSamples() container!
            if( !digit.GetSamples().empty() )
            {
                // We may need to remove the two last samples, if they are zeroes.
                for( int i=0; i<2; i++ )
                    if( !digit.GetSamples().empty() && digit.GetSamples().back()==1023 )
                        digit.GetSamples().pop_back();
                    else
                        break;

                assert(channel_info!=NULL);

                bool good_digit=true;
                // Now we check for the number of samples!
                if( digit.GetSamples().size()!=channel_info->samples )
                {
                    AddError(DaqError(DaqError::SADC_BAD_N_SAMPLES,DaqError::SOURCE_ID,GetSourceID(),
                                                                   DaqError::CHIP,adc.chip,
                                                                   DaqError::COUNTER,digit.GetSamples().size(),
                                                                   DaqError::COUNTER,channel_info->samples),opt);
                    good_digit = false;
                }

                // Check the channels order.
                if( digit.GetChannel()<=chan_prev )
                {
                    AddError(DaqError(DaqError::SADC_BAD_CHAN_N,DaqError::SOURCE_ID,GetSourceID(),
                                                                DaqError::CHIP,adc.chip),opt);
                    good_digit = false;
                }

                chan_prev=digit.GetChannel();

                // And create it!
                if( good_digit )
                    pre_digits.push_back(new Digit(digit));

                // Now we should clear the digit
                digit.GetIntegrals().clear();
                digit.GetSamples().clear();
                channel_info = NULL;
                
                // End of a channel digit pre-creation!
                // The work with a mapping will be done by the Decode() function call.
            }
            
            if( !data_end && expect_channel_info )
            {
                // OK, new digit!
                channel_info = reinterpret_cast<const ChannelInfo *>(start+i);
                digit.SetDataID(DataID(GetSourceID(),0,adc.chip,channel_info->channel));
                digit.GetIntegrals().push_back(channel_info->sum);
            }

            expect_channel_info = false;
            expect_data         = true;
        }

        // Check the data.
        if( !word_data && !word_channel_info )
            AddError(DaqError(DaqError::SADC_UNKNOWN_DATA,DaqError::SOURCE_ID,GetSourceID(),DaqError::COUNTER,adc.size),opt);
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Scan(DaqOption &opt)
{
    if( is_scaned )
        return;
    else
        is_scaned=true;

    if( !pre_digits.empty() )
        throw "ChipSADC::Scan(): an attempt to call Scan() without Clear()!";

    if( opt.GetSADCDataVer3().count(GetSourceID())>0 )
    {
        ScanDataVer3(opt);
        return;
    }
    
    Digit digit(DataID(0,0,0,0),DetID(""));
    int nwords=-1;

    const uint32 *p;
    for( p=GetDataStart(); p<GetDataEnd(); )
    {
        // We expect an ADC header
        const HeaderADC *h_adc = (HeaderADC*) p;

        const HeaderADC_compressed *h_adc_cmp = (HeaderADC_compressed*) p;
        CS::hufftree_t *hufftree = (CS::hufftree_t*)opt.GetChipSADC_HT();

        bool compressed = h_adc->GetMode()==HeaderADC::COMPRESSED;

        if( h_adc->GetMode()==HeaderADC::UNKNOWN )
        {
            AddError(DaqError(DaqError::SADC_UNKNOWN_MODE,DaqError::SOURCE_ID,GetSourceID()),opt);
            return;
        }

        // Now we check HeaderADC.

        if( compressed )
        {
            // Check last 8 bits of the event number.
            if( h_adc_cmp->event!=(255&event->GetHeader().GetEventNumberInBurst()) )
            {
                AddError(DaqError(DaqError::SADC_WRONG_EVENT, DaqError::SOURCE_ID,GetSourceID(),
                                                              DaqError::SL_EVENT_N,255&event->GetHeader().GetEventNumberInBurst(),
                                                              DaqError::EVENT_N,h_adc->event),opt);
              return;
            }
        }
        else
        {
            // Check last 12 bits of the event number.
            if( h_adc->event!=(4095&GetSLink().GetEventNumber()) )
            {
                AddError(DaqError(DaqError::SADC_WRONG_EVENT, DaqError::SOURCE_ID,GetSourceID(),
                                                              DaqError::SL_EVENT_N,4095&event->GetHeader().GetEventNumberInBurst(),
                                                              DaqError::EVENT_N,h_adc->event),opt);
                return;
            }
        }

        if( ((uint32(GetDataEnd()-p)<h_adc->size) && (!compressed))
            || ((uint32(GetDataEnd()-p)<h_adc_cmp->size) && compressed) )
        {
            AddError(DaqError(DaqError::SADC_H_ADC_BAD_SIZE,DaqError::SOURCE_ID,GetSourceID(),DaqError::COUNTER,h_adc->size,DaqError::VALUE,1),opt);
            return;
        }

        if( h_adc->size<1 || (p+h_adc->size)>GetDataEnd() )
        {
            AddError(DaqError(DaqError::SADC_H_ADC_BAD_SIZE,DaqError::SOURCE_ID,GetSourceID(),DaqError::COUNTER,h_adc->size,DaqError::VALUE,2),opt);
            return;
        }

        p++;    // go to the first word after the ADC header
        
        int chan_prev=-1;

        if( h_adc->GetMode()==HeaderADC::LATCH || h_adc->GetMode()==HeaderADC::ESPARSE || h_adc->GetMode()==HeaderADC::SPARSE )
        {
            if( opt.GetSADCDataVer2().count(GetSourceID())>0 )
            {
                ScanADCBlockFormat2(*h_adc,p,opt);
                p+=h_adc->size-1;
            }
            else
                // Decode ADC data block in the old format (and old code!)
                while( p<(((uint32*)h_adc)+h_adc->size) )
                {
                    if( h_adc->GetMode()==HeaderADC::LATCH || h_adc->GetMode()==HeaderADC::ESPARSE )
                    {
                        // This is latch all data or extended sparse mode

                        // We don't know how many there are samples.
                        if( nwords<0 && h_adc->size>1 )
                        {
                            // The task is to get the number of samples!
                            // For this we search for the data headers.

                            for( const uint32 *pp=p; pp<(((uint32*)h_adc)+h_adc->size); pp++ )
                            {
                                Data &d = *(Data*)pp;

                                if( d.id==2 )
                                  continue;

                                if( d.id!=1 )
                                {
                                    AddError(DaqError(DaqError::SADC_UNKNOWN_DATA,DaqError::SOURCE_ID,GetSourceID(),DaqError::COUNTER,h_adc->size),opt);
                                    return;
                                }

                                if( nwords==-1 )
                                    if( pp==p )
                                    {
                                        //printf("First data line is detected!\n");
                                        nwords=-2;
                                        continue;
                                    }
                                    else
                                    {
                                        AddError(DaqError(DaqError::SADC_NO_DATA_LINE,DaqError::SOURCE_ID,GetSourceID(),DaqError::COUNTER,h_adc->size),opt);
                                        return;
                                    }
                                else
                                {
                                    // at this moment this is not samples, but the number of words.
                                    nwords=pp-p;
                                    break;
                                }
                            }

                            if( nwords==-2 )
                            {
                                // Can not get number of samples: end of data rached!
                                // (There was only one data block, so one data header existed.)
                                // In that case we use adc block size.
                                nwords=h_adc->size-1;
                            }

                            // Must be bigger then 3 because of the sample integrals (3 words).                    
                            if( nwords<=3 )
                            {
                                AddError(DaqError(DaqError::SADC_SHORT_SAMPLES,DaqError::SOURCE_ID,GetSourceID(),DaqError::COUNTER,h_adc->size,DaqError::VALUE,nwords),opt);
                                return;
                            }
                        }

                        // Now it is time to read samples.
                        for( int i=0; i<nwords-3; i++,p++ )
                        {
                            Data &d = *(Data*)p;
                            digit.GetSamples().push_back(d.sample0);
                            digit.GetSamples().push_back(d.sample1);
                            digit.GetSamples().push_back(d.sample2);
                        }

                        // We may need to remove two last samples, if they are zeroes.
                        if( !digit.GetSamples().empty() && digit.GetSamples().back()==0 )
                        {
                            digit.GetSamples().pop_back();
                            if( !digit.GetSamples().empty() && digit.GetSamples().back()==0 )
                                digit.GetSamples().pop_back();
                        }
                    }

                    // We must have at least 3 words for channel data samples.
                    if( p+3>GetDataEnd() )
                    {
                        AddError(DaqError(DaqError::SADC_H_ADC_BAD_SIZE,DaqError::SOURCE_ID,GetSourceID(),
                                                                        DaqError::CHIP,h_adc->chip,
                                                                        DaqError::COUNTER,h_adc->size,
                                                                        DaqError::VALUE,3),opt);
                        return;
                    }

                    uint16 chan=1111, suppression=2;

                    // process sample integrals
                    for( int i=0; i<3; i++,p++ )
                    {
                        Integral &intg = *(Integral*)p;

                        if( intg.not_used_1!=0 || intg.not_used_2!=0 || intg.not_used_3!=2 )
                        {
                            AddError(DaqError(DaqError::SADC_BAD_INTEGRAL,DaqError::SOURCE_ID,GetSourceID(),
                                                                          DaqError::CHIP,h_adc->chip),opt);
                            //continue;
                        }

                        if( i>0 && (chan!=intg.channel || suppression!=intg.suppression) )
                            AddError(DaqError(DaqError::SADC_BAD_INTEGRALS,DaqError::SOURCE_ID,GetSourceID(),
                                                                          DaqError::CHIP,h_adc->chip),opt);

                        chan=intg.channel;
                        suppression=intg.suppression;

                        digit.GetIntegrals().push_back(intg.integral);
                    }

                    if( chan_prev!=-1 && chan<=chan_prev )
                        AddError(DaqError(DaqError::SADC_BAD_CHAN_N,DaqError::SOURCE_ID,GetSourceID(),
                                                                    DaqError::CHIP,h_adc->chip),opt);
                    chan_prev=chan;

                    // Finalize the digit
                    digit.SetDataID(DataID(GetSourceID(),0,h_adc->chip,chan));
                    digit.SetSuppression(suppression);
                    digit.SetOverflow(h_adc->overflow);

                    // And create it!
                    pre_digits.push_back(new Digit(digit));

                    // Now we should clear digit
                    digit.GetIntegrals().clear();
                    digit.GetSamples().clear();
                }
        }
        
        /* decoding for compressed sadc data (10 data words) */
        if( compressed )
        {
            uint32 *decoded_samples;
            int decoded_samples_sz = 0;
            unsigned int samples = 0;

            if (!hufftree)
                throw Exception("ChipSADC::Scan(): No huffman tree found! srcID=%d",GetSourceID());

            if( h_adc_cmp->size > 1)
            {
                int pos = 0;  
                
                decoded_samples_sz = sadc_huffman_decode(hufftree, p, h_adc_cmp->size - 1, &samples, &decoded_samples);
              
                if (decoded_samples_sz <= 0)
                    throw Exception("ChipSADC::Scan(): huffman decoding failed: srcID=%d",GetSourceID());
              
                /* now reconstruct original samples and store them in the digit */
                while(pos < decoded_samples_sz)
                {
                  uint32 prev_sample = 0;
                  int i, chan;
                  
                  chan = decoded_samples[pos++];
                  /* save first sample, because its needed to reconstruct all the others */
                  prev_sample = decoded_samples[pos++];
                  digit.GetSamples().push_back(prev_sample);
                  
                  /* reconstruct all the other samples */
                  for(i = 1; (unsigned int)i < samples; i++, pos++)
                  {
                    digit.GetSamples().push_back(( prev_sample + decoded_samples[pos])&4095);
                    prev_sample = (prev_sample + decoded_samples[pos])&4095;
                  }
                                                      
                  if( chan_prev!=-1 && chan<=chan_prev )
                    AddError(DaqError(DaqError::SADC_BAD_CHAN_N,DaqError::SOURCE_ID,GetSourceID(),
                                                                DaqError::CHIP,h_adc_cmp->chip),opt);
                  chan_prev=chan;
                  
                  // Finalize the digit
                  digit.SetDataID(DataID(GetSourceID(),0,h_adc_cmp->chip,chan));
                  //digit.SetSuppression(suppression);
                  digit.SetOverflow(h_adc_cmp->overflow);

                  // And create it!
                  pre_digits.push_back(new Digit(digit));
                  
                  // Now we should clear digit
                  digit.GetIntegrals().clear();
                  digit.GetSamples().clear();
                }
                p += h_adc_cmp->size - 1;
            }
        }
    }

    if( uint32(GetDataEnd()-p)!=0 )
        AddError(DaqError(DaqError::SADC_H_ADC_BAD_SIZE,DaqError::SOURCE_ID,GetSourceID(),DaqError::VALUE,4),opt);
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::ScanDataVer3(DaqOption &opt)
{
    bool debug=0;

    if( debug )
        Dump();

    for( const uint32 *p=GetDataStart(); p<GetDataEnd(); )
    {
        //printf("****************   size to scan %d\n",GetDataEnd()-p);
        assert(p<GetDataEnd());

    
        // p - points to an ADC header
        HeaderADC_v3 *adc = (HeaderADC_v3*)p;
        if( debug )
            adc->Print();

        if( adc->key==0 )
        {
            AddError(DaqError(DaqError::SADC_MODULE_MISSING, DaqError::SOURCE_ID,GetSourceID(),
                                                             DaqError::PORT,adc->hl_port,
                                                             DaqError::CHIP,adc->adc_id),       opt);
        }

        if( adc->key!=3 && !(adc->key==0 && adc->size==1 ) )
        {
            // Ooops... There is something which I don't know what it is...
            AddError(DaqError(DaqError::SADC_BAD_H_ADC,DaqError::SOURCE_ID,GetSourceID()),opt);
            break;
        }

        // Check that ADC block size does not go beyond the chip data block
        if( p+adc->size > GetDataEnd() )
        {
            AddError(DaqError(DaqError::SADC_H_ADC_BAD_SIZE,DaqError::SOURCE_ID,GetSourceID()),opt);
            break;
        }

        // OK, the ADC header looks OK, next time we can try to scan a next ADC header
        p += adc->size;
        //printf("next p:   end-p=%d\n",GetDataEnd()-p);

        //#warning "SADC-ADC header error bit is ignored (data version 2008)"
//         if( adc->error )
//         {
//             AddError(DaqError(DaqError::SADC_H_ADC_ERROR,DaqError::SOURCE_ID,GetSourceID()),opt);
//             continue;
//         }

        if( adc->event!=(4095&GetSLink().GetEventNumber()) )
        {
            AddError(DaqError(DaqError::SADC_WRONG_EVENT, DaqError::SOURCE_ID,GetSourceID(),
                                                          DaqError::CHIP,adc->adc_id,
                                                          DaqError::SL_EVENT_N,4095&event->GetHeader().GetEventNumberInBurst(),
                                                          DaqError::EVENT_N,adc->event),opt);
            break;
        }

//         if( adc->size==1 )
//         {
//             AddError(DaqError(DaqError::SADC_MODULE_MISSING, DaqError::SOURCE_ID,GetSourceID(),
//                                                           DaqError::PORT,adc->hl_port,
//                                                           DaqError::CHIP,adc->adc_id),opt);
//             continue;
//         }
// 
// 
        const DataHeader_v3 *data_header_next = (DataHeader_v3*)(((uint32*)adc)+1);
        for( int words_to_scan=adc->size-1; words_to_scan>0; )
        {
            const DataHeader_v3 *data_header=data_header_next;
            if( debug )
            {
                printf("pos=%zu   #=%d  %d\n",((uint32*)data_header)-GetDataStart(),words_to_scan,((uint32*)data_header)<GetDataEnd());
                data_header->Print("DataHeader: ");
            }

            if( data_header->id!=1 )
            {
                AddError(DaqError(DaqError::SADC_BAD_DATA_HEADER, DaqError::SOURCE_ID,GetSourceID(),DaqError::CHIP,adc->adc_id),opt);
                break;
            }

            if( debug )
                printf("GOOD DIGIT: hl_port=%d  adc_id=%d  chan_id=%d\n",adc->hl_port,adc->adc_id,data_header->chan_id);
            
            // Number of data words
            int data_words = (data_header->samples+1)/2, samples=data_header->samples;
            if( samples==0 || data_words==0 || data_words>words_to_scan )
            {
                AddError(DaqError(DaqError::SADC_DATA_OUTSIDE, DaqError::SOURCE_ID,GetSourceID()),opt);
                break;
            }
            words_to_scan -= data_words+1;
            data_header_next = (DataHeader_v3*)(((uint32*)data_header_next)+data_words+1);
            
            // Now we create a digit.
            DataID data_id( GetSourceID(), adc->hl_port, adc->adc_id, data_header->chan_id );
            Digit *digit = new Digit(data_id,DetID(""));
            pre_digits.push_back(digit);
            digit->GetIntegrals().push_back(data_header->sum);
            
            // Scan all data words (samples)
            for( int n=0; n<data_words; n++ )
            {
                const Data_v3 *d = (Data_v3*)(((uint32*)data_header)+n+1);
                uint16 id = d->id>>6;
                const bool line_ok =  (n==0) ? (id==1) : (id==2);
                if( !line_ok )
                {
                    AddError(DaqError(DaqError::SADC_BAD_DATA, DaqError::SOURCE_ID,GetSourceID(),DaqError::CHANNEL,data_header->chan_id),opt);
                    break;
                }

                assert(samples>0);
                digit->GetSamples().push_back(d->data1);
                samples--;
                if( samples>0 )
                {
                    digit->GetSamples().push_back(d->data2);
                    samples--;
                }
            }
            
            if( debug)
                digit->Print();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::HeaderADC_v3::Print(const char *prefix) const
{
    printf("%skey=%d hl_port=%d adc_id=%d error=%d size=%d event=%d\n",
            prefix,key,hl_port,adc_id,error,size,event);
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::DataHeader_v3::Print(const char *prefix) const
{
    printf("%sid=%d chan_id=%d zero=%d samples=%d sum=%d\n",
            prefix,id,chan_id,zero,samples,sum);
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Data_v3::Print(const char *prefix) const
{
    printf("%sid=%d data1=%d data2=%d\n",prefix,id,data1,data2);
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Decode(const Maps &maps,Digits &digits_list,DaqOption &opt)
{
    // Scan the chip data if it was not done yet.
    if( !is_scaned )
        Scan(opt);

    for( std::list<Digit*>::iterator it=pre_digits.begin(); it!=pre_digits.end(); it++ )
    {
        typedef Maps::const_iterator m_it; // Create a short name for map's iterator
        const pair<m_it,m_it> m_range = maps.equal_range((*it)->GetDataID()); // all maps with given data ID

        for( m_it c=m_range.first; c!=m_range.second; c++ )
        {
            const Digit *digit_map = dynamic_cast<Digit*>(c->second);
            if( digit_map==NULL )
                throw Exception("ChipSADC::Decode(): Internal error");
            const DigitCalib *digit_calib_map = dynamic_cast<const DigitCalib*>(digit_map);
            
            Digit *digit2 = NULL;
            DigitCalib *digit2_calib = NULL;
            if( digit_calib_map!=NULL )
            {
                digit2_calib = new DigitCalib(*digit_calib_map);
                digit2 = digit2_calib;
            }
            else
            {
                digit2 = new Digit(*digit_map);
            }
            
            *digit2 = **it;
            digit2->SetDetID(digit_map->GetDetID());
            digit2->SetX(digit_map->GetX());
            digit2->SetY(digit_map->GetY());

            assert(digit2->GetChip()==digit_map->GetChip());
            assert(digit2->GetChannel()==digit_map->GetChannel());

            digits_list.insert(pair<DetID,Digit*>(digit2->GetDetID(),digit2));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Clear(void)
{
    Chip::Clear();
    for( std::list<Digit*>::iterator it=pre_digits.begin(); it!=pre_digits.end(); it++ )
        delete *it;
    pre_digits.clear();
}

////////////////////////////////////////////////////////////////////////////////

vector<float> ChipSADC::Digit::GetNtupleData(void) const
{
    vector<float> v;
    v.push_back(GetX());
    v.push_back(GetY());

    // Push samples after the (x,y) coordinate.
    // This is 'feature' has been requested for PHAST SADC digits export.
    for( vector<uint16>::const_iterator it=GetSamples().begin(); it!=GetSamples().end(); it++ )
        v.push_back(*it);

    return v;
}

////////////////////////////////////////////////////////////////////////////////

ChipSADC::HeaderADC::Mode ChipSADC::HeaderADC::GetMode(void) const
{
    if( key==2 && sparse==0 )
        return LATCH;

    if( key==3 && sparse==1 )
        return ESPARSE;

    if( key==2 && sparse==1 )
        return SPARSE;
    
    if( key == 1)
        return COMPRESSED;

    return UNKNOWN;
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::HeaderADC::Print(const char *prefix) const
{
    printf("%s event=%4.4d size=%.4d sparse=%d ovfl=%d chip=%2d\n",
            prefix,event,size,sparse,overflow,chip);
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::HeaderADC_compressed::Print(const char *prefix) const
{
    printf("%s event=%4.4d huff_id=%d size=%.4d sparse=%d ovfl=%d chip=%2d\n",
            prefix,event,huff_id,size,sparse,overflow,chip);
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Data::Print(const char *prefix) const
{
    printf("%s samples=(%d,%d,%d)\n",prefix,sample0,sample1,sample2);
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Integral::Print(const char *prefix) const
{
    printf("%s integral=%4d suppression=%d channel=%2d\n",
            prefix,integral,suppression,channel);
}

////////////////////////////////////////////////////////////////////////////////

ChipSADC::Map::Map(const ObjectXML &o) :
  Chip::Map(o),
  port(0),
  chip(0),
  channel(0),
  x(0),
  y(0)
{
    if( GetName()!=o.GetName() || GetName()!="ChipSADC" )
        throw Exception("ChipSADC::Map::Map(): Internal error.");

    istringstream s(dec_line.c_str());

    string name;
    int chanN=-1;

    if( !IsOption("xy") )
        throw Exception("ChipSADC::Map::Map(): option \"xy\" is mandatory!");

    if( IsOption("FEM") )
    {
        s >> name >> det_ref >> source_id >> chip >> channel;
        chanN=1;
        chanF=chanL=channel;
        chanS=1;
        wireF=wireL=wireS=1;
        
        size_t pos=0;
        
        while(1)
        {
            pos=dec_line.find('(',pos);
            if( pos==string::npos )
                break;
            size_t end=dec_line.find(')',pos+1);
            if( end==string::npos )
                throw Exception("ChipSADC::Map::Map(): check the brackets () for \"FEM\" mapping!");
            
            
            string box=dec_line.substr(pos+1,end-pos-1);

            int x1,y1,x2,y2;
            if( 4!=sscanf(box.c_str(),"%d,%d,%d,%d",&x1,&y1,&x2,&y2) )
                throw Exception("ChipSADC::Map::Map(): bad Box description in \"FEM\" mapping: \"%s\"",box.c_str());

            boxes.push_back(DigitCalib::Box(x1,y1,x2,y2));
            
            pos=end+1;
        }
        
    }
    else
    {
        switch( GetVersion() )
        {
            case 1:
            {
                s >> name >> source_id >> chip >> channel >> x >> y;
                chanN=1;
                chanF=chanL=channel;
                chanS=1;
                wireF=wireL=x;
                wireS=1;

                break;
            }
            case 2:
            {
                uint16 chip_channel;
                s >> name >> source_id >> port >> chip_channel >> x >> y;
                channel=chip_channel&0xf; // last 4 bits
                chip = chip_channel>>4;   // ignore channel number.
                chanN=1;
                chanF=chanL=channel;
                chanS=1;
                wireF=wireL=x;
                wireS=1;

                break;
            }
            default:
                throw Exception("ChipSADC::Map::Map(): unknown version %d",GetVersion());
        }
    }

    if( s.fail() )
        throw Exception("ChipSADC::Map::Map(): bad format in line: %s",map_line.c_str());

    if( IsOption("data_format2") )
        data_fmt = "version2";
    if( IsOption("data_format3") )
        data_fmt = "version3";

    id=DetID(name);
    Check();
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Map::AddToMaps(Maps &maps,DaqOption &options) const
{
  const DataID data_id(GetSourceID(),port,chip,channel);

  if( maps.end()!=maps.find(data_id) && !IsMultiDigit() )
  {
      Print();
      Exception("ChipSADC::Map::AddToMaps(): map already exists").Print();
  }

  DigitCalib *digit_calib = NULL;
  Digit *digit = NULL;
  
  if( IsCalib() )
  {
      digit_calib = new DigitCalib(data_id,GetDetID());
      digit_calib->GetBoxes() = boxes;
      digit_calib->SetRefDet(det_ref);

      digit = digit_calib;
  }
  else
      digit = new Digit(data_id,GetDetID());
      
  
  digit->SetX(x);
  digit->SetY(y);
  maps.insert( pair<DataID,Digit*>(data_id,digit));
  maps.SetWires( GetDetID(), maps.GetWires(GetDetID())+1 );

  // Check the data format version number
  if( IsOption("data_format2") )
    options.GetSADCDataVer2().insert(GetSourceID());
  if( IsOption("data_format3") )
    options.GetSADCDataVer3().insert(GetSourceID());
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Map::Print(ostream &o,const string &prefix) const
{
  Chip::Map::Print(o,prefix);
  o<< prefix << "      ";
  char s[222];
  sprintf(s,"srcID=%d  port=%d  chip=%3d  channel=%3d   x=%2d  y=%2d\n",GetSourceID(),port,chip,channel,x,y);
  o << s;
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::Digit::Print(ostream &o,const string &prefix) const
{
    Chip::Digit::Print(o,prefix);
    o << prefix << "x=" << x << "   y=" << y << "\n";

    static_cast<ChipSADC::DataID>(GetDataID()).Print(prefix);

    o << prefix << (samples.empty()?"sparse mode,":"latch all mode,")
      << "    overflow=" << GetOverflow() << ",  suppresion=" << GetSuppression() << "\n";

    o << prefix << "Integrals:";
    for( vector<uint16>::const_iterator it=integrals.begin(); it!=integrals.end(); it++ )
        o << "   " << *it;
    o << "\n";

    if( !samples.empty() )
    {
        o << prefix << "Samples:";
        for( vector<uint16>::const_iterator it=samples.begin(); it!=samples.end(); it++ )
            o << "   " << *it;
        o << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::DigitCalib::Box::Print(ostream &o,const string &prefix) const
{
    o << prefix << "(" << x1 << "," << y1 << "," << x2 << "," << y2 << ")";
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::DigitCalib::Print(ostream &o,const string &prefix) const
{
    Chip::Digit::Print(o,prefix);

    static_cast<ChipSADC::DataID>(GetDataID()).Print(prefix);

    o << prefix << "References detector is \"" << det_ref << "\"\n";
    o << prefix << "Boxes:";
    for( vector<Box>::const_iterator b=boxes.begin(); b!=boxes.end(); b++ )
        b->Print(o," ");
    o << "\n";

    o << prefix << (GetSamples().empty()?"sparse mode,":"latch all mode,")
      << "    overflow=" << GetOverflow() << ",  suppresion=" << GetSuppression() << "\n";

    o << prefix << "Integrals:";
    for( vector<uint16>::const_iterator it=GetIntegrals().begin(); it!=GetIntegrals().end(); it++ )
        o << "   " << *it;
    o << "\n";

    if( !GetSamples().empty() )
    {
        o << prefix << "Samples:";
        for( vector<uint16>::const_iterator it=GetSamples().begin(); it!=GetSamples().end(); it++ )
            o << "   " << *it;
        o << "\n";
    }
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::DataID::Print(const string &prefix) const
{
    printf("%ssrcID=%d port=%d chip=%d channel=%d\n",
            prefix.c_str(),
            unsigned(u.s.src_id),unsigned(u.s.port),unsigned(u.s.chip),unsigned(u.s.chan));
}

////////////////////////////////////////////////////////////////////////////////

void ChipSADC::ChannelInfo::Print(const char *prefix) const
{
    printf("%schannel=%d  sum=%d  sample=%d\n",prefix,channel,sum,samples);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace
