#include "CsStrawTubesDetector-calib.h"
#include "TH1.h"
#include <cassert>
#include <stdexcept>
#include <sstream>

using namespace std;

namespace CS
{

float CsStrawTubesDetector_CalibXray::xray_correction(int channel,float y) const
{
    // First of all we have to find the channel
    map<int,vector<float> >::const_iterator ch=xray.find(channel);
    if( ch==xray.end() )
        throw std::logic_error("CsStrawTubesDetector_CalibXray::xray_correction(): channel was not found.");
    if( y<spacers_pos.front() )
        return ch->second.front();
    if( y>spacers_pos.back() )
        return ch->second.back();
    
    for( size_t i=1; i<spacers_pos.size(); i++ )
        if( spacers_pos[i-1]<=y && spacers_pos[i]>=y )
        {
            float dy=y-spacers_pos[i-1];
            float c = ch->second[i] - ch->second[i-1];
            return ch->second[i-1] + c * dy/(spacers_pos[i]-spacers_pos[i-1]);
        }
    throw std::logic_error("CsStrawTubesDetector_CalibXray::xray_correction(): internal problem");
}

void CsStrawTubesDetector_CalibXray::read(const string &data)
{
    istringstream is(data);

    spacers_pos=vector<float>(6,0);
    string line;
    char buf[22];

    getline(is,line);
    if( !is || 7!=sscanf(line.c_str(),"%s %g %g %g %g %g %g",buf,
           &spacers_pos[0],&spacers_pos[1],&spacers_pos[2],
           &spacers_pos[3],&spacers_pos[4],&spacers_pos[5]) )
        throw std::logic_error("Failed to read STRAW X-ray calibrations!");

    int ch=0;
    while( getline(is,line) )
    {
        float x[6];
        sscanf(line.c_str(),"%g %g %g %g %g %g",&x[0],&x[1],&x[2],&x[3],&x[4],&x[5]);
        xray[ch++] = vector<float>(x,x+6);
    }
}

void CsStrawTubesDetector_CalibXray::make_pictures(const string &stdc)
{
    char name[222];

    sprintf(name,"%s_Xcor_NEW",stdc.c_str());
    TH1F *h1 = new TH1F(name,"X-ray correction",222,0,222);

    for( map<int,vector<float> >::const_iterator ch=xray.begin(); ch!=xray.end(); ch++ )
    {
        assert(ch->second.size()==6);
        h1->Fill(ch->first,(ch->second[2]+ch->second[3])/2);
    }
}


void CsStrawTubesDetector_CalibT0::read(const string &data)
{
    istringstream is(data);

    string line;

    while( getline(is,line) )
    {
        if( line.length()==0 || line[0]=='#' )
            continue;

        int c,p;
        float T0;
        if( 3!=sscanf(line.c_str(),"%d %d %g",&c,&p,&T0) )
            throw std::logic_error("CsStrawTubesDetector_CalibT0::read(): bad line structure.");
        T0s[ChanPos(c,p)] = T0;
    }
}

float CsStrawTubesDetector_CalibT0::GetT0(int channel,int pos) const
{
    std::map<ChanPos,float>::const_iterator it=T0s.find(ChanPos(channel,pos));
    if( it==T0s.end() )
        throw std::logic_error("CsStrawTubesDetector_CalibT0::getT0() channel not found!");
    return it->second;
}

}
