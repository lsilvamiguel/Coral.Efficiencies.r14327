/*
 * GroupCamera.cpp
 *
 * This class maps at first the up and down data of one plane together to
 * generate time/energy difference plots for each szintilator.
 * Then, A and B ring informations are combined to generate the coincidence
 * plots. The hit information of Gandalf are organized in the 'Szinti' respectively
 * the 'ASzinti' class.
 *
 *
 *  Created on: Apr 18, 2012
 *      Author: gotzl
 */

#include "GroupCamera.h"
#include <sstream>
#include "TH1.h"
#include "TH2.h"


/*   	
 *	In HIGHRes mode, the 2D plots have much bigger bining. Due to a bug in root/coool
 *	it is then not possible to look at other groups than Camera.
 *	For monitoring, the extra granularity is not important.
 */
#define HIGH_RES	false


/*
 * calibration constants for the time of flight; each Block of two values
 * stands for combinations of one A szinti with two neighbouring Bs,
 * begining from A0 -> B0 and A0 -> B1: {B0,B1}
 */
/*
Valid before 106952
const double TOF_CALIB[24][2] = {	{-29.04,-26.70},	{-29.1,-29.46},		{-28.79,-28.3},		{-29.45,-29.58}, 	// A0 - A3
					{-30.92,-32.09},	{-31.97,-28.8},		{-28.59,-29.0},		{-28.24,-28.41},	// A4 - A7
					{-27.31,-27.48},	{-29.61,-28.42},	{-28.43,-27.49},	{-33.00,-31.05},		// A8 - A11
					{-26.52,-29.63},	{-29.93,-26.73},	{-28.27,-26.59},	{0,0},			// A12 - A15
					{-30.5,-30.36},		{-30.71,-43.60},	{-30.8,-16.21},		{-30.18,-32.63},	// A16 - A19
					{-29.32,-27.4},		{-27.49,-27.15},	{-26.8,-27.77},		{-28.3,0}		// A20 - A23
};
*/

/*
Valid after 106952
*/
      double TOF_CALIB[24][2] = {	{-27.81,-27.99},	{-29.19,-29.46},       	{-28.44,-28.67},       	{-29.20,-29.45}, 	// A0 - A3
					{-30.94,-32.26},	{-32.17,-28.97},       	{-29.40,-29.26},        {-29.58,-30.23},	// A4 - A7
					{-27.31,-27.48},	{-29.61,-28.42},	{-28.43,-27.49},	{-33.00,-31.05},	// A8 - A11
					{-26.52,-29.63},	{-29.93,-26.73},	{-28.27,-26.59},	{0,0},			// A12 - A15
					{-30.5,-30.36},		{-30.71,-43.60},	{-30.8,-16.21},		{-30.18,-32.63},	// A16 - A19
					{-29.32,-27.4},		{-27.49,-27.15},	{-26.8,-27.77},		{-28.3,-31.4}		// A20 - A23{B23,B0}
};

/*
 * Calibration constants for the timedifference in each slat, taken from laser such that
 * these differences are at zero
 * First iteration from laser run 106522
 */
/*
Valid before 106952
const double Adiff_CALIB[24] = {1.746,2.426,2.493,2.811,3.075,3.199,2.848,2.853,4.283,2.741,3.392,3.069,4.557,4.829,3.055,0,3.837,1.217,1.874,3.43,4.345,2.233,2.76,2.88};
const double Bdiff_CALIB[24] = {-0.673,-0.2665,1.281,3.807,-0.2368,-3.289,-1.397,-1.467,-1.185,2.805,-0.6711,1.03,-1.762,0,-1.855,0.3746,0.327,1.924,-1.55,-1.8,4.972,-0.6335,0.4,4.25};
*/

/*
Valid after 106952
*/
//                               A0     A1     A2     A3     A4     A5     A6     A7     A8     A9     A10    A11    A12    A13    A14    A15    A16    A17    A18    A19    A20    A21    A22    A23
      double Adiff_CALIB[24] = { 2.211, 2.861, 2.573, 2.935, 3.250, 3.200, 2.938, 2.906, 3.590, 2.883, 3.336, 3.213, 4.594, 3.964, 2.743, 0.000, 3.907, 1.138, 1.874, 3.430, 4.345, 2.207, 2.760, 3.016};
      double Bdiff_CALIB[24] = {-1.009,-0.250, 1.333, 3.728,-0.751,-3.111,-1.315,-1.738,-1.148, 2.805,-0.6711,1.03,-1.762,0,-1.855,0.3746,0.327,1.924,-1.55,-1.8,4.972,-0.6335,0.4,4.25};


/*
 * Calibration constants for the meantime in each slat
 * can be set to a reasonable value
 * First
 */
//                               A0     A1     A2     A3     A4     A5     A6     A7     A8     A9     A10    A11    A12    A13    A14    A15    A16    A17    A18    A19    A20    A21    A22    A23
      double Amean_CALIB[24] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};
      double Bmean_CALIB[24] = { 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000};


//Calibration constants det. by Robert leave commented, for crosschecking, Run 106522, constants seem consitent with run 106524
/*
//DiffT Offset is to put dt distr symmetric to 0
const double 
DiffT_OffA[24]={-1.747,-2.65,-2.493,-2.813,-3.25,-3.2,-2.737,0.,-4.30,-2.9,0.,-3.215,-4.266,-5.144,-2.65,0.,-3.82,-1.218,-1.847,0.,-4.379,-2.268,-2.580,-3.02};
const double 
DiffT_OffB[24]={0.673,0.267,-1.28,-3.806,0.729,3.256,1.474,1.257,1.184,-2.805,1.333,0.,1.762,0.,1.856,-0.375,0.,0.,+1.55,0.,-4.972,0.634,0.,-4.25};
//meanT Offset is to put meanTime distr symm to zero. Subtract AFTER div. by 2. +16 to all!
const double meanT_OffA[24]={2760.,2757.,2758.,2757.,2757.,2758.,2757.,0.,2758.,2756.,2756.,2756.,2757.,2758.,2756.,0.,0.,0.,2832.,0.,0.,0.,0.,0.};
const double meanT_OffB[24]={2789.,2786.,2787.,2728.,0.,0.,0.,0.,2785.,2785.,2785.,0.,2784.,0.,2786.,2786.,0.,0.,2861.,0.,2786.,2784.,2781.,2784.};
*/


//const double Z_CALIB = 20;


const double B_ZOFFSET = 42.5;

const float speedOfLightA = 13.;
const float speedOfLightB = 16.;

inline int Modulo(int value, int modulo)
{
  while(value < 0)
    value+=modulo;
  while(value >= modulo)
    value-=modulo;
  return value;
}


/*
 * Helper class to hold the information and the histos for one Szintilator.
 * In particular, the informations are the 2 Pmts (up/down) with possible multihit.
 */
class Szinti
{

	// time difference Up-Down
	TH1F* h1mult;
	// time difference Up-Down
	TH1F* h1dt;
	// time difference Up-Down for laser signals (large amplitudes)
	TH1F* h1dtlaser;
	// mean time of Up-Down
	TH1F* h1tmean;
	// time difference Up-Down
	TH1F* h1pos;
	// amplitude difference Up-Down
	TH1F* h1dE;
	// amplitude Up vs Down
	TH2F* h2E;
    // amplitude Up and down vs z
    TH2F* h2EUpVsZ; TH2F* h2EUpSinVsZ; TH2F* h2EDownVsZ;
	// 'Energy loss' vs time difference
	TH2F* h2EnergyLostVsTimeDifference;

public:


	enum Type
	{
		A=0,
		B
	};

	Szinti(int id,Type type,GroupCamera* group)
	{
		this->id=id;this->type=type;this->group=group;
		reset();

		if (group->isExpert()) {


			stringstream histname,histtitle;

			histname << group->GetName()<<"_coninc_mult_"<<(type==A?"A":"B")<< id;
			histtitle << group->GetName()<<" multiplicity of Up and Down coincidences"<<(type==A?"A":"B")<< id;
			h1mult = new TH1F(histname.str().c_str(),histtitle.str().c_str(),16,0,16);
			this->group->addHistogram(h1mult);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_dt_"<<(type==A?"A":"B")<< id;
			histtitle << group->GetName()<<" dt of Up and Down "<<(type==A?"A":"B")<< id;
			//A.F h1dt = new TH1F(histname.str().c_str(),histtitle.str().c_str(),2000,-80,80);
			h1dt = new TH1F(histname.str().c_str(),histtitle.str().c_str(),800,-50,50);
			h1dt->GetXaxis()->SetTitle("t_{up} - t_{down}");
			this->group->addHistogram(h1dt);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_dt_"<<(type==A?"A":"B")<< id << "_laser";
			histtitle << group->GetName()<<" dt of Up and Down "<<(type==A?"A":"B")<< id <<" (laser)";
			//A.F h1dt = new TH1F(histname.str().c_str(),histtitle.str().c_str(),2000,-80,80);
			h1dtlaser = new TH1F(histname.str().c_str(),histtitle.str().c_str(),800,-50,50);
			h1dtlaser->GetXaxis()->SetTitle("t_{up} - t_{down}");
			this->group->addHistogram(h1dtlaser);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_tmean_"<<(type==A?"A":"B")<< id;
			histtitle << group->GetName()<<" t mean of Up and Down "<<(type==A?"A":"B")<< id;
			h1tmean = new TH1F(histname.str().c_str(),histtitle.str().c_str(),800,-2500,-1500);
			h1tmean->GetXaxis()->SetTitle("(t_{up} + t_{down})*0.5");
			this->group->addHistogram(h1tmean);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_zpos_"<<(type==A?"A":"B")<< id;
			histtitle << group->GetName()<<" z Position in cm "<<(type==A?"A":"B")<< id;
			//A.F h1dt = new TH1F(histname.str().c_str(),histtitle.str().c_str(),2000,-80,80);
			h1pos = new TH1F(histname.str().c_str(),histtitle.str().c_str(),800,-400,400);
			h1pos->GetXaxis()->SetTitle("zpos in cm");
			this->group->addHistogram(h1pos);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_dE_"<<(type==A?"A":"B") << id;
			histtitle << group->GetName()<<" dE of Up and Down "<<(type==A?"A":"B") << id;
			h1dE = new TH1F(histname.str().c_str(),histtitle.str().c_str(),400,-1000,1000);
			h1dE->GetXaxis()->SetTitle("E_{up} - E_{down}");
			this->group->addHistogram(h1dE);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_E_Up_vs_Down_"<<(type==A?"A":"B") << id;
			histtitle << group->GetName()<<" E Up vs Down "<<(type==A?"A":"B") << id;
			h2E = new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:100),0,4000,(HIGH_RES?1000:100),0,4000);
			h2E->SetOption("colz");
			h2E->GetXaxis()->SetTitle("E_{down}");
			h2E->GetYaxis()->SetTitle("E_{up}");
			this->group->addHistogram(h2E);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_dE_vs_dt_"<<(type==A?"A":"B")<<id;
			histtitle << group->GetName()<<" dE vs dt of Up and Down "<<(type==A?"A":"B")<<id;
			h2EnergyLostVsTimeDifference = new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:100),-40,40,(HIGH_RES?1000:100),-600,600);
			h2EnergyLostVsTimeDifference->SetOption("colz");
			h2EnergyLostVsTimeDifference->GetYaxis()->SetTitle("E_{up} - E_{down}");
			h2EnergyLostVsTimeDifference->GetXaxis()->SetTitle("t_{up} - t_{down}");
			this->group->addHistogram(h2EnergyLostVsTimeDifference);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_E_Up_vs_z_"<<(type==A?"A":"B") << id;
			histtitle << group->GetName()<<" E Up vs z "<<(type==A?"A":"B") << id;
			h2EUpVsZ = new TH2F(histname.str().c_str(),histtitle.str().c_str(),100,-200,200,100,0,4000);
			h2EUpVsZ->SetOption("scat");
			h2EUpVsZ->GetXaxis()->SetTitle("z");
			h2EUpVsZ->GetYaxis()->SetTitle("E_{up}");
			this->group->addHistogram(h2EUpVsZ);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_E_Down_vs_z_"<<(type==A?"A":"B") << id;
			histtitle << group->GetName()<<" E Down vs z "<<(type==A?"A":"B") << id;
			h2EDownVsZ = new TH2F(histname.str().c_str(),histtitle.str().c_str(),100,-200,200,100,0,4000);
			h2EDownVsZ->SetOption("scat");
			h2EDownVsZ->GetXaxis()->SetTitle("z");
			h2EDownVsZ->GetYaxis()->SetTitle("E_{Down}");
			this->group->addHistogram(h2EDownVsZ);
		}
	}


	~Szinti() {
	}

	void reset() {
		hits_up=NULL;
		hits_down=NULL;
	}

	// pointer to the group
	GroupCamera* group;
	// channel (wire)
	int id;
	// name of the Group object
	Type type;
	// hits in upstream PMT
	const PlaneCamera::digits *hits_up;
	// hits in downstream PMT
	const PlaneCamera::digits *hits_down;

	// this function is called when all digits are present
	void FillHistos(){

		// fill the per channel plot. each slab has 2 bins, meaning 0 to 0.5 is A0 down, 0.5 to 1 is A0 up
	  if( type == B && id==18 && (hits_up==NULL || hits_down==NULL) ) {
	    //std::cout<<"id: "<<id<<"  hits_up: "<<hits_up
		//     <<"    hits_down: "<<hits_down<<std::endl;
	    //getchar();
	  }

		if ( hits_up!=NULL) {
			group->chan_multi->Fill(type==A?id:group->A_CHANNEL+id,hits_up->size());

			// multiplicity of up/down coincidences
			int multi=0;
			for(PlaneCamera::digits::const_iterator hit_up=hits_up->begin();hit_up!=hits_up->end();++hit_up)
			{
				group->chan_ampl->Fill(type==A?id:group->A_CHANNEL+id,(*hit_up)->getMaxAmplitude());
				group->chan_time->Fill(type==A?id:group->A_CHANNEL+id,(*hit_up)->GetTimeDecoded());
				group->chan_int->Fill(type==A?id:group->A_CHANNEL+id,(*hit_up)->getIntegral());
				group->chan_hrtime->Fill(type==A?id:group->A_CHANNEL+id,(*hit_up)->getHiResTime());
				
				if ( (*hit_up)->getOpMode()==CS::ChipGandalf::GADC_DEBUG
						|| (*hit_up)->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL)
					group->chan_ftime->Fill(type==A?id:group->A_CHANNEL+id,(*hit_up)->getFrameTime() );

				// only fill the up/down plots when we have data from both sides
				if ( hits_down==NULL) continue;

				for(PlaneCamera::digits::const_iterator hit_down=hits_down->begin();hit_down!=hits_down->end();++hit_down)
				{
				  //std::cout<<"type = "<<type<<"    id = "<<id<<std::endl;
				  if(type==B && id == 18) {
				    //std::cout<<"B16: tup = "<<(*hit_up)->GetTimeDecoded() <<"    tdn = "<<(*hit_down)->GetTimeDecoded()
				    //	     <<"    dt = "<<(*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded() << std::endl;
				    //printf("B16: tup = %f    tdn = %f    dt = %f\n",(*hit_up)->GetTimeDecoded(),(*hit_down)->GetTimeDecoded(),(*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded());
				    //printf("B16: tup = %f    tdn = %f    dt = %f\n",(*hit_up)->getTime(),(*hit_down)->getTime(),(*hit_up)->getTime() - (*hit_down)->getTime());
				    //std::cout<<"     tup coarse = "<<(*hit_up)->getCoarseTimeMSB()<<" "<<(*hit_up)->getCoarseTimeLSB()
				    //	     <<"  tup hires = "<<( (float)(*hit_up)->getHiResTime() )/1024<<std::endl;
				    //getchar();
				  }
					
					
					
					if (type==A) {

						double zpos = ((*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded()  + Adiff_CALIB[id] )*(speedOfLightA)/2.;
						double tmean = ((*hit_up)->GetTimeDecoded() + (*hit_down)->GetTimeDecoded() + Amean_CALIB[id])*0.5;

						group->chan_zpos->Fill(id,zpos);
						group->chan_tmean->Fill(id,tmean);

						if (group->isExpert()) {
							h2EUpVsZ->Fill( zpos ,((int)(*hit_up)->getMaxAmplitude()));
							h2EDownVsZ->Fill( zpos ,((int)(*hit_down)->getMaxAmplitude()));
							h1dE->Fill( ((int)(*hit_up)->getMaxAmplitude()) - ((int)(*hit_down)->getMaxAmplitude()) );
							h1dt->Fill( (*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded() + Adiff_CALIB[id] );
							if( (int)(*hit_up)->getMaxAmplitude()>1000 && (int)(*hit_down)->getMaxAmplitude()>1000  )
							  h1dtlaser->Fill( (*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded() + Adiff_CALIB[id] );
							h1pos->Fill( zpos );
							h1tmean->Fill( tmean );
						}
                    }
					else {

						double zpos = ((*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded()  + Bdiff_CALIB[id] )*(speedOfLightB)/2.;
						double tmean = ((*hit_up)->GetTimeDecoded() + (*hit_down)->GetTimeDecoded() + Bmean_CALIB[id])*0.5;

						group->chan_zpos->Fill(group-> A_CHANNEL+id,zpos);
						group->chan_tmean->Fill(group-> A_CHANNEL+id,tmean);
						
						if (group->isExpert()) {
							h2EUpVsZ->Fill( zpos ,((int)(*hit_up)->getMaxAmplitude()));
							h2EDownVsZ->Fill( zpos ,((int)(*hit_down)->getMaxAmplitude()));
							h1dE->Fill( ((int)(*hit_up)->getMaxAmplitude()) - ((int)(*hit_down)->getMaxAmplitude()) );
							h1dt->Fill( (*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded() + Bdiff_CALIB[id] );
							if( (int)(*hit_up)->getMaxAmplitude()>100 && (int)(*hit_down)->getMaxAmplitude()>100  )
							  h1dtlaser->Fill( (*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded() + Bdiff_CALIB[id] );
							h1tmean->Fill( (*hit_up)->GetTimeDecoded() + (*hit_down)->GetTimeDecoded() + Bmean_CALIB[id] );
							h1pos->Fill( zpos );
							h1tmean->Fill( tmean );
						}
					}
					if (group->isExpert()) {
						h2E->Fill( ((int)(*hit_up)->getMaxAmplitude()),((int)(*hit_down)->getMaxAmplitude()) );
						h2EnergyLostVsTimeDifference->Fill( (*hit_up)->GetTimeDecoded() - (*hit_down)->GetTimeDecoded(),
							((int)(*hit_up)->getMaxAmplitude()) - ((int)(*hit_down)->getMaxAmplitude()));
					}
					multi++;
				}
			}
			if (group->isExpert()) h1mult->Fill(multi);
		}
		if ( hits_down!=NULL) {
			group->chan_multi->Fill(0.5 + (type==A?id:group->A_CHANNEL+id),hits_down->size());
			for(PlaneCamera::digits::const_iterator hit_down=hits_down->begin();hit_down!=hits_down->end();++hit_down)
			{
				group->chan_ampl->Fill( 0.5 + (type==A?id:group->A_CHANNEL+id),(*hit_down)->getMaxAmplitude());
				group->chan_time->Fill( 0.5 + (type==A?id:group->A_CHANNEL+id),(*hit_down)->GetTimeDecoded());
				group->chan_int->Fill( 0.5 + (type==A?id:group->A_CHANNEL+id),(*hit_down)->getIntegral());
				group->chan_hrtime->Fill( 0.5 + (type==A?id:group->A_CHANNEL+id),(*hit_down)->getHiResTime());

				if ( (*hit_down)->getOpMode()==CS::ChipGandalf::GADC_DEBUG
						|| (*hit_down)->getOpMode()==CS::ChipGandalf::GADC_DEBUG_IL)
					group->chan_ftime->Fill( 0.5 + (type==A?id:group->A_CHANNEL+id),(*hit_down)->getFrameTime());
			}
		}
	}
};


/**
 * Extension of the Szinti class for an A-Ring Szintilator. This class has
 * 'connections' to neighboring B Szintilators. With them, the coincidence is checked.
 */
class ASzinti : public Szinti
{

	// vector holding the associated B Szintilators
	vector<Szinti*> connections;
	// tof plots for each combination
	map<int,TH1F*> h1Coinc;
	// tof plots for each combination
	map<int,TH1F*> h1Tof;
	// tof plots for each combination (laser signals)
	map<int,TH1F*> h1TofLaser;
	// energy plots for each combination
	map<int,TH2F*> h2EnergyAVsEnergyB;
	// tof vs energy in B plots for each combination
	map<int,TH2F*> h2EnergyBVsTof;
	// beta vs energy in B plots for each combination
	map<int,TH2F*> h2EnergyBVsBeta;
	// beta vs integral in B plots for each combination
	map<int,TH2F*> h2IntegralBVsBeta;
	// beta vs energy in A plots for each combination
	map<int,TH2F*> h2EnergyAVsBeta;
	// beta vs integral in A plots for each combination
	map<int,TH2F*> h2IntegralAVsBeta;
	// beta vs integral in B plots for each combination
	map<int,TH2F*> h2PosAVsPosB;
	// Eup * sin(theta) vs z in B plots for each combination
	map<int,TH2F*> h2EUpSinVsZ;

public:
	ASzinti(int id, GroupCamera* group) : Szinti(id,A,group) {}
	~ASzinti() {
		connections.clear();
		h1Tof.clear();
		h2EnergyAVsEnergyB.clear();
		h2EnergyBVsTof.clear();
		h2EnergyBVsBeta.clear();
		h2IntegralBVsBeta.clear();
		h2EnergyAVsBeta.clear();
		h2IntegralAVsBeta.clear();
		h2PosAVsPosB.clear();
		h2EUpSinVsZ.clear();
	}

	void FillHistos(){

	  Szinti::FillHistos();

		// if there is no hit at one side, we can't do the coinc plots
		if ( hits_up==NULL || hits_down==NULL ) return;

		// vector holding the time
		vector< GroupCamera::Hit* > hits;
		vector< GroupCamera::Hit* > hits_laser;
		// store coincidences in this szntilator
		for(PlaneCamera::digits::const_iterator hit_up=hits_up->begin();hit_up!=hits_up->end();++hit_up){
			for(PlaneCamera::digits::const_iterator hit_down=hits_down->begin();hit_down!=hits_down->end();++hit_down)
			{
				hits.push_back(new GroupCamera::Hit( id , *hit_up , *hit_down, Adiff_CALIB[id] , speedOfLightA ) );
							
				if( (int)(*hit_up)->getMaxAmplitude()>1000 && (int)(*hit_down)->getMaxAmplitude()>1000  ) {
				  hits_laser.push_back(new GroupCamera::Hit( id , (*hit_up) , (*hit_down), Adiff_CALIB[id] , speedOfLightA) );
				}
			}
		}


		// iterate through neighboring b szintilators and search for coincidences.

		// store the multiplicity for each coinc with B szinti
		map<int,int> multi;
		for(map<int,TH1F*>::iterator plot=h1Coinc.begin();plot!=h1Coinc.end();++plot)
			multi[(*plot).first] = 0;

		for(vector<Szinti*>::iterator b=connections.begin();b!=connections.end();++b) {

			if ( (*b)->hits_up==NULL || (*b)->hits_down==NULL ) continue;

			for(PlaneCamera::digits::const_iterator hit_up=(*b)->hits_up->begin();hit_up!=(*b)->hits_up->end();++hit_up){

				for(PlaneCamera::digits::const_iterator hit_down=(*b)->hits_down->begin();hit_down!=(*b)->hits_down->end();++hit_down){

					if (group->isExpert()) {
						

						for(vector< GroupCamera::Hit* >::iterator i=hits_laser.begin();i!=hits_laser.end();++i)
						{

						  if( (int)(*hit_up)->getMaxAmplitude()<100 || (int)(*hit_down)->getMaxAmplitude()<100  ) continue;
						  double tBup = (*hit_up)->GetTimeDecoded();
						  double tBdown = (*hit_down)->GetTimeDecoded();

						  //double tBup = (*hit_up)->getTime()*(*hit_up)->GetTimeUnit();
						  //double tBdown = (*hit_down)->getTime()*(*hit_down)->GetTimeUnit();

						  // TODO: is there an easear way to enumerate this calib constant ?
						  double tof = (tBup+tBdown)/2. - (*i)->mean + TOF_CALIB[id][ Modulo( ((*b)->id - id) , group->B_CHANNEL ) ];

						  if (h1TofLaser.find((*hit_up)->getChannel())!=h1TofLaser.end())
							h1TofLaser[(*hit_up)->getChannel()]->Fill( tof );
						  else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

						}
					}


					for(vector< GroupCamera::Hit* >::iterator a=hits.begin();a!=hits.end();++a)
					{
						
						GroupCamera::Hit *Bhit = new GroupCamera::Hit((*b)->id,(*hit_up),(*hit_down),Bdiff_CALIB[id],speedOfLightB,B_ZOFFSET);

						// TODO: is there an easear way to enumerate this calib constant ?
						double tof = Bhit->mean - (*a)->mean + TOF_CALIB[id][Modulo( (*b)->id - id , group->B_CHANNEL )];
                        double sinTheta=TMath::Sin(TMath::ATan((110.-25.)/(Bhit->zpos-(*a)->zpos)));
						double beta = TMath::Sqrt(pow((double)(110-25),(double)2)+pow((double)(Bhit->zpos-(*a)->zpos),(double)2))/tof/30.;

						group->protonCandidates.push_back(new GroupCamera::ProtonEvent((*a),Bhit,tof,beta));

						// print frame plots if condition is true and we have a debug event
						if (FRAME_PLOTS && beta==0.5) {
//							map<int,*PlaneGandalf>::iterator g_plane = group->planes[PlaneCamera::Aup]->g_planes.find( (*i)->chanUp );
//							if( g_plane != group->planes[PlaneCamera::Aup]->g_planes.end() ) {
//								TH1F *framePlot = (*g_plane).second->fFramePlotHist;
//								stringstream name;
//								name << "A"<< id << "_up_Frame_Plot_event_" << group->event->GetEventNumberInRun() <<" "<< multi[(*hit_up)->getChannel()];
//								TH1F *tmp = new TH1F(*framePlot);
//								tmp->SetTitle(name.str().c_str());
//								group->addHistogram(tmp);
//							}
//							if (group->planes[PlaneCamera::Adown]->g_planes.find( (*i)->chanDown )!=group->planes[PlaneCamera::Adown]->g_planes.end()
//									&& group->planes[PlaneCamera::Adown]->g_planes[(*i)->chanDown]->isDebug) {
//								TH1F *framePlot = group->planes[PlaneCamera::Adown]->g_planes[( (*i)->chanDown )]->fFramePlotHist;
//								stringstream name;
//								name << "A"<< id << "_down_Frame_Plot_event_" << group->event->GetEventNumberInRun()<<" "<< multi[(*hit_up)->getChannel()];
//								TH1F *tmp = new TH1F(*framePlot);
//								tmp->SetTitle(name.str().c_str());
//								group->addHistogram(tmp);
//							}
//							if (group->planes[PlaneCamera::Bup]->g_planes.find( (*hit_up)->getChannel() )!=group->planes[PlaneCamera::Bup]->g_planes.end()
//									&& group->planes[PlaneCamera::Bup]->g_planes[(*hit_up)->getChannel()]->isDebug) {
//								TH1F *framePlot = group->planes[PlaneCamera::Bup]->g_planes[( (*hit_up)->getChannel() )]->fFramePlotHist;
//								stringstream name;
//								name << "B"<< (*b)->id << "_up_Frame_Plot_event_" << group->event->GetEventNumberInRun()<<" "<< multi[(*hit_up)->getChannel()];
//								TH1F *tmp = new TH1F(*framePlot);
//								tmp->SetTitle(name.str().c_str());
//								group->addHistogram(tmp);
//							}
//							if (group->planes[PlaneCamera::Bdown]->g_planes.find( (*hit_down)->getChannel() )!=group->planes[PlaneCamera::Bdown]->g_planes.end()
//									&& group->planes[PlaneCamera::Bdown]->g_planes[(*hit_down)->getChannel()]->isDebug) {
//								TH1F *framePlot = group->planes[PlaneCamera::Bdown]->g_planes[( (*hit_down)->getChannel() )]->fFramePlotHist;
//								stringstream name;
//								name << "B"<< (*b)->id << "_down_Frame_Plot_event_" << group->event->GetEventNumberInRun()<<" "<< multi[(*hit_up)->getChannel()];
//								TH1F *tmp = new TH1F(*framePlot);
//								tmp->SetTitle(name.str().c_str());
//								group->addHistogram(tmp);
//							}
						}

						group->chan_tof->Fill( id+((*b)->id-id)*0.5, tof );

						if (group->isExpert()) {
							if (h1Tof.find((*hit_up)->getChannel())!=h1Tof.end())
								h1Tof[(*hit_up)->getChannel()]->Fill( tof );
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2EnergyAVsEnergyB.find((*hit_up)->getChannel())!=h2EnergyAVsEnergyB.end())
								h2EnergyAVsEnergyB[(*hit_up)->getChannel()]->Fill( Bhit->energy , (*a)->energy );
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2EnergyBVsTof.find((*hit_up)->getChannel())!=h2EnergyBVsTof.end())
								h2EnergyBVsTof[(*hit_up)->getChannel()]->Fill( tof , Bhit->energy ) ;
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2EnergyBVsBeta.find((*hit_up)->getChannel())!=h2EnergyBVsBeta.end() )
								h2EnergyBVsBeta[(*hit_up)->getChannel()]->Fill( beta , Bhit->energy ) ;
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2IntegralBVsBeta.find((*hit_up)->getChannel())!=h2IntegralBVsBeta.end())
								h2IntegralBVsBeta[(*hit_up)->getChannel()]->Fill( beta , Bhit->energyInt ) ;
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2EnergyAVsBeta.find((*hit_up)->getChannel())!=h2EnergyAVsBeta.end())
								h2EnergyAVsBeta[(*hit_up)->getChannel()]->Fill( beta , (*a)->energy ) ;
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2IntegralAVsBeta.find((*hit_up)->getChannel())!=h2IntegralAVsBeta.end())
								h2IntegralAVsBeta[(*hit_up)->getChannel()]->Fill( beta , (*a)->energyInt ) ;
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2PosAVsPosB.find((*hit_up)->getChannel())!=h2PosAVsPosB.end())
								h2PosAVsPosB[(*hit_up)->getChannel()]->Fill( Bhit->zpos , (*a)->zpos ) ;
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;

							if (h2EUpSinVsZ.find((*hit_up)->getChannel())!=h2EUpSinVsZ.end())
								h2EUpSinVsZ[(*hit_up)->getChannel()]->Fill( Bhit->zpos, (*hit_up)->getMaxAmplitude()*sinTheta ) ;
							else std::cout<<"GroupCamera::ERROR::there was a mixup in channel id:   A = "<<id<<";  B = " <<(*hit_up)->getChannel() << std::endl;
						}
						multi[(*hit_up)->getChannel()]++;
					}
				}
			}
		}
		if (group->isExpert())
			for(map<int,TH1F*>::iterator plot=h1Coinc.begin();plot!=h1Coinc.end();++plot)
				h1Coinc[(*plot).first]->Fill(multi[(*plot).first]);


	}

	/*
	 * add a connection to this szintilator. In case of the Prototype, each A szintialtor has
	 * 3 neighboring B szintilators it can make coincidences with. In case of Camera, there are 2.
	 */
	void addConnection(Szinti* bszinti)
	{

		if (group->isExpert()) {


			stringstream histname,histtitle;

			histname << group->GetName()<<"_coinc_mult_A" <<id<<" and B"<<bszinti->id;
			histtitle << group->GetName()<<" multiplicity of A" <<id<<" and B"<<bszinti->id << " coincidences";
			h1Coinc[bszinti->id] = new TH1F(histname.str().c_str(),histtitle.str().c_str(),16,0,16);
			group->addHistogram(h1Coinc[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_tof_A"<<id<<"_to_B"<<bszinti->id;
			histtitle << group->GetName()<<" Time of flight from A"<<id<<" to B"<<bszinti->id;
			h1Tof[bszinti->id]=(new TH1F(histname.str().c_str(),histtitle.str().c_str(),800,-120,120));
			h1Tof[bszinti->id]->GetXaxis()->SetTitle("t_{B} - t_{A}");
			group->addHistogram(h1Tof[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_tof_A"<<id<<"_to_B"<<bszinti->id<<"_laser";
			histtitle << group->GetName()<<" Time of flight form A"<<id<<" to B"<<bszinti->id << "(laser)";
			h1TofLaser[bszinti->id]=(new TH1F(histname.str().c_str(),histtitle.str().c_str(),800,-120,120));
			h1TofLaser[bszinti->id]->GetXaxis()->SetTitle("t_{B} - t_{A}");
			group->addHistogram(h1TofLaser[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName() << "_E_A"<<id<<"_vs_B"<<bszinti->id;
			histtitle << group->GetName() << " Energy deposit A"<<id<<" vs B"<<bszinti->id;
			h2EnergyAVsEnergyB[bszinti->id]=(new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:200),0,4000,(HIGH_RES?1000:200),0,4000));
			h2EnergyAVsEnergyB[bszinti->id]->SetOption("colz");
			h2EnergyAVsEnergyB[bszinti->id]->GetXaxis()->SetTitle("\\sqrt{E_{Aup}*E_{Adown}}");
			h2EnergyAVsEnergyB[bszinti->id]->GetYaxis()->SetTitle("\\sqrt{E_{Bup}*E_{Bdown}}");
			group->addHistogram(h2EnergyAVsEnergyB[bszinti->id]);
			histname.str("");histtitle.str("");
			histname << group->GetName() << "_E_B"<<bszinti->id<<"_vs_tof_to_A"<<id;
			histtitle << group->GetName() << " Energy deposit in B"<<bszinti->id<<" vs ToF to A "<<id;
			h2EnergyBVsTof[bszinti->id]=(new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:200),-50,50,(HIGH_RES?1000:200),0,4000));
			//h2EnergyBVsTof[bszinti->id]->SetOption("colz");
			h2EnergyBVsTof[bszinti->id]->GetXaxis()->SetTitle("t_{B} - t_{A}");
			h2EnergyBVsTof[bszinti->id]->GetYaxis()->SetTitle("\\sqrt{E_{Bup}*E_{Bdown}}");
			group->addHistogram(h2EnergyBVsTof[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName() << "_E_B"<<bszinti->id<<"_vs_beta_A"<<id;
			histtitle << group->GetName() << " Energy deposit in B"<<bszinti->id<<" vs Beta (A "<<id<<")";
			h2EnergyBVsBeta[bszinti->id]=(new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:100),0,1.2,(HIGH_RES?1000:100),0,3500));
			//h2EnergyBVsBeta[bszinti->id]->SetOption("colz");
			h2EnergyBVsBeta[bszinti->id]->GetXaxis()->SetTitle("\\beta");
			h2EnergyBVsBeta[bszinti->id]->GetYaxis()->SetTitle("\\sqrt{E_{Bup}*E_{Bdown}}");
			group->addHistogram(h2EnergyBVsBeta[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName() << "_Int_B"<<bszinti->id<<"_vs_beta_A"<<id;
			histtitle << group->GetName() << " Integral in B"<<bszinti->id<<" vs Beta (A "<<id<<")";
			h2IntegralBVsBeta[bszinti->id]=(new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:100),0,1.2,(HIGH_RES?1000:100),0,100000));
			//h2IntegralBVsBeta[bszinti->id]->SetOption("colz");
			h2IntegralBVsBeta[bszinti->id]->GetXaxis()->SetTitle("\\beta");
			h2IntegralBVsBeta[bszinti->id]->GetYaxis()->SetTitle("\\sqrt{Int_{Bup}*Int_{Bdown}}");
			group->addHistogram(h2IntegralBVsBeta[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName() << "_E_A"<<id<<"_vs_beta_B"<<bszinti->id;
			histtitle << group->GetName() << " Energy deposit in A"<<id <<" vs Beta B"<<bszinti->id;
			h2EnergyAVsBeta[bszinti->id]=(new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:100),0,1.2,(HIGH_RES?1000:100),0,3500));
			//h2EnergyAVsBeta[bszinti->id]->SetOption("colz");
			h2EnergyAVsBeta[bszinti->id]->GetXaxis()->SetTitle("\\beta");
			h2EnergyAVsBeta[bszinti->id]->GetYaxis()->SetTitle("\\sqrt{E_{Bup}*E_{Bdown}}");
			group->addHistogram(h2EnergyAVsBeta[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName() << "_Int_A"<<id<<"_vs_beta_B"<<bszinti->id;
			histtitle << group->GetName() << " Integral deposit in A"<<id <<" vs Beta B"<<bszinti->id;
			h2IntegralAVsBeta[bszinti->id]=(new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:100),0,1.2,(HIGH_RES?1000:100),0,100000));
			//h2IntegralAVsBeta[bszinti->id]->SetOption("colz");
			h2IntegralAVsBeta[bszinti->id]->GetXaxis()->SetTitle("\\beta");
			h2IntegralAVsBeta[bszinti->id]->GetYaxis()->SetTitle("\\sqrt{Int_{Bup}*Int_{Bdown}}");
			group->addHistogram(h2IntegralAVsBeta[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName() << "_zpos_A"<<id<<"_vs_zpos_B"<<bszinti->id;
			histtitle << group->GetName() << "z Position in A"<<id <<" vs z Position in B"<<bszinti->id;
			h2PosAVsPosB[bszinti->id]=(new TH2F(histname.str().c_str(),histtitle.str().c_str(),(HIGH_RES?1000:100),-300,300,(HIGH_RES?1000:100),-400,400));
			h2PosAVsPosB[bszinti->id]->SetOption("colz");
			h2PosAVsPosB[bszinti->id]->GetXaxis()->SetTitle("\\posB");
			h2PosAVsPosB[bszinti->id]->GetYaxis()->SetTitle("\\posA");
			group->addHistogram(h2PosAVsPosB[bszinti->id]);

			histname.str("");histtitle.str("");
			histname << group->GetName()<<"_E_Up_Sin_vs_z_A"<<id<<"_B"<<bszinti->id;
			histtitle << group->GetName()<<" E Up Sin theta vs z in A"<<id<<" and B"<<bszinti->id;
			h2EUpSinVsZ[bszinti->id] = new TH2F(histname.str().c_str(),histtitle.str().c_str(),100,-200,200,100,-400,400);
			h2EUpSinVsZ[bszinti->id]->SetOption("scat");
			h2EUpSinVsZ[bszinti->id]->GetXaxis()->SetTitle("z");
			h2EUpSinVsZ[bszinti->id]->GetYaxis()->SetTitle("E_{up} * sin (theta)");
			group->addHistogram(h2EUpSinVsZ[bszinti->id]);
		}

		connections.push_back(bszinti);
	}

};



GroupCamera::GroupCamera(const char* name): Group(name) {
	initialized=false;
}

void GroupCamera::Init(void)
{
	Group::Init();

	// get the planes
	for( vector<const Plane*>::const_iterator pp=fPlanes.begin(); pp!=fPlanes.end(); pp++ ) {
		const PlaneCamera *plane = dynamic_cast<const PlaneCamera *>(*pp);
		if (plane==NULL) continue;
		planes[ (int)plane->getType() ] = plane;
	}

	// check if we have found all 4 planes (A_up,A_down,B_up,B_down)
	// this is mandatory !
	for(int i=0;i<4;i++)
	{
		if (planes.find(i)==planes.end())
			return;//			throw "GroupCamera::Init::Could not find all planes";

	}

	// set the number of channels depending on the number of channels in the plane objects
	if (planes[PlaneCamera::Aup]->N_CHANNEL != planes[PlaneCamera::Adown]->N_CHANNEL )
		throw "GroupCamera::Init::'A' Planes have different channel size";
	A_CHANNEL=planes[PlaneCamera::Aup]->N_CHANNEL;

	if (planes[PlaneCamera::Bup]->N_CHANNEL != planes[PlaneCamera::Bdown]->N_CHANNEL )
		throw "GroupCamera::Init::'B' Planes have different channel size";
	B_CHANNEL=planes[PlaneCamera::Bup]->N_CHANNEL;

	// check channel size
//#ifdef PROTO
//	if (A_CHANNEL != 12) throw "GroupCamera::Init::'A' Planes number of channels differs from expected value (12)";
//#else
	if (A_CHANNEL != 24) throw "GroupCamera::Init::'A' Planes number of channels differs from expected value (24)";
//#endif
	if (B_CHANNEL != 24) throw "GroupCamera::Init::'B' Planes number of channels differs from expected value (24)";

#if USE_DATABASE == 1
	// read the calibrations
	ReadCalib();
#endif //USE_DATABASE

	// initialize histos
	stringstream name;
	name << fName << "_channelMultiplicity";
	chan_multi = new TH1F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL)*2,0,A_CHANNEL+B_CHANNEL);
	AddHistogram(chan_multi);

	name.str("");
	name << fName << "_channelAmplitude";
	chan_ampl = new TH2F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL)*2,0,A_CHANNEL+B_CHANNEL,(HIGH_RES?1000:200),0,4096);
	chan_ampl->SetOption("colz");
	chan_ampl->GetXaxis()->SetTitle("Channel");
	chan_ampl->GetYaxis()->SetTitle("ADC LSB");
	AddHistogram(chan_ampl);

	name.str("");
	name << fName << "_channelIntegral";
	chan_int = new TH2F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL)*2,0,A_CHANNEL+B_CHANNEL,(HIGH_RES?1000:200),0,60000);
	chan_int->SetOption("colz");
	AddHistogram(chan_int);

	name.str("");
	name << fName << "_channelTime";
	chan_time = new TH2F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL)*2,0,A_CHANNEL+B_CHANNEL,(HIGH_RES?1000:100),-2500,-1500);
	chan_time->SetOption("colz");
	chan_time->GetXaxis()->SetTitle("Channel");
	chan_time->GetYaxis()->SetTitle("ns");
	AddHistogram(chan_time);

	name.str("");
	name << fName << "_channelFrameTime";
	chan_ftime = new TH2F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL)*2,0,A_CHANNEL+B_CHANNEL,(HIGH_RES?1000:100),0,400);
	chan_ftime->SetOption("colz");
	chan_ftime->GetXaxis()->SetTitle("Channel");
	chan_ftime->GetYaxis()->SetTitle("n'th sample");
	AddHistogram(chan_ftime);

	name.str("");
	name << fName << "_channelHiResTime";
	chan_hrtime = new TH2F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL)*2,0,A_CHANNEL+B_CHANNEL,(HIGH_RES?1024:100),0,1024);
	chan_hrtime->SetOption("colz");
	chan_hrtime->GetXaxis()->SetTitle("Channel");
	chan_hrtime->GetYaxis()->SetTitle("HiRes LSB");
	AddHistogram(chan_hrtime);

	name.str("");
	name << fName << "_channelZpos";
	chan_zpos = new TH2F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL),0,A_CHANNEL+B_CHANNEL,(HIGH_RES?1000:100),-400,400);
	chan_zpos ->SetOption("colz");
	chan_zpos ->GetXaxis()->SetTitle("Channel");
	chan_zpos ->GetYaxis()->SetTitle("zpos in cm");
	AddHistogram(chan_zpos);

	name.str("");
	name << fName << "_channelTmean";
	chan_tmean = new TH2F(name.str().c_str(),name.str().c_str(),(A_CHANNEL+B_CHANNEL),0,A_CHANNEL+B_CHANNEL,(HIGH_RES?1000:100),-2500,-1500);
	chan_tmean ->SetOption("colz");
	chan_tmean ->GetXaxis()->SetTitle("Channel");
	chan_tmean ->GetYaxis()->SetTitle("(t_{A} + t_{B})*0.5");
	AddHistogram(chan_tmean);

	name.str("");
	name << fName << "_channelTOF";
	chan_tof = new TH2F(name.str().c_str(),name.str().c_str(),A_CHANNEL*2,0,A_CHANNEL,(HIGH_RES?1000:100),-120,120);
	chan_tof ->SetOption("colz");
	chan_tof ->GetXaxis()->SetTitle("Channel");
	chan_tof ->GetYaxis()->SetTitle("t_{B} - t_{A}");
	AddHistogram(chan_tof);

	name.str("");
	name << fName << "_number_of_tracks";
	track_mult = new TH1F(name.str().c_str(),name.str().c_str(),30,0,30);
	track_mult ->GetXaxis()->SetTitle("Tracks");
	AddHistogram(track_mult);

	// create szintilator objects holding the histograms
	for(int i=0;i<A_CHANNEL;i++)
		aszintis[i]=new ASzinti(i,this);

	for(int i=0;i<B_CHANNEL;i++)
	{
		bszintis[i]= new Szinti(i,Szinti::B,this);
		// create connections from A szinti to B szinti;
		// TODO: generate documentation for the labeling and the 'connections' I'm doing...

//#ifdef PROTO
//		if (i%2 == 0)
//			aszintis[i/2]->addConnection(bszintis[i]);
//		else {
//			aszintis[Modulo( (i-1)/2 ,A_CHANNEL)]->addConnection(bszintis[i]);
//			aszintis[Modulo( (i+1)/2 ,A_CHANNEL)]->addConnection(bszintis[i]);
//		}
//#else
		aszintis[Modulo(i-1,A_CHANNEL)]->addConnection(bszintis[i]);
		aszintis[i]->addConnection(bszintis[i]);
//#endif
	}

	initialized=true;
}


void GroupCamera::EndEvent(const CS::DaqEvent &event)
{
	if (!initialized) return;

	this->event = &event;

	Group::EndEvent();

	// reset all szinti objects
	for(map<int,ASzinti*>::iterator sz=aszintis.begin();sz!=aszintis.end();++sz)
		(*sz).second->reset();
	for(map<int,Szinti*>::iterator sz=bszintis.begin();sz!=bszintis.end();++sz)
		(*sz).second->reset();

	// link hits from planes to objects
	for(map<int, PlaneCamera::digits >::const_iterator chan=planes[PlaneCamera::Aup]->chanDigits.begin();
			chan!=planes[PlaneCamera::Aup]->chanDigits.end();++chan)
		aszintis[(*chan).first]->hits_up = &((*chan).second);

	for(map<int, PlaneCamera::digits >::const_iterator chan=planes[PlaneCamera::Bup]->chanDigits.begin();
			chan!=planes[PlaneCamera::Bup]->chanDigits.end();++chan)
		bszintis[(*chan).first]->hits_up = &((*chan).second);

	for(map<int, PlaneCamera::digits >::const_iterator chan=planes[PlaneCamera::Adown]->chanDigits.begin();
			chan!=planes[PlaneCamera::Adown]->chanDigits.end();++chan)
		aszintis[(*chan).first]->hits_down = &((*chan).second);

	for(map<int, PlaneCamera::digits >::const_iterator chan=planes[PlaneCamera::Bdown]->chanDigits.begin();
			chan!=planes[PlaneCamera::Bdown]->chanDigits.end();++chan)
		bszintis[(*chan).first]->hits_down = &((*chan).second);

	// fill histos and reconstruct proton candidates
	for(map<int,Szinti*>::iterator b=bszintis.begin();b!=bszintis.end();++b)
		(*b).second->FillHistos();
	for(map<int,ASzinti*>::iterator a=aszintis.begin();a!=aszintis.end();++a)
		(*a).second->FillHistos();



//	std::cout << "Found " << protonCandidates.size() << " possible Protons in event "<<event.GetEventNumberInRun() << std::endl;

	track_mult->Fill(protonCandidates.size());
	for(vector<ProtonEvent*>::iterator candidate=protonCandidates.begin();candidate != protonCandidates.end();++candidate) {
//		std::cout << " 			Aid: " << (*candidate)->Ahit->id << "     Bid: " << (*candidate)->Bhit->id << "   beta: << (*candidate)->beta << std::endl;
	}


	// clear old proton candidates at the end
	for(vector<ProtonEvent*>::iterator candidate=protonCandidates.begin();candidate != protonCandidates.end();++candidate)
		delete (*candidate);
	protonCandidates.clear();

}



///////////////////////////////////////////////////////////////////////////////
//////////////////////////// Calibration Database


#if USE_DATABASE == 1
void GroupCamera::ReadCalib()
{
  // if the calibTime is not set, we should not use calibrations
  if((planes[0]->calibTime)==0)
    return;

  tm * t;
  t = localtime ( &(planes[0]->calibTime) );

  std::cout<<"GroupCamera::ReadCalib() ==> "<<this->GetName()<<" reading calibrations !"<<std::endl;
  // read-in corresponding calibration constants
  try{
    ReadFromDataBase(calib_data,*t);

    if(calib_data.size() != (unsigned) B_CHANNEL) {
          std::cerr<<"Size of Calibration File is not correct ! Should be : "
      	    <<B_CHANNEL<<" Is "<<calib_data.size()<<" "
      	    <<t->tm_mday<<"."<<t->tm_mon+1<<"."<<t->tm_year+1900<<" "
            <<t->tm_hour<<":"<<t->tm_min<<":"<<t->tm_sec<<std::endl;
    }
    else {
      //let's put the calibration values into the arrays
      std::cerr<<"# chan  Amean     Bmean     Adiff     Bdiff     TOF_AiBi  TOF_AiBj (j=i+1)"<< std::endl;
      for(vector<CameraCalib>::iterator c=calib_data.begin();c!=calib_data.end();++c) {
        std::cerr<< (*c).ch<<": "<< (*c).Amean<<" "<< (*c).Bmean<<" "<< 
           (*c).Adiff<<" "<< (*c).Bdiff<<" "<< (*c).TOF_AiBi<<" "<< (*c).TOF_AiBj<<" "<< std::endl;
        
        Amean_CALIB[(*c).ch] = (*c).Amean;
        Bmean_CALIB[(*c).ch] = (*c).Bmean;
        Adiff_CALIB[(*c).ch] = (*c).Adiff;
        Bdiff_CALIB[(*c).ch] = (*c).Bdiff;
        TOF_CALIB[(*c).ch][0] = (*c).TOF_AiBi;
        TOF_CALIB[(*c).ch][1] = (*c).TOF_AiBj;

      }
    }

  }

  catch(CS::Exception& e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(const std::exception &e) {
    std::cerr<<e.what()<<std::endl;
  }
  catch(...) {
    std::cout<<"GroupCamera::ReadCalib() ==> "<<GetName()
	<<" calibrations, valid for ";
    std::cout<<t->tm_mday<<"."<<t->tm_mon+1<<"."<<t->tm_year+1900<<" "
	<<t->tm_hour<<":"<<t->tm_min<<":"<<t->tm_sec
	<<", not found in DB"<<std::endl;
  }

}

#endif //USE_DATABASE

