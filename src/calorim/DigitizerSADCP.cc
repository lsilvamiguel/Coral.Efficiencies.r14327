#include "DigitizerSADCP.h"
////////////////////////////////////////////////////////////////////////////////
DigitizerSADCP::DigitizerSADCP ( CsCalorimeter* parent, Type type )
              : DigitizerSADCBase(parent, type),
                mypulse_(NULL), plane_(-1) {}

// Destructor
DigitizerSADCP::~DigitizerSADCP  (void) {
  if ( fit_ofile_ )
    free(fit_ofile_);
};

void DigitizerSADCN::PrintType (){
  std::cout<<"PPPPPP\n";
}

////////////////////////////////////////////////////////////////////////////////
/*!
  \brief Simple function to extract signal using cfd algorithm
  \param sample {A baseline subtracted list of samples}
  \return {Pointer to the result. Null on no signal. (Has to be delted after use.)}
*/
CsDigitizerSADCP::cfd_result_t*
CsDigitizerSADCP::CFD(const CsDigitizerSADC::cfd_option_t &options,
                   const std::vector<unsigned int> &sample ){

  // Point to be returned initialized at first founfd signal
  cfd_result_t *res(NULL);

  // CFD implementation
  {

    // To remeber last slice
    float last_result(0);

    // Loop over valid samples defined by delay and sample.size()
    // Can be further restricrted by user seting smin and smax
    for( unsigned int i = max(options.delay, options.smin);
         i < min(sample.size(),(size_t)options.smax);
         i++) {

      const float this_result = sample[i] - sample[i - options.delay] * options.ampl;
      // CFD condition
      if( last_result > 0 && this_result <= 0 ) {
        // Threshold check
        if ( sample[i] > options.thr ) {

          // We have a signal so lets check if res exists
          if( !res ) {
            res = new cfd_result_t;
            res->ampl=0;
          }

          // Check for biggest signal
          if ( res->ampl < sample[i] ) {
            res->ampl = sample[i];
            res->time = i + ( this_result / (last_result - this_result) );
          }

        }
      }

      last_result = this_result;

    }

  }

  return res;


};

vector<short unsigned int> Correct_Background_Exp(const vector<short unsigned int> &_sample) {

/// \todo  This function only works for hardware based Baseline suptraction (at the moment)
/// \todo  TODO make it configurable in case of advantages

  // Time constant
  double tau[2];
  vector<short unsigned int> sample = _sample;


  // Baseline set by frontends
  short unsigned int baseline = 50;

  // Subtract baseline
  for( vector<short unsigned int>::iterator it = sample.begin();
       it != sample.end(); it++)
    *it -= baseline;

/// \todo  TODO replace tau extraction by fit or by parameter
/// \todo  TODO include error handling

  // Get tau for odd and even (should be identical except for errors)
  tau[0] = -2. / log( sample[2] / sample[0]);
  tau[1] = -2. / log( sample[3] / sample[1]);
  double Tau = ( tau[0] + tau[1] ) / 2.;

  // Correct samples again
  for( unsigned int i = 1;
       i < sample.size(); i++) {
    int corr = (int)round( sample[0] * exp( - i / Tau ) );
    sample[i] -= corr ;
  }
  return sample;

}

//////////////////////////////////////////////////////////
/// \todo  CsDigitizerSADC::PulseFit is not yet in production status!!!!!
const CsDigitizerSADC::result_t &
CsDigitizerSADC::PulseFit(const std::vector<uint16> &sample, const unsigned int icell) {

  assert( parent_->GetCalID() != 0 );


  const int id(parent_->GetCalID()),
    x(parent_->GetColumnOfCell( icell )),
    y(parent_->GetRowOfCell( icell ));

  ClearResult();

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,00)

  // Initilize spline
  if( !PulseManager::Instance()->GetPulse(id,0) ) {

    {
      list<string> spline_name;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "SPLINE_FILE", spline_name ) ) {
        unsigned int i = 0;
        for( list<string>::iterator it = spline_name.begin(); it != spline_name.end(); ++it) {
          try {
            cout << "CsDigitizerSADC::PulseFit: Register " << *it
                 << " as " <<  parent_->GetName() << ":" << i << endl;
            MyPulse *mypulse = new MyPulse(it->c_str());
            if( mypulse )
              PulseManager::Instance()->Insert(id, i, mypulse );
            else
              throw -66;
            i++;
          } catch ( int &err ) {
            /*
            for( vector< MyPulse*>::iterator it=parent_->mypulse_.begin();
                 it!=parent_->mypulse_.end(); ++it) {
              if (*it) {
                delete *it;
                *it = NULL;
              }
            }
            parent_->mypulse_.clear();
            */
            throw Reco::Exception("%s in %s at %i : Catched %i while init mypulse_\n",
                                  __FILE__, __FUNCTION__, __LINE__, err);
          }
        }
      } else {
        throw Reco::Exception("%s in %s at %i : Uncomplete configuration : No SPLINE_FILE_NAME for %s\n",
                              __FILE__, __FUNCTION__, __LINE__, parent_->GetName().c_str());
      }

    }

  }


  // Get other configurations
  //cout << DAQDetName << endl;
  if( plane_ == -1 ) {
    if( id == 1 ) // ECAL01
      plane_ = atoi(parent_->GetName().c_str());
    else if ( id == 2 ) { //ECAL02
      float xcell = x - 30.7, ycell = y - 23.4;
      if ( ( xcell > -23 && xcell < -15 && ycell >- 5  && ycell < 5 )
           ||( xcell > -15 && xcell <-11 && ycell > -8  && ycell < 8 )
           ||( xcell > -11 && xcell < 13 && ycell > -11 && ycell < 12 )
           ||( xcell > 13 && xcell < 17 && ycell > -8  && ycell < 8 ) )
        plane_ = 0;
      else if ( xcell > 17 && xcell < 25 && ycell > -5 && ycell < 5 )
        plane_ = 1;
      else
        plane_ = 2;
    } else
      plane_ = 0;
  }
  if( !mypulse_ ) {
    mypulse_ = PulseManager::Instance()->GetPulse(id, plane_);
    {
      string value;
      if( CsOpt::Instance()->getOpt( "CALDEBUG", "FIT_OFILE", value ) )
        {
          asprintf( &fit_ofile_, "%s", value.c_str());
        } else {
          fit_ofile_ = NULL;
        }
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_BASE", value ) )
        {
          fit_base_ = atof(value.c_str());
        } else {
          fit_base_ = 0.;
        }
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_RANGE_TO_FRONT", value ) )
        {
          fit_range_to_front_ = -atof(value.c_str());
        }
      else
        fit_range_to_front_ = -10.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_RANGE_TO_BACK", value ) )
        {
          fit_range_to_back_ = atof(value.c_str());
        }
      else
        fit_range_to_back_ = 10.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_CUT_SLOPE", value ) )
        {
          fit_cut_slope_ = atof(value.c_str());
        }
      else
        fit_cut_slope_ = 1000.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_CUT_OFFSET", value ) )
        {
          fit_cut_offset_ = atof(value.c_str());
        }
      else
        fit_cut_offset_ = 1000.;
       if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_CUT_DEBUG", value ) )
        {
          fit_cut_debug_ = atof(value.c_str());
        }
      else
        fit_cut_debug_ = 0.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_SAMPLE_CERR", value ) )
        {
          fit_sample_cerr_ = atof(value.c_str());
        }
      else
        fit_sample_cerr_ = 7.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "FIT_SAMPLE_DERR", value ) )
        {
          fit_sample_derr_ = atof(value.c_str());
        }
      else
        fit_sample_derr_ = 0.;
      if( CsOpt::Instance()->getOpt( parent_->GetName(), "AC_CUT_MIN", value ) )
        {
          ac_cut_min_ = atof(value.c_str());
        }
      else
        ac_cut_min_ = 0.;

    }


  }


  if( !mypulse_  )
    throw Reco::Exception("%s in %s at %i : Uncomplete configuration : %s: No SPLINE for plane %i (%s)\n",
                          __FILE__, __FUNCTION__, __LINE__, parent_->GetName().c_str(), plane_, DAQDetName.c_str());

  TF1 mypulse_func_("f",mypulse_,0.,32.,3,"MyPulse");


  // Only continue with non empty sample
  if( !sample.size() )
    return result_;

  unsigned short Sample[sample.size()];

  double base[2] = {fit_base_,fit_base_};
  if ( fit_base_ == 0. ) {
    unsigned int use_base_diff = 6; // has to be even
    {
      unsigned int nbase = 0;
      for ( unsigned int i = 0; i < sample.size() && i < use_base_diff; i+=2 ) {
        base[0] += sample[i];
        base[1] += sample[i+1];
        nbase++;
      }
      base[0] /= nbase;
      base[1] /= nbase;
    }
  }

  // Generate histogram
  TGraphAsymmErrors g( sample.size() );

  // Get starting parameters for fit and fill histogram
  list<int> localMax;
  {
    int max = 0;
    for ( unsigned int i = 0; i < sample.size(); i++ ) {
      Sample[i] = sample[i];
      if ( max < 0 && sample[i] >= sample[i] - 1 )
        max = i;
      else if ( sample[max] < sample[i] )
        max = i;
      else if ( ( (sample[max]-base[max&1]) * 0.8 > (sample[i]-base[i&1]) ) ||
                ( i + 1 ==  sample.size() ) ) {
        localMax.push_back(max);
        max = -1;
      }

      g.SetPoint(i, i, sample[i] - base[i&1] );
      double ampl_err[2] = {
        (fit_sample_cerr_+fit_sample_derr_*sqrt(sample[i] - base[i&1])),
        (fit_sample_cerr_+fit_sample_derr_*sqrt(sample[i] - base[i&1])),
      };
      if( sample[i] == overflow_amplitude_ )
        ampl_err[1] = 4. * overflow_amplitude_;
      g.SetPointError(i, 0.3 / 12.6,  0.3 / 12.6, ampl_err[0], ampl_err[1]);

    }
    if( max > 0 && (unsigned int)max < sample.size() -1 )
      localMax.push_back(max);
  }


  if ( !localMax.size() )
    return result_;

  unsigned int Max( 1 << 31);
  for ( list<int>::iterator it = localMax.begin();
        it !=  localMax.end(); ++it ) {

    int use(0), max(*it);
    if ( ac_cut_min_ == -1 ||
         sample[max] - base[max&1] <= ac_cut_min_ ) {
      use = 1;
    } else {
      int pre(0), post(0);
      float sighalf =  (sample[max]-base[max&1])* 0.5;
      while ( max - pre >= 0 &&
              sample[max - pre]-base[(max - pre)&1]  > sighalf )
        pre++;
      while ( (unsigned int)max + post < sample.size() &&
              sample[max + post]-base[(max + post)&1]  > sighalf ) post++;
      if (  pre + post - 1 > 5  ) {
        use = 1;
      }
    }
    if ( use ) {


      if ( Max == (unsigned int)(1 << 31) )
        Max = *it;

      else if ( sample[Max] < sample[*it] )
        Max = *it;

    }


  }

  if ( Max == (unsigned int)(1 << 31) )
    return result_;


  // Prepare function
  mypulse_func_.SetParameters(Max , sample[Max] - base[Max&1], 0. );

  // Execute fit
  g.Fit(&mypulse_func_, "Q", "",
        max(0. , Max - fit_range_to_front_) ,
        min( (double)sample.size(), Max + fit_range_to_back_));

  TF1 *fpulse = g.GetFunction("f");
  double time = fpulse->GetParameter(0);
  double timeErr = fpulse->GetParError(0);
  double ampl = fpulse->GetParameter(1);
  double amplErr = fpulse->GetParError(1);
  double fbase =  fpulse->GetParameter(2);
  double fbaseErr = fpulse->GetParError(2);
  double chi2 = fpulse->GetChisquare();
  double ndf = fpulse->GetNDF();
  double prob = fpulse->GetProb();
  unsigned int run = CsEvent::Instance()->getRunNumber();
  unsigned int eventNb = CsEvent::Instance()->getEventNumberInRun();


  if (time < 9 || time > 22 ) return result_;

  if ( fit_cut_debug_ && fit_cut_debug_ <= fabs( ampl / (sample[Max] - (base[Max&1] + fbase) ) - 1 ) ) return result_;

/// \todo  TODO: Add checks to assert success of fit

  time *= sadc_clock_;
  timeErr *= sadc_clock_;

  time -= 27.8;


  //Analysis output to file

  if ( fit_ofile_ ) {
    static FILE *file = NULL;
    if (!file) {
      char* fname;
      asprintf( &fname, "%s.%f.%f",
                fit_ofile_, -fit_range_to_front_, fit_range_to_back_);
      cout << "Opening file '" << fname << "'\t";
      cout.flush();
      file = fopen( fname , "w");
      free( fname );
      cout << "Done"<< endl;
      assert( file );
    }

    double diff = ampl / (sample[Max] - (base[Max%2] + fbase) ) - 1;

    unsigned int size = sample.size();

    fwrite(&run, sizeof(unsigned int), 1, file);
    fwrite(&eventNb, sizeof(unsigned int), 1, file);
    fwrite(&id, sizeof(int), 1, file);
    fwrite(&plane_, sizeof(int), 1, file);
    fwrite(&x, sizeof(int), 1, file);
    fwrite(&y, sizeof(int), 1, file);
    fwrite(&diff, sizeof(double), 1 , file);
    fwrite(&ampl, sizeof(double), 1 , file);
    fwrite(&chi2, sizeof(double), 1 , file);
    fwrite(&ndf, sizeof(double), 1 , file);
    fwrite(&prob, sizeof(double), 1 , file);
    fwrite(&time, sizeof(double),1 , file);
    fwrite(base, sizeof(double), 2, file);
    fwrite(&fbase, sizeof(double), 1, file);
    fwrite(&size, sizeof(unsigned int), 1, file);
    fwrite(Sample, sizeof(unsigned short), size, file);

  }

  if( ampl && fit_cut_slope_ != 1000. &&  fit_cut_offset_ != 1000. ) {
    if( (chi2 / ndf / ampl * 200.) > ( fit_cut_slope_ * ampl + fit_cut_offset_ ) )
      return result_;
  }

  // Fill result structure


  result_.ampl = ampl;
  result_.time = time;
  result_.base = fbase+base[0];
  result_.base_diff = base[0] - base[1];
  result_.dampl = amplErr;
  result_.dtime = timeErr;
  result_.dbase = fbaseErr;
  result_.chi2 = chi2;
  result_.ndf = ndf;
  result_.ready = 1;

  return result_;

#else
  throw Reco::Exception("CsDigitizerSADC::FitPulse Requires ROOT-Version >= 5.16/00! \n");
  return result_;
#endif
};
