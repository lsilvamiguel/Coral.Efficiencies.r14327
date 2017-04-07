/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/CalorimeterOptions.cc,v $
   $Date: 2011/02/04 17:59:18 $
   $Revision: 1.17 $
   -------------------------------------------------------------------------
*/

// --- Internal files ----
#include "Calorimeter.h"

#include <fstream>
#include <sstream>

using namespace std;

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

Calorimeter::Options::Options() {
    options_ok                            = false;

    debug                                 = -1;
    calo_type                             = ElectroMagnetic;
    particle_default                      = CalorimeterParticle::GAMMA;
    mixed_blocks                          = false;
    is_real_data                          = true;

    readout_sparsified_                   = true;
    delta_spars_mode_                     = 3.;

    calibration_flag                      = false;
    calib_histo_units_in_gev              = false;
    calibration_energy                    = 10.;
    calibration_energy_min                = 0;
    calibration_energy_max                = 0;
    calibration_min                       = 0;
    calibration_max                       = 0;
    calib_histo_size                      = 1000;
    calib_histo_range                     = 5.;
    new_calib_stat_min                    = 10.;
    recalc_from_primary_calib             = false;
    recalc_from_primary_time_calib        = false;
    led_histo_size                        = 500;
    hist_led_max                          = 5000.;
    calibration_uncertainty               = 0;
    level_particle_calo_association_histo = 1;
    ngam_max4histo                        = 50;

    default_calibration                   = 1.0;
    default_time0_calibration             = 0.;
    enable_time_calibration               = false;
    calibration_time_gate                 = 5.0;
    scale_calibration                     = 1.0;
    tolerate2keepold_calib                = 0.0000001;
    scale_calibration_by_led              = false;
    use_ProbEnT                           = false;
    use_cell_calib_correction             = false;
    use_led_ref_calib_correction          = false;
    use_led_ref_calib_correction_in_spills= false;
    scale_cells_threshols                 = 1.0;

    vertex_position[0]                    = 0.0;
    vertex_position[1]                    = 0.0;
    vertex_position[2]                    = 0.0;

    position_correction[0]                = 0.0;
    position_correction[1]                = 0.0;
    position_correction[2]                = 0.0;

    position_offset[0]                    = 0.0;
    position_offset[1]                    = 0.0;
    position_offset[2]                    = 0.0;

    store_noise_info                      = false;
    store_rndm_noise_info                 = false;
    store_gamma_noise_info                = false;
    store_noise_histo                     = false;
    update_ecut_using_prob_cell_noise     = false;
    prob_cell_noise                       = 0.01;
    prob_gamma_noise                      = 0.001;
    stat_min_cell_noise                   = 1.0;
    stat_min_gamma_noise                  = 1.0;

    scale_old_calib_info                  = false;

    base_energy_noise                     = 2.0;

    store_leds_all                        = false;
    store_in_spill_leds_all             = false;
    led_statistic_inspect                 = 0.;
    reset_after_led_inspect               = true;

    ecut_extend_list                      = 10.0;
    readout_term                          = 0;     // reseting readout term; loosy solution

    use_preshower_correction              = false;
    preshower_length                      = 0.;

    store_back_data                       = false;

    fortran_bad_cell_check                = 0;
    fortran_fast_reconstruction           = 0;
    fortran_hbook                         = 0;
    fortran_variable_thresholds           = 0;
    ecal2_time_0                          = 0.;
    ecal2_time_width                      = 5.;

    flag_fast_online                      = false;
    print_general_info                    = false;
    print_unmaped_data                    = true;
    print_bad_data                        = true;
    do_reconstruction                     = true;
    calib_event_reconstruction            = false;

//////////////////////////////////////////////////
    ecal0_nonlinear_calib_correction      = false;
//////////////////////////////////////////////////

    correct_for_digitization              = true;
    use_time_in_reconstruction            = true;
    combined_reconstruction               = false;
    kolosov_reconstruction                = true;
    fortran_reconstruction                = false;
    repeat_reconstruction                 = -1;
    reco_cluster_search                   = false;
    reco_cluster_search_only              = false;
    reco_use_hint_particles_info          = false;
    add_time_to_reco_particles            = false;
    recover_bad_cells                     = false;
    use_simple_fit                        = true;
    fitting_type                          = Simple;
    add_search                            = false;
    hicut_add_search_reco2                = 50.;
    hicut_add_search                      = 100000.;
    ecut_add_search                       = 100000.;
    fill_histos                           = true;
    fill_only_reco_histos_in_reco_test    = true;
    histos_level                          = HISTLVL_NORMAL;
    make_profiles                         = false;
    fill_raw_histos                       = false;
    test_histo_raw_info                   = true;
    test_histo_raw_info_level             = 0;
    fill_internal_histos                  = false;
    fill_fit_info_histo                   = false;
    fill_internal_correlations            = false;
    fill_external_correlations            = false;
    fill_calib_spill_correlations_histo   = false;
    hist_energy_max                       = 200;
    hist_time_min                         = -50.;
    hist_time_max                         = 50.;
    hist_time_nbins_1d                    = 200;
    hist_time_nbins_2d                    = 50;
    hist_energy_nbins_1d                  = 1000;
    hist_energy_nbins_2d                  = 50;

    need_cells_XY_histo                   = false;
    particle_energy_threshold             = 0.3;
    cell_energy_threshold                 = 0.1;
//  energy_leak_for_SimpleFit             = 0.06; // Former default value for several years
    energy_leak_for_SimpleFit             = 0.00;  // Default was changed Thu Aug  6 15:01:02 CEST 2009 this require some changes in coral options
    distmax_for_cells_search              = 0;
    tolerance_for_nearby_cells            = 0.1;
    check_front_surface_position          = true;
    init_xy_structure                     = true;
    cluster_search_cell_energy_threshold  = 0.02;
    cluster_search_cell_ampl_deviation    = 2;
    show_all_cells                        = true;

    // options for ReconstructionCombined
    correct_longitudinal_leakage          = true;
    min_shower_distance                   = 1.;
    cell_adc_threshold                    = 0.;
    use_combined_calibration              = false;
    combined_calibration                  = "";

    monitor_histo_show_fit                = false;
    print_bad_cells_info                  = false;
    monitor_db_calib_level                = 2;
    monitor_db_calib_min_value            = 0.001;
    monitor_db_calib_max_value            = 0.05;
    monitor_db_led_level                  = 2;
    monitor_db_led_min_value              = 10.;
    monitor_db_led_max_value              = 3000.;
    monitor_db_porog_level                = 2;
    monitor_db_porog_min_value            = 0.;
    monitor_db_porog_max_value            = 2.0;

    inspect_led_Ncells_no_signal          = 0;
    inspect_led_Ncells_small_signal       = 0;
    inspect_led_Fcells_good_signal        = 0.9999999;

    fit_led_histo_at_end_of_job           = false;
    debug_reconstruction                  = false;
    debug_reconstruction_test             = false;
    debug_reconstruction_combined         = false;
    ecut4mgg_histo                        = 0.;

    sort_particles                        = false;
    blackbox_chk_file                     = NULL;


    for( int i=0; i<50; i++ )
    {
        double_options_int[i] = false;
        bool_options_int[i]   = false;
    }

    // options for Monte Carlo construction
    mc_default_calibration                = 1.0;
    mc_make_real_digitization             = false;

    mc_make_fiadc_digitization            = false;
    mc_fiadc_sparce_delta                 = 3;

    mc_data_extraction_threshold          = 0.01;

    add_shower_fluctuations_for_FMC       = true;

    mc_smear_response                     = true;
    mc_smear_constant                     = true;
    mc_smear_stochastic                   = true;
    mc_smear_readout                      = true;
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Options::SetRecoOptions( const string &opt_ )
{
  if( options_ok  )
    throw Exception("SetRecoOptions is closed already!");

  if( opt_.length()==0 )
    return;

  string opt(opt_);  // create copy which is non-const
  istringstream oo_orig(opt.c_str());

  for( size_t i=0; i<opt.length(); i++ )
    opt[i] = toupper(opt[i]);

  istringstream oo(opt.c_str());

  string o;
  while( oo>>o )
  {
    string o_orig;

    oo_orig >> o_orig;
    string value_orig="";

    size_t n=o.find('=');
    string name="";
    string value="";
    string name_index="";
    string sam_index="";
    size_t nind=0;
    int index_opt =-1;

    if( n==string::npos || n==0 || n==o.length() )
      goto BAD_OPTION;
    name=string(o,0,n);
    value=string(o,n+1);

    value_orig=string(o_orig,n+1);

    nind=name.find('[');
    if( !(nind==string::npos || nind == 0 || nind == name.length()) )
    {
      name_index=string(name,0,nind);
      sam_index=string(name,nind+1,name.length()-nind-2);
      index_opt = atoi( sam_index.c_str() );
    }

    if( name=="DEBUG" )
    {
      debug = atoi(value.c_str());
    }
    else
    if( name=="READOUT_SPARSIFIED" )
    {
      if( value=="YES" )
        readout_sparsified_=true;
      else
      if( value=="NO" )
        readout_sparsified_=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="DELTA_SPARS_MODE" )
    {
      delta_spars_mode_ = atof(value.c_str());
    }
    else
    if( name=="CALIB" )
    {
      if( value=="YES" )
        calibration_flag=true;
      else
      if( value=="NO" )
        calibration_flag=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="SCALE_CALIBRATION" )
    {
      scale_calibration = atof(value.c_str());
    }
    else
    if( name=="TOLERATE2KEEPOLD_CALIB" )
    {
      tolerate2keepold_calib = atof(value.c_str());
    }
    else
    if( name=="SCALE_CELLS_THRESHOLDS" )
    {
      scale_cells_threshols = atof(value.c_str());
    }
    else
    if( name=="CORRECT_POSITION_X" )
    {
      position_correction[0] = atof(value.c_str());
    }
    else
    if( name=="CORRECT_POSITION_Y" )
    {
      position_correction[1] = atof(value.c_str());
    }
    else
    if( name=="CORRECT_POSITION_Z" )
    {
      position_correction[2] = atof(value.c_str());
    }
    else
    if( name=="OFFSET_POSITION_X" )
    {
      position_offset[0] = atof(value.c_str());
    }
    else
    if( name=="OFFSET_POSITION_Y" )
    {
      position_offset[1] = atof(value.c_str());
    }
    else
    if( name=="OFFSET_POSITION_Z" )
    {
      position_offset[2] = atof(value.c_str());
    }
    else
    if( name=="VERTEX_POSITION_X" )
    {
      vertex_position[0] = atof(value.c_str());
    }
    else
    if( name=="VERTEX_POSITION_Y" )
    {
      vertex_position[1] = atof(value.c_str());
    }
    else
    if( name=="VERTEX_POSITION_Z" )
    {
      vertex_position[2] = atof(value.c_str());
    }
    else
    if( name=="USE_PROBENT" )
    {
      if( value=="YES" )
        use_ProbEnT=true;
      else
      if( value=="NO" )
        use_ProbEnT=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="SCALE_CALIBRATION_BY_LED" )
    {
      if( value=="YES" )
        scale_calibration_by_led=true;
      else
      if( value=="NO" )
        scale_calibration_by_led=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="USE_LED_REF_CALIB_CORRECTION" )
    {
      if( value=="YES" )
        use_led_ref_calib_correction=true;
      else
      if( value=="NO" )
        use_led_ref_calib_correction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="USE_LED_REF_CALIB_CORRECTION_IN_SPILLS" )
    {
      if( value=="YES" )
        use_led_ref_calib_correction_in_spills=true;
      else
      if( value=="NO" )
        use_led_ref_calib_correction_in_spills=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="USE_EDEPCORR" || name=="EDEPCORR" )
    {
      if( value=="YES" )
        edep_corr = "_EdepCorr";
      else
      if( value=="NO" )
        edep_corr = "";
      else
	edep_corr = value_orig;
    }
    else
    if( name=="USE_TISDEPCORR" || name=="TISDEPCORR" )
    {
      if( value=="YES" )
        tisdep_corr = "_TiSdepCorr";
      else
      if( value=="NO" )
        tisdep_corr = "";
      else
        tisdep_corr = value_orig;
    }
    else
    if( name=="USE_CELL_CALIB_CORRECTION" )
    {
      if( value=="YES" )
        use_cell_calib_correction=true;
      else
      if( value=="NO" )
        use_cell_calib_correction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="CALIB_ENERGY" )
    {
      calibration_energy = atof(value.c_str());
    }
    else
    if( name=="CALIB_ENERGY_MIN" )
    {
      calibration_energy_min = atof(value.c_str());
    }
    else
    if( name=="CALIB_ENERGY_MAX" )
    {
      calibration_energy_max = atof(value.c_str());
    }
    else
    if( name=="CALIB_UNCERTAINTY" )
    {
      calibration_uncertainty = atof(value.c_str());
    }
    else
    if( name=="CALIB_MIN" )
    {
      calibration_min = atof(value.c_str());
    }
    else
    if( name=="CALIB_MAX" )
    {
      calibration_max = atof(value.c_str());
    }
    else
    if( name=="NEW_CALIB_STAT_MIN" )
    {
      new_calib_stat_min = atof(value.c_str());
    }
    else
    if( name=="DEFAULT_CALIBRATION" )
    {
      default_calibration = atof(value.c_str());
    }
    else
    if( name=="MC_DEFAULT_CALIBRATION" )
    {
      mc_default_calibration = atof(value.c_str());
    }
    else
    if( name=="DEFAULT_TIME0_CALIBRATION" )
    {
      default_time0_calibration = atof(value.c_str());
    }
    else
    if( name=="ENABLE_TIME_CALIB" )
    {
      if( value=="YES" )
        enable_time_calibration=true;
      else
      if( value=="NO" )
        enable_time_calibration=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="RECALC_PRIM_CALIB" )
    {
      if( value=="YES" )
        recalc_from_primary_calib=true;
      else
      if( value=="NO" )
        recalc_from_primary_calib=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="RECALC_PRIM_TIME_CALIB" )
    {
      if( value=="YES" )
        recalc_from_primary_time_calib=true;
      else
      if( value=="NO" )
        recalc_from_primary_time_calib=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="STORE_NOISE_INFO" )
    {
      if( value=="YES" )
        store_noise_info=true;
      else
      if( value=="NO" )
        store_noise_info=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="SCALE_OLD_CALIB_INFO" )
    {
      if( value=="YES" )
        scale_old_calib_info=true;
      else
      if( value=="NO" )
        scale_old_calib_info=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="STORE_RNDM_NOISE_INFO" )
    {
      if( value=="YES" )
        store_rndm_noise_info=true;
      else
      if( value=="NO" )
        store_rndm_noise_info=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="STORE_GAMMA_NOISE_INFO" )
    {
      if( value=="YES" )
        store_gamma_noise_info=true;
      else
      if( value=="NO" )
        store_gamma_noise_info=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="STORE_NOISE_HISTO" )
    {
      if( value=="YES" )
        store_noise_histo=true;
      else
      if( value=="NO" )
        store_noise_histo=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="USE_PROB_NOISE" )
    {
      if( value=="YES" )
        update_ecut_using_prob_cell_noise=true;
      else
      if( value=="NO" )
        update_ecut_using_prob_cell_noise=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="PROB_CELL_NOISE" )
    {
      prob_cell_noise = atof( value.c_str() );
    }
    else
    if( name=="PROB_GAMMA_NOISE" )
    {
      prob_gamma_noise = atof( value.c_str() );
    }
    else
    if( name=="STAT_MIN_CELL_NOISE" )
    {
      stat_min_cell_noise = atof( value.c_str() );
    }
    else
    if( name=="STAT_MIN_GAMMA_NOISE" )
    {
      stat_min_gamma_noise = atof( value.c_str() );
    }
    else
    if( name=="BASE_ENERGY_NOISE" )
    {
      base_energy_noise = atof( value.c_str() );
    }
    else
    if( name=="STAT_LED_INSPECT" )
    {
      led_statistic_inspect = atof( value.c_str() );
    }
    else
    if( name=="RESET_LED_INSPECT" )
    {
      if( value=="YES" )
        reset_after_led_inspect=true;
      else
      if( value=="NO" )
        reset_after_led_inspect=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="STORE_LEDS_ALL" )
    {
      if( value=="YES" )
        store_leds_all=true;
      else
      if( value=="NO" )
        store_leds_all=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="STORE_IN_SPILL_LEDS_ALL" )
    {
      if( value=="YES" )
        store_in_spill_leds_all=true;
      else
      if( value=="NO" )
        store_in_spill_leds_all=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="ADD_SMEARING" )
    {
      if( value=="YES" )
        mc_smear_response=true;
      else
      if( value=="NO" )
        mc_smear_response=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="ADD_CONSTANT_FLUCTUATIONS" )
    {
      if( value=="YES" )
        mc_smear_constant=true;
      else
      if( value=="NO" )
        mc_smear_constant=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="ADD_STOCHASTIC_FLUCTUATIONS" )
    {
      if( value=="YES" )
        mc_smear_stochastic=true;
      else
      if( value=="NO" )
        mc_smear_stochastic=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="ADD_READOUT_FLUCTUATIONS" )
    {
      if( value=="YES" )
        mc_smear_readout=true;
      else
      if( value=="NO" )
        mc_smear_readout=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="ADD_SHOWER_FLUCTUATIONS_FOR_FMC" )
    {
      if( value=="YES" )
        add_shower_fluctuations_for_FMC=true;
      else
      if( value=="NO" )
        add_shower_fluctuations_for_FMC=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="RESET_READOUT_TERM" )
    {
      readout_term = atof(value.c_str());
    }
    else
    if( name=="MC_DATA_EXTRACTION_THRESHOLD" )
    {
      mc_data_extraction_threshold = atof(value.c_str());
    }
    else
    if( name=="ECUT_EXTEND_LIST" )
    {
      ecut_extend_list = atof(value.c_str());
    }
    else
    if( name=="MC_MAKE_REAL_DIGITIZATION" )
    {
      if( value=="YES" )
        mc_make_real_digitization=true;
      else
      if( value=="NO" )
        mc_make_real_digitization=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="MC_MAKE_FIADC_DIGITIZATION" )
    {
      if( value=="YES" )
        mc_make_fiadc_digitization=true;
      else
      if( value=="NO" )
        mc_make_fiadc_digitization=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="MC_FIADC_SPARCE_DELTA" )
    {
      mc_fiadc_sparce_delta = atof(value.c_str());
    }
   else
    if( name=="TIME_GATE_CALIB" )
    {
      calibration_time_gate = atof(value.c_str());
    }
    else
    if( name=="CLUSTER_SEARCH_CELL_ENERGY_THRESHOLD" )
    {
      cluster_search_cell_energy_threshold = atof(value.c_str());
    }
    else
    if( name=="CLUSTER_SEARCH_CELL_AMPL_DEVIATION" )
    {
      cluster_search_cell_ampl_deviation = atof(value.c_str());
    }
    else
    if( name=="FIT" )
    {
      if( value=="SIMPLE" )
      {
        use_simple_fit=true;
	fitting_type = Simple;
      }
      else
      if( value=="NORMAL" )
      {
        use_simple_fit=false;
	fitting_type = Normal;
      }
      else
      if( value=="NOFIT" )
      {
        use_simple_fit=false;
	fitting_type = NoFit;
      }
      else
      {
        goto BAD_OPTION;
      }
    }
    else
    if( name=="ADD_SEARCH" )
    {
      if( value=="YES" )
        add_search=true;
      else
      if( value=="NO" )
        add_search=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="HICUT_ADD_SEARCH" )
    {
      hicut_add_search = atof(value.c_str());
    }
    else
    if( name=="HICUT_ADD_SEARCH_RECO2" )
    {
      hicut_add_search_reco2 = atof(value.c_str());
    }
    else
    if( name=="ECUT_ADD_SEARCH" )
    {
      ecut_add_search = atof(value.c_str());
    }
    else
    if( name=="HIST" )
    {
      if( value=="YES" )
      {
        fill_histos=true;
      }
      else
      if( value=="NO" )
      {
        fill_histos=false;
      }
      else
        goto BAD_OPTION;
    }
    else
    if( name=="FILL_ONLY_RECO_HIST_IN_RECO_TEST" )
    {
      if( value=="YES" )
      {
        fill_only_reco_histos_in_reco_test=true;
      }
      else
      if( value=="NO" )
      {
        fill_only_reco_histos_in_reco_test=false;
      }
      else
        goto BAD_OPTION;
    }
    else
    if( name=="HIST_LEVEL")
    {
      if( value=="NORMAL" )
      {
        histos_level=Options::HISTLVL_NORMAL;
      }
      else
      if( value=="VERBOSE" )
      {
        histos_level=Options::HISTLVL_VERBOSE;
      }
      else
      if( value=="DEBUG" )
      {
        histos_level=Options::HISTLVL_DEBUG;
      }
      else
        goto BAD_OPTION;
    }
    else
    if( name=="MAKE_PROFILES" )
    {
      if( value=="YES" )
        make_profiles=true;
      else
      if( value=="NO" )
        make_profiles=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="RAW_HIST" )
    {
      if( value=="YES" )
        fill_raw_histos=true;
      else
      if( value=="NO" )
        fill_raw_histos=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="TEST_HISTO_RAW_INFO" )
    {
      if( value=="YES" )
        test_histo_raw_info=true;
      else
      if( value=="NO" )
        test_histo_raw_info=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="TEST_HISTO_RAW_INFO_LEVEL" )
    {
      test_histo_raw_info_level=atoi(value.c_str());;
    }
    else
    if( name=="FILL_INTERNAL_HIST" )
    {
      if( value=="YES" )
        fill_internal_histos=true;
      else
      if( value=="NO" )
        fill_internal_histos=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="FILL_FIT_INFO_HIST" )
    {
      if( value=="YES" )
        fill_fit_info_histo=true;
      else
      if( value=="NO" )
        fill_fit_info_histo=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="FILL_INTERNAL_CORRELATIONS" )
    {
      if( value=="YES" )
        fill_internal_correlations=true;
      else
      if( value=="NO" )
        fill_internal_correlations=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="FILL_EXTERNAL_CORRELATIONS" )
    {
      if( value=="YES" )
        fill_external_correlations=true;
      else
      if( value=="NO" )
        fill_external_correlations=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="FILL_CALIB_SPILL_CORRELATIONS_HISTO" )
    {
      if( value=="YES" )
        fill_calib_spill_correlations_histo=true;
      else
      if( value=="NO" )
        fill_calib_spill_correlations_histo=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="HIST_CELLS_DXDY" )
    {
      if( value=="YES" )
        need_cells_XY_histo=true;
      else
      if( value=="NO" )
        need_cells_XY_histo=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="FLAG_FAST_ONLINE" )
    {
      if( value=="YES" )
        flag_fast_online=true;
      else
      if( value=="NO" )
        flag_fast_online=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="HIST_ENERGY_MAX" )
    {
      hist_energy_max = atof(value.c_str());
    }
    else
    if( name=="HIST_TIME_MIN" )
    {
      hist_time_min = atof(value.c_str());
    }
    else
    if( name=="HIST_TIME_MAX" )
    {
      hist_time_max = atof(value.c_str());
    }
    else
    if( name=="LED_HISTO_SIZE" )
    {
       led_histo_size=atoi(value.c_str());
    }
    else
    if( name=="NGAM_MAX_FOR_HISTO" )
    {
       ngam_max4histo=atoi(value.c_str());
    }
     else
    if( name=="LEVEL_ASSOCIATION_HISTO" )
    {
       level_particle_calo_association_histo=atoi(value.c_str());
    }
   else
    if( name=="HIST_LED_MAX" )
    {
      hist_led_max = atof(value.c_str());
    }
    else
    if( name=="PRINT_GENERAL_INFO" )
    {
      if( value=="YES" )
        print_general_info=true;
      else
      if( value=="NO" )
        print_general_info=false;
      else
        goto BAD_OPTION;
    }
     else
    if( name=="PRINT_UNMAPPED_DATA" )
    {
      if( value=="YES" )
        print_unmaped_data=true;
      else
      if( value=="NO" )
        print_unmaped_data=false;
      else
        goto BAD_OPTION;
    }
   else
    if( name=="PRINT_BAD_DATA" )
    {
      if( value=="YES" )
        print_bad_data=true;
      else
      if( value=="NO" )
        print_bad_data=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="PARTICLE_DEFAULT" )
    {
      particle_default = CalorimeterParticle::ParticleID(atoi(value.c_str()));
    }
    else
    if( name=="PARTICLE_ENERGY_THRESHOLD" )
    {
      particle_energy_threshold = atof(value.c_str());
    }
    else
    if( name=="ENERGY_LEAK_FOR_SIMPLE_FIT" )
    {
      energy_leak_for_SimpleFit = atof(value.c_str());
    }
    else
    if( name=="ENERGY_CUT_FOR_MGG_HISTO" )
    {
      ecut4mgg_histo = atof(value.c_str());
    }
    else
    if( name=="DISTMAX_FOR_CELLS_SEARCH" )
    {
      distmax_for_cells_search = atof(value.c_str());
    }
    else
    if( name=="MIN_SHOWER_DISTANCE" )
    {
      min_shower_distance = atof(value.c_str());
    }
    else
    if( name=="CORRECT_LONGITUDINAL_LEAKAGE")
    {
      if( value=="TRUE")
        correct_longitudinal_leakage=true;
      else
      if( value=="FALSE")
        correct_longitudinal_leakage=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="CELL_ADC_THRESHOLD" )
    {
      cell_adc_threshold = atof(value.c_str());
    }
    else
    if( name=="USE_RECO_COMBINED_CALIB" )
    {
      if( value=="TRUE")
        use_combined_calibration=true;
      else
      if( value=="FALSE")
        use_combined_calibration=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="TOLERANCE_FOR_NEARBY_CELLS" )
    {
      tolerance_for_nearby_cells = atof(value.c_str());
    }
    else
    if( name=="CELL_ENERGY_THRESHOLD" )
    {
      cell_energy_threshold = atof(value.c_str());
    }
    else
    if( name=="USE_PRESHOWER_CORRECTION" )
    {
      if( value=="YES" )
        use_preshower_correction=true;
      else
      if( value=="NO" )
        use_preshower_correction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="PRE_SHOWER_LENGTH" )
    {
      preshower_length = atof(value.c_str());
    }
    else
    if( name=="DO_RECONSTRUCTION" )
    {
      if( value=="YES" )
        do_reconstruction=true;
      else
      if( value=="NO" )
        do_reconstruction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="CALIB_EVENT_RECONSTRUCTION" )
    {
      if( value=="YES" )
        calib_event_reconstruction=true;
      else
      if( value=="NO" )
        calib_event_reconstruction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="ECAL0_NONLINEAR_CALIB_CORRECTION" )
    {
      if( value=="YES" )
        ecal0_nonlinear_calib_correction=true;
      else
      if( value=="NO" )
        ecal0_nonlinear_calib_correction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="CORRECT_FOR_DIGITIZATION" )
    {
      if( value=="YES" )
        correct_for_digitization=true;
      else
      if( value=="NO" )
        correct_for_digitization=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="REPEAT_RECONSTRUCTION" )
    {
       repeat_reconstruction=atoi(value.c_str());
    }
    else
    if( name=="RECONSTRUCTION" )
    {
      combined_reconstruction=false;
      fortran_reconstruction=false;
      kolosov_reconstruction=false;
      if( value=="COMBINED" )
        combined_reconstruction=true;
      else if( value=="FORTRAN" )
        fortran_reconstruction=true;
      else if( value=="KOLOSOV" )
        kolosov_reconstruction=true;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="USE_TIME_IN_RECONSTRUCTION" )
    {
      if( value=="YES" )
        use_time_in_reconstruction=true;
      else
      if( value=="NO" )
        use_time_in_reconstruction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="STORE_BACK_DATA" )
    {
      if( value=="YES" )
        store_back_data=true;
      else
      if( value=="NO" )
        store_back_data=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="FAST_RECONSTRUCTION" )
    {
      if      ( value=="YES" ) fortran_fast_reconstruction = 1;
      else if ( value=="NO"  ) fortran_fast_reconstruction = 0;
      else goto BAD_OPTION;
    }
    else
    if( name=="HBOOK" )
    {
      if      ( value=="YES" ) fortran_hbook = 1;
      else if ( value=="NO"  ) fortran_hbook = 0;
      else goto BAD_OPTION;
    }
     else
    if( name=="VARIABLE_THRESHOLDS" )
    {
      if      ( value=="NO"   ) fortran_variable_thresholds = 0;
      else if ( value=="YES"  ) fortran_variable_thresholds = 1;
      else if ( value=="READ" ) fortran_variable_thresholds = 2;
      else goto BAD_OPTION;
    }
    else
    if( name=="BAD_CELL_CHECK" )
    {
      if      ( value=="YES" ) fortran_bad_cell_check = 1;
      else if ( value=="NO"  ) fortran_bad_cell_check = 0;
      else goto BAD_OPTION;
    }
    else
    if( name=="ECAL2_TIME_0" )
    {
      ecal2_time_0= atof( value.c_str() );
    }
    else
    if( name=="ECAL2_TIME_WIDTH" )
    {
      ecal2_time_width= atof( value.c_str() );
    }
    else
    if( name=="CALIB_HISTO_SIZE" )
    {
       calib_histo_size=atoi(value.c_str());
    }
    else
    if( name=="TIME1D_HISTO_SIZE" )
    {
       hist_time_nbins_1d=atoi(value.c_str());
    }
    else
    if( name=="TIME2D_HISTO_SIZE" )
    {
       hist_time_nbins_2d=atoi(value.c_str());
    }
    else
    if( name=="ENERGY1D_HISTO_SIZE" )
    {
       hist_energy_nbins_1d=atoi(value.c_str());
    }
    else
    if( name=="ENERGY2D_HISTO_SIZE" )
    {
       hist_energy_nbins_2d=atoi(value.c_str());
    }
    else
    if( name=="CALIB_HISTO_RANGE" )
    {
      calib_histo_range = atof(value.c_str());
    }
    else
    if( name=="CLUSTER_SEARCH" )
    {
      if( value=="YES" )
        reco_cluster_search=true;
      else
      if( value=="NO" )
        reco_cluster_search=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="CLUSTER_SEARCH_ONLY" )
    {
      if( value=="YES" )
        reco_cluster_search_only=true;
      else
      if( value=="NO" )
        reco_cluster_search_only=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="RECOVER_BAD_CELLS" )
    {
      if( value=="YES" )
        recover_bad_cells=true;
      else
      if( value=="NO" )
        recover_bad_cells=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="USE_HINT_PARTICLES" )
    {
      if( value=="YES" )
        reco_use_hint_particles_info=true;
      else
      if( value=="NO" )
        reco_use_hint_particles_info=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="ADD_TIME_TO_RECO_PARTICLES" )
    {
      if( value=="YES" )
        add_time_to_reco_particles=true;
      else
      if( value=="NO" )
        add_time_to_reco_particles=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="CHECK_FRONT_SURFACE_POSITION" )
    {
      if( value=="YES" )
        check_front_surface_position=true;
      else
      if( value=="NO" )
        check_front_surface_position=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="INIT_XY_STRUCTURE" )
    {
      if( value=="YES" )
        init_xy_structure=true;
      else
      if( value=="NO" )
        init_xy_structure=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="GUI_SHOW_ALL_CELLS" )
    {
      if( value=="YES" )
        show_all_cells=true;
      else
      if( value=="NO" )
        show_all_cells=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="PRINT_BAD_CELLS_INFO" )
    {
      if( value=="YES" )
        print_bad_cells_info=true;
      else
      if( value=="NO" )
        print_bad_cells_info=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="MONITOR_HISTO_SHOW_FIT" )
    {
      if( value=="YES" )
        monitor_histo_show_fit=true;
      else
      if( value=="NO" )
        monitor_histo_show_fit=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="MONITOR_DB_CALIB_LEVEL" )
    {
      monitor_db_calib_level = (int)atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_CALIB_MIN" )
    {
      monitor_db_calib_min_value = atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_CALIB_MAX" )
    {
      monitor_db_calib_max_value = atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_LED_LEVEL" )
    {
      monitor_db_led_level = (int)atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_LED_MIN" )
    {
      monitor_db_led_min_value = atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_LED_MAX" )
    {
      monitor_db_led_max_value = atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_POROG_LEVEL" )
    {
      monitor_db_porog_level = (int)atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_POROG_MIN" )
    {
      monitor_db_porog_min_value = atof(value.c_str());
    }
    else
    if( name=="MONITOR_DB_POROG_MAX" )
    {
      monitor_db_porog_max_value = atof(value.c_str());
    }
    else
    if( name=="INSPECT_LED_NCELLS_NO_SIGNAL" )
    {
      inspect_led_Ncells_no_signal = (int)atof(value.c_str());
    }
    else
    if( name=="INSPECT_LED_NCELLS_SMALL_SIGNAL" )
    {
      inspect_led_Ncells_small_signal = (int)atof(value.c_str());
    }
    else
    if( name=="INSPECT_LED_FCELLS_GOOD_SIGNAL" )
    {
      inspect_led_Fcells_good_signal = atof(value.c_str());
    }
    else
    if( name=="FIT_LED_HISTO_AT_END_OF_JOB" )
    {
      if( value=="YES" )
        fit_led_histo_at_end_of_job=true;
      else
      if( value=="NO" )
        fit_led_histo_at_end_of_job=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="DEBUG_RECONSTRUCTION" )
    {
      if( value=="YES" )
        debug_reconstruction=true;
      else
      if( value=="NO" )
        debug_reconstruction=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="DEBUG_RECONSTRUCTION_TEST" )
    {
      if( value=="YES" )
        debug_reconstruction_test=true;
      else
      if( value=="NO" )
        debug_reconstruction_test=false;
      else
        goto BAD_OPTION;
    }
    else
    if( name=="DEBUG_RECONSTRUCTION_COMBINED" )
    {
      if( value=="YES" )
        debug_reconstruction_combined=true;
      else
      if( value=="NO" )
        debug_reconstruction_combined=false;
      else
        goto BAD_OPTION;
    }
    else
    if( index_opt >=0 && name_index=="SET_MISC_DOUBLE_OPTION" )
    {
      double dv = atof(value.c_str());
      SetMiscDoubleOption( index_opt, dv );
    }
    else
    if( index_opt >=0 && name_index=="SET_MISC_BOOL_OPTION" )
    {
      if( value=="YES" )
        SetMiscBoolOption( index_opt, true );
      else
      if( value=="NO" )
        SetMiscBoolOption( index_opt, false );
      else
        goto BAD_OPTION;
    }
    else
    if( index_opt >=0 && name_index=="SET_MISC_TIME_OPTION" )
    {
      tm t;
      char c;
      sscanf(value.c_str(),"%d%c%d%c%d%c%d%c%d%c%d",
          &t.tm_year,   &c,
          &t.tm_mon,    &c,
          &t.tm_mday,   &c,
          &t.tm_hour,   &c,
          &t.tm_min,    &c,
          &t.tm_sec        );
      t.tm_year -= 1900;
      t.tm_mon  -= 1;
      time_t topt = mktime( &t);
      SetMiscTimeOption( index_opt, topt );
    }
    else
    if( name=="SET_OUTPUT_INFO" )
    {
      for( int i=0; i< (int)output_info_names.size(); i++ )
      {
	if( output_info_names[i] == value_orig)
	{
          cerr <<" OUTPUT_INFO name already exist " << i << " ! " << value_orig << " We are more strict now ! Please Clean-up your options " << endl;
	  goto BAD_OPTION;
	}
      }
      output_info_names.push_back( value_orig );
    }
    else
    if( name=="SET_INPUT_INFO" )
    {
      for( int i=0; i< (int)input_info_names.size(); i++ )
      {
	if( input_info_names[i] == value_orig)
	{
          cerr <<" INPUT_INFO name already exist " << i << " ! " << value_orig << " We are more strict now ! Please Clean-up your options " << endl;
	  goto BAD_OPTION;
	}
      }
      input_info_names.push_back( value_orig );
    }
    else
      if( name=="SORT_PARTICLES" )
        if( value=="YES" )
          sort_particles=true;
        else if( value=="NO" )
          sort_particles=false;
        else
          goto BAD_OPTION;
    else
      if( name=="BLACKBOX_CHK_FILE" )
        if( value!="" )
          blackbox_chk_file = new ofstream( value.c_str() );
        else
          blackbox_chk_file = NULL;
    else
    {
      throw Exception("Calorimeter::Options::SetRecoOptions(): Unknown option %s", o.c_str());
    }

  }

  return;

 BAD_OPTION:
  throw Exception("Calorimeter::Options::SetRecoOptions(): Bad option %s", o.c_str());
}

////////////////////////////////////////////////////////////////////////////////

double Calorimeter::Options::GetMiscDoubleOption(int iopt) const
{
  assert( 0 <= iopt && iopt < 50 );

  if( !double_options_int[iopt] )
    throw Exception( "Calorimeter::Options::GetMiscDoubleOption(%i):  This option has "
		     "never been initialized!", iopt );

  return misc_double_options[iopt];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Options::SetMiscDoubleOption(int iopt, double opt)
{
  assert( 0 <= iopt && iopt < 50 );

  double_options_int[iopt] = true;
  misc_double_options[iopt] = opt;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::Options::ProbeMiscDoubleOption(int iopt) const
{
  assert( 0 <= iopt && iopt < 50 );
  return double_options_int[iopt];
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::Options::GetMiscBoolOption(int iopt) const
{
  assert( 0 <= iopt && iopt < 50 );

  if( !bool_options_int[iopt] )
    throw Exception( "Calorimeter::Options::GetMiscBoolOption(%i):  This option has "
		     "never been initialized!", iopt );
  return misc_bool_options[iopt];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Options::SetMiscBoolOption(int iopt, bool opt)
{
  assert( 0 <= iopt && iopt < 50 );

  bool_options_int[iopt] = true;
  misc_bool_options[iopt] = opt;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::Options::ProbeMiscBoolOption(int iopt) const
{
  assert( 0 <= iopt && iopt < 50 );
  return bool_options_int[iopt];
}

////////////////////////////////////////////////////////////////////////////////

time_t Calorimeter::Options::GetMiscTimeOption(int iopt) const
{
  assert( 0 <= iopt && iopt < 50 );
  if( !time_options_int[iopt] )
    throw Exception( "Calorimeter::Options::GetMiscTimeOption(%i):  This option has "
		     "never been initialized!", iopt );
  return misc_time_options[iopt];
}

////////////////////////////////////////////////////////////////////////////////

void Calorimeter::Options::SetMiscTimeOption(int iopt, time_t opt)
{
  assert( 0 <= iopt && iopt < 50 );

  time_options_int[iopt] = true;
  misc_time_options[iopt] = opt;
}

////////////////////////////////////////////////////////////////////////////////

bool Calorimeter::Options::ProbeMiscTimeOption(int iopt) const
{
  assert( 0 <= iopt && iopt < 50 );
  return time_options_int[iopt];
}

////////////////////////////////////////////////////////////////////////////////


} // namespace Reco
