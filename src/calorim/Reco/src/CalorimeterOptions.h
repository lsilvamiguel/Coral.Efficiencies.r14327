#ifndef inside_Calorimeter_h
#  error Only for inclusion from inside Calorimeter.h!
#endif


class Options
{
 public:
  // =========================================================================
  // Enumeration for histogram level
  // =========================================================================
  enum HistLevel{
   HISTLVL_NORMAL,
   HISTLVL_VERBOSE,
   HISTLVL_DEBUG
  };

  //! Constructor
  Options();

  /*! \brief Parsing of RecoOptions

    RecoOptions are given as a whitespace-separated list of key=value
    pairs.  No space is allowed around the = sign.  Everything is
    case-insensitive.

    Some of the available options are:
    \item debug=#            Set debug level. 0 - no debug
    \item calib=yes/no       Set calibration procedure ON or OFF
    \item calib_energy=#     Set calibration energy in GeV
    \item fit=simple/normal  Set fit function to simple or normal
    \item hist=yes/no        Set filling histograms ON or OFF
    \item particle_energy_threshold=# set particles energy threshhold
    \item cluster_search=yes/no
    \item use_hint_particles=yes/no

    Read the source for more options and to find out the default values.
  */
  void SetRecoOptions( const std::string &opt_ );



  bool Ok   ( void ) const { return options_ok; }
  void SetOK( void )       { options_ok = true;  }
  void More( void )        { options_ok = false; }
  void ResetOutputInfoNames ( void ) { output_info_names.clear(); }
  void ResetInputInfoNames ( void )  { input_info_names.clear(); }

  /// \return true if option was initialized
  bool  		ProbeMiscDoubleOption	  (int iopt) const;
  /// \return external bool otion iopt
  double		GetMiscDoubleOption   (int iopt) const;
  /// \return external bool otion iopt
  void  		SetMiscDoubleOption   (int iopt,double opt);
  /// \return true if option was initialized
  bool  		ProbeMiscTimeOption	(int iopt) const;
  /// \return external double option iopt
  time_t		GetMiscTimeOption   (int iopt) const;
  /// \return external bool option iopt
  void  		SetMiscTimeOption   (int iopt,time_t opt);
  /// \return true if option was initialized
  bool  		ProbeMiscBoolOption	(int iopt) const;
  /// \return external bool otion iopt
  bool  		GetMiscBoolOption     (int iopt) const;
  /// \return external bool otion iopt
  void  		SetMiscBoolOption     (int iopt,bool opt);

  // ==========================================
  //  Attributes, data
  // ==========================================

  int                   debug;

  /// The calorimeter's type.
  CalorimeterType       calo_type;

  /// Default particle type.
  CalorimeterParticle::ParticleID particle_default;

  /// Have calo blocks of different sizes? This is not an option actually.
  bool                  mixed_blocks;

  /// Options for data processing.
  bool                  is_real_data;
  /// Options for readout electronics etc.
  bool                  readout_sparsified_;
  double                delta_spars_mode_;

  // Options for calibration.
  bool                  calibration_flag;
  bool                  calib_histo_units_in_gev;
  double                calibration_energy;
  double                calibration_energy_min;
  double                calibration_energy_max;
  double                calibration_min;
  double                calibration_max;
  double                default_calibration;
  double                default_time0_calibration;
  bool                  enable_time_calibration;
  double                calibration_time_gate;
  double                scale_calibration;
  double                tolerate2keepold_calib;
  bool                  scale_calibration_by_led;
  bool                  use_ProbEnT;
  bool                  use_cell_calib_correction;
  bool                  use_led_ref_calib_correction;
  bool                  use_led_ref_calib_correction_in_spills;

/// energy-dependent per-cell calibration corrections
  std::string           edep_corr;
  /// time-in-spill-dependent per-cell calibration corrections
  std::string           tisdep_corr;
  double                vertex_position    [3];
  double                position_correction[3];
  double                scale_cells_threshols;
  int                   calib_histo_size;
  double                calib_histo_range;
  double                new_calib_stat_min;
  bool                  recalc_from_primary_calib;
  bool                  recalc_from_primary_time_calib;
  double                position_offset[3]; // offsets for position counters

  double                calibration_uncertainty;

  // Options for data preparations (NOISE etc.)
  bool                  store_noise_info;
  bool                  store_rndm_noise_info;
  bool                  store_gamma_noise_info;
  bool                  store_noise_histo;
  bool                  update_ecut_using_prob_cell_noise;
  double                prob_cell_noise;
  double                prob_gamma_noise;
  double                stat_min_cell_noise;
  double                stat_min_gamma_noise;

  double                base_energy_noise;

  bool                  store_leds_all;
  bool                  store_in_spill_leds_all;
  bool                  leds_all_tolerate_bad_preformance;
  double                led_statistic_inspect;
  bool                  reset_after_led_inspect;

  bool                  scale_old_calib_info;

  double                ecut_extend_list;

  double                readout_term;    // reseting readout term; loosy solution

  /// Options for on-line, if true skip proper calculation of energy in cells
  bool                  flag_fast_online;

  bool                  print_general_info;
  bool                  print_unmaped_data;
  bool                  print_bad_data;

  /// Preshower options
  bool                  use_preshower_correction;
  double                preshower_length;

  /// General options for reconstruction
  bool                  do_reconstruction;
  bool                  combined_reconstruction;    // Future standard reconstruction
  bool                  kolosov_reconstruction;     // Previous standard reconstruction
  bool                  fortran_reconstruction;     // Alternative reconstruction of Anatoly Lednev
  bool                  use_time_in_reconstruction; // Reconstruction of Anatoly Lednev can use time info
  bool                  calib_event_reconstruction; // call the reconstruction for calibration events

  bool                  ecal0_nonlinear_calib_correction;

  int                   repeat_reconstruction;
  bool                  reco_cluster_search;
  bool                  reco_cluster_search_only;
  bool                  reco_use_hint_particles_info;
  bool                  recover_bad_cells;
  bool                  add_time_to_reco_particles;
  bool                  correct_for_digitization;

  /// Specific options for reconstruction
  bool                  use_simple_fit;
  FitMethod             fitting_type;
  bool                  add_search;
  double                hicut_add_search;
  double                hicut_add_search_reco2;
  double                ecut_add_search;

  /// Do we want to fill test-histograms?
  bool                  fill_histos;
  bool                  fill_only_reco_histos_in_reco_test;
  HistLevel             histos_level;
  bool                  make_profiles;
  bool                  fill_raw_histos;
  bool                  test_histo_raw_info;
  int                   test_histo_raw_info_level;
  bool                  fill_fit_info_histo;
  bool                  fill_internal_histos;
  bool                  fill_internal_correlations;
  bool                  fill_external_correlations;
  bool                  fill_calib_spill_correlations_histo;
  double                hist_energy_max;
  double                hist_time_min;
  double                hist_time_max;
  int                   hist_time_nbins_1d;
  int                   hist_time_nbins_2d;
  int                   hist_energy_nbins_1d;
  int                   hist_energy_nbins_2d;
  bool                  need_cells_XY_histo;

  int                   led_histo_size;
  double                hist_led_max;
  int                   level_particle_calo_association_histo;
  int                   ngam_max4histo;

  /// Options for reconstruction parameters
  double                particle_energy_threshold;
  double                cell_energy_threshold;

  double                cluster_search_cell_energy_threshold;
  double                cluster_search_cell_ampl_deviation;

  /// Options for specific reconstruction metods
  double                energy_leak_for_SimpleFit;
  double                distmax_for_cells_search;

  /// Options for "fortran" reconstruction
  int                   fortran_bad_cell_check;
  int                   fortran_fast_reconstruction;
  int                   fortran_hbook;
  int                   fortran_variable_thresholds;
  double		ecal2_time_0;
  double		ecal2_time_width;
  std::vector<double>   fortran_showerprof_a;
  std::vector<double>   fortran_showerprof_b;

  /// Options for "combined" reconstruction
  bool                  correct_longitudinal_leakage;
  double                min_shower_distance;
  double                cell_adc_threshold;
  bool                  use_combined_calibration;
  std::string           combined_calibration;
  std::vector<double>   param_energy_error;
  std::vector<double>   param_time_error;
  std::vector<int>      combined_allowed_showers;

  /// Store RawData from CalorimeterParticle, option is used during calibration
  bool                   store_back_data;

  /// Options specific for initial geometry settings
  double                tolerance_for_nearby_cells;
  bool                   check_front_surface_position;
  bool                   init_xy_structure;

  /// GUI options
  bool                  show_all_cells;

  /// Monitoring options
  bool                  monitor_histo_show_fit;
  bool                  print_bad_cells_info;
  int                   monitor_db_calib_level;
  double                monitor_db_calib_min_value;
  double                monitor_db_calib_max_value;
  int                   monitor_db_led_level;
  double                monitor_db_led_min_value;
  double                monitor_db_led_max_value;
  int                   monitor_db_porog_level;
  double                monitor_db_porog_min_value;
  double                monitor_db_porog_max_value;

  /// Data inspection options
  int                   inspect_led_Ncells_no_signal;
  int                   inspect_led_Ncells_small_signal;
  double                inspect_led_Fcells_good_signal;


  /// Options for data processing
  bool                  fit_led_histo_at_end_of_job;

  /// Options for debugging.
  bool                  debug_reconstruction;
  bool                  debug_reconstruction_test;
  bool                  debug_reconstruction_combined;

  /// More specific options
  double                ecut4mgg_histo;

  /// Misc options (normaly for external usage of information related to Calorimeter)
  double                misc_double_options[50];
  bool                  double_options_int[50];
  bool                  misc_bool_options[50];
  bool                  bool_options_int[50];
  time_t                misc_time_options[50];
  bool                  time_options_int[50];

  // Option to sort particles
  bool                  sort_particles;
  std::ostream*         blackbox_chk_file;

  int                   digitizer_method;
  int                   digitizer_method_led;


  std::vector< std::string>    output_info_names;
  std::vector< std::string>    input_info_names;



  // =========================================================================
  // stuff from former OptionsMC class
  // =========================================================================

  double            mc_default_calibration;
  /// Fast-Monte-Carlo options
  bool              add_shower_fluctuations_for_FMC;
  double            mc_data_extraction_threshold;

  bool              mc_make_real_digitization;
  bool              mc_make_fiadc_digitization;
  double            mc_fiadc_sparce_delta;
  /// Monte-Carlo options for data smearing
  bool              mc_smear_response;
  bool              mc_smear_constant;
  bool              mc_smear_stochastic;
  bool              mc_smear_readout;


 private:
  bool              options_ok;

};

