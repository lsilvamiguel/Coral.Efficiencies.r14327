#ifndef inside_Calorimeter_h
#  error Only for inclusion from inside Calorimeter.h!
#endif

    class AlignmentInTime
    {
      public:
        AlignmentInTime ( void ) {Clear(); Configure(true,true,true);}
        void Clear ( void ) {map_time2position_.clear();}
        void Configure ( bool use_relative_counter_x, bool use_relative_counter_y, bool interpolate ) {
                                       use_relative_counter_x_= use_relative_counter_x;
                                       use_relative_counter_y_= use_relative_counter_y;
                                       interpolate_ = interpolate; }
        std::pair<bool,std::vector< double > > GetPosition ( const time_t  &time ) const;
      public:
      bool use_relative_counter_x_;
      bool use_relative_counter_y_;
      bool interpolate_;
      std::vector< std::pair< time_t, std::vector< double> >  >          map_time2position_;
    };
