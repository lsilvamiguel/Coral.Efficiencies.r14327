// User functions b.g. 9/10/2000

// Called after initialization
void CoralUserInit( void );

// Called after event readout
void CoralUserEvent( void );

// Called before end of program
void CoralUserEnd( void );

void  FillDeltaSADC( const CS::Chip::Digits& digits);
// std::pair <int,int>  MiniMaxSADC(  const std::vector<uint16>& samples );
