
#ifndef TWatches_h
#define TWatches_h

#include <iostream>
#include <iomanip>
#include <cstring>   // for strlen()

float SysTime(float&);

/*!
  \brief Stopwatch
  Multi (20) channel watches for code timing
 
*/

class TWatches {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;


private:

  signed char   state  [21]; // 0 - ready    1 - started     -1 - stoped
  float  fStarts[21];
  float  fStops [21];
  double fSum   [21];
  
public:
  //! Constructor
  TWatches()   
  {
    NobjCreated++;
    for (int i=0; i < 21; i++) {
      state[i]= 0;
      fSum[i] = 0;
    }
  };

  //! Destructor
  ~TWatches()   
  {
    NobjDestructed++;
  }

  //! Reset all stopwatches
  void Reset() {
    for (int i=0; i < 21; i++) {
      state[i]= 0;
      fSum[i] = 0;
    }
  };
  
  //! Start stopwatch # i
  void Start(int i) {
    if(i > 20) {
      std::cout<<"Stop Watch # " <<i<<"   Max. possible = 20"<<std::endl;
      return;
    }
    if(state[i] == 1){
      std::cout<<"Stop Watch # "<<i<<"  already started "<<std::endl;
      return;
    } else {
      state[i] = 1;
      SysTime(fStarts[i]);
    }
    return;
  };

  //! Stop stopwatch # i
  float Stop(int i){
    if(i > 20) {
      std::cout<<"Stop Watch # " <<i<<"   Max. possible = 20"<<std::endl;
      return(0);
    }
    if(state[i] == 1){
      state[i] = -1;
      SysTime(fStops[i]);
      fSum[i] += fStops[i]-fStarts[i];
      return(fStops[i]-fStarts[i]);
    } else {
      std::cout<<"Stop Watch # "<<i<<"  has not been started "<<std::endl;
      return(0);
    }
  };

  //! Returns measured time in sec for stapwatch # i
  float Time(int i){
    if(i > 20) {
      std::cout<<"Stop Watch # " <<i<<"   Max. possible = 20"<<std::endl;
      return(0);
    }
    if(state[i] == -1){
      return(fStops[i]-fStarts[i]);
    } else {
      std::cout<<"Stop Watch # "<<i<<"  has not been started/stoped"<<std::endl;
      return(0);
    }
  };

 //! Returns sum of measured times in sec for stapwatch # i
  float SumTime(int i){
    if(i > 20) {
      std::cout<<"Stop Watch # " <<i<<"   Max. possible = 20"<<std::endl;
      return(0);
    }
    if(state[i] == -1){
      return(fSum[i]);
    } else {
      std::cout<<"Stop Watch # "<<i<<"  has not been started/stoped"<<std::endl;
      return(0);
    }
  };

  /*!
    Print measured time for stopwatch #i with some comment (optional).
    With no arguments, print time for all stopwatches
  */  
  void Print(int i=0, const char* comment = ""){
    if(i > 20) {
      std::cout<<"Stop Watch # " <<i<<"   Max. possible = 20"<<std::endl;
      return;
    }
    if(i == 0){
      for(int j =0; j < 20; j++){
	if(state[j] == -1) {
	  std::cout<<"Stop Watch # "<<j<<" :  "<<std::setprecision(5)<<std::setw(9)<<(fStops[j]-fStarts[j])<<" sec. ";
	  std::cout<<std::endl;
	}
      }
    } else {
      if(state[i] == -1) {
	if(strlen(comment) == 0) {
	  std::cout<<"Stop Watch # "<<i<<" :  "<<std::setprecision(5)<<std::setw(9)<<(fStops[i]-fStarts[i])<<" sec. ";
	  std::cout<<std::endl;
	} else {
	  std::cout<<comment<<" "<<std::setprecision(5)<<std::setw(9)<<(fStops[i]-fStarts[i])<<" sec. "<<std::endl;
	} 
      }
    }
  };

};
#endif







