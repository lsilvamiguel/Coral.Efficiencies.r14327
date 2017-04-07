
#include <iostream>
#include "Variable.h"

ClassImp(Variable);


void Variable::MyDump() {
  std::cout<<"-------------------------------"<<std::endl;
  std::cout<<fName<<" "<<fNvalues<<" "<<fMin<<" "<<fMax<<" "<<std::endl;
  for(int i=0;i<fNvalues;i++) 
    std::cout<<fValues[i]<<" ";
  std::cout<<std::endl;
}

