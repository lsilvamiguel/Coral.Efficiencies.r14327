#ifndef CondDbHandlerError_h
#define CondDbHandlerError_h
#include <iostream>
#include <string>
#include <exception>

class CondDbHandlerError : public std::exception{
public:
  CondDbHandlerError(string str) : _str(str){};
  void Print(){cout << "Exceptin : ";cout<< _str<< endl;}
  void Abort(){cout << "Exceptin : ";cout<< _str<< endl;abort();}
  void Exit(){cout << "Exceptin : ";cout<< _str<< endl;exit(1);}
private:
  string _str;
};
#endif



