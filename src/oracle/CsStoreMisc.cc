#include "CsStoreMisc.h"
#include <cstdio>
#include <vector>

using namespace std;

string itostr(int i) {
   const int buflen = 100;
   char	buff[buflen];
   snprintf(buff, buflen, "%d", i);
   return buff;
}

string utostr(unsigned u) {
   const int buflen = 100;
   char	buff[buflen];
   snprintf(buff, buflen, "%u", u);
   return buff;
}

void AddValueToVector(vector<uint8>& buffer, const void* value, size_t size)
{
  const uint8* p = (const uint8*)value;
  for(size_t i = 0; i < size; i++)
    buffer.push_back(p[i]);
}

