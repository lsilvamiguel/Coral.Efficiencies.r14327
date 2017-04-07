#undef  __STRICT_ANSI__
#include "CoralUser.h"
#include "CsOraStore.h"

using namespace std;

static CsOraStore* store;

bool CoralUserSetup(int argc, char* argv[])
{
  return true;
}

bool CoralUserInit()
{
  store = CsOraStore::Instance();
  return true;
}

bool CoralUserEvent()
{
  if(store)
    {
      unsigned eventSize = store->getRawBufferLength();
      uint8* buffer = store->rawBuffer();
      cout << "Event size: " << eventSize << endl;
      for(int i = 0; i < 10; i++)
	cout << "byte[" << i << "] = " << (int)buffer[i] << endl;
    }
  return true;
}

bool CoralUserEnd()
{
  return true;
}
