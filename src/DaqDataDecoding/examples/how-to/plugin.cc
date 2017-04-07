#include "DaqEventsManager.h"

using namespace CS;

extern "C"
{
void ddd_call(const CS::DaqEventsManager &m);
}

void ddd_call(const CS::DaqEventsManager &m)
{
    printf("Plugin: %zu digits.\n",m.GetEventDigits().size());
}
