#include "MMan_2triplets.h"

#include "TROOT.h"

extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("simpleMP","DATE ROOT Monitor Module",initfuncs);

int main() {
  MMan mm("/home/data2002/root",21510,"",0,2,0.9);
  mm.Loop(20000,0);
}

