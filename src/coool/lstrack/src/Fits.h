#include "TH1.h"

//func = "gauspol0" or "gaus"
void FitCTime(TH1F *h, float& sigma, float& mean, float min=-500,
	      float max=500, float minbg=-500,float maxbg=-200);

void FitResiduals(TH1F *h, float& sigma, float& mean);
