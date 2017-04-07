#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,22,00)

#include <cmath>

#include "ToyCalorimeter.h"

#include <CalorimeterParticle.h>
#include <CellDataRaw.h>
#include <CellType.h>
#include <ReconstructionCombined.h>
#include <ShowerProfileLednev.h>

#include <TFile.h>
#include <TRandom3.h>

///////////////////////////////////////////////////////////////////////////
// global things
const unsigned int randomSeed(0);
const unsigned int maxNumber(10);

double F(const double& x, const double& y) {
    double result(0.);

    double showerProfileA[] = { 1.191, -0.191 };
    double showerProfileB[] = { 6.664, 43.777 };

    for (unsigned int i=0; i<2; i++)
        result += showerProfileA[i]
                  * (atan(x/showerProfileB[i])
                     + atan(y/showerProfileB[i])
                     + atan2(x*y, showerProfileB[i]*sqrt(showerProfileB[i]*showerProfileB[i]+x*x+y*y))
                    );
                
    result /= (2.*M_PI);
    result += 0.25;
    
    return result;
}

double G(const double& x, const double& y) {
    double sum(0.);

    sum+=F(x+38.3/2., y+38.3/2.);
    sum-=F(x+38.3/2., y-38.3/2.);
    sum-=F(x-38.3/2., y+38.3/2.);
    sum+=F(x-38.3/2., y-38.3/2.);

    return sum;
}

///////////////////////////////////////////////////////////////////////////
void pass00(TDirectory* dir) {
    std::cout << "pass 00: plot profiles" << std::endl;
    std::cout << "    processing ." << std::flush;

    Reco::CellType cellType("toyCellType", 37.3, 37.3, 45.,
                            38.3, 38.3,
                            2.59, 42., 0.014725,
                            0.065, 0.02, 0.01,
                            0.8686, 0.,
                            "GAMS");
    std::vector<double> params;
    params.push_back(1.191);
    params.push_back(-0.191);
    params.push_back(6.664);
    params.push_back(43.777);
    Reco::ShowerProfile* showerProfile = new Reco::ShowerProfileLednev(0, cellType, params);

    TH2D* profile = new TH2D("profile", "profile", 800, -160., 160., 800, -160., 160.);
    profile->SetOption("COLZ");
    TH2D* derivX  = new TH2D("derivX",  "derivX",  800, -160., 160., 800, -160., 160.);
    derivX ->SetOption("COLZ");
    TH2D* derivY  = new TH2D("derivY",  "derivY",  800, -160., 160., 800, -160., 160.);
    derivY ->SetOption("COLZ");

    TH2D* ratio1x1 = new TH2D("ratio1x1", "ratio1x1", 800, -160., 160., 800, -160., 160.);
    ratio1x1->SetMaximum(1.);
    ratio1x1->SetMinimum(0.);
    ratio1x1->SetOption("COLZ");
    TH2D* ratio3x3 = new TH2D("ratio3x3", "ratio3x3", 800, -160., 160., 800, -160., 160.);
    ratio3x3->SetMaximum(1.);
    ratio3x3->SetMinimum(0.);
    ratio3x3->SetOption("COLZ");
    TH2D* ratio5x5 = new TH2D("ratio5x5", "ratio5x5", 800, -160., 160., 800, -160., 160.);
    ratio5x5->SetMaximum(1.);
    ratio5x5->SetMinimum(0.);
    ratio5x5->SetOption("COLZ");
    TH2D* ratio7x7 = new TH2D("ratio7x7", "ratio7x7", 800, -160., 160., 800, -160., 160.);
    ratio7x7->SetMaximum(1.);
    ratio7x7->SetMinimum(0.);
    ratio7x7->SetOption("COLZ");

    for (int i=1; i<=profile->GetNbinsX(); i++) {
        if ( ((i+1)*50%(profile->GetNbinsX()/50))==0 )
            std::cout << "." << std::flush;
        double x = profile->GetXaxis()->GetBinCenter(i);
        for (int j=1; j<=profile->GetNbinsY(); j++) {
            double y = profile->GetYaxis()->GetBinCenter(j);
            profile->Fill(x, y, showerProfile->GetEnergyInCell(x, y));
            derivX ->Fill(x, y, showerProfile->GetDeDxInCell(x, y));
            derivY ->Fill(x, y, showerProfile->GetDeDyInCell(x, y));

            double ratio = showerProfile->GetEnergyInCell(remainder(x, 38.3), remainder(y, 38.3));
            ratio1x1->Fill(x, y, ratio);

            for (int k=-1; k<=1; k++) {
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+k*38.3, remainder(y, 38.3)-  38.3);
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+k*38.3, remainder(y, 38.3)+  38.3);
            }
            for (int k=0; k<=0; k++) {
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)-  38.3, remainder(y, 38.3)+k*38.3);
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+  38.3, remainder(y, 38.3)+k*38.3);
            }
            ratio3x3->Fill(x, y, ratio);

            for (int k=-2; k<=2; k++) {
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+k*38.3, remainder(y, 38.3)-2*38.3);
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+k*38.3, remainder(y, 38.3)+2*38.3);
            }
            for (int k=-1; k<=1; k++) {
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)-2*38.3, remainder(y, 38.3)+k*38.3);
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+2*38.3, remainder(y, 38.3)+k*38.3);
            }
            ratio5x5->Fill(x, y, ratio);

            for (int k=-3; k<=3; k++) {
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+k*38.3, remainder(y, 38.3)-3*38.3);
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+k*38.3, remainder(y, 38.3)+3*38.3);
            }
            for (int k=-2; k<=2; k++) {
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)-3*38.3, remainder(y, 38.3)+k*38.3);
                ratio += showerProfile->GetEnergyInCell(remainder(x, 38.3)+3*38.3, remainder(y, 38.3)+k*38.3);
            }
            ratio7x7->Fill(x, y, ratio);
        }
    }

    profile->Write(); delete profile;
    derivX ->Write(); delete derivX;
    derivY ->Write(); delete derivY;

    std::cout << "." << std::endl;

    delete showerProfile;
}

///////////////////////////////////////////////////////////////////////////
// pass 0
// test derivatives
void pass10(TDirectory* dir) {
    std::cout << "pass 10: test derivatives" << std::endl;
    std::cout << "    processing " << std::flush;

    TRandom3 zufall;
    zufall.SetSeed(randomSeed);

    ToyCalorimeter* toyCalorimeter = new ToyCalorimeter();
    toyCalorimeter->SetHistogramsBaseDir(dir->mkdir("Reco"));

    double maxDiff_f(0.);
    double maxDiff_g_d_f(0.);
    double maxDiff_g_d_g(0.);
    double maxDiff_g_f_g(0.);

    double maxDiff_g_check[4]  = { 0., 0., 0., 0. };
    double maxDiff_g_Rcheck[4] = { 0., 0., 0., 0. };
    for (unsigned int c=0; c<maxNumber; c++) {
        if ( ((c+1)*50%maxNumber)==0 )
            std::cout << "." << std::flush;
        double e = zufall.Uniform(0.2, 200.);
        double x = zufall.Uniform(-20., 20.);
        double y = zufall.Uniform(-20., 20.);
        double t = zufall.Uniform(-50., 50.);

        std::vector<Reco::ReconstructionCombined::FitInfo> cells;
        for (int i=-7; i<=7; i++) {
            for (int j=-7; j<=7; j++) {
                double eCell = e*G(j*38.3-x, i*38.3-y);
                cells.push_back(Reco::ReconstructionCombined::FitInfo(toyCalorimeter, Reco::CellDataRaw(15*(i+7) + (j+7), eCell, t)));
                cells.back().SetEnergyErr( sqrt(0.15*0.15*eCell + (0.015*eCell)*(0.015*eCell) + 0.05*0.05));
                cells.back().SetTimeErr( sqrt((0.82/eCell) + (1.71/(eCell*eCell)) + 0.34) );
            }
        }

        std::vector<Reco::ReconstructionCombined::FitShower> showers;
        showers.push_back(Reco::ReconstructionCombined::FitShower(e, x, y, 0., 0., 0., 0., 0.));
        Reco::ReconstructionCombined::FitInput fitInput(toyCalorimeter, cells, 4*showers.size());

        for (unsigned int i=0; i<1000; i++) {
            double par[4];
            par[0] = zufall.Uniform(0.2, 200.);
            par[1] = zufall.Uniform(-20., 20.);
            par[2] = zufall.Uniform(-20., 20.);
            par[3] = zufall.Uniform(-50., 50.);

            double fEval = fitInput.DoEval(par);

            double gDeri[4];
            for (unsigned int j=0; j<4; j++)
                gDeri[j] = fitInput.DoDerivative(par, j);
                
            double gGrad[4];
            fitInput.Gradient(par, gGrad);

            double fFdF, gFdF[4];
            fitInput.FdF(par, fFdF, gFdF);

            double gCheck[4];
            for (unsigned int j=0; j<4; j++) {
                double t1Par[4], t2Par[4];
                for (unsigned int k=0; k<4; k++) {
                    t1Par[k] = par[k];
                    t2Par[k] = par[k];
                }
                const double diff(1.e-4);
                t1Par[j] -= diff;
                t2Par[j] += diff;
                
                double v1 = fitInput.DoEval(t1Par);
                double v2 = fitInput.DoEval(t2Par);
                gCheck[j] = (v2-v1)/(2.*diff);
            }


            if (fabs(fFdF-fEval)>maxDiff_f)
                maxDiff_f=fabs(fFdF-fEval);

            for (unsigned int j=0; j<4; j++) {
                if (fabs(gDeri[j]-gFdF[j])>maxDiff_g_d_f)
                    maxDiff_g_d_f=fabs(gDeri[j]-gFdF[j]);
                if (fabs(gDeri[j]-gGrad[j])>maxDiff_g_d_g)
                    maxDiff_g_d_g=fabs(gDeri[j]-gGrad[j]);
                if (fabs(gFdF[j]-gGrad[j])>maxDiff_g_f_g)
                    maxDiff_g_f_g=fabs(gFdF[j]-gGrad[j]);

                if (fabs(gFdF[j]-gCheck[j])>maxDiff_g_check[j]) {
                    maxDiff_g_check[j]=fabs(gFdF[j]-gCheck[j]);
                    std::cout << i << ", " << j << ": " << maxDiff_g_check[j] << "(" << gFdF[j] << "," << gCheck[j] << ")" << std::endl;
                }
                if ((fabs(gFdF[j]-gCheck[j])/gFdF[j])>maxDiff_g_Rcheck[j])
                    maxDiff_g_Rcheck[j]=fabs(gFdF[j]-gCheck[j])/gFdF[j];
            }
        }
    }

    dir->cd();

    delete toyCalorimeter;

    std::cout << " finished." << std::endl;

    std::cout << "maximum difference for function value ('DoEval' and 'FdF'):                 " << maxDiff_f << std::endl;
    std::cout << "maximum difference for derivative value ('DoDerivative' and 'FdF'):         " << maxDiff_g_d_f << std::endl;
    std::cout << "maximum difference for derivative value ('DoDerivative' and 'Gradient'):    " << maxDiff_g_d_g << std::endl;
    std::cout << "maximum difference for derivative value ('FdF' and 'Gradient'):             " << maxDiff_g_f_g << std::endl << std::endl;
    std::cout << "maximum difference for derivative between from values and function:         " << maxDiff_g_check[0] << std::endl;
    std::cout << "                                                                            " << maxDiff_g_check[1] << std::endl;
    std::cout << "                                                                            " << maxDiff_g_check[2] << std::endl;
    std::cout << "                                                                            " << maxDiff_g_check[3] << std::endl;
    std::cout << "maximum difference for derivative between from values and function (ratio): " << maxDiff_g_Rcheck[0] << std::endl;
    std::cout << "                                                                            " << maxDiff_g_Rcheck[1] << std::endl;
    std::cout << "                                                                            " << maxDiff_g_Rcheck[2] << std::endl;
    std::cout << "                                                                            " << maxDiff_g_Rcheck[3] << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// pass 1
// one cluster, no cell thresholds
void pass20(TDirectory* dir) {
    std::cout << "pass 20: one cluster, no cell thresholds" << std::endl;
    std::cout << "    processing " << std::flush;

    TRandom3 zufall;
    zufall.SetSeed(randomSeed);

    ToyCalorimeter* toyCalorimeter = new ToyCalorimeter();
    toyCalorimeter->SetHistogramsBaseDir(dir->mkdir("Reco"));

    TH1D* hCount = new TH1D("hCount", "number of reconstructed showers",  10, -0.5, 9.5);
    TH1D* hResE  = new TH1D("hResE",  "input - shower energy",           200, -20., 20.);
    TH1D* hResX  = new TH1D("hResX",  "input - shower position in X",    200, -20., 20.);
    TH1D* hResY  = new TH1D("hResY",  "input - shower position in Y",    200, -20., 20.);
    for (unsigned int c=0; c<maxNumber; c++) {
        if ( ((c+1)*50%maxNumber)==0 )
            std::cout << "." << std::flush;
        double x = zufall.Uniform(-20., 20.);
        double y = zufall.Uniform(-20., 20.);
        double e = zufall.Uniform(0.2, 200.);

        std::vector<Reco::CellDataRaw> cells;
        for (int i=-7; i<=7; i++) {
            for (int j=-7; j<=7; j++) {
                cells.push_back(Reco::CellDataRaw(15*(i+7) + (j+7), e*G(j*38.3-x, i*38.3-y), 0.));
                cells.back().SetTimeErr(1e20);
            }
        }

        toyCalorimeter->StoreData(cells);
        toyCalorimeter->Reconstruction();
        std::vector<Reco::CalorimeterParticle> showers = toyCalorimeter->GetCalorimeterParticles();

        hCount->Fill(showers.size());
        for (std::vector<Reco::CalorimeterParticle>::const_iterator it=showers.begin(); it!=showers.end(); it++) {
            hResE->Fill(e - it->GetE());
            hResX->Fill(x - it->GetX());
            hResY->Fill(y - it->GetY());
        }
    }

    dir->cd();
    hCount->Write(); delete hCount;
    hResE->Write();  delete hResE;
    hResX->Write();  delete hResX;
    hResY->Write();  delete hResY;

    delete toyCalorimeter;

    std::cout << " finished." << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// pass 2
// one cluster, with cell thresholds (0.2 GeV)
void pass30(TDirectory* dir) {
    std::cout << "pass 30: one cluster, with cell thresholds (0.2 GeV)" << std::endl;
    std::cout << "    processing " << std::flush;

    TRandom3 zufall;
    zufall.SetSeed(randomSeed);

    ToyCalorimeter* toyCalorimeter = new ToyCalorimeter(0.2);
    toyCalorimeter->SetHistogramsBaseDir(dir->mkdir("Reco"));

    TH1D* hCount = new TH1D("hCount", "number of reconstructed showers",  10, -0.5, 9.5);
    TH1D* hResE  = new TH1D("hResE",  "input - shower energy",           200, -20., 20.);
    TH1D* hResX  = new TH1D("hResX",  "input - shower position in X",    200, -20., 20.);
    TH1D* hResY  = new TH1D("hResY",  "input - shower position in Y",    200, -20., 20.);
    for (unsigned int c=0; c<maxNumber; c++) {
        if ( ((c+1)*50%maxNumber)==0 )
            std::cout << "." << std::flush;
        double x = zufall.Uniform(-20., 20.);
        double y = zufall.Uniform(-20., 20.);
        double e = zufall.Uniform(0.2, 200.);

        std::vector<Reco::CellDataRaw> cells;
        for (int i=-7; i<=7; i++) {
            for (int j=-7; j<=7; j++) {
                cells.push_back(Reco::CellDataRaw(15*(i+7) + (j+7), e*G(j*38.3-x, i*38.3-y), 0.));
                cells.back().SetTimeErr(1e20);
            }
        }

        toyCalorimeter->StoreData(cells);
        toyCalorimeter->Reconstruction();
        std::vector<Reco::CalorimeterParticle> showers = toyCalorimeter->GetCalorimeterParticles();

        hCount->Fill(showers.size());
        for (std::vector<Reco::CalorimeterParticle>::const_iterator it=showers.begin(); it!=showers.end(); it++) {
            hResE->Fill(e - it->GetE());
            hResX->Fill(x - it->GetX());
            hResY->Fill(y - it->GetY());
        }
    }

    dir->cd();
    hCount->Write(); delete hCount;
    hResE->Write();  delete hResE;
    hResX->Write();  delete hResX;
    hResY->Write();  delete hResY;

    delete toyCalorimeter;

    std::cout << " finished." << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// pass 3
// one cluster, with cell thresholds (0.2 GeV), with noise
void pass40(TDirectory* dir) {
    std::cout << "pass 40: one cluster, with cell thresholds (0.2 GeV), with noise" << std::endl;
    std::cout << "    processing " << std::flush;

    TRandom3 zufall;
    zufall.SetSeed(randomSeed);

    ToyCalorimeter* toyCalorimeter = new ToyCalorimeter(0.2);
    toyCalorimeter->SetHistogramsBaseDir(dir->mkdir("Reco"));

    TH1D* hCount = new TH1D("hCount", "number of reconstructed showers",  10, -0.5, 9.5);
    TH1D* hResE  = new TH1D("hResE",  "input - shower energy",           200, -20., 20.);
    TH1D* hResX  = new TH1D("hResX",  "input - shower position in X",    200, -20., 20.);
    TH1D* hResY  = new TH1D("hResY",  "input - shower position in Y",    200, -20., 20.);
    for (unsigned int c=0; c<maxNumber; c++) {
        if ( ((c+1)*50%maxNumber)==0 )
            std::cout << "." << std::flush;
        double x = zufall.Uniform(-20., 20.);
        double y = zufall.Uniform(-20., 20.);
        double e = zufall.Uniform(0.2, 200.);

        std::vector<Reco::CellDataRaw> cells;
        for (int i=-7; i<=7; i++) {
            for (int j=-7; j<=7; j++) {
                cells.push_back(Reco::CellDataRaw(15*(i+7) + (j+7), e*G(j*38.3-x, i*38.3-y) + zufall.Gaus(0, 0.07), 0.));
                cells.back().SetTimeErr(1e20);
            }
        }

        toyCalorimeter->StoreData(cells);
        toyCalorimeter->Reconstruction();
        std::vector<Reco::CalorimeterParticle> showers = toyCalorimeter->GetCalorimeterParticles();

        hCount->Fill(showers.size());
        for (std::vector<Reco::CalorimeterParticle>::const_iterator it=showers.begin(); it!=showers.end(); it++) {
            hResE->Fill(e - it->GetE());
            hResX->Fill(x - it->GetX());
            hResY->Fill(y - it->GetY());
        }
    }

    dir->cd();
    hCount->Write(); delete hCount;
    hResE->Write();  delete hResE;
    hResX->Write();  delete hResX;
    hResY->Write();  delete hResY;

    delete toyCalorimeter;

    std::cout << " finished." << std::endl;
}

int main(int argc, char* argv[]) {
    TFile out("hist.root", "RECREATE");
    
    pass00(out.mkdir("pass00"));
    pass10(out.mkdir("pass10"));
    pass20(out.mkdir("pass20"));
    pass30(out.mkdir("pass30"));
    pass40(out.mkdir("pass40"));

    out.Write();
    out.Close();

    return 0;
}

#else

#include <iostream>

int main(int argc, char* argv[]) {
    std::cerr << "Need at least ROOT 5.22 to be able to run this test!" << std::endl;
    return 1;
}

#endif
