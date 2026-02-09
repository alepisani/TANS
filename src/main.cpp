#include "../include/const.h"
#include "../include/event.h"
#include "../include/point.h"
#include "../include/particle.h"
#include <iostream>
#include "TApplication.h" 
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoTube.h"
#include "TColor.h"
#include "TVirtualGeoTrack.h" 
#include "TGeoTrack.h"
#include "TPolyLine3D.h"
#include <thread>
#include "TSystem.h" 
#include "TRandom3.h"
#include "TF1.h"
#include "TCanvas.h"
using namespace std;


int main(int argc, char **argv) {

    TApplication app("app", &argc, argv);

    int seed = 0;
    gRandom->SetSeed(seed);

    event ev;
    ev.RunFullSimulation();
    
    //ev.eventsimulation();
    //ev.display_event();


    std::atomic<bool> shouldExit(false);
    
    std::thread inputThread([&shouldExit]() {
        std::cin.get();
        shouldExit = true;
    });
    
    while (!shouldExit) {
        gSystem->ProcessEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    
    inputThread.join();
    

    /*TF1 *f_mult = new TF1("f_mult", "[0]*ROOT::Math::gamma_pdf(x, [1], [2])", 0, 150);

    // Esempio parametri per LHC (valori indicativi)
    double n_mean = 30.0; 
    double k = 2.0; 
    f_mult->SetParameters(1, k, n_mean/k);
    TCanvas *c_mult = new TCanvas("c_mult", "Multiplicity Distribution", 800, 600);
    c_mult->cd();
    f_mult->Draw();
    c_mult->SaveAs("../data/multiplicity_distribution.png");*/
    
    /*TF1 *f_eta = new TF1("f_eta", "[0]*(1+[2]*x*x) / (cosh([1]*x)*cosh([1]*x))", -2.5, 2.5);
    f_eta->SetParameters(1.0, 0.4, 0.18); // Esempio di parametri

    TCanvas *c1 = new TCanvas("c1", "Distribution of eta", 800, 600);
    c1->cd();
    f_eta->Draw();
    c1->SaveAs("../data/eta_distribution.png");*/

    return 0;

}