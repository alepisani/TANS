#include "../include/const.h"
#include "../include/event.h"
#include "../include/point.h"
#include "../include/particle.h"
#include "../include/reconstruction.h"
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
    //reconstruction reco;
    ev.RunFullSimulation();
    
    //ev.eventsimulation();
    //ev.display_event();

    //double z_reco = reco.reco_z(ev);
    //cout << "Z ricostruito: " << z_reco << " mm" << endl;
    





    
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
    
    return 0;

}