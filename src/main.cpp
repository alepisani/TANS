#include "../include/const.h"
#include "../include/evento.h"
#include "../include/point.h"
#include "../include/tracklet.h"
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
using namespace std;


int main(int argc, char **argv) {

    TApplication app("app", &argc, argv);

    int seed = 0;
    gRandom->SetSeed(seed);

    evento ev;
    ev.setmultiplicity();
    ev.generate_vertex();
    ev.event();
    for(int i = 0; i < ev.points_L1.size(); i++){
        cout << "L1" << endl;
        cout << ev.points_L1[i] << endl;
    }
    for(int i = 0; i < ev.points_L2.size(); i++){
        cout << "L2" << endl;
        cout << ev.points_L2[i] << endl;
    }
    ev.smearing();
    for(int i = 0; i < ev.points_L1.size(); i++){
        cout << "L1s" << endl;
        cout << ev.points_L1[i] << endl;
    }
    for(int i = 0; i < ev.points_L2.size(); i++){
        cout << "L2s" << endl;
        cout << ev.points_L2[i] << endl;
    }
    ev.display_event();





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