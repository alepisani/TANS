#include "../include/const.h"
#include "../include/evento.h"
#include "../include/traccia.h"
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
using namespace std;


int main(int argc, char **argv) {

    TApplication app("app", &argc, argv);

    tracklet trkl;
    cout << trkl << endl;
    trkl.generate_eta();
    trkl.generate_phi();
    trkl.generate_theta();
    trkl.generate_vertex();
    cout << trkl << endl;

    evento ev;
    ev.display_event();

    std::atomic<bool> shouldExit(false);
    
    std::thread inputThread([&shouldExit]() {
        std::cin.get();
        shouldExit = true;
    });
    
    // Loop custom invece di app.Run()
    while (!shouldExit) {
        gSystem->ProcessEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    
    inputThread.join();
    return 0;

    return 0;
}