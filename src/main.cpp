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
#include <TStopwatch.h>
using namespace std;

int main() {

    int seed = 0;
    gRandom->SetSeed(seed);

    TStopwatch t;
    t.Start();
    
    event ev;
    ev.RunFullSimulation();
    
    t.Stop();

    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "| time to process all the event and simulation \n";
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "| Real time: " << t.RealTime() << " s\n";
    cout << "| CPU time: " << t.CpuTime()  << " s\n";
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << boolalpha << "| multiple scattering: " << multiple_scattering_on << endl;
    cout << "| distribution used: ";
    if(get_data_from_kinem){
        cout << " from ./data/kinem.root" << endl;
    }
    else {cout << " uniform distribution" << endl;}
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;

    cout << endl;
    
    return 0;

}