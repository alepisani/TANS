#include "../include/const.h"
#include "../include/event.h"
#include "../include/point.h"
#include "../include/particle.h"
#include "../include/reconstruction.h"
#include "../include/simulation.h"
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
    
    simulation simu;
    simu.sim(); 
    //reconstruction
    //analysis

    /**     
     * reconstruction()
     *      leggi il ttree di simulation
     *      for(getentries)
     *          create vector di hit e tracklet
     *          reserve sui vector con multiplicity
     *          for(hitL12)
     *              riempi vector
     *              tracklet
     *              taglio in phi
     *              histogramma di z_rec, prendi la moda e i due bin vicini
     *              cosa intendo per tracce non ricostruite?
     *                  se ho ricostruito < 1 (2?) tracce 
     *                  se hist ha solo 1 entries per ogni bin
     *                  come tratti picchi multipli nell'istogramma? 
     *                      fai scan manuale e consideri qui bin sopra una certa soglia (correlata con il picco massimo)
     *              trovato il massimo fai la running window sui valori di z nel bin della moda (e quelli adiacenti)
     *          riempi ttree
     * 
     * analysis()
     *      legge ttree di ricostruzione e simulazione e fa i plot del caso
     */

    //event ev;
    //ev.RunFullSimulation();
    
    t.Stop();

    cout << endl;
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "| time to process all the " << nEvents << " event and simulation \n";
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