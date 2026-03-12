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

    cout << endl;
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "| MonteCarlo simulation: " << endl; 
    cout << "+-----------------------------------------------+" << endl;
    cout << "| nEvents = " << nEvents << endl;
    cout << boolalpha << "| multiple scattering: " << multiple_scattering_on << endl;
    cout << "| distribution used: ";
    if(get_data_from_kinem){
        cout << " from ./data/kinem.root" << endl;
    }
    else {cout << " uniform distribution" << endl;}
    cout << "| selection on delta phi = " << delta_phi << " Rad" << endl;
    cout << "| running window width = " << half_window * 2 << " mm" << endl;
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << endl;

    TStopwatch t;
    t.Start();
    
    simulation simu;
    simu.sim(); 
    cout << endl; cout << endl; 
    
    reconstruction reco;
    reco.reco();
    //analysis

    /**     
     * reconstruction()
     *      voglio che apra MC/simualtion e crei reconstrucion con i dati che mi interessano
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
     *      voglio che apre MC/reconstrucion e scrive i plot in MC/analysis
     *      se non ricostruisco qualcosa per qualche motivo mi faccio restituire un dummy value
     *      il dummy value lo uso per stimare efficeinza
     */

    /**
     * AGGIUNGI LA VERTIà MONTECARLO NEL FILE simulation.root (?) oppure
     * apri due file contemporaneamente e prendi info 
     * add z_true in simulation.root
     */


    
    t.Stop();

    cout << endl;
    cout << endl;
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "| time to process all the " << nEvents << " event and simulation \n";
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << "| Real time: " << t.RealTime() << " s\n";
    cout << "| CPU time: " << t.CpuTime()  << " s\n";
    cout << "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+" << endl;
    cout << endl;
    
    return 0;

}