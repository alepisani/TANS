#include "../include/const.h"
#include "../include/evento.h"
#include "../include/traccia.h"
#include "../include/point.h"
#include "../include/tracklet.h"
#include <iostream>
#include "TApplication.h" // NECESSARIO per app stand-alone
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoTube.h"
#include "TColor.h"
#include "TVirtualGeoTrack.h" // NECESSARIO per GetTrack
#include "TGeoTrack.h"
#include "TPolyLine3D.h"
using namespace std;

/*int main(int argc = 0, char** argv = nullptr) {

    TApplication app("app", &argc, argv);

    //inizializzo ambiente
    TGeoManager *geom = new TGeoManager("Hierarchy", "Esempio Tracker");
    
    // definizione materiali
    TGeoMaterial *mat = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMedium *med = new TGeoMedium("Vacuum", 1, mat);

    // TOP Volume
    TGeoVolume *top = geom->MakeBox("TOP", med, 100, 100, 100);
    geom->SetTopVolume(top);
    
    // Geometria (Beam pipe e Layer)
    TGeoVolume *beampipe = geom->MakeTube("beam_pipe", med, beam_pipe_radius, beam_pipe_radius + beam_pipe_thickness, beam_pipe_lenght / 2.);
    beampipe->SetLineColor(kGreen);
    top->AddNode(beampipe, 1);

    //layer1 interno
    TGeoVolume *layer1 = geom->MakeTube("Interno", med, layer1_radius, layer1_radius + layer1_thickness, layer1_lenght / 2.);
    layer1->SetLineColor(kRed);
    top->AddNode(layer1, 1);

    //layer2 esterno
    TGeoVolume *layer2 = geom->MakeTube("Esterno", med, layer2_radius, layer2_radius + 1, layer2_lenght / 2.);
    layer2->SetLineColor(kBlue);
    top->AddNode(layer2, 1);

    geom->CloseGeometry();

    //creazione delle tracce
    int n_punti = 2;
    TPolyLine3D *linea = new TPolyLine3D(n_punti);
    TPolyLine3D *beam = new TPolyLine3D(n_punti);
    
    //SetPoint(indice, x, y, z)
    linea->SetPoint(0, 0, 0, 0);
    linea->SetPoint(1, 50, 50, 50);
    beam->SetPoint(0, 0, 0, -beam_pipe_lenght / 2.);
    beam->SetPoint(1, 0, 0, +beam_pipe_lenght / 2.);

    //grafica
    linea->SetLineColor(kYellow);
    linea->SetLineWidth(4);
    beam->SetLineColor(kBlue);
    beam->SetLineWidth(5);
    
    // Disegno
    // Usiamo un Canvas per assicurarci che il contesto grafico sia inizializzato
    //TCanvas *c1 = new TCanvas("c1", "Viewer", 800, 800);
    top->Draw("ogl");
    
    // Disegna le tracce sulla geometria esistente
    linea->Draw("same");
    beam->Draw("same");


    

    app.Run(); 

    return 0;
    //commento per  edere se fuznione

}*/

int main() {

    tracklet trkl;
    cout << trkl << endl;
    trkl.generate_eta();
    trkl.generate_phi();
    trkl.generate_theta();
    trkl.generate_vertex();
    cout << trkl << endl;



    return 0;
}