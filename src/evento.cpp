#include "../include/evento.h"
#include "../include/const.h"
#include "TApplication.h" 
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMaterial.h"
#include "TGeoTube.h"
#include "TColor.h"
#include "TVirtualGeoTrack.h" 
#include "TGeoTrack.h"
#include "TPolyLine3D.h"
#include <iostream>
using namespace std;

//evento::evento() : TObject(), x(0), y(0), z(0), molteplicita(0), rand() {
    //costruttore di default
//}

//evento::evento() {}

//evento::evento(double x, double y, double z, int molteplicita) : TObject(), x(x), y(y), z(z), molteplicita(molteplicita), rand() {
    //costruttore con parametri
//}



void evento::setMolteplicita() {
    molteplicita = static_cast<int>(rand.Uniform(1, 50));
    //manca distribuzione assegnata
    //molteplicita = 10;
}

void evento::display_event(){
    
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
    top->Draw("ogl");
    
    // Disegna le tracce sulla geometria esistente
    linea->Draw("same");
    beam->Draw("same");

}







std::ostream &operator<<(std::ostream &output, const evento & ev) {
    output << "x del vertice: " << ev.x << endl;
    output << "y del vertice: " << ev.y << endl;
    output << "z del vertice: " << ev.z << endl;
    output << "molteplicita: " << ev.molteplicita << endl;
    return output;
}