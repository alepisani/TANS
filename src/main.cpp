#include "../include/const.h"
#include <iostream>
#include "TGeoManager.h"   // Il "regista" della geometria
#include "TGeoVolume.h"    // Per la gestione dei volumi 
#include "TGeoMaterial.h"  // Per definire materiali e medie (indici di rifrazione, X0, ecc.)
#include "TGeoTube.h"      // Per la forma specifica del cilindro (TGeoShape)
#include "TColor.h"        // Se vuoi usare costanti di colore come kBlue o kRed
using namespace std;

int main(){

    // 1. Inizializzazione del Manager e del volume contenitore (Mondo)
    TGeoManager *geom = new TGeoManager("Hierarchy", "Esempio Intercapedine");
    TGeoMedium *med = 0; // Usiamo il medium nullo per semplicità di disegno
    TGeoVolume *top = geom->MakeBox("TOP", med, 50, 50, 50);
    geom->SetTopVolume(top);

    // BEAM PIPE
    TGeoVolume *beampipe = geom->MakeTube("beam pipe", med, beam_pipe_radius, beam_pipe_radius + beam_pipe_thickness, beam_pipe_lenght / 2);
    beampipe->SetLineColor(kGreen);
    beampipe->SetFillColor(kGreen);

    // 2. Creazione del Cilindro Esterno (R=20, dz=30)
    TGeoVolume *extCyl = geom->MakeTube("Esterno", med, layer2_radius, layer2_radius + 1, layer2_lenght / 2);
    extCyl->SetLineColor(kBlue);
    extCyl->SetFillColor(kBlue);
    //extCyl->SetTransparency(50); 

    // 3. Creazione del Cilindro Interno (R=10, dz=30)
    TGeoVolume *intCyl = geom->MakeTube("Interno", med, layer1_radius, layer1_radius + layer1_thickness, layer2_lenght / 2);
    intCyl->SetLineColor(kRed);
    intCyl->SetFillColor(kRed);

    // 4. Costruzione della gerarchia
    // Posizioniamo il cilindro esterno nel Top, e quello interno dentro l'esterno
    top->AddNode(beampipe, 1);
    top->AddNode(extCyl, 1);
    top->AddNode(intCyl, 1); 

    // 5. Chiusura e Disegno
    geom->CloseGeometry();
    top->Draw("ogl");

    return 0;
}