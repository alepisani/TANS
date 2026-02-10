#include "../include/particle.h"
#include "../include/point.h"
#include "../include/const.h"
#include "../include/event.h"
#include "../include/reconstruction.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace std;

ClassImp(reconstruction);


/*void reconstruction::reco_z(event ev){
    
TH1D *hist_z_vtx = new TH1D("hist_reco", "Ricostruzione Vertice Z; z_{cand} [mm]; Conteggi", 200, -beam_pipe_lenght, beam_pipe_lenght);

    double z1, z2, r1, r2, phi1, phi2, z_cand;

    cout << ev.points_L1.size() << endl;
    for(size_t i = 0; i < ev.points_L1.size(); ++i) {
        
        z1 = ev.points_L1[i].get_z();
        r1 = layer1_radius; // Usiamo il raggio nominale per stabilità
        phi1 = ev.points_L1[i].get_phi();

        for(size_t j = 0; j < ev.points_L2.size(); ++j) {
            
            z2 = ev.points_L2[j].get_z();
            r2 = layer2_radius;
            phi2 = ev.points_L2[j].get_phi();

            // Calcolo del Delta Phi con correzione periodicità
            double dphi = std::abs(phi1 - phi2);
            if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;

            // Taglio sui "Tracklets" (combinatoria)
            if(dphi <= delta_phi) {
                // Formula corretta dell'intersezione con l'asse r=0
                z_cand = (r2 * z1 - r1 * z2) / (r2 - r1);
                hist_z_vtx->Fill(z_cand);
            }
        }
}

// 4. Gestione Grafica
TCanvas *c_reco = new TCanvas("c_reco", "Canvas Ricostruzione", 800, 600);
hist_z_vtx->SetLineColor(kRed);
hist_z_vtx->Draw();

// Fondamentale per vedere il risultato in un'applicazione compilata:
c_reco->Update();
c_reco->Modified();

}*/

double reconstruction::reco_z(event ev){
    
TH1D *hist_z_vtx = new TH1D("hist_reco", "Ricostruzione Vertice Z; z_{cand} [mm]; Conteggi", 200, -beam_pipe_lenght, beam_pipe_lenght);

    double z1, z2, r1, r2, phi1, phi2, z_cand;

    cout << ev.points_L1.size() << endl;
    for(size_t i = 0; i < ev.points_L1.size(); ++i) {
        
        z1 = ev.points_L1[i].get_z();
        r1 = layer1_radius; // Usiamo il raggio nominale per stabilità
        phi1 = ev.points_L1[i].get_phi();

        for(size_t j = 0; j < ev.points_L2.size(); ++j) {
            
            z2 = ev.points_L2[j].get_z();
            r2 = layer2_radius;
            phi2 = ev.points_L2[j].get_phi();

            // Calcolo del Delta Phi con correzione periodicità
            double dphi = std::abs(phi1 - phi2);
            if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;

            // Taglio sui "Tracklets" (combinatoria)
            if(dphi <= delta_phi) {
                // Formula corretta dell'intersezione con l'asse r=0
                z_cand = (r2 * z1 - r1 * z2) / (r2 - r1);
                hist_z_vtx->Fill(z_cand);
            }
        }
}

// 4. Gestione Grafica
TCanvas *c_reco = new TCanvas("c_reco", "Canvas Ricostruzione", 800, 600);
hist_z_vtx->SetLineColor(kRed);
hist_z_vtx->Draw();

// Fondamentale per vedere il risultato in un'applicazione compilata:
c_reco->Update();
c_reco->Modified();

double z_reco = hist_z_vtx->GetBinCenter(hist_z_vtx->GetMaximumBin()); //per ottenere la moda dell'istogramma
return z_reco;

}




/*
void reconstruction::reco_z(event& ev) {
    
// 1. Setup dell'istogramma
// Usiamo un nome univoco o lo cancelliamo se esiste già per evitare memory leak in loop
TH1D *hist_z_vtx = new TH1D("hist_reco", "Ricostruzione Vertice Z; z_{cand} [cm]; Conteggi", 
200, -beam_pipe_lenght, beam_pipe_lenght);

// 2. Variabili di appoggio
double z1, z2, r1, r2, phi1, phi2, z_cand;

// 3. Accesso ai membri dell'evento passato per riferimento
// Nota: assumo che points_L1 e points_L2 siano public o abbiano getter
const auto& hits1 = ev.points_L1; 
const auto& hits2 = ev.points_L2;

cout << hits1.size() << endl;
for(size_t i = 0; i < hits1.size(); ++i) {
    
z1 = hits1[i].get_z();
r1 = layer1_radius; // Usiamo il raggio nominale per stabilità
phi1 = hits1[i].get_phi();

for(size_t j = 0; j < hits2.size(); ++j) {
    
z2 = hits2[j].get_z();
r2 = layer2_radius;
phi2 = hits2[j].get_phi();

// Calcolo del Delta Phi con correzione periodicità [0, pi]
double dphi = std::abs(phi1 - phi2);
if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;

// Taglio sui "Tracklets" (combinatoria)
if(dphi <= delta_phi) {
    // Formula corretta dell'intersezione con l'asse r=0
    z_cand = (r2 * z1 - r1 * z2) / (r2 - r1);
    hist_z_vtx->Fill(z_cand);
}
}
}

// 4. Gestione Grafica
TCanvas *c_reco = new TCanvas("c_reco", "Canvas Ricostruzione", 800, 600);
hist_z_vtx->SetLineColor(kRed);
hist_z_vtx->Draw();

// Fondamentale per vedere il risultato in un'applicazione compilata:
c_reco->Update();
c_reco->Modified();

}

*/


