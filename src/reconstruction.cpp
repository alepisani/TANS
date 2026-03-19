#include "../include/particle.h"
#include "../include/point.h"
#include "../include/const.h"
#include "../include/event.h"
#include "../include/reconstruction.h"
#include "../include/simulation.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

ClassImp(reconstruction);

void reconstruction::reset(){

    HitL1.clear();
    HitL2.clear();
    z_cand.clear();

}

double reconstruction::running_window(){

    int maxfrequency = 1;
    double vertex = 999.; //initialised to a dummy value
    int size_zcand = this->get_zcand().size();

    /** 
     * if only one tracklet is recostructed we are goint to take
     * the intersection of the tracklet with the z-axis as reconstructed vertex
     */
    //if(size_zcand == 1) return this->get_zcand()[0];

    //creation of the window to iterate on each value
    for(int i = 0; i < size_zcand; i++){
        
        double window_in = z_cand[i] - half_window; 
        double window_out = z_cand[i] + half_window;
        
        int frequency = 0;        
        double running_sum = 0.;
        
        //check the window frequency 
        for(int j = 0; j < size_zcand; j++){ 

            if(z_cand[j] >= window_in && z_cand[j] <= window_out){
                
                frequency++;
                running_sum += z_cand[j];

            }

        }

        //find the maximum and consider the vertex as the mean of the window
        if(frequency > maxfrequency){
            
            maxfrequency = frequency;
            vertex = running_sum / frequency;  

        }
    }
    
    return vertex;
    
}

void reconstruction::reco(){

    TFile input_file("../data/simulation.root");

    simulation sim;
    
    TTree *tree_input = (TTree*)input_file.Get("MC");
    
    TFile output_file("../data/reconstruction.root", "RECREATE");
    output_file.SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kLZ4);
    output_file.SetCompressionLevel(1);
    
    TTree *tree_output = new TTree("Reco", "reconstruction");
    
    TClonesArray *hitsL1 = nullptr;
    TClonesArray *hitsL2 = nullptr;
    double z_rec;
    tree_input->SetBranchAddress("HitsL1", &hitsL1);
    tree_input->SetBranchAddress("HitsL2", &hitsL2);
    tree_output->Branch("Z_Reco", &z_rec, "Z_Reco/D");

    cout << "reconstruction in progress: " << endl;

    for(int nEvent = 0; nEvent < tree_input->GetEntries() ; nEvent++){
        
        tree_input->GetEvent(nEvent);
        
        this->reset();
        
        for(int i = 0; i < hitsL1->GetEntries(); i++){

            point* p1 = (point*)hitsL1->At(i);
            HitL1.push_back(*p1);

        }

        for(int i = 0; i < hitsL2->GetEntries(); i++){

            point* p2 = (point*)hitsL2->At(i);
            HitL2.push_back(*p2);

        }

        double z1, z2, r1, r2, phi1, phi2;
        int size_L1 = HitL1.size();
        int size_L2 = HitL2.size();

        //for each event create all the z candidates
        for(int i = 0; i < size_L1; i++){

            z1 = HitL1[i].get_z();
            r1 = HitL1[i].get_R();
            phi1 = HitL1[i].get_phi();

            for(int j = 0; j < size_L2; j++){

                z2 = HitL2[j].get_z();
                r2 = HitL2[j].get_R();
                phi2 = HitL2[j].get_phi();

                //taking into accout the periodicity
                double dphi = std::abs(phi1 - phi2);  
                if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;
                
                //selection on phi and geometry
                if(dphi <= delta_phi) {

                    z_rec = (r2 * z1 - r1 * z2) / (r2 - r1); 

                    z_cand.push_back(z_rec);
                
                }

            }

        }

        //sort on the z_cand vector
        std::sort(z_cand.begin(), z_cand.end());
        z_rec = this->running_window();

        //if(!(z_rec == 999)) 
        tree_output->Fill();

        //loading bar
        if (nEvent % 1000 == 0 || nEvent == nEvents - 1) {
            sim.printProgressBar(nEvent + 1, nEvents, 30);
        }

    }
    
    tree_output->Write();
    output_file.Close();
    input_file.Close();
    
}

void reconstruction::analysis(){

    TFile simulation_file("../data/simulation.root");
    TFile reconstruction_file("../data/reconstruction.root");
    
    TFile output_file("../data/analysis.root", "RECREATE");
    output_file.SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kLZ4);
    output_file.SetCompressionLevel(1);

    TH1D* hResidui = new TH1D("hResidui", "Residui; (z_{gen} - z_{reco}) [#mum]; Conteggi", 500, -1000, 1000);
    TH2D* hResVsMult = new TH2D("hResVsMult", "Residui vs Mult; Molteplicita'; Residuo [#mum]", 50, 1, 50, 200, -1000, 1000);
    TH2D* hResVsZtrue = new TH2D("hResVsZtrue", "Residui vs Ztrue; z [mm]; Residuo [#mum]", 50, -150, +150, 200, -1000, 1000);

    TH1F* hEvents_z = new TH1F("hEvents_z", "Events vs Z_{true}; z [mm]; Counts", 50, -200, +200);
    TH1F* hSuccess_z = new TH1F("hSuccess_z", "Successful reconstructions vs Z_{true}; z [mm]; Counts", 50, -200, +200);

    TH1F* hEvents_mult = new TH1F("hEvents_mult", "Events vs Multiplicity; Multiplicity; Counts", 50, 0, 50);
    TH1F* hSuccess_mult = new TH1F("hSuccess_mult", "Successful reconstructions vs Multiplicity; Multiplicity; Counts", 50, 0, 50);

    // Get the trees
    TTree *tree_simulation = (TTree*)simulation_file.Get("MC");
    TTree *tree_reconstruction = (TTree*)reconstruction_file.Get("Reco");
    
    double z_reco, z_real, residual;
    TClonesArray *vertex = nullptr;
    int multiplicity;
    
    // Set branch addresses
    tree_reconstruction->SetBranchAddress("Z_Reco", &z_reco);
    tree_simulation->SetBranchAddress("Vertex", &vertex);
    tree_simulation->SetBranchAddress("Multiplicity", &multiplicity);
    
    int nEvents = tree_simulation->GetEntries();

    if(nEvents != tree_reconstruction->GetEntries()){

        std::cerr << "Warning, entries on simulation and reconstruction does not correspond \n";

    }

    cout << endl;
    cout << "analysis in progress: " << endl;

    for(int nEvent = 0; nEvent < nEvents; nEvent++){

        tree_simulation->GetEvent(nEvent);
        tree_reconstruction->GetEvent(nEvent);
        
        // Extract the vertex z coordinate
        point* vtx = (point*)vertex->At(0);
        z_real = vtx->get_z();
        
        // Fill efficiency histograms
        hEvents_z->Fill(z_real);
        hEvents_mult->Fill(multiplicity);
        
        // Calculate residual
        if(z_reco != 999){

            residual = (z_real - z_reco) * 1000; //um
            hResidui->Fill(residual);
            hResVsMult->Fill(multiplicity, residual);
            hResVsZtrue->Fill(z_real, residual);
            hSuccess_z->Fill(z_real);
            hSuccess_mult->Fill(multiplicity);

        }
        
        //loading bar
        if (nEvent % 1000 == 0 || nEvent == nEvents - 1) {
            simulation sim;
            sim.printProgressBar(nEvent + 1, nEvents, 30);
        }

    }

    // ── helper function to find the stdev of the distribution ─────────────────────
    auto fitAndGetSigma = [](TH1D* proj, TH1D* hOut, int bin) {
        if(proj->GetEntries() <= 10) return;

        double mean  = proj->GetMean();
        double sigma = proj->GetStdDev();
        if(sigma <= 0) return;

        // gaussian fit in the range (mean +- 2sigma)
        TF1* gaus = new TF1("gaus_fit", "gaus", mean - 2*sigma, mean + 2*sigma);
        int fitStatus = proj->Fit(gaus, "RQ0"); // R=range, Q=quiet, 0=no draw

        if(fitStatus == 0 && gaus->GetParameter(2) > 0){
            hOut->SetBinContent(bin, gaus->GetParameter(2)); 
            hOut->SetBinError(bin,   gaus->GetParError(2));  
        }

        delete gaus;
    };

    // ── Resolution vs Multiplicity ──────────────────────────────────────────────
    int nbins_mult = hResVsMult->GetNbinsX();
    TH1D* hResVsMult_1D = new TH1D("hResVsMult_1D",
        "Resolution vs Multiplicity; Multiplicity; #sigma [#mum]",
        nbins_mult, 0, 50);

    for(int i = 1; i <= nbins_mult; i++){
        TH1D* proj = hResVsMult->ProjectionY(Form("proj_mult_%d", i), i, i);
        fitAndGetSigma(proj, hResVsMult_1D, i);
        delete proj;
    }

    // ── Resolution vs Z ─────────────────────────────────────────────────────────
    int nbins_z = hResVsZtrue->GetNbinsX();
    TH1D* hResVsZTrue_1D = new TH1D("hResVsZTrue_1D",
        "Resolution vs Z_{true}; z [mm]; #sigma [#mum]",
        nbins_z, -150, 150);

    for(int i = 1; i <= nbins_z; i++){
        TH1D* proj = hResVsZtrue->ProjectionY(Form("proj_z_%d", i), i, i);
        fitAndGetSigma(proj, hResVsZTrue_1D, i);
        delete proj;
    }
    
    cout << endl;
    
    // Create TEfficiency objects
    TEfficiency* effZ = new TEfficiency(*hSuccess_z, *hEvents_z);
    effZ->SetName("Eff_vs_z");
    effZ->SetTitle("Efficiency vs Z_{true}; z [mm]; Efficiency");
    
    TEfficiency* effMult = new TEfficiency(*hSuccess_mult, *hEvents_mult);
    effMult->SetName("Eff_vs_Mult");
    effMult->SetTitle("Efficiency vs Multiplicity; Multiplicity; Efficiency");
    
    // Write histogram to output file
    hResidui->Write();
    hResVsMult->Write();
    hResVsZtrue->Write();
    effZ->Write();
    effMult->Write();
    hResVsMult_1D->Write();
    hResVsZTrue_1D->Write();

    output_file.Close();
    simulation_file.Close();
    reconstruction_file.Close();
    
}