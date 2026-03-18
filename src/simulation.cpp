#include "../include/simulation.h"
#include "../include/const.h"
#include "../include/point.h"
#include "../include/event.h"
#include "TRandom3.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

simulation::simulation():TObject(){}

void simulation::sim(){
    
    TH1D* hist_eta = nullptr;
    TH1I* hist_mult = nullptr;
    TFile* hist_kinem = nullptr;
    if(get_data_from_kinem) {
        hist_kinem = TFile::Open("../data/kinem.root");
        if(hist_kinem && !hist_kinem->IsZombie()) {
            hist_eta = (TH1D*)hist_kinem->Get("heta2");
            hist_mult = (TH1I*)hist_kinem->Get("hm");
            if(!hist_eta || !hist_mult) {
                std::cerr << "Warning: histograms not found in the given file\n";
            }
        } else {
            std::cerr << "Warning: could not open the given file\n";
        }
    }

    TFile hfile("../data/simulation.root","RECREATE");

    /**
     * this is another algorithm to store data.
     * it is way faster (on 1000000 events is 3 times faster)
     * as downside we have larger .root file
     * compressionlevel = 1, faster and doesnt care about file sizing
     * compressionlevel = 9, slower but optimizes the file sizing
     * 
     * atm i prefer 1 to maximise time in the simulation
     * (will be faster in reading)
     */
    hfile.SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kLZ4);
    hfile.SetCompressionLevel(1);

    TTree *tree = new TTree("MC","simulation");
    
    TClonesArray *ptrVTX = new TClonesArray("point",100);
    TClonesArray *ptrhitsL1 = new TClonesArray("point",100);
    TClonesArray *ptrhitsL2 = new TClonesArray("point",100);
    TClonesArray &VTX = *ptrVTX;
    TClonesArray &hitsL1 = *ptrhitsL1;
    TClonesArray &hitsL2 = *ptrhitsL2;
    int multiplicity;
    
    /**
     * bigger buffer to make it faster
     * not a problem in this case, few branches. with larger branches this 
     * lead to ram overusage
     */
    tree->Branch("Vertex",&ptrVTX, 256000);
    tree->Branch("HitsL1",&ptrhitsL1, 256000);
    tree->Branch("HitsL2",&ptrhitsL2, 256000);
    tree->Branch("Multiplicity", &multiplicity, "Multiplicity/I");
    
    event ev;
    
    cout << "simulation in progress: " << endl;

    for(int i = 0; i < nEvents; i++){
        
        ev.reset();
        ev.single_event(ptrVTX, ptrhitsL1, ptrhitsL2, hist_mult, hist_eta);
        multiplicity = ev.get_multiplicity();
        
        tree->Fill();
        
        ptrVTX->Clear();
        ptrhitsL1->Clear();
        ptrhitsL2->Clear();

        //loading bar
        if (i % 1000 == 0 || i == nEvents - 1) {
            printProgressBar(i + 1, nEvents, 30);
        }

    }
    
    //hfile.Write();
    tree->Write("", TObject::kOverwrite);
    hfile.Close();

}

void simulation::printProgressBar(int current, int total, int barWidth = 30) {
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    std::ostringstream oss;
    oss << "\r[";

    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            oss << "█";
        else
            oss << " ";
    }

    oss << "] ";
    oss << std::setw(3) << int(progress * 100.0) << "% ";
    oss << "(events=" << current << "/" << total << ")";

    std::string output = oss.str();
    size_t terminal_width = 80;
    if (output.size() < terminal_width)
        output += std::string(terminal_width - output.size(), ' ');

    std::cout << output << std::flush;
}