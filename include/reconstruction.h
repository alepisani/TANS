#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include <iostream>
using namespace std;
#include "TObject.h"


class reconstruction : public TObject {

    public:

    //void reco_z(event);
    double reco_z(event*, class TH1D*);
    double running_window();

    std::vector <double> z_candidates;  
    std::vector <double> vertex_candidate;

    ClassDef(reconstruction,1);

};



#endif