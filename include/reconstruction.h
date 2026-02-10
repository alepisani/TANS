#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H
#include <iostream>
using namespace std;
#include "TObject.h"


class reconstruction : public TObject {

    public:

    //void reco_z(event);
    double reco_z(event);


    /**
     * combinatoriale su tuttle tracklet
     * calcolo z
     * fill hist
     * find max
     * 
     * recTolerance:
     * deltaPhi:   0.01    # [rad] maximum angular difference between two points to be used in reconstruction
     * zBinWidth:  0.5     # [cm]  bin width used to reconstruct z coordinate of the vertex
     * 
     * 
     */
    


    private:

    ClassDef(reconstruction,1);

};



#endif