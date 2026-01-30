#include "../include/tracklet.h"
#include "../include/point.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
using namespace std;

ClassImp(tracklet);

TRandom3 rnd(0);

// COSTRUTTORE
tracklet::tracklet() : point(), theta(0), eta(0), phi(0) {}

// DISTRUTTORE
tracklet::~tracklet() {}

void tracklet::generate_theta() {
    do {
        theta = rnd.Uniform(0,M_PI); 
    } while (abs(-log(tan(theta/2))) < 1.);
}

void tracklet::generate_phi(){
    phi=rnd.Uniform(0, 2*M_PI); 
}

void tracklet::generate_eta(){
    eta=-log(tan(theta/2));
}

void tracklet::generate_vertex(){
    double x = rnd.Gaus(0,0.1);
    double y = rnd.Gaus(0,0.1);
    double z = rnd.Gaus(0,53.);
    point.set_point(x,y,z);
}

std::ostream &operator<<(std::ostream &output, const tracklet &trkl){
    output << "Tracklet info: " << endl;
    output << "Point: " << trkl.point << endl;
    output << "Theta: " << trkl.theta << endl;
    output << "Phi: " << trkl.phi << endl;
    output << "Eta: " << trkl.eta << endl;
    return output;
}

