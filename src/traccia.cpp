#include "../include/traccia.h"
#include <iostream>
#include <cmath>
using namespace std;

void traccia::setTheta(){
    //distribuzione assegnata??
    do {
        theta = rand.Uniform(0,M_PI); 
    } while (abs(-log(tan(theta/2))) < 1.);
}

void traccia::setPhi(){
    phi=rand.Uniform(0, 2*M_PI); //in teoria sarebbe distribuzione assegnata?
}

void traccia::setEta(){
    eta=-log(tan(theta/2));
}

std::ostream &operator<<(std::ostream &output, const traccia &track){
    output << "Theta: " << track.theta << endl;
    output << "Phi: " << track.phi << endl;
    output << "Eta: " << track.eta << endl;
    return output;
}