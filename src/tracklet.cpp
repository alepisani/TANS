#include "../include/tracklet.h"
#include "../include/point.h"
#include "../include/const.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
using namespace std;

ClassImp(tracklet);


// COSTRUTTORE
//tracklet::tracklet() : point(), point(), theta(0), eta(0), phi(0) {}
tracklet::tracklet() {}
// DISTRUTTORE
tracklet::~tracklet() {}


void tracklet::generate_theta(){
    theta = 2*atan(exp(-eta));
    /*do {
        theta = gRandom->Uniform(0,M_PI); 
    } while (abs(-log(tan(theta/2))) >= 1.);*/
}

//AL CONTRARIO! GENERO ETA SU UNA FUNZIONE CHE MI DA LUI, CON ACCETTANZA IN BASE A DOVE STA IL VERTICE

void tracklet::generate_phi(){
    phi=gRandom->Uniform(0, 2 * M_PI); 
}

void tracklet::generate_eta(){
    //eta=-log(tan(theta/2));
    eta = gRandom->Uniform(-1,1); //mi serve la distribuzione assegnata in eta

}

void tracklet::set_points(point point0, point point1){

    point_int = point0;
    point_ext = point1;

}

void tracklet::set_point_ext(point point){

    point_ext = point;

}

void tracklet::set_point_int(point point){

    point_int = point;

}

void tracklet::set_theta(double t){

    theta = t;

}

void tracklet::set_phi(double p){

    phi = p;

}

point tracklet::find_intersection(double radius) {
    
    double theta = this->get_theta();
    double phi = this->get_phi();
    double c1 = sin(theta) * cos(phi);
    double c2 = sin(theta) * sin(phi);
    double c3 = cos(theta);
    double x0 = this->point_int.get_x();
    double y0 = this->point_int.get_y();
    double z0 = this->point_int.get_z();

    double delta = (x0 * c1 + y0 * c2) * (x0 * c1 + y0 * c2) - (c1 * c1 + c2 * c2) * (x0 * x0 + y0 * y0 - radius * radius);
    double t = (-(x0 * c1 + y0 * c2)+sqrt(delta))/(c1 * c1 + c2 * c2);

    double x_ext = x0 + c1 * t;
    double y_ext = y0 + c2 * t;
    double z_ext = z0 + c3 * t;

    point p;
    p.set_point(x_ext, y_ext, z_ext);

    //make point_ext
    return p;
}

double tracklet::multiple_scattering(int Z, double X0, double thickness) {
    double p=700; //MeV/c 
    double beta=1.;
    double theta0=13.6/(beta*p)*Z*sqrt(thickness/(sin(theta)*X0))*(1+0.038*log(thickness/X0)); //theta nel piano
    double theta_rms=theta0*sqrt(2); //rms nello spazio
    double theta_ms=gRandom->Gaus(0, theta_rms);
    //return theta_ms;

    double thetaMS = gRandom->Gaus(0, theta_rms);
    double phiMS = 2. * M_PI * gRandom->Rndm();
    double theta = this->get_theta();
    double phi = this->get_phi();

    // versor after MS in reference system of the particle
    double vec[3] = {cos(phiMS)*sin(thetaMS), sin(phiMS)*sin(thetaMS), cos(thetaMS)};

    //rotate the vector in the lab system
    double rotMat[3][3] =  {{-sin(phi), -cos(phi)*cos(theta),   cos(phi)*sin(theta)},
                            {cos(phi),  -sin(phi)*cos(theta),   sin(phi)*sin(theta)},
                            {0.,        sin(theta),             cos(theta)}};
    
    double vecp[3];
    for(int i=0; i<3; i++)  vecp[i] = vec[i];

    for(int i=0; i<3; i++)
    {
        vec[i] = 0;
        for(int j=0; j<3; j++)  vec[i] += rotMat[i][j] * vecp[j];
    }
    
    double newPhi = atan(vec[1]/vec[0]);
    double newTheta = vec[2];

    return newTheta;


}


void tracklet::rotate (double theta_p, double phi_p) {
    double mr[3][3];

    mr[0][0] = sin(phi);
    mr[1][0] = cos(phi);
    mr[2][0] = 0.;
    mr[0][1] = cos(phi)*cos(theta);
    mr[1][1] = cos(theta)*sin(phi);
    mr[2][1] = sin(theta);
    mr[0][2] = sin(theta)*cos(phi);
    mr[1][2] = sin(theta)*sin(phi);
    mr[2][2] = cos(theta);

    double cd[3]; //coseni direttori
    double cdp[3];  
    cdp[0] = sin(theta_p)*cos(phi_p);
    cdp[1] = sin(theta_p)*sin(phi_p);
    cdp[2] = cos(theta_p);

    for (int i=0; i<3; i++) {
        cd[i]=0.;
        for (int j=0; j<3; j++){
            cd[i]+=mr[i][j]*cdp[j];
        }
    }
    
    theta=acos(cd[2]);
    phi=atan2(cd[1], cd[0]);

}

std::ostream &operator<<(std::ostream &output, const tracklet &trkl){
    output << "Tracklet info: " << endl;
    output << "Point_int (mm): " << trkl.point_int << endl;
    output << "Point_ext (mm): " << trkl.point_ext << endl;
    output << "Theta (rad): " << trkl.theta << endl;
    output << "Phi (rad): " << trkl.phi << endl;
    output << "Eta: " << trkl.eta << endl;
    return output;
}

