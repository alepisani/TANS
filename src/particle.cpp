#include "../include/particle.h"
#include "../include/point.h"
#include "../include/const.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"


using namespace std;

ClassImp(particle);

particle::particle(const point& pt, double p, double eta_in): TObject(), pt(pt), phi(p), eta(eta_in) {

    theta = 2.0 * atan(exp(-eta));

}

void particle::generate_theta(TH1D* hist_eta) {

    /**
     * generate eta and then theta either from uniform or data/kinem.root
     */
    
    double eta_min = -2.;
    double eta_max = +2.;

    if (!get_data_from_kinem) {

        eta = gRandom->Uniform(eta_min, eta_max);

    } 
    else{

        eta = hist_eta->GetRandom();

    }

    theta = 2.0 * atan(exp(-eta));
}

void particle::generate_phi(){

    phi = gRandom->Uniform(0, 2 * M_PI); 

}

void particle::set_theta(double t){

    theta = t;

}

void particle::set_phi(double p){

    phi = p;

}

void particle::set_point(const point& p){

    pt = p;

}

void particle::find_intersection(double radius){

    /**
     * this function find the intersection between the particle trajectory and the detector layers
     * take in input the radius where there is intersection
     * 
     * trajectory of the particle: 
     * x = x0 + t * sin(theta) * cos(phi)
     * y = y0 + t * sin(theta) * sin(phi)
     * z = z0 + t * cos(theta)    
     * 
    */

    double theta = this->get_theta();
    double phi = this->get_phi();
    double c1 = sin(theta) * cos(phi);
    double c2 = sin(theta) * sin(phi);
    double c3 = cos(theta);
    double x0 = this->pt.get_x();
    double y0 = this->pt.get_y();
    double z0 = this->pt.get_z();

    double delta = (x0 * c1 + y0 * c2) * (x0 * c1 + y0 * c2) - (c1 * c1 + c2 * c2) * (x0 * x0 + y0 * y0 - radius * radius);
    double t = (-(x0 * c1 + y0 * c2)+sqrt(delta))/(c1 * c1 + c2 * c2);

    double x_ext = x0 + c1 * t;
    double y_ext = y0 + c2 * t;
    double z_ext = z0 + c3 * t;
    class point pt(x_ext, y_ext, z_ext);

    this->set_point(pt);

}

void particle::rotate(double theta_p, double phi_p) {

    /**
     * rotate the particle's trajectory by a given angle (after multiple scattering)
     */

    //rotation matrix
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

    //direction cosines of the particle's trajectory, before and after the rotation
    double cd[3]; 
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
    
    //theta and phi in the lab reference frame, after the rotation
    theta=acos(cd[2]);
    phi=atan2(cd[1], cd[0]);

}

void particle::multiple_scattering(int Z, double X0, double thickness) {

    /**
     * multiple scattering of the particle, and consequent rotation of the trajectory
     */
    
    double p = 700; // MeV/c 
    double beta = 1.;
    double path = thickness/sin(theta);
    double theta0 = 13.6/(beta*p)*Z*sqrt(path/X0)*(1+0.038*log(path/X0)); //theta in plane
    double theta_rms = theta0*sqrt(2); //rms in 3D space
    double theta_ms = gRandom->Gaus(0, theta_rms);
    double phi_ms = gRandom->Uniform(0, 2 * M_PI);

    this->rotate(theta_ms, phi_ms);

}

std::ostream &operator<<(std::ostream &output, const particle &prtl){
    output << "+------------------------------particle" << endl;
    output << "| Point (mm): " << prtl.pt << endl;
    output << "| Theta (rad): " << prtl.theta << endl;
    output << "| Phi (rad): " << prtl.phi << endl;
    output << "+--------------------------------------" << endl;
    return output;
}

