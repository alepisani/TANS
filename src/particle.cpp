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

void particle::generate_theta(TH1D* hist_eta){

    if(!get_data_from_kinem){

        //Generate theta, considering detector's geometric acceptance
        //considering the asymmetric position of the vertex and safety margins
        
        double x_vtx = pt.get_x();
        double y_vtx = pt.get_y();
        double z_vtx = pt.get_z();
        
        // vertex radial distance from the beamline
        double r_vtx = sqrt(x_vtx*x_vtx + y_vtx*y_vtx);
        
        // distances from the vertex to the end of the layers in both directions
        double z_forward = layer1_lenght * 0.5 - z_vtx;   
        double z_backward = layer1_lenght * 0.5 + z_vtx;  
        
        // available radial distance from the vertex
        double r_available = (layer1_radius - r_vtx);
        
        // maximum acceptance angles in both directions
        double theta_fwd = atan(r_available / z_forward);       
        double theta_bwd = M_PI - atan(r_available / z_backward); 
        
        // conversion to pseudorapidity limits
        double eta_max = -log(tan(theta_fwd / 2.0));      
        double eta_min = -log(tan(theta_bwd / 2.0));      
        
        eta = gRandom->Uniform(eta_min, eta_max);
        
        theta = 2.0 * atan(exp(-eta));
    }

    if(get_data_from_kinem){
        
        eta = hist_eta->GetRandom();
        
        theta = 2.0 * atan(exp(-eta));
        
    }

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
    double theta0 = 13.6/(beta*p)*Z*sqrt(thickness/(sin(theta)*X0))*(1+0.038*log(thickness/X0)); //theta in plane
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

