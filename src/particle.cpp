#include "../include/particle.h"
#include "../include/point.h"
#include "../include/const.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
using namespace std;

ClassImp(particle);

particle::particle(point pt, double t, double p):TObject(), pt(pt), theta(t), phi(p){

    eta = -log(tan(theta/2));

}

void particle::generate_theta(){

    // Generazione con accettanza geometrica del detector
    // considerando la posizione asimmetrica del vertice e margini di sicurezza
    const double z_det_max = 135.0;      // mm (metà della lunghezza 270mm)
    const double r_det_max = 40.0;       // mm (raggio layer1)
    const double safety_factor = 0.95;   // Fattore di sicurezza per margine conservativo
    
    // Posizione attuale della particella (dal vertice)
    double x_vtx = pt.get_x();
    double y_vtx = pt.get_y();
    double z_vtx = pt.get_z();
    
    // Raggio attuale del vertice nel piano trasversale
    double r_vtx = sqrt(x_vtx*x_vtx + y_vtx*y_vtx);
    
    // Calcolo distanze asimmetriche verso +z e -z
    double z_forward = z_det_max - z_vtx;   // Distanza verso +z (forward)
    double z_backward = z_det_max + z_vtx;  // Distanza verso -z (backward)
    
    // Raggio disponibile dal vertice (con fattore di sicurezza)
    double r_available = (r_det_max - r_vtx) * safety_factor;
    if (r_available < 0.5) r_available = 0.5;  // Protezione minima
    
    // Angoli massimi di accettanza nelle due direzioni
    double theta_fwd = atan(r_available / z_forward);       // Angolo massimo verso +z
    double theta_bwd = M_PI - atan(r_available / z_backward); // Angolo massimo verso -z
    
    // Conversione a pseudorapidità
    double eta_max = -log(tan(theta_fwd / 2.0));      // Limite massimo (forward)
    double eta_min = -log(tan(theta_bwd / 2.0));      // Limite minimo (backward)
    
    // Genera eta uniformemente entro i limiti asimmetrici
    eta = gRandom->Uniform(eta_min, eta_max);
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

void particle::set_point(point p){

    pt = p;

}

void particle::find_intersection(double radius){

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
    class point pt;
    pt.set_point(x_ext, y_ext, z_ext);

    this->set_point(pt);

}

void particle::rotate(double theta_p, double phi_p) {
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

void particle::multiple_scattering(int Z, double X0, double thickness) {
    double p=700; //MeV/c 
    double beta=1.;
    double theta0=13.6/(beta*p)*Z*sqrt(thickness/(sin(theta)*X0))*(1+0.038*log(thickness/X0)); //theta nel piano
    double theta_rms=theta0*sqrt(2); //rms nello spazio
    double theta_ms = gRandom->Gaus(0, theta_rms);
    double phi_ms = gRandom->Uniform(0, 2*M_PI);

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
