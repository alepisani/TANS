#include "../include/point.h"
#include <cmath>
#include <iostream>
#include <cmath>

using namespace std;


ClassImp(point);

// COSTRUTTORE
point::point() : x(0), y(0), z(0) {}

// DISTRUTTORE
point::~point() {}

void point::set_point(double a, double b, double c){
    x = a;
    y = b;
    z = c;
}

void point::generate_theta(){
    theta = 2*atan(exp(-eta));
}

void point::generate_phi(){
    phi=gRandom->Uniform(0, 2 * M_PI); 
}

void point::generate_eta(){
    //eta=-log(tan(theta/2));
    eta = gRandom->Uniform(-1,1); //mi serve la distribuzione assegnata in eta

}

void point::set_theta(double t){
    theta=t;
}
void point::set_phi(double p){
    phi = p;
}

void point::set_points(point point0, point point1){

    point_int = point0;
    point_ext = point1;

}

void point::set_point_ext(point point){

    point_ext = point;

}

void point::set_point_int(point point){

    point_int = point;

}

void point::set_cilindrical(){

    R = sqrt(x*x + y*y);
    phi = atan2(y, x);

}

void point::update_coordinates(){

    x = R * cos(phi);
    y = R * sin(phi);

}

point point::extend_segment(double theta, double phi, double lenght){

    double c1 = sin(theta) * cos(phi);
    double c2 = sin(theta) * sin(phi);
    double c3 = cos(theta);
    double x0 = this->get_x();
    double y0 = this->get_y();
    double z0 = this->get_z();

    double delta = (x0 * c1 + y0 * c2) * (x0 * c1 + y0 * c2) - (c1 * c1 + c2 * c2) * (x0 * x0 + y0 * y0 - lenght * lenght);
    double t = (-(x0 * c1 + y0 * c2)+sqrt(delta))/(c1 * c1 + c2 * c2);

    double x = x0 + c1 * t;
    double y = y0 + c2 * t;
    double z = z0 + c3 * t;

    this->set_point(x, y, z);
    return *this;

}

double point::multiple_scattering(int Z, double X0, double thickness) {
    double p=700; //MeV/c 
    double beta=1.;
    double theta0=13.6/(beta*p)*Z*sqrt(thickness/(sin(theta)*X0))*(1+0.038*log(thickness/X0)); //theta nel piano
    double theta_rms=theta0*sqrt(2); //rms nello spazio
    double theta_ms=gRandom->Gaus(0, theta_rms);
    return theta_ms;

}


void point::rotate (double theta_p, double phi_p) {
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

std::ostream &operator<<(std::ostream &output, const point &point){
    output << "point coordinate (mm): (" << point.x << ", " 
                                    << point.y << ", " 
                                    << point.z << ")" << endl;
    output << "Point_int (mm): " << point.point_int << endl;
    output << "Point_ext (mm): " << point.point_ext << endl;
    output << "Theta (rad): " << point.theta << endl;
    output << "Phi (rad): " << point.phi << endl;
    output << "Eta: " << point.eta << endl;
    return output;
}