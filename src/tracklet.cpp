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

void tracklet::generate_theta() {
    do {
        theta = gRandom->Uniform(0,M_PI); 
    } while (abs(-log(tan(theta/2))) >= 1.);
}

void tracklet::generate_phi(){
    phi=gRandom->Uniform(-M_PI, +M_PI); 
}

void tracklet::generate_eta(){
    eta=-log(tan(theta/2));
}

void tracklet::set_points(point point0, point point1){

    point_int = point0;
    point_ext = point1;

}

void tracklet::set_point_ext(point point){

    point_ext = point;

}

void tracklet::find_beampipe_intersection() {
    
    double theta = this->get_theta();
    double phi = this->get_phi();
    // Direzione del raggio (versore)
    double dx = sin(theta) * sin(phi);  // componente x
    double dy = cos(theta);              // componente y (theta=0 -> y massimo)
    double dz = sin(theta) * cos(phi);  // componente z (phi=0 -> z massimo)
    
    // Punto iniziale
    double x0 = this->point_int.get_x();
    double y0 = this->point_int.get_y();
    double z0 = this->point_int.get_z();
    
    // Equazione parametrica della retta: P(t) = P0 + t*d
    // Equazione del cilindro (asse lungo y): x^2 + z^2 = R^2
    // Sostituendo: (x0 + t*dx)^2 + (z0 + t*dz)^2 = R^2
    
    // Coefficienti dell'equazione quadratica: a*t^2 + b*t + c = 0
    double a = dx*dx + dz*dz;
    double b = 2.0 * (x0*dx + z0*dz);
    double c = x0*x0 + z0*z0 - beam_pipe_radius*beam_pipe_radius;
    
    double discriminant = b*b - 4.0*a*c;
    
    if (discriminant < 0) {
        // Nessuna intersezione
        std::cerr << "Nessuna intersezione con il cilindro!" << std::endl;
        //return punto(0, 0, 0);
    }
    
    // Due soluzioni possibili
    double t1 = (-b - sqrt(discriminant)) / (2.0*a);
    double t2 = (-b + sqrt(discriminant)) / (2.0*a);
    
    // Scegli la soluzione positiva più piccola (primo punto di intersezione)
    double t;
    if (t1 > 0 && t2 > 0) {
        t = std::min(t1, t2);
    } else if (t1 > 0) {
        t = t1;
    } else if (t2 > 0) {
        t = t2;
    } else {
        // Entrambe negative: il punto è oltre il cilindro nella direzione opposta
        std::cerr << "Intersezione nella direzione opposta!" << std::endl;
        //return punto(0, 0, 0);
    }
    
    // Calcola il punto di intersezione
    double x_int = x0 + t * dx;
    double y_int = y0 + t * dy;
    double z_int = z0 + t * dz;
    point p;
    p.set_point(x_int, y_int, z_int);

    //make point_ext
    this->point_ext = p;
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

