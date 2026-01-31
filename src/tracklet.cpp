#include "../include/tracklet.h"
#include "../include/point.h"
#include "../include/const.h"
#include "TRandom3.h"
#include <iostream>
#include <cmath>
using namespace std;

ClassImp(tracklet);

TRandom3 rnd(0);

// COSTRUTTORE
//tracklet::tracklet() : point(), point(), theta(0), eta(0), phi(0) {}
tracklet::tracklet() {}
// DISTRUTTORE
tracklet::~tracklet() {}

void tracklet::generate_theta() {
    do {
        theta = rnd.Uniform(0,M_PI); 
    } while (abs(-log(tan(theta/2))) >= 1.);
}

void tracklet::generate_phi(){
    phi=rnd.Uniform(-M_PI, +M_PI); 
}

void tracklet::generate_eta(){
    eta=-log(tan(theta/2));
}

void tracklet::set_points(point point0, point point1){

    point_int = point0;
    point_ext = point1;

}

void tracklet::find_beampipe_intersection(double x, double y, double z, double theta, double phi) {
    // Vettore direzione della retta
    double dx = sin(theta) * sin(phi);
    double dy = cos(theta);
    double dz = sin(theta) * cos(phi);
    
    // Intersezione con cilindro infinito (asse y)
    // (x + t*dx)² + (z + t*dz)² = R²
    
    double R = beam_pipe_radius;
    
    double a = dx*dx + dz*dz;
    double b = 2.0 * (x*dx + z*dz);
    double c = x*x + z*z - R*R;
    
    double discriminant = b*b - 4*a*c;
    
    if (discriminant < 0) {
        cout << "ERRORE" << endl; // Nessuna intersezione con il cilindro
    }
    
    double t1 = (-b - sqrt(discriminant)) / (2*a);
    double t2 = (-b + sqrt(discriminant)) / (2*a);
    
    // Scegli la prima intersezione positiva
    double t = (t1 > 0) ? t1 : t2;
    
    if (t <= 0) {
        cout << "ERRORE" << endl; // Intersezione dietro al punto di partenza
    }
    
    // Calcola il punto di intersezione
    double x1 = x + t * dx;
    double y1 = y + t * dy;
    double z1 = z + t * dz;
    // Verifica che il punto sia dentro i limiti della beam pipe
    double half_length = beam_pipe_lenght / 2.0;
    
    if (y1 > -half_length || y1 < half_length) {
        
        point_ext.set_point(x1, y1, z1);
        
    }

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

