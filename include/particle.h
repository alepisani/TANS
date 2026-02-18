#ifndef PARTICLE_H
#define PARTICLE_H
#include <iostream>
#include "TObject.h"
#include "point.h"
#include "TH1.h"
using namespace std;

class particle : public TObject {

    public:
    particle():TObject(), pt(), theta(0.), phi(0.), eta(0.){}
    particle(const point&, double, double);
    
    void generate_theta(TH1D* hist_eta = nullptr);
    void generate_phi();
    
    double get_theta() const {return theta;}
    double get_phi() const {return phi;}
    point& get_point() {return pt;}
    
    void set_theta(double);
    void set_phi(double);
    void set_point(const point&);

    void find_intersection(double);
    void multiple_scattering(int Z, double X0, double t);
    void rotate (double theta_p, double phi_p);

    friend std::ostream &operator<<(std::ostream &output, const particle &particle);


    private:
        point pt;
        double theta;
        double phi;
        double eta;  

    ClassDef(particle, 1);
};

#endif