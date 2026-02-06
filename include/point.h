#ifndef POINT_H
#define POINT_H
#include <iostream>
using namespace std;
#include "TObject.h"

class point : public TObject {

    public:
    point();              
    virtual ~point();     
    void set_point(double, double, double);
    void generate_theta();
    void generate_phi();
    void generate_eta();
    double get_theta() const {return theta;}
    double get_phi() const {return phi;}
    double get_eta() const {return eta;}
    void set_theta(double);
    void set_phi(double);

    //auto get_point() {return point;}
    double get_x() const {return x;}
    double get_y() const {return y;}
    double get_z() const {return z;}
    double get_R() const {return R;}
    //double get_phi() const {return phi;}
    void set_cilindrical();
    //void set_phi(double);
    void update_coordinates();
    point extend_segment(double, double, double);

    auto get_point_int() const {return point_int;}
    auto get_point_ext() const {return point_ext;}
    void set_points(point, point);
    void set_point_ext(point);
    void set_point_int(point);
    point find_intersection(double);
    double multiple_scattering(int Z, double X0, double t);
    void rotate (double theta_p, double phi_p);
  
    friend std::ostream &operator<<(std::ostream &output, const point &point);

    private:
        //cordinate planari
        double x;
        double y;
        double z;  

        double theta; 
        double phi; 
        double eta; 
        
        //coordinate cilindriche
        double R;
        // z coincinde con coordinate planari
        //double phi;

    ClassDef(point, 1);
};

#endif