#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"
#include "./particle.h"
#include "./point.h"

class event : public TObject {

public:
   event():TObject(), multiplicity(0), pnt(), vertex(){}
   event(int, point, point);
   
   int get_multiplicity() const {return multiplicity;}
   const point& get_point() const {return pnt;}
   const point& get_vertex() const {return vertex;}
   
   void setmultiplicity(TH1I* hist_mult = nullptr);
   void set_vertex(const point& vtx);

   void RunFullSimulation();
   void printProgressBar(int, int, int);

   friend std::ostream &operator<<(std::ostream &output, const event &ev);

   std::vector<point> points_BP;        //will contain all points on bp
   std::vector<point> points_L1;        //will contain all points on l1
   std::vector<point> points_L2;        //will contain all points on l2
     

  
private:
   int multiplicity;
   point pnt;
   point vertex;

   ClassDef(event,1);
};





#endif
