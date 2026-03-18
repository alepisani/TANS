#ifndef SIMULATION_H
#define SIMULATION_H

#include "TObject.h"


class simulation : public TObject {

public:
    simulation();
    
    void sim();
    void printProgressBar(int, int, int);

   ClassDef(simulation,1);
};


#endif