#include "../include/evento.h"
#include <iostream>
using namespace std;

//evento::evento() : TObject(), x(0), y(0), z(0), molteplicita(0), rand() {
    //costruttore di default
//}

//evento::evento() {}

//evento::evento(double x, double y, double z, int molteplicita) : TObject(), x(x), y(y), z(z), molteplicita(molteplicita), rand() {
    //costruttore con parametri
//}



void evento::setMolteplicita() {
    molteplicita = static_cast<int>(rand.Uniform(1, 50));
    //manca distribuzione assegnata
    //molteplicita = 10;
}

std::ostream &operator<<(std::ostream &output, const evento & ev) {
    output << "x del vertice: " << ev.x << endl;
    output << "y del vertice: " << ev.y << endl;
    output << "z del vertice: " << ev.z << endl;
    output << "molteplicita: " << ev.molteplicita << endl;
    return output;
}