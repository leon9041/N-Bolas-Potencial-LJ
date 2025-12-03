#ifndef CAJA_H
#define CAJA_H

#include <vector>
#include <fstream>
#include "Bola.h"

class Caja {
public:
    double W, H;
    std::vector<Bola> bolas;

    // Parámetros LJ
    double sigma;
    double epsilon;
    double r_cut;       
    double U_rcut;      

    // Variables para Presión (Virial)
    double virial_total; // Acumulador del virial para el paso actual

    Caja(double W_, double H_);
    void setParametrosLJ(double s, double e, double rc);
    
    // Cambiamos a inicialización en grilla
    void inicializarGrilla(int N, double vmax);
    void pasoTemporal(double dt);
    void calcularFuerzas();
    void guardarEstado(std::ofstream& out, double t) const;
    
    double energiaPotencial() const;
    double energiaTotal() const;
    
    // Nuevo método para obtener presión instantánea
    double obtenerPresion(double area) const;
};

#endif