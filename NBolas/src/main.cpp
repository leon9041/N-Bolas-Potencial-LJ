/**
 * @file     main.cpp
 * @brief    simulacion de gas con potencial de Lennard-Jones y Velocity-Verlet
 * @author   Angie Gomez, Leonardo Tovar
 * @date     29/10/25
 * @version  1.0
 * @license  owner
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <cstdlib>
#include <cmath>
#include "Caja.h"
#include <clocale>

int main() {
    std::setlocale(LC_ALL, "C");
    std::filesystem::create_directories("results");

    int N;
    double epsilon, sigma, vmax;

    std::cout << ">>> Simulacion Lennard-Jones (Fluido Infinito) <<<\n";
    std::cout << "1. Ingrese N (ej. 100): "; std::cin >> N;
    std::cout << "2. Ingrese Epsilon (ej. 1.0): "; std::cin >> epsilon;
    std::cout << "3. Ingrese Sigma (ej. 0.1): "; std::cin >> sigma;
    std::cout << "4. Ingrese vmax (ej. 0.01): "; std::cin >> vmax;
    // --- Parámetros ---
    double W = 1.0, H = 1.0;            
    double dt = 0.000005;         
    double T_total = 1.0;     
    
    double r_cut = 3.0 * sigma; 
    
    // Inicialización
    Caja caja(W, H);
    caja.setParametrosLJ(sigma, epsilon, r_cut); 
    
    // USAR NUEVA INICIALIZACIÓN EN GRILLA
    caja.inicializarGrilla(N, vmax);
    
    // Calcular fuerzas iniciales
    caja.calcularFuerzas();

    std::ofstream out("results/salida.dat");
    out << std::fixed << std::setprecision(6);

    // Archivo de datos (Ahora incluye PRESION)
    std::ofstream datFile("results/datos_graficas.dat");
    datFile << "tiempo E_cin E_pot E_tot Presion\n";

    double t = 0.0;
    int steps = static_cast<int>(T_total / dt);

    std::cout << "\nIniciando (" << steps << " pasos)...\n";

    for (int s = 0; s < steps; ++s) {
        if (s % 50 == 0) caja.guardarEstado(out, t);
        
        caja.pasoTemporal(dt);
        t += dt;

        if (s % 50 == 0) {
            double E_tot = caja.energiaTotal();
            double E_pot = caja.energiaPotencial();
            double E_cin = E_tot - E_pot;
            
            // CALCULAR PRESIÓN
            double Presion = caja.obtenerPresion(W * H);

            datFile << t << " " << E_cin << " " << E_pot << " " << E_tot << " " << Presion << "\n";
            
            if (s % 2000 == 0) {
                std::cout << "Paso " << s << "/" << steps 
                          << " P: " << std::setprecision(3) << Presion 
                          << " Etot: " << E_tot << "\n";
            }
        }
    }

    out.close();
    datFile.close();

    std::cout << "\nGenerando graficas...\n";
    int ret = std::system("python scripts/graficas_animacion.py");
    if (ret != 0) std::system("python3 scripts/graficas_animacion.py");

    return 0;

}
