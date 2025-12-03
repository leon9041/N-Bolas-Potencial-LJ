/**
 * @file    caja.cpp
 * @brief   Implementacion de POO para Dinámica Molecular con LJ (PBC + Virial)
 */

#include "Caja.h"
#include <random>
#include <iostream>
#include <cmath>

Caja::Caja(double W_, double H_) : W(W_), H(H_) {
    sigma = 0.0;
    epsilon = 0.0;
    r_cut = 0.0;
    U_rcut = 0.0;
    virial_total = 0.0;
}

// ---------------- Inicializar en Grilla (Evita solapamientos) ----------------
// ---------------- Inicializar en Grilla (Evita solapamientos) ----------------
void Caja::inicializarGrilla(int N, double vmax) {
    bolas.clear();
    
    // Calcular cuántas partículas caben por fila (raíz cuadrada)
    int n_side = std::ceil(std::sqrt(N)); 
    double d_grid = W / n_side; // Distancia entre partículas
    
    // Inicializar velocidades aleatorias
    std::mt19937 rng(45); // Semilla fija
    std::uniform_real_distribution<double> uv(-vmax, vmax);

    int count = 0;
    for (int i = 0; i < n_side; ++i) {
        for (int j = 0; j < n_side; ++j) {
            if (count >= N) break;
            
            // Posición centrada en la celda de la grilla
            double x = (i + 0.5) * d_grid;
            double y = (j + 0.5) * d_grid;
            
            // Radio visual (no afecta física LJ, solo para pintar)
            double r_vis = sigma * 0.5; 
            
            bolas.push_back(Bola(count, x, y, uv(rng), uv(rng), r_vis, 1.0));
            count++;
        }
    }
    std::cout << "Inicializadas " << count << " particulas en grilla.\n";
}

// ---------------- Un paso de simulación (Velocity-Verlet con PBC) ----------------
void Caja::pasoTemporal(double dt) {
    // FASE 1: Paso completo de Posición y medio de Velocidad
    for (auto& b : bolas) {
        b.pos.x += b.vel.x * dt + 0.5 * b.acc_old.x * dt * dt;
        b.pos.y += b.vel.y * dt + 0.5 * b.acc_old.y * dt * dt;
        
        b.vel.x += 0.5 * b.acc_old.x * dt;
        b.vel.y += 0.5 * b.acc_old.y * dt;
        
        // --- APLICAR PBC (Envolver partículas) ---
        if (b.pos.x < 0) b.pos.x += W;
        else if (b.pos.x >= W) b.pos.x -= W;
        
        if (b.pos.y < 0) b.pos.y += H;
        else if (b.pos.y >= H) b.pos.y -= H;
    }

    // FASE 2: Cálculo de Fuerzas (y Virial)
    calcularFuerzas(); 

    // FASE 3: Completar la Velocidad
    for (auto& b : bolas) {
        b.vel.x += 0.5 * b.acc.x * dt;
        b.vel.y += 0.5 * b.acc.y * dt;
    }
}

// ---------------- Calcular Fuerzas LJ + PBC + Virial ----------------
void Caja::calcularFuerzas() {
    // 1. Reiniciar aceleraciones y virial
    virial_total = 0.0;
    for (auto& b : bolas) {
        b.acc_old = b.acc; 
        b.acc = Vec2(0.0, 0.0); 
    }

    double r_cut_sq = r_cut * r_cut;
    double half_W = 0.5 * W;
    double half_H = 0.5 * H;
    double sigma6_base = sigma * sigma * sigma * sigma * sigma * sigma; // sigma^6

    size_t N = bolas.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            Vec2 r_vec = bolas[i].pos - bolas[j].pos;

            // --- Distancia Mínima de Imagen (PBC) ---
            if (r_vec.x > half_W) r_vec.x -= W;
            else if (r_vec.x < -half_W) r_vec.x += W;
            if (r_vec.y > half_H) r_vec.y -= H;
            else if (r_vec.y < -half_H) r_vec.y += H;

            double r2 = r_vec.norm2();
            
            // Evitar división por cero (por si acaso)
            if (r2 < 1e-10) r2 = 1e-10;

            if (r2 < r_cut_sq) {
                // Cálculo optimizado (sin sqrt innecesarias)
                double r2_inv = 1.0 / r2;
                double r6_inv = r2_inv * r2_inv * r2_inv;
                
                double s6_r6 = sigma6_base * r6_inv; // (sigma/r)^6
                
                // F(r)/r = (24*epsilon/r^2) * [ 2*(sigma/r)^12 - (sigma/r)^6 ]
                double factor = (24.0 * epsilon * r2_inv) * (2.0 * (s6_r6 * s6_r6) - s6_r6);
                
                Vec2 F = r_vec * factor;

                // Ley de Newton
                bolas[i].acc.x += F.x / bolas[i].masa;
                bolas[i].acc.y += F.y / bolas[i].masa;
                bolas[j].acc.x -= F.x / bolas[j].masa;
                bolas[j].acc.y -= F.y / bolas[j].masa;
                
                // --- CÁLCULO DEL VIRIAL (Clausius) ---
                // Virial par = F_ij dot r_ij
                virial_total += F.dot(r_vec);
            }
        }
    }
}

// ---------------- Calcular Presión Instantánea ----------------
double Caja::obtenerPresion(double area) const {
    // P = (N*k_B*T + Virial/dimensión) / Area
    // Aquí usamos la forma mecánica: P = (2*E_cin + Virial) / (D * Area)
    // D = 2 (dimensiones)
    
    double E_cin_total = 0.0;
    for(const auto& b : bolas) E_cin_total += b.energiaCin();
    
    // Presión = (Sum(m*v^2) + Sum(F*r)) / (2 * Area)
    // Nota: 2*E_cin = Sum(m*v^2)
    return (2.0 * E_cin_total + virial_total) / (2.0 * area);
}

// ---------------- Energía Potencial con PBC ----------------
double Caja::energiaPotencial() const {
    double E_pot = 0.0;
    size_t N = bolas.size();
    double r_cut_sq = r_cut * r_cut;
    double half_W = 0.5 * W;
    double half_H = 0.5 * H;
    double sigma6_base = sigma * sigma * sigma * sigma * sigma * sigma;

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            Vec2 r_vec = bolas[i].pos - bolas[j].pos;
            
            // PBC
            if (r_vec.x > half_W) r_vec.x -= W;
            else if (r_vec.x < -half_W) r_vec.x += W;
            if (r_vec.y > half_H) r_vec.y -= H;
            else if (r_vec.y < -half_H) r_vec.y += H;
            
            double r2 = r_vec.norm2();
            
            if (r2 < r_cut_sq) {
                double r2_inv = 1.0 / r2;
                double r6_inv = r2_inv * r2_inv * r2_inv;
                double s6_r6 = sigma6_base * r6_inv; 
                
                // V = 4*epsilon * [ (sigma/r)^12 - (sigma/r)^6 ]
                double U_LJ = 4.0 * epsilon * ((s6_r6 * s6_r6) - s6_r6);
                
                E_pot += (U_LJ - U_rcut);
            }
        }
    }
    return E_pot;
}

// ... Resto de funciones (guardarEstado, setParametrosLJ, energiaTotal) igual ...
void Caja::guardarEstado(std::ofstream& out, double t) const {
    out << t;
    for (const auto& b : bolas) {
        out << " " << b.pos.x << " " << b.pos.y << " " << b.vel.x << " " << b.vel.y;
    }
    out << "\n";
}

double Caja::energiaTotal() const {
    double E_cin = 0;
    for (const auto& b : bolas) E_cin += b.energiaCin();
    return E_cin + energiaPotencial();
}

void Caja::setParametrosLJ(double s, double e, double rc) {
    sigma = s;
    epsilon = e;
    r_cut = rc;
    if (r_cut > 0) {
        double r2 = r_cut*r_cut;
        double r2_inv = 1.0/r2;
        double r6_inv = r2_inv*r2_inv*r2_inv;
        double s6 = sigma*sigma*sigma*sigma*sigma*sigma;
        double s6_r6 = s6 * r6_inv;
        U_rcut = 4.0 * epsilon * (s6_r6*s6_r6 - s6_r6);
    } else U_rcut = 0.0;
}