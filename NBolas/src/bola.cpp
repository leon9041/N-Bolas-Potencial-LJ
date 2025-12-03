/**
 * @file     bola.cpp
 * @brief    simulacion de Bola en caja con Velocity Verlet
 * @author   Angie Gomez, Leonardo Tovar
 * @date     02/12/25
 * @version  1.0
 * @license  owner
 */

#include "Bola.h"
#include <iostream>

// iniciar la partícula 
Bola::Bola(int id_, double x, double y, double vx, double vy, double r, double m)
    : pos(x, y), vel(vx, vy), acc(0.0, 0.0), acc_old(0.0, 0.0), radio(r), masa(m), id(id_) {}

// Paso medio de velocidad (Velocity-Verlet)
void Bola::vel_half_step(const Vec2& acc, double dt) {
    vel.x += 0.5 * acc.x * dt;
    vel.y += 0.5 * acc.y * dt;
}

// Paso completo de posición
void Bola::pos_full_step(double dt) {
    pos.x += vel.x * dt;
    pos.y += vel.y * dt;
}

// Completar la velocidad (segunda mitad del paso)
void Bola::vel_full_step(const Vec2& acc, double dt) {
    vel.x += 0.5 * acc.x * dt;
    vel.y += 0.5 * acc.y * dt;
}

// ---------------- Rebotar con paredes ----------------
void Bola::rebotePared(double W, double H) {
    // Rebote elástico PERFECTO con paredes
    bool rebotó = false;
    
    if (pos.x - radio < 0.0) {
        pos.x = radio;  // Posición exacta en la pared
        vel.x = -vel.x;
        rebotó = true;
    } else if (pos.x + radio > W) {
        pos.x = W - radio;
        vel.x = -vel.x;
        rebotó = true;
    }

    if (pos.y - radio < 0.0) {
        pos.y = radio;
        vel.y = -vel.y;
        rebotó = true;
    } else if (pos.y + radio > H) {
        pos.y = H - radio;
        vel.y = -vel.y;
        rebotó = true;
    }
    
    // Si rebotó, asegurar que no quede fuera de bounds
    if (rebotó) {
        if (pos.x < radio) pos.x = radio;
        if (pos.x > W - radio) pos.x = W - radio;
        if (pos.y < radio) pos.y = radio;
        if (pos.y > H - radio) pos.y = H - radio;
    }
}
// ---------------- Detectar si dos bolas se tocan ----------------
bool Bola::estaSolapadaCon(const Bola& otra) const {
    Vec2 d = pos - otra.pos;
    return d.norm() < (radio + otra.radio);
}

// ---------------- Resolver colisión entre dos bolas iguales ----------------
void Bola::resolverColision(Bola& otra) {
    Vec2 d = pos - otra.pos;
    double dist = d.norm();

    if (dist == 0.0) return; // Evita división entre cero

    Vec2 n = d / dist; // vector unitario normal entre bolas

    // velocidad relativa
    Vec2 vr = vel - otra.vel;
    double vn = vr.dot(n);

    // si se alejan, no hacer nada
    if (vn >= 0) return;

    // Intercambio de componentes normales (colisión elástica igual masa)
    double vi_n = vel.dot(n);
    double vj_n = otra.vel.dot(n);
    double dui = (vj_n - vi_n);
    vel.x += dui * n.x;
    vel.y += dui * n.y;

    double duj = (vi_n - vj_n);
    otra.vel.x += duj * n.x;
    otra.vel.y += duj * n.y;

    // Separar las bolas si se solapan
    double overlap = (radio + otra.radio) - dist;
    if (overlap > 0) {
        double shift = overlap / 2.0;
        pos.x += n.x * shift;
        pos.y += n.y * shift;
        otra.pos.x -= n.x * shift;
        otra.pos.y -= n.y * shift;
    }
}

// ---------------- Energía cinética ----------------
double Bola::energiaCin() const {
    return 0.5 * masa * (vel.x * vel.x + vel.y * vel.y);
}

// ---------------- Imprimir datos ----------------
void Bola::print(std::ostream& os) const {
    os << pos.x << " " << pos.y << " " << vel.x << " " << vel.y;
}

// rebote para medir presión
bool Bola::reboteParedConImpulso(double W, double H, double& impulso) {
    impulso = 0.0;      // reinicia el impulso acumulado
    bool choco = false; // indica si hubo al menos un rebote

    // Rebote con paredes verticales (izquierda y derecha)
    if (pos.x - radio < 0.0) {
        pos.x = radio;
        impulso += 2.0 * masa * std::abs(vel.x);
        vel.x = -vel.x;
        choco = true;
    }
    else if (pos.x + radio > W) {
        pos.x = W - radio;
        impulso += 2.0 * masa * std::abs(vel.x);
        vel.x = -vel.x;
        choco = true;
    }

    // Rebote con paredes horizontales (abajo y arriba)
    if (pos.y - radio < 0.0) {
        pos.y = radio;
        impulso += 2.0 * masa * std::abs(vel.y);
        vel.y = -vel.y;
        choco = true;
    }
    else if (pos.y + radio > H) {
        pos.y = H - radio;
        impulso += 2.0 * masa * std::abs(vel.y);
        vel.y = -vel.y;
        choco = true;
    }

    return choco;
}


