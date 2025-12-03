#ifndef BOLA_H
#define BOLA_H
#include <cmath>
#include <ostream>

struct Vec2 {
    double x, y;
    Vec2() : x(0), y(0) {}
    Vec2(double x_, double y_) : x(x_), y(y_) {}
    Vec2 operator+(const Vec2& o) const { return Vec2(x + o.x, y + o.y); }
    Vec2 operator-(const Vec2& o) const { return Vec2(x - o.x, y - o.y); }
    Vec2 operator*(double s) const { return Vec2(x * s, y * s); }
    Vec2 operator/(double s) const { return Vec2(x / s, y / s); }
    double dot(const Vec2& o) const { return x * o.x + y * o.y; }
    double norm2() const { return x * x + y * y;  }
    double norm() const { return std::sqrt(norm2()); }
};

class Bola {
public:
    Vec2 pos;
    Vec2 vel;
    Vec2 acc;       // Nueva aceleración del paso actual
    Vec2 acc_old;   // Aceleración del paso anterior
    double radio;
    double masa;
    int id;
    
    Bola(int id_, double x, double y, double vx, double vy, double r = 0.01, double m = 1.0);
    void vel_half_step(const Vec2& acc, double dt);
    void pos_full_step(double dt);
    void vel_full_step(const Vec2& acc, double dt);

    void rebotePared(double W, double H);     
    bool estaSolapadaCon(const Bola& otra) const;
    void resolverColision(Bola& otra);
    double energiaCin() const;
    void print(std::ostream& os) const;
    //versión para contar impulsos (presión)
    bool reboteParedConImpulso(double W, double H, double& impulso);
    void calcular_fuerza_y_aceleracion(const Bola& otra, double sigma, double epsilon);
};

#endif
