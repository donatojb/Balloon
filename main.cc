#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>
#include <cmath>
#include <chrono>
#include "matplotlibcpp.h"

//g++ main.cc -std=c++11 -I/usr/include/python3.8/ -lpython3.8 -lpthread -lutil -ldl -Xlinker -export-dynamic


using namespace std;
namespace plt = matplotlibcpp;

typedef vector<double> V;
typedef vector<V> VV;

const double kb = 1.3806e-23;

void init_v(double T, double m, V &vx, V &vy) {
    
    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count());
    double sd = sqrt(kb*T/m);
    normal_distribution<double> distribution(0.0, sd);
    for (double &vxi : vx) vxi = distribution(generator);
    for (double &vyi : vy) vyi = distribution(generator);
        
}

void init_p(double l, V &x, V &y) {
    
    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> distribution(0.0, l);
    for (double &xi : x) xi = distribution(generator);
    for (double &yi : y) yi = distribution(generator);
        
}

void force(double A, double B, V &x, V &y, V &fx, V&fy) {

    fx = V(fx.size(), 0); 
    fy = V(fy.size(), 0); 
    
    for (int i = 0; i < (int)x.size(); ++i) 
        for (int j = 0; j < (int)x.size(); ++j) 
            if (j != i) {
                //calculem força de i sobre j
                double r2 = pow(x[j]-x[i],2) + pow(y[j]-y[i],2);
                double f = A/pow(r2,7) - B/pow(r2,4);
                fx[j] += f*(x[j]-x[i]);
                fy[j] += f*(y[j]-y[i]);
            }
}

void primera_iter(double dt, double m, VV &p, V &v, V &f) {
    for (int i = 0; i < (int)p[0].size(); ++i) 
        p[0][i] = p[1][i] + v[i]*dt + 0.5*(f[i]/m)*dt*dt;

}

int main() {
    
    
    
    
    //paramtres de les particules i la caixa
    int N = 100;
    double m = 10e-27;
    double sig = 2.576e-10;
    double eps = 10.2*kb;
    double T = 300;
    double lx = 1e-7;
    double ly = 1e-7;
    
    //parametres de simulacio
    double dt = 10e-13;
    double A = 4*12*eps*pow(sig,12);
    double B = 4*6*eps*pow(sig,6);
    
    //vectors de informacio
    VV x = VV(2, V(N));
    VV y = VV(2, V(N));
    V vx = V(N);
    V vy = V(N);
    V fx = V(N);
    V fy = V(N);
    
    //inicialitzar velocitats
    init_v(T, m, vx, vy);
    
    //inicialitzar posicions
    init_p(lx, x[1], y[1]);
    
    //inicialitzar forces
    force(A, B, x[1], y[1], fx, fy);
    
    //fer primera iteracio de les posicions sense verlet
    primera_iter(dt, m, x, vx, fx);
    primera_iter(dt, m, y, vy, fy);

    plt::figure(1); 
       
    plt::Plot plot("primera");
    
    plt::plot(x[0], y[0], "or"); 
   // plt::xlim(0.0,lx);
    //plt::ylim(0.0,ly);
    plt::show();
    

}
