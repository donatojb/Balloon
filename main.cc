#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>
#include <cmath>
#include <chrono>
#include "matplotlibcpp.h"

//g++ main.cc -std=c++11 -I/usr/include/python3.8/ -lpython3.8 

//g++ main.cc -std=c++11 -I/home/donato/.local/lib/python3.8/site-packages/numpy/core/include -I/home/donato/.local/lib/python3.8/site-packages/tensorflow/include/external/local_config_python/python_include/ -lpython3.8 -lpthread -lutil -ldl -Xlinker -export-dynamic

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

void force(double m, double A, double B, V &x, V &y, V &fx, V&fy) {

    fx = V(fx.size(), 0); 
    fy = V(fy.size(), 0); 
    
    for (int i = 0; i < (int)x.size(); ++i) 
        for (int j = 0; j < (int)x.size(); ++j) 
            if (j != i) {
                //calculem força de i sobre j
                double r2 = pow(x[j]-x[i],2) + pow(y[j]-y[i],2);
                double f = A/pow(r2,7) - B/pow(r2,4);
                fx[j] += f*(x[j]-x[i]);
                fy[j] += f*(y[j]-y[i]) - 9.81*m*1e10;
            }
}

void primera_iter(double dt, double m, VV &p, V &v, V &f) {
    for (int i = 0; i < (int)p[0].size(); ++i) 
        p[0][i] = p[1][i] + v[i]*dt + 0.5*(f[i]/m)*dt*dt;
}

void next_iter(double l, double m, double dt, V &f, VV &p, V &v) {

    for (int i = 0; i < (int)p[0].size(); ++i) {
        
        
        
        // formula de l'algorisme de verlett
        double aux = p[0][i];
        double aux2 = p[1][i];
        p[0][i] = 2*aux - p[1][i] + f[i]/m *dt*dt;
        p[1][i] = aux;
        
        // condicions de vora
        if (p[0][i] < 0) {
            p[1][i] *= -1;
            p[0][i] *= -1;
        }
        else if (p[0][i] > l) {
            p[1][i] = 2*l - p[1][i];
            p[0][i] = 2*l - p[0][i];
        
        // calcul de la velocitat
        v[i] = (p[0][i] - aux2)/2*dt;
    }
}

double kineticE (double m, V &vx, V &vy) {

    double E = 0;
    for (int i = 0; i < (int)vx.size(); ++i) {
        
        E += 0.5 * m * (vx[i]*vx[i] + vy[i]*vy[i]);
    }
    return E;

}

int main() {
    
    
    
    
    //paramtres de les particules i la caixa
    int Niter = 100000;
    int N = 1000;
    //double m = 10e-27;
    double m =   6.6335209e-26;
    double sig = 2.576e-10;
     
    //double sig = 2.576e-10;
    double sig = 3.405e-10;
    //double eps = 10.2*kb;
    double eps = 119.8*kb;
    double T = 300;
    double lx = 1e-7;
    double ly = 1e-7;
    
    //parametres de simulacio
    double dt = 1e-13;
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
    force(m, A, B, x[1], y[1], fx, fy);
    
    //fer primera iteracio de les posicions sense verlet
    primera_iter(dt, m, x, vx, fx);
    primera_iter(dt, m, y, vy, fy);
       
       
    plt::title("Simulació");
    plt::Plot plot("Heli", "or");
    plt::legend();
    plt::xlim(0.0, lx);
    plt::ylim(0.0, ly);
    
    plot.update(x[1], y[1]);
    plt::pause(0.1);
    
    for (int k = 0; k < Niter; ++k) {
        plot.update(x[0], y[0]);
        plt::pause(0.005);
        
        force(m, A, B, x[0], y[0], fx, fy);
        next_iter(lx, m, dt, fx, x, vx);
        next_iter(ly, m, dt, fy, y, vy);
        
        double E = kineticE(m, vx, vy);
        cout << E << endl;
    }
    
    
   
    

}
