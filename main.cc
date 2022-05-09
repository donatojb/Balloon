#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>

using namespace std;

typedef vector<double> V;
typedef vector<V> VV;

const double k = 1.3806e-23;

void init_v(double T, double m, V &v) {
    
    default_random_engine generator;
    normal_distribution<double> distribution(0.0, sqrt(k*T/m));
    for (double &vi : v) vi = distribution(generator);
        
}

void init_p(double l, V &p) {
    
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0, l);
    for (double &pi : p) pi = distribution(generator);
        
}

void force(double A, double B, V &x, V &y, V &fx, V&fy) {

    fx = V(fx.size(), 0); 
    fy = V(fy.size(), 0); 
    
    for (int i = 0; i < x.size(); ++i) 
        for (int j = 0; j < x.size(); ++j) 
            if (j != i) {
            
                double r = sqrt(x[i]
            
            
            }
        

}

int main() {
    
    
    
    
    //paramtres de les particules i la caixa
    int N = 100;
    double m = 10e-27;
    double sig = 10e-10;
    double eps = 10e-10;
    double T = 300;
    double lx = 1e-8;
    double ly = 1e-8;
    
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
    init_v(T, m, vx);
    init_v(T, m, vy);
    
    //inicialitzar posicions
    init_p(lx, x[1]);
    init_p(ly, y[1]);
    
    //inicialitzar forces
    
    
    
}
    
