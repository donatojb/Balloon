#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>
#include <cmath>
#include <chrono>

using namespace std;

typedef vector<int> VE;
typedef vector<double> V;
typedef vector<V> VV;
typedef vector<VV> VVV;

const double kb = 1.3806e-23;
const double pi = M_PI;

   
//paramtres de les particules i la caixa
int Niter = 5000000;
VE N = {100, 100};
V m = {6.6335209e-26, 6.6335209e-25};
VV sig = {{2.576e-10, 2.997e-10}, {2.997e-10, 3.418e-10}};
VV eps = {{10.2*kb, 34.8*kb}, {34.8*kb, 119*kb}}; 
double T = 300;
double lx = 3*1e-8;
double ly = 3*1e-8;
double rmax = 10*1e-10;
    
//parametres de simulacio
double dt = 1e-15;
VV A = VV(N.size(), V(N.size()));
VV B = VV(N.size(), V(N.size()));

    
//vectors de informacio
VVV x = VVV(N.size());
VVV y = VVV(N.size());
VV vx = VV(N.size());
VV vy = VV(N.size());
VV fx = VV(N.size());
VV fy = VV(N.size());




void escriure(V &v, bool salt){
	for (double vi : v) cout << vi << " ";
	if (salt) cout << endl;
}

//a diu si iteracio anterior o actual
void visualitza(int a) {
	for (int k = 0; k < (int)N.size(); ++k) 
		escriure(x[k][a], false);
	cout << endl;
	for (int k = 0; k < (int)N.size(); ++k)	
		escriure(y[k][a], false);
	cout << endl;

}


//part numero indicador del tipus particula
void init_p() {
    
    for (int k = 0; k < N.size(); ++k) {
        default_random_engine generator;
        generator.seed(chrono::system_clock::now().time_since_epoch().count());
        uniform_real_distribution<double> distributionx(0.0, lx);
        for (double &xi : x[k][1]) xi = distributionx(generator);
        uniform_real_distribution<double> distributiony(0.0, ly);
        for (double &yi : y[k][1]) yi = distributiony(generator);
    }
        
}

void init_v() {
    for (int k = 0; k < N.size(); ++k) {
        default_random_engine generator;
        generator.seed(chrono::system_clock::now().time_since_epoch().count());
        double sd = sqrt(kb*T/m[k]);
        normal_distribution<double> distribution(0.0, sd);
        for (double &vxi : vx[k]) vxi = distribution(generator);
        for (double &vyi : vy[k]) vyi = distribution(generator);
    }
}

void init_param() {

    for (int k = 0; k < N.size(); ++k) {
        x[k] = VV(2, V(N[k]));
        y[k] = VV(2, V(N[k]));
        vx[k] = V(N[k]);
        vy[k] = V(N[k]);
        fx[k] = V(N[k]);
        fy[k] = V(N[k]);
    }
    
    for (int k = 0; k < N.size(); ++k)
        for (int j = 0; j < N.size(); ++j) {
            A[k][j] = 4*12*eps[k][j]*pow(sig[k][j],12);
            B[k][j] = 4*6*eps[k][j]*pow(sig[k][j],6);
    }
    

}


//a diu si utilitzar la x actual o l'anterior
void force(int a) {
    
    for (int k = 0; k < (int)N.size(); ++k) {
        
        fx[k] = V(fx[k].size(), 0); 
        fy[k] = V(fy[k].size(), 0); 
        
        for (int h = 0; h < (int)N.size(); ++h) {
            
            for (int i = 0; i < N[h]; ++i) {
                for (int j = 0; j < N[k]; ++j) {
                    
                    fy[k][j] += -9.81e10*m[k];
                    if (j != i or h != k) {
                        //calculem força de i sobre j- lennard jones
                        double r2 = pow(x[k][a][j]-x[h][a][i],2) + pow(y[k][a][j]-y[h][a][i],2);
                        if(r2<rmax*rmax){
                            double f = A[k][h]/pow(r2,7) - B[k][h]/pow(r2,4);
                            fx[k][j] += f*(x[k][a][j]-x[h][a][i]);
                            fy[k][j] += f*(y[k][a][j]-y[h][a][i]);
                        }
                    }
                }
            }
        }
    }
}

void primera_iter() {
    
    force(1);
    
    for (int k = 0; k < (int)N.size(); ++k)
        for (int i = 0; i < N[k]; ++i) {
            x[k][0][i] = x[k][1][i] + vx[k][i]*dt + 0.5*(fx[k][i]/m[k])*dt*dt;
            y[k][0][i] = y[k][1][i] + vy[k][i]*dt + 0.5*(fy[k][i]/m[k])*dt*dt;
        }
}

void next_iteri(double l, double m, V &f, VV &p, V &v) {
   
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
        }
        // calcul de la velocitat
        v[i] = (p[0][i] - aux2)/2*dt;
    }
    
}

void next_iter() {

    force(0);
    for (int i = 0; i < (int)N.size(); ++i) {
        next_iteri(lx, m[i], fx[i], x[i], vx[i]);
        next_iteri(ly, m[i], fy[i], y[i], vy[i]);
    }
}



double kineticE () {
    double E = 0;
    for (int j = 0; j < N.size(); ++j) {
        for (int i = 0; i < (int)vx[j].size(); ++i) {
            
            E += 0.5 * m[j] * (vx[j][i]*vx[j][i] + vy[j][i]*vy[j][i]);
        }
     }
     return E;
}



int main() {
    
    //inicialitzar parametres
    init_param();
    
    //inicialitzar velocitats
    init_v();
    
    //inicialitzar posicions
    init_p();
    
    //fer primera iteracio de les posicions sense verlet
    primera_iter();
    visualitza(1);
    
	
    for (int k = 0; k < Niter; ++k) {
        next_iter();
        if (k%100 == 0)  visualitza(0);
    }
    
    

}
