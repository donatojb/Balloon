#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>
#include <cmath>
#include <chrono>
#include <fstream>

using namespace std;

typedef vector<int> VE;
typedef vector<double> V;
typedef vector<V> VV;
typedef vector<VV> VVV;

//constants
const double kb = 1.3806e-23;
const double pi = M_PI;

   
//paramtres de les particules 
VE N = {10, 10};
V m = {6.6335209e-26, 6.6335209e-25};
VV sig = {{2.576e-10, 2.997e-10}, {2.997e-10, 3.418e-10}};
VV eps = {{10.2*kb, 34.8*kb}, {34.8*kb, 119*kb}}; 

//parametres d'ambient
double T = 300;
double g = 0;
double lx = 3*1e-8;
double ly = 3*1e-8;
    
//parametres de simulacio
int Niter = 300;
double dt = 1e-15;
double rmax = 10*1e-10;

    
//vectors de informacio
VVV x = VVV(N.size());
VVV y = VVV(N.size());
VV vx = VV(N.size());
VV vy = VV(N.size());
VV fx = VV(N.size());
VV fy = VV(N.size());
VV A = VV(N.size(), V(N.size()));
VV B = VV(N.size(), V(N.size()));

//fitxer de sortida
fstream fout;



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
        
        for (int i = 0; i < N[k]; ++i) {
            //força de la gravetat
            fy[k][i] += g*m[k];

            // iterem per la resta de particules (j != i (del mateix tipus) o l!=k (de tipus diferents))
            for (int l = 0; l < (int)N.size(); ++l) { 
                for (int j = 0; j < N[l]; ++j) {
                    if (j != i or l != k) {

                        //calculem força de j sobre i amb lennard jones
                        double r2 = pow(x[k][a][j]-x[l][a][i],2) + pow(y[k][a][j]-y[l][a][i],2);
                        if(r2<rmax*rmax){
                            double f = A[k][l]/pow(r2,7) - B[k][l]/pow(r2,4);
                            fx[k][i] += f*(x[k][a][i]-x[l][a][j]);
                            fy[k][i] += f*(y[k][a][i]-y[l][a][j]);
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
    for (int k = 0; k < (int)N.size(); ++k) {
        next_iteri(lx, m[k], fx[k], x[k], vx[k]);
        next_iteri(ly, m[k], fy[k], y[k], vy[k]);
    }
}


//k es el tipus de particula
double kineticE (int k) {
    double KE = 0;
    for (int i = 0; i < N[k]; ++i) {
            KE += 0.5 * m[k] * (vx[k][i]*vx[k][i] + vy[k][i]*vy[k][i]);
     }
    return KE;
}

double potentialE () {

    double PE = 0;

    for (int k = 0; k < (int)N.size(); ++k) {

        for (int i = 0; i < N[k]; ++i) {

            // iterem per la resta de particules (j != i (del mateix tipus) o l!=k (de tipus diferents))
            for (int l = 0; l < (int)N.size(); ++l) { 
                for (int j = 0; j < N[l]; ++j) {
                    if (j != i or l != k) {

                        //calculem el potencial de j sobre i amb lennard jones
                        double r2 = pow(x[k][0][j]-x[l][0][i],2) + pow(y[k][0][j]-y[l][0][i],2);
                        if(r2 < rmax*rmax) {
                            PE += (A[k][l]/12)/pow(r2,6) - (B[k][l]/6)/pow(r2,3);
                        }
                    }
                }
            }
        }
    }

    return PE/2;
}




void guardar(int it) {
	fout << it << endl;
    fout << potentialE() << endl;
    for (int k = 0; k < (int)N.size(); ++k) {
        fout << kineticE(k) << endl;

        bool first = true;
        for (int i = 0; i < N[k]; ++i) {
            if (first) first = false;
            else fout << " ";
            fout << x[k][0][i];
        }
        fout << endl;

        first = true;
        for (int i = 0; i < N[k]; ++i) {
            if (first) first = false;
            else fout << " ";
            fout << y[k][0][i];
        }
        fout << endl;
    }
}

void init_fout() {

    fout.open("in.fo", ios::out);
    fout << dt << endl;
    fout << lx << endl << ly << endl;
    fout << (int)N.size() << endl;
    for (int k = 0; k < (int)N.size(); ++k) fout << N[k] << endl;

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

    //inicialitzar el fitxer de sortida
    init_fout();

    for (int it = 0; it < Niter; ++it) {
        next_iter();
        if (it%100 == 0)  guardar(it);
    }
    
    

}
