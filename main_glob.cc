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

//valors amb els que crec que es conserva energia cinetica:
//Argo, l = 5e-8, dt= 1e-15, T=300
   
//paramtres de les particules i la caixa
int Niter = 30000;
VE N = {200, 100};
V m = {1e-26, 6.6335209e-25};//6.6335209e-26     
VV sig = {{2.576e-10, 2.997e-10}, {2.997e-10, 3.418e-10}};
VV eps = {{10.2*kb, 34.8*kb}, {34.8*kb, 119*kb}}; 
double T = 300;
double lx = 3*1e-8;
double ly = 3*1e-8;
double rmax = 10*1e-10;
    
//parametres de simulacio
double dt = 1e-17;
VV A = VV(N.size(), V(N.size()));
VV B = VV(N.size(), V(N.size()));

    
//vectors de informacio
VVV x = VVV(N.size());
VVV y = VVV(N.size());
VV vx = VV(N.size());
VV vy = VV(N.size());
VV fx = VV(N.size());
VV fy = VV(N.size());

//param globus
double Rad = lx/3; //radi
double p_ox = lx/2; // centre
double p_oy = ly/2;
double d_o; // distancia repos molla
double k = 400; //fixar valor molla


void escriure(V &v, bool salt){
	for (double vi : v) cout << vi << " ";
	if (salt) cout << endl;
}

//a diu si iteracio anterior o actual
void visualitza(int a) {
	for (int i = 0; i < (int)N.size(); ++i) 
		escriure(x[i][a], false);
	cout << endl;
	for (int i = 0; i < (int)N.size(); ++i)	
		escriure(y[i][a], false);
	cout << endl;

}

void init_v() {
    for (int i = 0; i < N.size(); ++i) {
        default_random_engine generator;
        generator.seed(chrono::system_clock::now().time_since_epoch().count());
        double sd = sqrt(kb*T/m[i]);
        normal_distribution<double> distribution(0.0, sd);
        for (double &vxi : vx[i]) vxi = distribution(generator);
        for (double &vyi : vy[i]) vyi = distribution(generator);
    }
        
}

//part numero indicador del tipus particula
void init_p(int part) {
    
    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> distribution(0.0, lx);
    for (double &xi : x[part][1]) xi = distribution(generator);
    for (double &yi : y[part][1]) yi = distribution(generator);
        
}


void init_pglob(int part) {
    for (int i = 0; i < N[part]; ++i) {
        x[part][1][i] = p_ox + Rad*cos((2*pi*i)/N[part]);
        y[part][1][i] = p_oy + Rad*sin((2*pi*i)/N[part]);
    }
    double dx = x[part][1][0]-x[part][1][1];
    double dy = y[part][1][0]-y[part][1][1];
    d_o = sqrt(dx*dx + dy*dy);
    
}


void init_pcercle(int part) {
	double R = Rad*(0.9);
	default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < N[part]; ++i) {
		double r = R*sqrt(distribution(generator));
		double theta = distribution(generator)*2*pi;
        x[part][1][i] = p_ox + r*cos(theta);
        y[part][1][i] = p_oy + r*sin(theta);
    }
	
}


//a diu si utilitzar la x actual o l'anterior
void force(int a) {
    for (int l = 0; l < (int)N.size(); ++l) {
        fx[l] = V(fx[l].size(), 0); 
        fy[l] = V(fy[l].size(), 0); 
        for (int h = 0; h < (int)N.size(); ++h) {
            for (int i = 0; i < (int)x[h][a].size(); ++i) 
                for (int j = 0; j < (int)x[l][a].size(); ++j) {
                    if (l == h and l == 1) {
                        //força elastica entre particules contigues
                        if (abs(i-j) == 1 or abs(i-j) == N[1]-1) {
                            double dx = x[l][a][j]-x[h][a][i];
                            double dy = y[l][a][j]-y[h][a][i];
                            double r = sqrt(dx*dx+dy*dy);
                            fx[l][j] += -0.5*k*(dx-d_o*(dx/r));
                            fy[l][j] += -0.5*k*(dy-d_o*(dy/r));                                
                        }
                    }
                    else if (j != i or h != l) {
                        //calculem força de i sobre j- lennard jones
                        double r2 = pow(x[l][a][j]-x[h][a][i],2) + pow(y[l][a][j]-y[h][a][i],2);
                        if(r2<rmax*rmax){
                            double f = A[l][h]/pow(r2,7) - B[l][h]/pow(r2,4);
                            fx[l][j] += f*(x[l][a][j]-x[h][a][i]);
                            fy[l][j] += f*(y[l][a][j]-y[h][a][i]);
                        }
                    }
                }
        }
    }
}

void primera_iter() {
    for (int j = 0; j < (int)N.size(); ++j)
        for (int i = 0; i < (int)x[j][0].size(); ++i) {
            x[j][0][i] = x[j][1][i] + vx[j][i]*dt + 0.5*(fx[j][i]/m[j])*dt*dt;
            y[j][0][i] = y[j][1][i] + vy[j][i]*dt + 0.5*(fy[j][i]/m[j])*dt*dt;
        }
}

void next_iter(double l, double m, V &f, VV &p, V &v) {
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
    for (int i = 0; i < N.size(); ++i) {
        x[i] = VV(2, V(N[i]));
        y[i] = VV(2, V(N[i]));
        vx[i] = V(N[i]);
        vy[i] = V(N[i]);
        fx[i] = V(N[i]);
        fy[i] = V(N[i]);
    }
    
    for (int i = 0; i < N.size(); ++i)
        for (int j = 0; j < N.size(); ++j) {
            A[i][j] = 4*12*eps[i][j]*pow(sig[i][j],12);
            B[i][j] = 4*6*eps[i][j]*pow(sig[i][j],6);
    }
    
    
    //inicialitzar velocitats
    init_v();
    
    //inicialitzar posicions
    
    //init_p(0);
    init_pcercle(0);
    init_pglob(1);
    
    visualitza(1);
    //inicialitzar forces
    force(1);
    
    //fer primera iteracio de les posicions sense verlet
    primera_iter();
	
	
    for (int k = 0; k < Niter; ++k) {
        
        force(0);
        for (int i = 0; i < (int)N.size(); ++i) {
            next_iter(lx, m[i], fx[i], x[i], vx[i]);
            next_iter(ly, m[i], fy[i], y[i], vy[i]);
        }
        if (k%100 == 0)  visualitza(0);
        
    }
    
    

}
