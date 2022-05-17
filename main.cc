#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>
#include <cmath>
//#include <omp.h>
#include <chrono>

using namespace std;

typedef vector<int> VE;
typedef vector<VE> VVE;
typedef vector<double> V;
typedef vector<V> VV;
typedef vector<VV> VVV;

const double kb = 1.3806e-23;
const double pi = M_PI;


int seed = chrono::system_clock::now().time_since_epoch().count();
//int seed = 89;
   
//paramtres de les particules i la caixa
int Niter = 1000000;
int part_globo = 2;
VE N = {100, 1000, 100};
V m = {1e-26, 2e-26,  6.6335209e-25};//6.6335209e-26     
VV sig = {{2.576e-10, 2.576e-10, 2.997e-9},{2.576e-10, 2.576e-10, 2.997e-9}, {2.997e-9, 2.997e-9, 3.418e-10}};
VV eps = {{10.2*kb, 10.2*kb, 34.8*kb},{10.2*kb, 10.2*kb, 34.8*kb},  {34.8*kb, 34.8*kb, 119*kb}}; 
double T = 10000; 
double lx = 3*1e-8;
double ly = 5*lx;
double rmax = 10*1e-10;
    
//parametres de simulacio
double dt = 1e-17;
int Nq = 500;
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
double k = 5000;//20000 //fixar valor molla


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
        generator.seed(seed);
        double sd = sqrt(kb*T/m[i]);
        normal_distribution<double> distribution(0.0, sd);
        for (double &vxi : vx[i]) vxi = distribution(generator);
        for (double &vyi : vy[i]) vyi = distribution(generator);
    }
        
}

//part numero indicador del tipus particula
void init_p(int part) {
    
    default_random_engine generator;
    generator.seed(seed);
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
    d_o = 0.75*sqrt(dx*dx + dy*dy);
    
}


void init_pcercle(int part) {
	double R = Rad*(0.9);
	default_random_engine generator;
    generator.seed(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < N[part]; ++i) {
		double r;
        if (i%2 == 0) r = R*sqrt(distribution(generator));
		else r = Rad*1.1 + (lx/2-Rad)*sqrt(distribution(generator));
        double theta = distribution(generator)*2*pi;
        x[part][1][i] = p_ox + r*cos(theta);
        y[part][1][i] = p_oy + r*sin(theta);
    }
	
}

bool no_colocar(double r, string on) {
    if (on == "dins") return r > 0.9*Rad;
    else if (on == "fora") return r < 1.5*Rad;
    else return (r < 1.1*Rad and r > 0.9*Rad);
}
void init_pquadriculat(int part, string on) {
    default_random_engine generator;
    generator.seed(seed);
    VVE Ocupat(Nq, VE(Nq, 0));
    Ocupat[0][0] = 1;
    uniform_int_distribution<int> distribution(1, Nq-1);
    for (int i = 0; i < N[part]; ++i) {
        int xi = 0;
        int yi = 0;
        double r = Rad;
        while (Ocupat[xi][yi] or no_colocar(r, on)) {
            xi = distribution(generator);
            yi = distribution(generator);
            double xd = lx*double(xi)/Nq;
            double yd = ly*double(yi)/Nq;
            r = sqrt((xd-p_ox)*(xd-p_ox) + (yd-p_oy)*(yd-p_oy));
        }
        Ocupat[xi][yi] = 1;
        x[part][1][i] = lx*double(xi)/Nq;
        y[part][1][i] = ly*double(yi)/Nq;
    }
}

//a diu si utilitzar la x actual o l'anterior
void force(int a) {
    //#pragma omp parallel
    //#pragma omp single
    for (int l = 0; l < (int)N.size(); ++l) {
        fx[l] = V(fx[l].size(), 0); 
        fy[l] = V(fy[l].size(), 0); 
        for (int j = 0; j < N[l]; ++j) {
            //#pragma omp task
            for (int h = 0; h < (int)N.size(); ++h)
                for (int i = 0; i < N[h]; ++i) {
                    if (l == h and l == part_globo) {
                        //força elastica entre particules contigues
                        if (abs(i-j) == 1 or abs(i-j) == N[part_globo]-1) {
                            double dx = x[l][a][j]-x[h][a][i];
                            double dy = y[l][a][j]-y[h][a][i];
                            double r = sqrt(dx*dx+dy*dy);
                            fx[l][j] += -k*(dx-d_o*(dx/r));
                            fy[l][j] += -k*(dy-d_o*(dy/r));                                
                        }
                    }
                    else if (j != i or h != l) {
                        //calculem força de i sobre j- lennard jones
                        double r2 = pow(x[l][a][j]-x[h][a][i],2) + pow(y[l][a][j]-y[h][a][i],2);
                        if(r2<rmax*rmax){
                            double f = A[l][h]/pow(r2,7) - B[l][h]/pow(r2,4);
                            fx[l][j] += f*(x[l][a][j]-x[h][a][i]);
                            fy[l][j] += f*(y[l][a][j]-y[h][a][i])- 9.81*m[j];
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
    //init_pcercle(0);
    init_pquadriculat(0, "dins");
    init_pquadriculat(1, "fora");
    init_pglob(part_globo);
    
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
        if (k%10 == 0)  visualitza(0);
        
    }
    
    

}
