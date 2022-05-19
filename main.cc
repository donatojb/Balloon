#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>

using namespace std;

typedef vector<string> VS;
typedef vector<int> VE;
typedef vector<VE> VVE;
typedef vector<double> V;
typedef vector<V> VV;
typedef vector<VV> VVV;

//constants
const double kb = 1.3806e-23;
const double pi = M_PI;

   
//parametres d'ambient
double T;
double g;
double lx;
double ly;
   
//paramtres de les particules 
VE N;
V m;
VV sig;
VV eps; 
VS loc;

   
//parameteres del globus
int glob; 
double Rg; //radi
double p_ox; //centre 
double p_oy;
double d_o; // distancia repos molla
double Kg ;//20000 //fixar valor molla
double elong;  //elongacio inicial

    
//parametres de simulacio
int Niter;
double dt;
double rmax;
int Nq; //longitud de la quadricula
int seed = chrono::system_clock::now().time_since_epoch().count();



//vectors de informacio
VVV x;
VVV y;
VV vx;
VV vy;
VV fx;
VV fy;
VV A; 
VV B;

//fitxer de sortida
fstream fout;


void read_param() {
    
    //parameteres d'ambient
    cin >> T;
    cout << T;
    cin >> g;
    cout << g;
    cin >> lx >> ly;
    
    //parametres de particula
    int Npart; cin >> Npart;
    N = VE(Npart);
    m = V(Npart);
    sig = VV(Npart, V(Npart));
    eps = VV(Npart, V(Npart));
    loc = VS(Npart);
    
    for (int k = 0; k < Npart; ++k) {
        cin >> N[k] >> m[k] >> loc[k];
        for (int l = 0; l < Npart; ++l) {
            cin >> sig[k][l] >> eps[k][l];
            eps[k][l] *= kb;
        }
    }
    
    //parametres de globus
    cin >> glob >> Rg >> p_ox >> p_oy >> Kg >> elong;
    
    //parametres de simulacio
    cin >> Niter >> dt >> rmax >> Nq >> seed;
    if (seed == 0) seed = chrono::system_clock::now().time_since_epoch().count();
    

}




void init_param() {
    
    x = VVV(N.size());
    y = VVV(N.size());
    vx = VV(N.size());
    vy = VV(N.size());
    fx = VV(N.size());
    fy = VV(N.size());
    A = VV(N.size(), V(N.size()));
    B = VV(N.size(), V(N.size()));

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


void init_pglob(int part) {
    for (int i = 0; i < N[part]; ++i) {
        x[part][1][i] = p_ox + Rg*cos((2*pi*i)/N[part]);
        y[part][1][i] = p_oy + Rg*sin((2*pi*i)/N[part]);
    }
    
    double dx = x[part][1][0]-x[part][1][1];
    double dy = y[part][1][0]-y[part][1][1];
    d_o = (1/elong)*sqrt(dx*dx + dy*dy);
}


bool no_colocar(double r, string on) {
    if (on == "dins") return r > 0.9*Rg;
    else if (on == "fora") return r < 1.5*Rg;
    else return (r < 1.1*Rg and r > 0.9*Rg);
}

void init_pk(int part, string on) {

    if (on == "glob") {
        
        init_pglob(part);
        return;
    }
    
    default_random_engine generator;
    generator.seed(seed);
    VVE Ocupat(Nq, VE(Nq, 0));
    cout << Nq << endl;
    Ocupat[0][0] = 1;
    uniform_int_distribution<int> distribution(1, Nq-1);
    for (int i = 0; i < N[part]; ++i) {
        int xi = 0;
        int yi = 0;
        double r = Rg;
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


void init_p() {
    
    for (int k = 0; k < N.size(); ++k) 
        init_pk(k, loc[k]);
        
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



//a diu si utilitzar la x actual o l'anterior
void force(int a) {
    
    for (int k = 0; k < (int)N.size(); ++k) {
        
        fx[k] = V(N[k], 0); 
        fy[k] = V(N[k], 0); 
        
        for (int i = 0; i < N[k]; ++i) {
            //força de la gravetat
            fy[k][i] += g*m[k];

            // iterem per la resta de particules (j != i (del mateix tipus) o l!=k (de tipus diferents))
            for (int l = 0; l < (int)N.size(); ++l) { 
                for (int j = 0; j < N[l]; ++j) {
                    
                    //força entre particules del globo
                    if (l == k and k == glob ) {
                        if (abs(i-j) == 1 or abs(i-j) == N[glob]-1) {
                            double dx = x[k][a][i]-x[l][a][j];
                            double dy = y[k][a][i]-y[l][a][j];
                            double r = sqrt(dx*dx+dy*dy);
                            fx[k][i] += -Kg*(dx-d_o*(dx/r));
                            fy[k][i] += -Kg*(dy-d_o*(dy/r));     
                        }                      
                    }
                    // força entre les demes particules
                    else if (j != i or l != k) {

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
                    
                    //potencial entre les particules del globo
                    if (l == k and k == glob ) {
                        if (abs(i-j) == 1 or abs(i-j) == N[glob]-1) {
                            double dx = x[k][0][i]-x[l][0][j];
                            double dy = y[k][0][i]-y[l][0][j];
                            double r2 = dx*dx+dy*dy;
                            PE += 0.5*Kg*r2;
                        }                      
                    }
                    // potencial entre les demes particules
                    else if (j != i or l != k) {
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
    
    //llegir input
    read_param();
    
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
        if (it%1 == 0)  guardar(it);
        if (it%100 == 0) cout << it << endl;
    }
    
    

}
