#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <math.h>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream> 

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
int Npart;
VE N;
V m;
VV sig;
VV eps; 
VV rmax;
VS loc;

   
//parameteres del globus
int glob; 
double Rg; //radi
double p_ox; //centre 
double p_oy;
double d_o; // distancia repos molla
double Kg ; //fixar valor molla
double elong; //elongacio inicial

    
//parametres de simulacio
int Niter;
double dt;
int Nq; //longitud de la quadricula
int fs;
int seed;


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
ifstream fin;

void read_double(double &var) {

    string line;
    getline(fin, line);
    stringstream ss;
    ss << line;
    string aux;
    ss >> aux >> var;
}

void read_int(int &var) {

    string line;
    getline(fin, line);
    stringstream ss;
    ss << line;
    string aux;
    ss >> aux >> var;
}

void read_part(int k) {

    string line;
    getline(fin, line);
    stringstream ss;
    ss << line;
    string aux;
    
    ss >> aux >> N[k] >> m[k] >> loc[k];
    
    for (int l = 0; l < Npart; ++l) 
        ss >> sig[k][l];

    for (int l = 0; l < Npart; ++l) {
        ss >> eps[k][l];
        eps[k][l] *= kb;
    }
    
    for (int l = 0; l < Npart; ++l) 
        ss >> rmax[k][l];
        
}

void read_param(string infile) {
    
    fin.open(infile, std::ios_base::in);
    if (!fin.is_open()) {
        cerr << "No es pot obrir el fitxer - '"
             << infile << "'" << endl;
        exit(EXIT_FAILURE);
    }
    string line;
    
    //parameteres d'ambient
    read_double(T);
    read_double(g);
    read_double(lx);
    read_double(ly);
    
    getline(fin, line);
    read_int(Npart);
    getline(fin,line); getline(fin,line);
    
    
    //parametres de particula
    N = VE(Npart);
    m = V(Npart);
    sig = VV(Npart, V(Npart));
    eps = VV(Npart, V(Npart));
    rmax = VV(Npart, V(Npart));
    loc = VS(Npart);
    
    for (int k = 0; k < Npart; ++k)
        read_part(k);
    getline(fin, line);
    
    //parametres de globus
    read_int(glob);
    read_double(Rg);
    read_double(p_ox);
    read_double(p_oy);
    read_double(Kg);
    read_double(elong);
    getline(fin, line);
    
    //parametres de simulacio
    read_int(Niter);
    read_double(dt);
    read_int(Nq);
    read_int(fs);
    read_int(seed);
    if (seed == 0) seed = chrono::system_clock::now().time_since_epoch().count();    

}




void init_param() {
    
    x = VVV(Npart);
    y = VVV(Npart);
    vx = VV(Npart);
    vy = VV(Npart);
    fx = VV(Npart);
    fy = VV(Npart);
    A = VV(Npart, V(Npart));
    B = VV(Npart, V(Npart));

    for (int k = 0; k < Npart; ++k) {
        x[k] = VV(2, V(N[k]));
        y[k] = VV(2, V(N[k]));
        vx[k] = V(N[k]);
        vy[k] = V(N[k]);
        fx[k] = V(N[k]);
        fy[k] = V(N[k]);
    }
    
    for (int k = 0; k < Npart; ++k)
        for (int j = 0; j < Npart; ++j) {
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
    else if (on == "fora") return r < 1.1*Rg;
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
    
    for (int k = 0; k < Npart; ++k) 
        init_pk(k, loc[k]);
        
}


void init_v() {
    for (int k = 0; k < Npart; ++k) {
        if (k != glob) {
            default_random_engine generator;
            generator.seed(seed);
            double sd = sqrt(kb*T/m[k]);
            normal_distribution<double> distribution(0.0, sd);
            for (double &vxi : vx[k]) vxi = distribution(generator);
            for (double &vyi : vy[k]) vyi = distribution(generator);
        }
    }
}


void init_fout(string outfile) {

    fout.open(outfile, ios::out);
    fout << dt << endl;
    fout << lx << endl << ly << endl;
    fout << Npart << endl;
    for (int k = 0; k < Npart; ++k) fout << N[k] << endl;

}


//a diu si utilitzar la x actual o l'anterior
void force(int a) {
    
    for (int k = 0; k < Npart; ++k) {
        
        fx[k] = V(N[k], 0); 
        fy[k] = V(N[k], 0); 
        
        for (int i = 0; i < N[k]; ++i) {
            //for??a de la gravetat
            fy[k][i] += g*m[k];

            // iterem per la resta de particules (j != i (del mateix tipus) o l!=k (de tipus diferents))
            for (int l = 0; l < Npart; ++l) { 
                for (int j = 0; j < N[l]; ++j) {
                    
                    //for??a entre particules del globo
                    if (l == k and k == glob ) {
                        if (abs(i-j) == 1 or abs(i-j) == N[glob]-1) {
                            double dx = x[k][a][i]-x[l][a][j];
                            double dy = y[k][a][i]-y[l][a][j];
                            double r = sqrt(dx*dx+dy*dy);
                            fx[k][i] += -Kg*((r-d_o)*(dx/r));
                            fy[k][i] += -Kg*((r-d_o)*(dy/r));     
                        }                      
                    }
                    // for??a entre les demes particules
                    else if (j != i or l != k) {

                        //calculem for??a de j sobre i amb lennard jones
                        double r2 = pow(x[k][a][i] - x[l][a][j],2) + pow(y[k][a][i] - y[l][a][j],2);
                        if(r2 < rmax[k][l]*rmax[k][l]){
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
    
    for (int k = 0; k < Npart; ++k) {
        for (int i = 0; i < N[k]; ++i) {
            x[k][0][i] = x[k][1][i] + vx[k][i]*dt + 0.5*(fx[k][i]/m[k])*dt*dt;
            y[k][0][i] = y[k][1][i] + vy[k][i]*dt + 0.5*(fy[k][i]/m[k])*dt*dt;
        }
    }

}

void next_iteri(double l, double m, V &f, VV &p, V &v) {
   
    for (int i = 0; i < (int)p[0].size(); ++i) {
        
        // formula de l'algorisme de verlett
        double aux = p[0][i];
        double aux2 = p[1][i];
        
        
        p[0][i] = 2*aux - aux2 + f[i]/m *dt*dt;
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
        v[i] = (p[0][i] - aux2)/(2*dt);

    }
    
}

void next_iter() {

    force(0);
    for (int k = 0; k < Npart; ++k) {
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

    for (int k = 0; k < Npart; ++k) {

        for (int i = 0; i < N[k]; ++i) {

            // iterem per la resta de particules (j != i (del mateix tipus) o l!=k (de tipus diferents))
            for (int l = 0; l < Npart; ++l) { 
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
                        double r2 = pow(x[k][0][i]-x[l][0][j],2) + pow(y[k][0][i]-y[l][0][j],2);
                        if(r2 < rmax[k][l]*rmax[k][l]) {
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
    for (int k = 0; k < Npart; ++k) {
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
        
        first = true;
        for (int i = 0; i < N[k]; ++i) {
            if (first) first = false;
            else fout << " ";
            fout << vx[k][i];
        }
        fout << endl;
        
        first = true;
        for (int i = 0; i < N[k]; ++i) {
            if (first) first = false;
            else fout << " ";
            fout << vy[k][i];
        }
        fout << endl;
    }
    
}



int main(int argc, char *argv[])  {
    
    //llegir input
    read_param(argv[1]);
    
    
    //inicialitzar parametres
    init_param();
    
    //inicialitzar velocitats
    init_v();
    //inicialitzar posicions
    init_p();
    
    //fer primera iteracio de les posicions sense verlet
    primera_iter();

    //inicialitzar el fitxer de sortida
    init_fout(argv[2]);

    for (int it = 0; it < Niter; ++it) {
        next_iter();
        if (it%fs == 0)  guardar(it);
        if (it%fs == 0) cout << it << endl;
    }
    
    

}
