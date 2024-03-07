#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
using namespace std;
//CASE OF A CYLINDER WITH nf FINS
//Declaration of global variables

//FISICAL PARAMETERS
double Rint = 0; // Internal radius of the fin
double Rext = 10; // External radius of the fin
double ef = 10E-3; // Fin's vertical length [m]
double eb = 0; // Distance between fins [m]
double Text = 10; // Temp at exterior [Cº]
double T_wall = 100; // Temp of the wall r = Rint [Cº]
//double lambda_f[2] = [200,10]; // lambda of the material (homogeneous case) [a_nT^n,...a_1T,a_0] (function of T)
double lambda_f = 100;
double alfa_exterior = 0; //alpha in the length of the fin with outter fluid
double alfa_extrem = 0; //alpha in the end of fin
double nf = 5; //Number of fins
//NUMRICAL PARAMETERS
const int N = 5; //Number of control volumes
double delta = 1E-6; //Convergence criteria
double Tstart = 30; //Initial temp map (supposed equal in all geometry) [Cº]

// DEFFINE VECTOR LENGTH AND OPERATION VARIABLES
double T[N+2] = {0};
double rw[N+1] = {0}; //fronteres entre volumns de control i extrems
double rp[N+2] = {0};
double a_E[N+2] = {0};
double a_W[N+2] = {0};
double a_P[N+2] = {0};
double b_P[N+2] = {0};
double Ap[N] = {0};
double S[N+2] = {0};

double dpw = 0; // malla uniforme
double dpe = 0;
double deltaR = 0;




class vec { //Vector deffinition class
    public:

    vec(){
        deltaR = (Rext-Rint)/N;
        dpw = deltaR/2; // per malla uniforme
        dpe = deltaR/2;

        for (int i = 0; i<N+1;i++){ //Definició de el vector r_w
            rw[i] = Rint + (i)*deltaR;
            //cout << rw[i] << "\n";
        }
        rp[0] = Rint; //Definició components inicials del vector r_w
        rp[N+1] = Rext;

        for (int i = 0; i<N+1;i++){ // definició de les arees laterals (parets entre Volumns de control)
            S[i] = rw[i] * 2 *M_PI * ef ;
            //cout << S[i] << "\n";
        }
        for (int i = 0; i<N;i++){ //definició area superior i inferior dels volumns de control
            Ap[i] = 2 * M_PI * pow(rw[i+1],2)-pow(rw[i],2);
            cout << Ap[i] << "\n";
        }
    }
};

class Mapa_inicial { //Initial temp map deffinition class
public:

    Mapa_inicial(){ //Definim el mapa inicial de temperatures
        for(int i = 0;i<(N+2);i++){
            T[i] = Tstart;
        }
    }
};

class Coefs_discrtitzacio { //Coef deffinition class
public:

    Coefs_discrtitzacio(){
        //Primer node
        a_W[0] = 0;
        a_E[0] = 0;
        a_P[0] = 1;
        b_P[0] = T_wall;
        //Ultim node
        a_W[N+1] = lambda_f/dpw;
        a_E[N+1] = 0;
        a_P[N+1] = a_W[N+1] + alfa_extrem;
        b_P[N+1] = alfa_extrem*Text;
        //Nodes intermitjos
        for (int i = 1; i< (N+1);i++){
            a_W[i] = lambda_f * S[i-1] /dpw;
            a_E[i] = lambda_f * S[i] /dpe;
            b_P[i] = alfa_exterior * Text * Ap[i-1];
            a_P[i] = a_E[i] + a_W[i] + alfa_exterior * Ap[i-1];
        }
    }
};

class Solver { //Initial temp map deffinition class
public:

    void gauss_seidel(){ //Definim el mapa inicial de temperatures

        }
    }
    void TDMA (){ }
};


int main() {
    //for(int i = 0; i<sizeof(T);i++) {

    vec  v;
    Mapa_inicial map;
    Coefs_discrtitzacio coef;
    Solver solver;


    return 0;
}
