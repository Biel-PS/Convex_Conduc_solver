#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
using namespace std;
//CASE OF A CYLINDER WITH nf FINS
//Declaration of global variables

//FISICAL PARAMETERS
const double Rint = 15E-3; // Internal radius of the fin
const double Rext = 75E-3; // External radius of the fin
const double ef = 4E-3; // Fin's vertical length [m]
const double eb = 5E-3; // Distance between fins [m]
const double Text = 25; // Temp at exterior [Cº]
const double T_wall = 72.13; // Temp of the wall r = Rint [Cº]
//double lambda_f[2] = [200,10]; // lambda of the material (homogeneous case) [a_nT^n,...a_1T,a_0] (function of T)
const double lambda_f = 58;
const double alfa_exterior = 110; //alpha in the length of the fin with outter fluid
const double alfa_extrem = 0; //alpha in the end of fin
const double nf = 163; //Number of fins
//NUMRICAL PARAMETERS
const int N = 300; //Number of control volumes
const double delta = 1E-10; //Convergence criteria
const double Tstart = 10; //Initial temp map (supposed equal in all geometry) [Cº]


// DEFFINE VECTOR LENGTH AND OPERATION VARIABLES
double T[N+2] = {0};
double T_n[N+2] = {0};
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
double Q_f = 0;
double Q_b = 0;
double Q_base_fin = 0;




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
            Ap[i] = 2 * M_PI * (pow(rw[i+1],2)-pow(rw[i],2));
            //cout << Ap[i] << "\n";
        }
    }
};

class Mapa_inicial { //Initial temp map deffinition class
public:

    Mapa_inicial(){ //Definim el mapa inicial de temperatures
        for(int i = 0;i<(N+2);i++){
            T[i] = Tstart;
            T_n[i] = T[i] + 10; // Ens assegurem que per la primera iteració no es compleix el criteri d'acceptació
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
        int cont = 0;
        while (true){
            if(norm(T_n,T,delta)){
                cout << cont << "\n";
                break;
            }
            else{
                for (int i = 0; i<(N+2); i++){
                    T[i] = T_n[i];
                }
                T_n[0] = (a_E[0]*T[1] + b_P [0])/a_P[0];
                for (int i = 1; i<N+1;i++){
                    T_n[i] = (a_E[i]*T[i+1] + a_W[i]*T[i-1] + b_P [i])/a_P[i];
                }
                T_n[N+1] = (a_W[N+1]*T[N] + b_P [N+1])/a_P[N+1];
            }
        cont++;}


    }

    void TDMA (){ }
    void Qfin (){
        Q_f = 0;
        for (int i = 1 ; i<N+1;i++){
            Q_f += alfa_exterior*(T[i]-Text)*Ap[i-1];
        }
        Q_f += alfa_extrem*(T[N+1] - Text)*2*M_PI *Rext * ef;
       // Q_f *= nf;
        Q_base_fin  = - lambda_f* ((T[1]-T[0])/dpe) * 2 * M_PI * Rint * ef;

    }
    void Qbody (){
        Q_b  = alfa_exterior*(T_wall - Text)*2 * M_PI * Rint * eb;
    }
    static bool norm (double T_1[N+2],double T_2[N+2],double delta){
        bool done = true;
        for (int i = 0; i< N+2;i++){
             if (abs(T_1[i]-T_2[i]) > delta){
                 done = false;
                 break;
             }
        }
        return sqrt(done);
    }
};


int main() {
    //for(int i = 0; i<sizeof(T);i++) {

    vec  v;
    Mapa_inicial map;
    Coefs_discrtitzacio coef;
    Solver solver;
    solver.gauss_seidel();
    for (int i = 0; i<(N+2); i++){
        cout << T[i] << "\n";
    }
    solver.Qfin();
    solver.Qbody();

    cout <<" Q_fin: " << Q_f <<"\n Q_body: "<<Q_b <<"\n Q_base_fin: "<<Q_base_fin;


    return 0;
}
