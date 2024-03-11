#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;
//CASE OF A CYLINDER WITH nf FINS
//Declaration of global variables

//FISICAL PARAMETERS
const double Rint = 1; // Internal radius of the fin
const double Rext = 3; // External radius of the fin
const double ef = 4E-3; // Fin's vertical length [m]
const double eb = 5E-3; // Distance between fins [m]
const double Text = 25; // Temp at exterior [Cº]
const double T_wall = 200; // Temp of the wall r = Rint [Cº]
//double lambda_f[2] = [200,10]; // lambda of the material (homogeneous case) [a_nT^n,...a_1T,a_0] (function of T)
const double lambda_f = 100;
const double alfa_exterior = 200; //alpha in the length of the fin with outter fluid
const double alfa_extrem = 50; //alpha in the end of fin
const double nf = 163; //Number of fins
//NUMRICAL PARAMETERS
const int N = 500; //Number of control volumes
const double delta = 1E-5; //Convergence criteria
const double Tstart = 100; //Initial temp map (supposed equal in all geometry) [Cº]


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
double S[N+1] = {0};


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
           // cout << rw[i] << "\n";
        }
        rp[0] = rw[0] ;
        rp[N+1] = rw[N];
        for (int i = 1; i<N+1;i++){ //Definició de el vector r_w
            rp[i] = (rw[i] + rw[i-1])/2;
        }
       /* for(double & i : rp){
            cout << i << "\n";
        }*/

        for (int i = 0; i<N+1;i++){ // definició de les arees laterals (parets entre Volumns de control)
            S[i] = rw[i] * 2 *M_PI * ef ;
          // cout << S[i] << "\n";
        }
        for (int i = 0; i<N;i++){ //definició area superior i inferior dels volumns de control
            Ap[i] = 2 * M_PI * (pow(rw[i+1],2.0)-pow(rw[i],2.0));
           // cout << Ap[i] << "\n";
        }
    }
};

class Mapa_inicial { //Initial temp map deffinition class
public:

    Mapa_inicial(){ //Definim el mapa inicial de temperatures
        for(double & i : T){
            i = Tstart;
            //cout << i <<"\n";
          //  T_n[i] = T[i] + 10; // Ens assegurem que per la primera iteració no es compleix el criteri d'acceptació
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

    static void gauss_seidel(){ //Definim el mapa inicial de temperatures
        int cont = 0;
        while (true){
            if(norma(T_n,T)){
                //cout << cont << "\n";
                break;
            }
            else{
                for (int i = 0; i<(N+2); i++){
                    T[i] = T_n[i];
                }
                T_n[0] = (a_E[0]*T[1] + b_P [0])/a_P[0];
                T_n[N+1] = (a_W[N+1]*T[N] + b_P[N+1])/a_P[N+1];
                for (int i = 1; i< N+1; i++){
                    T_n[i] = (a_E[i]*T[i+1] + a_W[i]*T[i-1] + b_P[i])/a_P[i];
                }
            }
            cont++;
        }
    }

    void TDMA (){ }
    static void Qfin (){
        Q_f = 0;

        for (int i = 1 ; i<N+1;i++){
            Q_f += alfa_exterior*(T_n[i]-Text)*Ap[i-1];
        }
        Q_f += alfa_extrem * (T_n[N+1] - Text) * 2 * M_PI * Rext * ef;
        // Q_f *= nf;
        Q_base_fin  = - lambda_f*((T_n[1]-T_n[0])/dpe) * (2 * M_PI * Rint * ef);

    }
    static void Qbody (){
        Q_b  = alfa_exterior*(T_wall - Text)* 2 * M_PI * Rint * eb;
    }
    static bool norma (double T_1[N+2],double T_2[N+2]){
        bool done = true;
        double T_cache = 0;
        for (int i = 0; i< N+2;i++){
            T_cache +=pow((T_1[i]-T_2[i]),2.0);
            }
        cout << sqrt(T_cache) << "\n";
        if (sqrt(T_cache) >= delta){
            done = false;
        }
        return done;
    }
};


int main() {
    //for(int i = 0; i<sizeof(T);i++) {

    vec  v;
    Mapa_inicial map;
    Coefs_discrtitzacio coef;
    Solver solver;
    Solver::gauss_seidel();

    for (double i : T){
     // cout << i << "\n";
    }

    Solver::Qfin();
    Solver::Qbody();

    cout <<" Q_fin: " << Q_f <<"\n Q_body: "<<Q_b <<"\n Q_base_fin: "<<Q_base_fin << "\n Error : "<<abs(Q_base_fin - Q_f) ;


    return 0;
}