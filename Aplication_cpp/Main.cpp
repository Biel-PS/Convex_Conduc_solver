#include <iostream>
using namespace std;
//CASE OF A CYLINDER WITH nf FINS
//Declaration of global variables

//FISICAL PARAMETERS
double Rint = 0; // Internal radius of the fin
double Rext = 0; // External radius of the fin
double ef = 0; // Fin's vertical length [m]
double eb = 0; // Distance between fins [m]
double Text = 10; // Temp at exterior [Cº]
double T_wall = 100; // Temp of the wall r = Rint [Cº]
//double lambda_f[2] = [200,10]; // lambda of the material (homogeneous case) [a_nT^n,...a_1T,a_0] (function of T)
double alfa_exterior = 0; //alpha in the length of the fin with outter fluid
double alfa_extrem = 0; //alpha in the end of fin
double nf = 5; //Number of fins
//NUMRICAL PARAMETERS
const int N = 10; //Number of control volumes
double delta = 1E-6; //Convergence criteria
double Tstart = 30; //Initial temp map (supposed equal in all geometry) [Cº]
double T[N+2] = {0};

class vec { //Vector deffinition class
    public:
    vec(int N){
        for (int i = 0; i<sizeof(T);i++) {
            T[i] = 0;
        }
    }
};



int main() {
    //for(int i = 0; i<sizeof(T);i++) {
        cout << T[3];
   // }



    return 0;
}