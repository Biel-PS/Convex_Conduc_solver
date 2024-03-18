//define _USE_MATH_DEFINES;

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;
// CASE FOR A 2D SECTION (EXAMPLE FROM LESON) by BIEL PUJADAS SURIOL, started: 15/03/2024
// This case its being dessigned as a square that can be discretized and where the temperatures at the top, bottom, left
// and right are all known data.

//NUMRICAL PARAMETERS
    //Mesh parameters
    const int N = 10; //Number of divisions in vertical axis (rows)
    const int M = 10; // Number of divisions in horizontal axis (columns)
    const int Num_materials = 2; // Specify the number of materials in the grid
    //NOTE: the boundary_material_coordinates only works for rectangular distributions of material within the control surface.
    const double boundary_material_coordinates[Num_materials][4] = {{0,int(N+1)/2,0,(M+1)},{int((N+1)/2) +1,(N+1),0,M+1}} ; // Specify the boundary of each material as [material][0 --> start and 1--> finish row, 2-->start and 3-->finish column]

    //NOTE: THE START AND END ROW/COLUMN IN THIS VECTOR IS CONSIDERED OF THE MATERIAL
    //THAT IS SPECIFIED, ROW 1 TO 2 MEANS THAT NODES IN ROW 1 AND 2 ARE OF THE SELECTED MATERIAL!!!Ç

    //REMEMBER: there are N+1 rows and M+1 columns because of the wall nodes (zero is the beginning)!!

    const double rho[Num_materials] = {0}; //rho is suposed constant with T and t
    const double Cp[Num_materials] = {0}; //Cp is suposed constant with T and t

    //Temporal and convergence parameters
    const double delta_convergence = 1E-9; //Convergence criteria
    const double delta_t = 0.001; // Time increment [s]
    const double t_init = 0; // initial time [s]
    const double t_end = 10; // end time [s]
    const double beta = 0.5; // =0 explícit, =0.5 Charles-Nicholson, = 1 Implícit
    const double relaxation = 1;

//FISICAL PARAMETERS
    //Geometrical parameters
    const double H = 2; // Height of the square plane
    const double L = 2; // Length of the square plane
    const double W = 2; // Profundity of the square plane in case a 3d case with 2d heat transfer is wanted

    //External convection temperatures in the boundary
    const double Tnorth = 200; // Temperature in the north wall
    const double Tsouth = 100;// Temperature in the south wall
    const double Teast = 50; // Temperature in the east wall
    const double Twest = 80; // Temeprature in the west wall
    //External convection constants in the boundary
    const double alfa_n = 0; //north wall
    const double alfa_s = 0; //south wall
    const double alfa_e = 100; //east wall
    const double alfa_w = 100; //west wall

    //Initial temperature map
    //const double initial_temp_map[] = []; // IN CASE THE INITIAL MAP IS DISCRETIZED, IT WOULD BE PUT IN HERE
    const double Tstart = 100; // INITIAL MAP IN CASE ALL THE NODES AT THE SAME TEMP AT T = 0;

    //Deffinition of the lambda parameter in heach material
    const int max_lambda_degree = 0; // degree of the max polynomial lambda (if its a constant: = 0), if its a line = 1, etc
    const double lambda_f[Num_materials][max_lambda_degree + 1] = {0}; //lambda of each material [material][a_0,a_1*T,a_2*T^2,...] (function of T)

//Control variables
    bool excep = false; //if this variable becomes true at any point, the program stops.

// DEFFINE VECTOR LENGTH AND OPERATION VARIABLES
//NOTE: THE ORIGIN OF COORDINATES IS CONSIDERED IN THE INTERSECCION OF THE W AND S WALLS (BOTTOM LEFT VERTICE)!!!0
    //Postion of heach control volume (ONLY central node)
    double x_p [N][M] = {0};
    double y_p [N][M] = {0};
    int Material_matrix[N+2][M+2] = {0}; //assign the index of every material to each node

    //Position of every node
    double x_all [N+2][M+2] = {0};
    double y_all [N+2][M+2] = {0};

    //Temperature and needed coefficients
    double T[2][N+2][M+2] = {0}; //Temperature [0 --> past delta time, 1--> actual delta time][row][column]
    double a[N+2][M+2] = {0};
    double bp[N+2][M+2] = {0};

    //Surfaces and distance between nodes/surfaces (uniform mesh)
    double S_h = 0;
    double S_v = 0;
    double dpv = 0; // distance between vertical nodes and surfaces
    double dph = 0; // distance between horizontal nodes and surfaces

//METHODS THAT WILL BE USED

static void vec_geometric_deff (){ //initial vector deffinition method
    //Compute the coordinates of every node

    dph = (L/(2*M));
    dpv = (H/(2*N));

    for (int i = 0; i<N+2; i++){
        x_all[i][0] = 0;
        x_all[i][1] = x_all[i][0] + dph;
        for (int j = 2; j<M+1; j++){
            x_all[i][j] = x_all[i][j-1] + 2*dph;
        }
        x_all[i][M+1] = L;
    }
    for (int j = 0; j<M+2; j++){
        y_all[0][j] = 0;
        y_all[1][j] = y_all[0][j] + dpv;
        for (int i = 2; i<N+1; i++){
            y_all[i][j] = y_all[i-1][j] + 2*dpv;
        }
        y_all[N+1][j] = H;
    }
    /*  for (int i = 0; i<N+2; i++){
       for (int j = 0; j<M+2; j++){
           cout << "X: " << x_all[i][j] << "      Y: " <<y_all[i][j] <<"\n";
       }
   }*/

    //Compute the surface of every n,s and w,e between control volumes.
    S_h = 2*dph*W;
    S_v = 2*dpv*W;
   // cout << S_h << " " << S_v;

   //Define the empty material matrix

    for (int i = 0;i<N+2 ; i++){
        for (int j = 0; j<M+2 ; j++){
            Material_matrix[i][j] = -1;
        }
    }

   //Compute the material matrix, where every node is assigned a material
   for (int k = 0; k < Num_materials; k++ ) {
       for (int i = 0; i < N + 2; i++) {
           for (int j = 0; j < M + 2; j++) {

               if(boundary_material_coordinates[k][0] <= i && boundary_material_coordinates[k][1] >= i &&
               boundary_material_coordinates[k][2] <= j && boundary_material_coordinates[k][3] >= j)
               {
                   Material_matrix[i][j] = k;
               }
           }
       }
   }
   //Print the matrix of material (JUST FOR DEBUGGING)
    /*for (int i = N+1; i>=0; i--){
        for (int j = 0; j<M+2; j++){
            cout << Material_matrix[i][j] << " ";
        }
        cout << "\n";
    }*/

   //CONTROL THAT ALL THE MATERIAL TYPES HAVE BEEN ASSIGNED
   for (int i = 0;i<N+2 ; i++){
       if(excep){
           break;
       }
       for (int j = 0; j<M+2 ; j++){
            if(Material_matrix[i][j] == -1){
                cout << "boundary_material_coordinates has an empty spot, the boundaries may be wrongly defined >:( \n";
                excep = true;
                break;
            }
       }
   }




}

static void Calc_coeff(){

}

static double harmonic_mean (double x, double y,double dxy,double dx1,double dy1){
    return dxy/(dx1/x + dy1/y); //method that calculates the harmonic mean of two values.
}

/*
class Mapa_inicial { //Initial temp map deffinition class
public:

    Mapa_inicial(){ //Definim el mapa inicial de temperatures
        for(double & i : T_n){
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

    static void TDMA (){
        P[0] = a_E[0]/a_P[0];
        R[0] = b_P[0]/a_P[0];
        for (int i = 1; i<N+1;i++){
            P[i] = (a_E[i])/(a_P[i] - a_W[i] * P[i-1]);
            R[i] = (b_P[i] + a_W[i] * R[i-1])/(a_P[i] - a_W[i]*P[i-1]);
        }
        P[N+1] = 0;
        R[N+1] = (b_P[N+1] + a_W[N+1] * R[N])/(a_P[N+1] - a_W[N+1]*P[N]);

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
                T_n[0] = P[0]*T[1] + R[0];
                for (int i = 1; i<N+1;i++){
                    T_n[i] = P[i]*T[i+1] + R[i];
                }
                T_n[N+1] = R[N+1];
            }
            cont++;
        }
    }
    static void Qfin (){
        Q_f = 0;
        Q_base_fin = 0;

        for (int i = 1 ; i<N+1;i++){
            Q_f += alfa_exterior*(T_n[i]-Text)*Ap[i-1];
        }
        Q_f += alfa_extrem * (T_n[N+1] - Text) * 2 * M_PI * Rext * ef;
        // Q_f *= nf;
        Q_base_fin  = - lambda_f*((T_n[1]-T_n[0])/dpe) * (2 * M_PI * Rint * ef) ;


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
        //cout << sqrt(T_cache) << "\n";
        if (sqrt(T_cache) > delta){
            done = false;
        }
        return done;
    }
};
*/

int main() {
    //for(int i = 0; i<sizeof(T);i++) {
/*
    vec  v;
    Mapa_inicial map;
    Coefs_discrtitzacio coef;
    Solver solver;
    Solver::TDMA();

    for (double i : T){
        cout << i << "\n";
    }

    Solver::Qfin();
    Solver::Qbody();

    cout <<" Q_fin: " << Q_f <<"\n Q_body: "<<Q_b <<"\n Q_base_fin: "<<Q_base_fin << "\n Error : "<<abs(Q_base_fin - Q_f) ;
    */
    vec_geometric_deff();

    return 0;
}