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
    const int M = 10; // Number of divisions in horizontal axis (0umns)
    const int Num_materials = 2; // Specify the number of materials in the grid
    //NOTE: the boundary_material_coordinates only works for rectangular distributions of material within the control surface.
    const double boundary_material_coordinates[Num_materials][4] = {{0,int(N+1)/2,0,(M+1)},{int((N+1)/2) +1,(N+1),0,M+1}} ; // Specify the boundary of each material as [material][0 --> start and 1--> finish row, 2-->start and 3-->finish 0umn]

    //NOTE: THE START AND END ROW/0UMN IN THIS VECTOR IS CONSIDERED OF THE MATERIAL
    //THAT IS SPECIFIED, ROW 1 TO 2 MEANS THAT NODES IN ROW 1 AND 2 ARE OF THE SELECTED MATERIAL!!!Ç

    //REMEMBER: there are N+1 rows and M+1 0umns because of the wall nodes (zero is the beginning)!!

    const double rho[Num_materials] = {100,200}; //rho is suposed constant with T and t
    const double Cp[Num_materials] = {50,80}; //Cp is suposed constant with T and t
    const double qv_p[Num_materials] = {10e6,50e6};
    
    //Temporal and convergence parameters
    const double delta_convergence = 1E-9; //Convergence criteria
    const double delta_t = 0.001; // Time increment [s]
    const double t_init = 0; // initial time [s]
    const double t_end = 10; // end time [s]
    const double Beta = 0.5; // =0 explícit, =0.5 Charles-Nicholson, = 1 Implícit
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
    const double lambda_f[Num_materials][max_lambda_degree + 1] = {{50},{80}}; //lambda of each material [material][a_0,a_1*T,a_2*T^2,...] (function of T)

//Control variables
    bool excep = false; //if this variable becomes true at any point, the program stops.
    bool first_time_Tmap = true;
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
    double T[2][N+2][M+2] = {0}; //Temperature [0 --> past delta time, 1--> actual delta time][row][0umn]
    double T_estimada[2][N+2][M+2] = {0};
    double Q_p[2][N+2][M+2] = {0};    
   
    
    /* double a[N+2][M+2] = {0};
    double bp[N+2][M+2] = {0};*/

    //Surfaces and distance between nodes/surfaces (uniform mesh)
    double S_h = 0;
    double S_v = 0;
    double V_p = 0;
    double dpv = 0; // distance between vertical nodes and surfaces
    double dph = 0; // distance between horizontal nodes and surfaces

//METHODS THAT WILL BE USED

double harmonic_mean (double x, double y,double dxy,double dx1,double dy1){
    return dxy/(dx1/x + dy1/y); //method that calculates the harmonic mean of two values.
}

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
    V_p = dpv*dph*W;
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
                // set for the correct range the value of the index of the corresponding material
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
                    cout << "ERROR 01 :boundary_material_coordinates has an empty spot, the boundaries may be wrongly defined ! \n";
                excep = true;
                break;
            }
       }
   }
}

static void Mapa_inicial (){
    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < M + 2; j++) {
            T[0][N][M] = Tstart;
        }
    }
}
static void Mapa_estimat (){
    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < M + 2; j++) {
            T_estimada[1][N][M] = T[0][N][M];
            Q_p[0][i][j] = Q_p[1][i][j];
        }
    }
}

static void solver_gauss_seidel (){
    // Define method variables
    double dpv_half = dpv/2,dph_half = dph/2;

    // Defne the coefficients
    double a_N = 0, a_E = 0, a_S = 0, a_W = 0, a_p = 0,b_p = 0;
    //lambdas at surface between nodes
    double lam_n = 0, lam_e = 0, lam_s = 0, lam_w = 0;
    //lambdas at the contiguous nodes
    double lam_p = 0,lam_E = 0, lam_W = 0, lam_N = 0, lam_S = 0;
    
    //Compute coefficient for boundary nodes at EAST AND WEST sides
    for (int i = 1; i < N + 1; i++) {
            //Compute lambda at the WEST BOUNDARY
            lam_p = 0,lam_E = 0, lam_W = 0, lam_N = 0, lam_S = 0;
            for (int k = 0; k< max_lambda_degree+1; k++){
                lam_p += lambda_f[Material_matrix[i][0]][k] + pow(T[1][i][0],k);
                lam_E += lambda_f[Material_matrix[i][1]][k] + pow(T[1][i][1],k);
            }
            lam_e = harmonic_mean(lam_p,lam_E,dph,dph_half,dph_half);

            a_E = lam_e/dph;
            a_p = a_E + alfa_w;
            b_p = alfa_w * Twest;

            //We compute the temperature
            T[1][i][0] = (a_E* T[1][i][1] +  b_p)/a_p;

            //COMPUTE LAMNDA AT THE EAST BOUNDARY

        lam_p = 0,lam_E = 0, lam_W = 0, lam_N = 0, lam_S = 0;
        for (int k = 0; k< max_lambda_degree+1; k++){
            lam_p += lambda_f[Material_matrix[i][M+1]][k] + pow(T[1][i][M+1],k);
            lam_W += lambda_f[Material_matrix[i][M]][k] + pow(T[1][i][M],k);
        }
        lam_w = harmonic_mean(lam_p,lam_W,dph,dph_half,dph_half);

        a_W = lam_w/dph;
        a_p = a_W + alfa_e;
        b_p = alfa_e * Teast;

        //We compute the temperature
        T[1][i][M+1] = (a_W* T[1][i][M] +  b_p)/a_p;
    }
    
    // Compute the coefficients at the NORTH and South side

    for (int j = 1; j < M+1; j++) {
        //Compute lambda at the NORTH BOUNDARY
        lam_p = 0,lam_E = 0, lam_W = 0, lam_N = 0, lam_S = 0;
        for (int k = 0; k< max_lambda_degree+1; k++){
            lam_p += lambda_f[Material_matrix[M+1][j]][k] + pow(T[1][M+1][j],k);
            lam_S += lambda_f[Material_matrix[M][j]][k] + pow(T[1][M][j],k);
        }
        lam_s = harmonic_mean(lam_p,lam_S,dpv,dpv_half,dpv_half);

        a_S = lam_s/dpv;
        a_p = a_S + alfa_n;
        b_p = alfa_n * Tnorth;

        //We compute the temperature
        T[1][M+1][j] = (a_S* T[1][M][j] +  b_p)/a_p;

        //COMPUTE LAMNDA AT THE SOUTH BOUNDARY

        lam_p = 0,lam_E = 0, lam_W = 0, lam_N = 0, lam_S = 0;
        for (int k = 0; k< max_lambda_degree+1; k++){
            lam_p += lambda_f[Material_matrix[0][j]][k] + pow(T[1][0][j],k);
            lam_N += lambda_f[Material_matrix[1][j]][k] + pow(T[1][1][j],k);
        }
        lam_n = harmonic_mean(lam_p,lam_N,dpv,dph_half,dph_half);

        a_N = lam_n/dpv;
        a_p = a_N + alfa_s;
        b_p = alfa_s * Tsouth;

        //We compute the temperature
        T[1][0][j] = (a_N* T[1][1][j] +  b_p)/a_p;
    }
    
    // Compute coefficients for the internal nodes
    for (int i = 1; i < N + 1; i++) {
        for (int j = 1; j < M + 1; j++) {
            //Compute lambda at the node
            lam_p = 0,lam_E = 0, lam_W = 0, lam_N = 0, lam_S = 0;
            for (int k = 0; k< max_lambda_degree+1; k++){
              lam_p += lambda_f[Material_matrix[i][j]][k] + pow(T[1][i][j],k);
              lam_E += lambda_f[Material_matrix[i][j+1]][k] + pow(T[1][i][j+1],k);
              lam_N += lambda_f[Material_matrix[i+1][j]][k] + pow(T[1][i+1][j],k);
              lam_W += lambda_f[Material_matrix[i][j-1]][k] + pow(T[1][i][j-1],k);
              lam_S += lambda_f[Material_matrix[i-1][j]][k] + pow(T[1][i-1][j],k);
            }
            lam_n = harmonic_mean(lam_p,lam_N,dpv,dpv_half,dpv_half);
            lam_e = harmonic_mean(lam_p,lam_E,dph,dph_half,dph_half);
            lam_s = harmonic_mean(lam_p,lam_N,dpv,dpv_half,dpv_half);
            lam_w = harmonic_mean(lam_p,lam_W,dph,dph_half,dph_half);

            a_E = Beta * lam_e * S_h/dph;
            a_S = Beta * lam_s * S_v/dpv;
            a_W = Beta * lam_w * S_h/dph;
            a_N = Beta * lam_n * S_v/dpv;
            
            a_p = a_E + a_W +a_N + a_S + rho[Material_matrix[i][j]] * V_p * Cp[Material_matrix[i][j]] / delta_t;
            b_p = qv_p[Material_matrix[i][j]] * V_p + rho[Material_matrix[i][j]] * Cp[Material_matrix[i][j]] * T[1][i][j] / delta_t + (1+Beta) * Q_p[0][i][j];
            Q_p[1][i][j] = ((-1*(T[0][i][j-1] - T[0][i][j]))*a_W - (-1*(T[0][i][j+1] - T[0][i][j]))*a_E +
                    (-1*(T[0][i-1][j] - T[0][i][j])) * a_S + ((-1*(T[0][i+1][j] - T[0][i][j]))*a_N ))/Beta + qv_p[Material_matrix[i][j]] * V_p;
            
            //We compute the temperature
            
            T[1][i][j] = (a_E* T[1][i][j+1] + a_W* T[1][i][j-1] + a_S * T[1][i-1][j] +  a_N*T[1][i+1][j]+ b_p)/a_p;
        }
    }
}


static void Calc_coeff(){

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

    vec_geometric_deff();
    Mapa_inicial();
    Mapa_estimat();
    solver_gauss_seidel();

    return 0;
}