// CASE FOR A 2D HEAT TRANSFER PROBLEM WITH 4 MATERIALS by BIEL PUJADAS SURIOL, started: 15/03/2024

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
using namespace std;
// This case its being dessigned as a square that can be discretized and where the temperatures at the top, bottom, left
// and right are all known data.

//NUMRICAL PARAMETERS
//Mesh parameters
const int N =120; //Number of divisions in vertical axis (rows)
const int M =120;// Number of divisions in horizontal axis (Columns)
const int Num_materials = 4; // Specify the number of materials in the grid

const double p_1[2] = {0.50,0.40};
const double p_2[2] = {0.50,0.70};
const double p_3[2] = {1.10,0.80};

//Geometrical parameters
const double H = 0.8; // Height of the square plane
const double L = 1.10; // Length of the square plane
const double W = 1; // Profundity of the square plane in case a 3d case with 2d heat transfer is wanted

//REMEMBER: there are N+1 rows and M+1 columns because of the wall nodes (zero is the beginning)!!


const double rho[Num_materials] = {1500,1600,1900,2500}; //rho is suposed constant with T and t
const double Cp[Num_materials] = {750,770,810,930}; //Cp is suposed constant with T and t
const double qv_p[Num_materials] = {0,0,0,0};

//Temporal and convergence parameters
const double delta_convergence = 1E-8; //Convergence criteria
const double delta_t = 0.1; // Time increment [s]
const double t_init = 0; // initial time [s]
double t_actual = t_init;
const double t_end = 5000; // end time [s]
const double Beta =0.5; // =0 explícit, =0.5 Charles-Nicholson, = 1 Implícit
const double relaxation = 1.05;
const int ratio_print_t = round(50/delta_t);
//FISICAL PARAMETERS

//External convection temperatures in the boundary
const double Tnorth = 0; // Temperature in the north wall
const double Tsouth = 23;// Temperature in the south wall
const double Teast = 0; // Temperature in the east wall
const double Twest = 33; // Temeprature in the west wall
//External convection constants in the boundary
const double alfa_n = 0; //north wall
//const double alfa_s = 0; //south wall
const double alfa_e = 0; //east wall
const double alfa_w = 9; //west wall

const double q_w = 54.55; //[W/m2]

//Initial temperature map
//const double initial_temp_map[] = []; // IN CASE THE INITIAL MAP IS DISCRETIZED, IT WOULD BE PUT IN HERE
const double Tstart = 8.0; // INITIAL MAP IN CASE ALL THE NODES AT THE SAME TEMP AT T = 0;

//Deffinition of the lambda parameter in heach material
const int max_lambda_degree = 0; // degree of the max polynomial lambda (if its a constant: = 0), if its a line = 1, etc
const double lambda_f[Num_materials][max_lambda_degree + 1] = {{170},{140},{200},{140}}; //lambda of each material [material][a_0,a_1*T,a_2*T^2,...] (function of T)


double lambda_vector [N+2][M+2][4] = {0};
//for each wall 0--> north, 1--> east, 2--> west, 3--> south
//Control variables
bool excep = false; //if this variable becomes true at any point, the program stops.

// DEFFINE VECTOR LENGTH AND OPERATION VARIABLES
//NOTE: THE ORIGIN OF COORDINATES IS CONSIDERED IN THE INTERSECCION OF THE W AND S WALLS (BOTTOM LEFT VERTICE)!!!0


int Material_matrix[N+2][M+2] = {-1}; //assign the index of every material to each node

//Position of every node
double x_all [N+2][M+2] = {0};
double y_all [N+2][M+2] = {0};

//Temperature and needed coefficients
double T[2][N+2][M+2] = {0}; //Temperature [0 --> past delta time, 1--> actual delta time][row][0umn]
double T_cache[2][N+2][M+2] = {0};

//Surfaces and distance between nodes/surfaces (uniform mesh)
double S_h = 0;
double S_v = 0;
double V_p = 0;
const double dpv = (H/(N)); // distance between vertical nodes and surfaces
const double dph = (L/(M)); // distance between horizontal nodes and surfaces
const double dph_half = dph/2;
const double dpv_half = dpv/2;

double dn = 0;
double ds = 0;
double de = 0;
double dw = 0;

double Heat_parets [4][int(t_end/delta_t)] = {0};//for each wall 0--> north, 1--> east, 2--> west, 3--> south


//METHODS THAT WILL BE USED

double harmonic_mean (double x, double y,double dxy,double dx1,double dy1){
    return dxy/(dx1/x + dy1/y); //method that calculates the harmonic mean of two values.
}

static bool norma (double T_1[2][N+2][M+2],double T_2[2][N+2][M+2]){
    bool done = true;
    double T_mem = 0.0;
    for (int i = 0; i< N+2;i++){
        for (int j = 0; j<M+2; j++){
            if (T_mem <  abs(T_1[0][i][j]-T_2[1][i][j])){
                T_mem = abs(T_1[0][i][j]-T_2[1][i][j]);
            }
        }
    }
    if (T_mem > delta_convergence){
        done = false;
    }
    return done;
}



static void vec_geometric_deff (){ //initial vector deffinition method
    //Compute the coordinates of every node

    for (int i = 0; i<N+2; i++){
        x_all[i][0] = 0;
        x_all[i][1] = x_all[i][0] + dph_half;
        for (int j = 2; j<M+1; j++){
            x_all[i][j] = x_all[i][j-1] + dph;
        }
        x_all[i][M+1] = L;
    }
    for (int j = 0; j<M+2; j++){
        y_all[0][j] = 0;
        y_all[1][j] = y_all[0][j] + dpv_half;
        for (int i = 2; i<N+1; i++){
            y_all[i][j] = y_all[i-1][j] + dpv;
        }
        y_all[N+1][j] = H;
    }
    //Compute the surface of every n,s and w,e between control volumes.
    S_h = dph*W;
    S_v = dpv*W;
    V_p = dpv*dph*W;

    //Define the empty material matrix
    for (int i = 0;i<N+2 ; i++){
        for (int j = 0; j<M+2 ; j++){
            Material_matrix[i][j] = -1;
        }
    }
    //Compute the material matrix, where every node is assigned a material
    for (int i = 1; i < N + 1; i++) {
        for (int j = 1; j < M+1; j++) {
            if (x_all[i][j] <= p_1[0] && y_all[i][j] < p_1[1]) {
                Material_matrix[i][j] = 0;
            }
            else if (x_all[i][j] < p_1[0] && y_all[i][j] >= p_1[1]) {
                Material_matrix[i][j] = 2;
            }
            else if (x_all[i][j] > p_1[0] && y_all[i][j] <= p_2[1]) {
                Material_matrix[i][j] = 1;
            }
            else if (x_all[i][j] > p_1[0] && y_all[i][j] > p_2[1]) {
                Material_matrix[i][j] = 3;
            }
        }
    }

    for(int i = 0; i<N+2;i++){
        Material_matrix[i][0] = Material_matrix[i][1];
        Material_matrix[i][M+1] = Material_matrix[i][M];
    }
    for(int j = 0; j<M+2;j++){
        Material_matrix[0][j] = Material_matrix[1][j];
        Material_matrix[N+1][j] = Material_matrix[N][j];
    }

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

    //CALCULATE THE LAMBDA IN EACH BOUNDARY SURFACE
    //lambdas at the contiguous nodes
    double lam_p = 0,lam_E = 0, lam_W = 0, lam_N = 0, lam_S = 0;

    for (int i = 1; i < N + 1; i++) {
        for (int j = 1; j < M + 1; j++) {
            //At the internal control volumes
            dn = y_all[i+1][j]-y_all[i][j];
            ds = y_all[i][j]-y_all[i-1][j];
            de = x_all[i][j+1]-x_all[i][j];
            dw = x_all[i][j]-x_all[i][j-1];

            lam_p = lambda_f[Material_matrix[i][j]][0];
            lam_E = lambda_f[Material_matrix[i][j + 1]][0];
            lam_N = lambda_f[Material_matrix[i + 1][j]][0];
            lam_W = lambda_f[Material_matrix[i][j - 1]][0];
            lam_S = lambda_f[Material_matrix[i - 1][j]][0];

            lambda_vector[i][j][0] = harmonic_mean(lam_p, lam_N, dn, dn/2, dn/2);
            lambda_vector[i][j][1] = harmonic_mean(lam_p, lam_E, de, de/2, de/2);
            lambda_vector[i][j][3] = harmonic_mean(lam_p, lam_S, ds, ds/2, ds/2);
            lambda_vector[i][j][2] = harmonic_mean(lam_p, lam_W, dw, dw/2, dw/2);

        }
    }
    // Set the lambndas at the boundary
    for(int i = 1; i<N+1;i++){
        //west boundary
        lam_E = lambda_f[Material_matrix[i][0]][0];
        lambda_vector[i][0][1] = lam_E;

        //east boundary
        lam_W = lambda_f[Material_matrix[i][M+1]][0];
        lambda_vector[i][M+1][2] = lam_W;
    }
    for(int j = 1; j<M+1;j++){
        //North boundary
        lam_S = lambda_f[Material_matrix[N+1][j]][0];
        lambda_vector[N+1][j][3] = lam_S;

        //south boundary
        lam_N = lambda_f[Material_matrix[0][j]][0];
        lambda_vector[0][j][0] = lam_N;
    }

}

static void Next_delta_t (){
    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < M + 2; j++) {
            T[0][i][j] = T[1][i][j];
        }
    }

}
static void Mapa_inicial (){
    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < M + 2; j++) {
            T[0][i][j] = Tstart;
        }
    }
}
static void Mapa_estimat (){
    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < M + 2; j++) {
            T[1][i][j] = T[0][i][j];
        }
    }

}

static void cache (){
    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < M + 2; j++) {
            T_cache[0][i][j] = T[1][i][j];
        }
    }
}

static void solver_gauss_seidel () {
    // Defne the coefficients
    double a_N = 0, a_E = 0, a_S = 0, a_W = 0, a_p = 0,b_p = 0;
    //lambdas at surface between nodes
    double lam_n = 0, lam_e = 0, lam_s = 0, lam_w = 0;



    //Compute coefficient for boundary nodes at EAST AND WEST sides
    for (int i = 1; i < N + 1; i++) {
        lam_e = lambda_vector[i][0][1];
        de = x_all[i][1]-x_all[i][0];

        a_E = lam_e/(de);
        a_p = a_E + alfa_w;
        b_p = alfa_w * Twest;

        //COMPUTE TEMPERATURE AT THE WEST BOUNDARY
        T[1][i][0] = (a_E* T[1][i][1] +  b_p)/a_p;

        //COMPUTE TEMPERATURE AT THE EAST BOUNDARY
        T[1][i][M+1] = 8 + 0.005*t_actual;
    }

    // Compute the coefficients at the NORTH and South side

    for (int j = 1; j < M+1; j++) {
        //Compute lambda at the NORTH BOUNDARY
        lam_s = lambda_vector[N+1][j][3];
        lam_n = lambda_vector[0][j][0];

        a_S = lam_s/(y_all[N+1][0]-y_all[N][0]);
        a_p = a_S;

        b_p = q_w;
        //b_p = 0;
        //We compute the temperature
        T[1][N+1][j] = (a_S* T[1][N][j] +  b_p)/a_p;
       // T[1][N+1][j] = T[1][N][j]; //FOR ADIABATIC WALL

        //COMPUTE LAMNDA AT THE SOUTH BOUNDARY
//        a_N = lam_n/(y_all[1][0]-y_all[0][0]);
//        a_p = a_N;
//        b_p = 0;

         T[1][0][j] = Tsouth; //FOR ISOTHERMAL WALL
         //T[1][0][j] = (a_N* T[1][1][j] +  b_p)/a_p; //FOR ADIABATIC WALL

    }
// Compute the vertex
    for (int i = 0; i<4; i++) {
        switch (i) {
            case 0:                              // bottom left vertice
                T [1][0][0]  = T [1][0][1];
            case 1:                              //bottom right vertice
                T [1][0][M+1] =T [1][0][M];
            case 2:                              // top left vertice
                T [1][N+1][0] = T[1][N+1][1];
            case 3:                              // top right vertice
                T [1][N+1][M+1] = T [1][N+1][M];
        }
    }
    // Compute coefficients for the internal nodes
    for (int i = 1; i < N + 1; i++) {
        for (int j = 1; j < M + 1; j++) {

            lam_n = lambda_vector[i][j][0];
            lam_s = lambda_vector[i][j][3];
            lam_w = lambda_vector[i][j][2];
            lam_e = lambda_vector[i][j][1];

            dn = y_all[i+1][j]-y_all[i][j];
            ds = y_all[i][j]-y_all[i-1][j];
            de = x_all[i][j+1]-x_all[i][j];
            dw = x_all[i][j]-x_all[i][j-1];

            a_E = lam_e * S_v/de;
            a_S = lam_s * S_h/ds;
            a_W = lam_w * S_v/dw;
            a_N = lam_n * S_h/dn;

            a_p = (a_E + a_W +a_N + a_S)*Beta + rho[Material_matrix[i][j]] * V_p * Cp[Material_matrix[i][j]] / delta_t;



            b_p = qv_p[Material_matrix[i][j]]*Beta * V_p + V_p *rho[Material_matrix[i][j]] * Cp[Material_matrix[i][j]] * T[0][i][j] / delta_t + (1-Beta) *((T[0][i][j-1] - T[0][i][j])*a_W +(T[0][i][j+1] - T[0][i][j])*a_E   + (T[0][i-1][j] - T[0][i][j]) * a_S   + (T[0][i+1][j] - T[0][i][j])*a_N + qv_p[Material_matrix[i][j]] * V_p);

            T[1][i][j] = T_cache[0][i][j] + relaxation*((Beta*a_E* T[1][i][j+1] + Beta*a_W* T[1][i][j-1] + Beta*a_S * T[1][i-1][j] +  Beta*a_N*T[1][i+1][j]+ b_p)/a_p - T_cache[0][i][j]);
        }

    }
}
static void Balanç(int cont){
    double total_heat_flux = 0;
    for(int i = 1; i<N+1;i++){
        //west boundary
        Heat_parets[2][cont] += alfa_w*(T[1][i][0] - Twest)*W*p_3[1]/N;
        //east boundary
        Heat_parets[1][cont] += -lambda_vector[i][M+1][2]*((T[1][i][M+1]-T[1][i][M])/(x_all[i][M+1] - x_all[i][M]))*(p_3[1]/N)*W;
    }
    for(int j = 1; j<M+1;j++){
        //North boundary
        Heat_parets[0][cont] += -lambda_vector[N+1][j][3]*((T[1][N+1][j]-T[1][N][j])/(y_all[N+1][j]-y_all[N][j-1]))*(p_3[0]/M)*W;
        //south boundary
        Heat_parets[3][cont] += -lambda_vector[0][j][0]*((T[1][1][j]-T[1][0][j])/(y_all[1][j]-y_all[0][j-1]))*(p_3[0]/M)*W;
    }
}

int main() {
    clock_t tStart = clock();
    bool first = true;
    int cont = 0;

    double Q_est_total = 0;
    double Q_west_total = 0;
    double Q_nord_total = 0;
    double Q_south_total = 0;

    vec_geometric_deff();
    Mapa_inicial();
    Mapa_estimat();
    string s;
    while (t_actual <= (t_end+delta_t)){
        while((!norma(T_cache,T) || first)&& !excep){
            first = false;
            cache(); // start the cache that will go into norm
            solver_gauss_seidel();
        }

        first = true;

        if(cont%ratio_print_t== 0) {
            cout << "Process at: " << round(t_actual / (t_end + delta_t) * 100 ) << " %  \n";
        }
        Balanç(cont);

        cont++;
        Next_delta_t();
        t_actual += delta_t;
    }
    for (int i = 0; i< int(t_end/delta_t);i++){
        Q_nord_total += Heat_parets[0][i]*delta_t;
        Q_est_total += Heat_parets[1][i]*delta_t;
        Q_west_total += Heat_parets[2][i]*delta_t;
        Q_south_total += Heat_parets[3][i]*delta_t;
    }
    cout << "Calor total nord: " << Q_nord_total*1e-6 << " MJ\n";
    cout << "Calor total est: " << Q_est_total*1e-6 << " MJ\n";
    cout << "Calor total oest: " << Q_west_total*1e-6 << " MJ\n";
    cout << "Calor total sud: " << Q_south_total*1e-6 << " MJ\n";
    // Keep all the data in .csv files to later do the post-processing
    std::ofstream Tempfile;
    std::ofstream pos_file_x;
    std::ofstream pos_file_y;
    std::ofstream Data;

    Tempfile.open ("Temp_map.csv");
    pos_file_x.open ("posx_map.csv");
    pos_file_y.open ("posy_map.csv");
    Data.open("readme.txt"); //This file is used to identify what parameters where used in the simulation

    for (int i = N+1; i>=0; i--){
        for (int j = 0; j<M+2; j++){
            Tempfile << T[1][i][j] << ",";
            pos_file_x << x_all[i][j] << ",";
            pos_file_y << y_all[i][j] << ",";
        }
        Tempfile << "\n";
        pos_file_x << "\n";
        pos_file_y << "\n";
    }
    double Time_exec =(double)(clock()-tStart)/CLOCKS_PER_SEC;
    cout << "Time taken to execute: " << Time_exec << "s";
    Data << "Beta = " << Beta << "\n" << "Delta t = " << delta_t <<" s\n";
    Data << "Delta de convergència = " << delta_convergence << "\n";
    Data << "Volums de control eix vertical = " << N <<"\n";
    Data << "Volums de control eix horitzontal = " << M <<"\n";
    Data << "Parametre de relaxació = " << relaxation <<"\n";
    Data << "Temps d'execució = " << Time_exec << " s \n\n";

    Data << "Calor total nord: " << Q_nord_total*1e-6 << " MJ\n";
    Data << "Calor total est: " << Q_est_total*1e-6 << " MJ\n";
    Data << "Calor total oest: " << Q_west_total*1e-6 << " MJ\n";
    Data << "Calor total sud: " << Q_south_total*1e-6 << " MJ\n\n";

    Data << "Temperatura mapa inicial = " << Tstart << " ºC \n";

    Tempfile.close();
    pos_file_x.close();
    pos_file_y.close();

    return 0;
}