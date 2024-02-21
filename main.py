import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

#geometric_par --> is a numerical array witch contains the geometric data of any case
# if spehere  (S)--> [polar [deg],azimutal [deg],radius [m]]
# if cylinder (C)--> [0, height [m], radius [m]]
# if plane block (P) --> [z length [m], y length [m], x length [m]]
# for the plane block, q is 1D in x direction

geometric_par = [1,1,15E-3]

type = 'P' #S for sphere, C for cylinder, P for plane block (type must be char!)

#geometric_temp --> is a numerical array witch contains the Temperature data of all any case
# All cases --> [innerT [k], OuterT [K]]
# NOTE: innerT is considered to be the yz plane in x = 0 in case of P, center axis in C and center point on S
# OuterT is set at R surfaces in S and C, and the plane zy in x = x lenght in case of type = P.

geometric_temp = [273.15+30,273.15+80]

qv = 1E7 #internal generated heat.
# print(qv)

lambda_transport = [0.5] #Deffine the number of lambdas from closest to x = 0 (or r = 0) to furthest (one for each iteration).
alfa_convex = [30,200] #Deffine the number of convexion coeff from closest to x = 0 (or r = 0) to furthest.
start_run = True
Position = 0;

q_v,lam,C1,C2 = sp.symbols ('q_v lambda C_1 C_2')
# if there is no alfa between boundaries set it as 0.
class Plane_block:
    def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex,qv,Position,previous_iter_T,pos_queue):
        self.Geometry = geometric_par
        self.Temp = geometric_temp
        self.Lambda = lambda_transport[Position]
        self.Alfa = alfa_convex
        self.previous_iter_T = previous_iter_T
        self.pos_queue = pos_queue
        self.qv = qv
        self.position = Position
        # Definition of transfer equation
        x = sp.symbols ('x')
        self.T = -q_v/(2*lam)*x**2 + C1*x + C2
        self.q_x = q_v*x-lam*C1
    def Solve_for_unique_body (self):
        tw1,tw2,x = sp.symbols('T_{w1} T_{w2} x')

        c1 = (tw2 - tw1) / (self.Geometry[2]) + (self.qv * self.Geometry[2]) / (2 * self.Lambda)


        q_ca = self.Alfa[0]*(self.Temp[0]-tw1)
        q_a = - self.Lambda*c1

        q_cb = self.Alfa[1] * (tw2 - self.Temp[1])
        q_b = self.qv * self.Geometry[2] - self.Lambda * c1

        tw1_sol = sp.solve((q_ca-q_a),tw1)

        q_b = q_b.subs(tw1,tw1_sol[0])

        tw2_sol = sp.solve((q_cb-q_b),tw2)

        # sp.pprint(tw2_sol)

        tw1_sol = tw1_sol[0].subs(tw2,tw2_sol[0])
        # print(c1)
        c1 = c1.subs([ (tw2,tw2_sol[0]), (tw1, tw1_sol)])
        c2 = tw1_sol-273.15
        # print(c1)
        self.T = self.T.subs([(q_v,self.qv), (lam,self.Lambda), (C1,c1), (C2,c2)])
        self.q_x = self.q_x.subs([(q_v, self.qv), (lam, self.Lambda), (C1, c1)])
        return self.T, self.q_x
    def Plot_T_and_q_x (self,step):
        x = sp.symbols('x')
        pos_Temp_max = 0
        first = True
        num_positions = int((self.Geometry[2])/step + 1)
        x_ = np.linspace(0,self.Geometry[2],num_positions)
        Temp = np.zeros(num_positions)
        q_x = np.zeros(num_positions)
        cache = 0;
        for i in range(0,num_positions):
            Temp[i] = self.T.subs(x,x_[i])
            q_x[i] = self.q_x.subs(x,x_[i])
            if cache != 0 and Temp[i-1] > Temp[i] and first == True:
                pos_Temp_max = i-1
                first = False
            cache = i


        print('Tmax is ' + str(round(Temp[pos_Temp_max],2)) + ' cº and its found in position ' + str(x_[pos_Temp_max]*1E3) + ' mm')
        print('T at x=0 is T_0= ' + str(round(Temp[0],2)) + ' cº and T at x= ' + str(x_[cache]*1E3) + 'mm is T_e = ' + str(round(Temp[cache],2)) + 'ºc')
        print('q at x=0 is q_0=' + str(round(q_x[0],2)) + ' w/m2 and q at x= ' + str(x_[cache]*1E3) + 'mm is q_e = ' + str(round(q_x[cache],2)) + 'w/m2')
        plt.figure()
        plt.title("Temperature in function of material lenght")
        plt.plot(x_*1E3,Temp,linewidth=2, markersize=12,color = 'red')
        plt.xlabel("Material lenght [mm]")
        plt.ylabel("Temperature [cº]")
        plt.grid()
        plt.figure()
        plt.title("Heat transfer in function of material lenght")
        plt.plot(x_*1E3,q_x)
        plt.xlabel("Material lenght [mm]")
        plt.ylabel("Heat [w/m2]")
        plt.grid()
        plt.show()

    def Solve_for_multiple (self):
        next_T = sp.symbols('T_next')
        tw1, tw2, x = sp.symbols('T_{w1} T_{w2} x')
        if self.pos_queue== 0:
            c1 = (tw2 - tw1) / (self.Geometry[2]) + (self.qv * self.Geometry[2]) / (2 * self.Lambda)
            q_ca = self.Alfa[0] * (self.Temp[0] - tw1)
            q_a = - self.Lambda * c1






        q_cb = self.Alfa[1] * (tw2 - next_T)
        q_b = self.qv * self.Geometry[2] - self.Lambda * c1

        tw1_sol = sp.solve((q_ca - q_a), tw1)

        q_b = q_b.subs(tw1, tw1_sol[0])

        tw2_sol = sp.solve((q_cb - q_b), tw2)

        # sp.pprint(tw2_sol)

        tw1_sol = tw1_sol[0].subs(tw2, tw2_sol[0])
        # print(c1)
        c1 = c1.subs([(tw2, tw2_sol[0]), (tw1, tw1_sol)])
        c2 = tw1_sol - 273.15
        # print(c1)
        self.T = self.T.subs([(q_v, self.qv), (lam, self.Lambda), (C1, c1), (C2, c2)])
        self.q_x = self.q_x.subs([(q_v, self.qv), (lam, self.Lambda), (C1, c1)])
        return self.T, self.q_x


#
# class Sphere:
#     def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex):
#         self.Geometry = geometric_par
#         self.Temp = geometric_temp
#         self.Lambda = lambda_transport
#         self.Alfa = alfa_convex
#
#         r = sp.symbols ('r')
#         T = -q_v/(6*lam)*r**2 + (c1/r) + c2
#         q_r = (1/3)*q_v*r+lam*c1/(r**2)
#
#
#
#
# class Cylinder:
#     def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex,previous_iter_q,first):
#         self.Geometry = geometric_par
#         self.Temp = geometric_temp
#         self.Lambda = lambda_transport
#         self.Alfa = alfa_convex
#         self.previous_iter_q = previous_iter_q
#         self.first = first
#         r = sp.symbols ('r')
#         T = -q_v/(4*lam)*r**2 + c1*sp.log(r) + c2
#         q_r = (1/2)*q_v*r-lam*c1/(r)
#         qw
#     def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex,qv,Position):
#         self.Geometry = geometric_par
#         self.Temp = geometric_temp
#         self.Lambda = lambda_transport[Position]
#         self.Alfa = alfa_convex
#         self.qv = qv
#         # Definition of transfer equation
#         r = sp.symbols('r')
#         self.T = -q_v / (4 * lam) * r ** 2 + c1 * sp.log(r) + c2
#         self.q_r = (1 / 2) * q_v * r - lam * c1 / (r)
#
#     def Solve_for_unique_body (self):
#         tw1,tw2,x = sp.symbols('T_{w1} T_{w2} x')
#
#         c1 = (tw2 - tw1) / (self.Geometry[2]) + (self.qv * self.Geometry[2]) / (2 * self.Lambda)
#
#
#         q_ca = self.Alfa[0]*(self.Temp[0]-tw1)
#         q_a = - self.Lambda*c1
#
#         q_cb = self.Alfa[1] * (tw2 - self.Temp[1])
#         q_b = self.qv * self.Geometry[2] - self.Lambda * c1
#
#         tw1_sol = sp.solve((q_ca-q_a),tw1)
#
#         q_b = q_b.subs(tw1,tw1_sol[0])
#
#         tw2_sol = sp.solve((q_cb-q_b),tw2)
#
#         # sp.pprint(tw2_sol)
#
#         tw1_sol = tw1_sol[0].subs(tw2,tw2_sol[0])
#         # print(c1)
#         c1 = c1.subs([ (tw2,tw2_sol[0]), (tw1, tw1_sol)])
#         c2 = tw1_sol-273.15
#         # print(c1)
#         self.T = self.T.subs([(q_v,self.qv), (lam,self.Lambda), (C1,c1), (C2,c2)])
#         self.q_x = self.q_x.subs([(q_v, self.qv), (lam, self.Lambda), (C1, c1)])
#         return self.T, self.q_x
#
#
#
#     def Plot_T_and_q_x (self,step):
#         x = sp.symbols('x')
#         pos_Temp_max = 0
#         first = True
#         num_positions = int((self.Geometry[2])/step + 1)
#         x_ = np.linspace(0,self.Geometry[2],num_positions)
#         Temp = np.zeros(num_positions)
#         q_x = np.zeros(num_positions)
#         cache = 0;
#         for i in range(0,num_positions):
#             Temp[i] = self.T.subs(x,x_[i])
#             q_x[i] = self.q_x.subs(x,x_[i])
#             if cache != 0 and Temp[i-1] > Temp[i] and first == True:
#                 pos_Temp_max = i-1
#                 first = False
#             cache = i
#
#
#         print('Tmax is ' + str(Temp[pos_Temp_max]) + ' cº and its found in position ' + str(x_[pos_Temp_max]*1E3) + ' mm')
#         print('T at x=0 is T_0= ' + str(Temp[0]) + ' cº and T at x= ' + str(x_[cache]*1E3) + 'mm is T_e = ' + str(Temp[cache]) + 'ºc')
#         print('q at x=0 is q_0=' + str(q_x[0]) + ' w/m2 and q at x= ' + str(x_[cache]*1E3) + 'mm is q_e = ' + str(q_x[cache]) + 'w/m2')
#         plt.figure()
#         plt.title("Temperature in function of material lenght")
#         plt.plot(x_*1E3,Temp,linewidth=2, markersize=12,color = 'red')
#         plt.xlabel("Material lenght [mm]")
#         plt.ylabel("Temperature [cº]")
#         plt.grid()
#         plt.figure()
#         plt.title("Heat transfer in function of material lenght")
#         plt.plot(x_*1E3,q_x)
#         plt.xlabel("Material lenght [mm]")
#         plt.ylabel("Heat [w/m2]")
#         plt.grid()
#         plt.show()


match type:
    case 'P':


        pos_queue = 0 #pos_queue 0 if first, 1, 2,...
        pros_queue = 0 #pos_queue 0 if first, 1, 2,...
        previous_iter_T = 0
        number_of_bodyes_to_study = 5
        Control_vector = np.zeros(number_of_bodyes_to_study)


        Body = Plane_block (geometric_par,geometric_temp,lambda_transport,alfa_convex,qv,Position,previous_iter_T,pos_queue)
        [T,q_x] = Body.Solve_for_unique_body()
        # sp.pprint(T)
        Body.Plot_T_and_q_x(1E-4)



    #
    # case 'C':
    #     Body = Cylinder(geometric_par, geometric_temp, lambda_transport, alfa_convex)
    # case 'S':
    #     Body = Sphere(geometric_par, geometric_temp, lambda_transport, alfa_convex)
     case _:
        print("Please, input a valid type value")
        start_run = False

pos_queue = 0  # pos_queue 0 if first, 1, 2,...
previous_iter_T = 0
number_of_bodyes_to_study = 5
Control_vector = np.zeros(number_of_bodyes_to_study)

for i in range(0,number_of_bodyes_to_study):
