import numpy as np
import sympy as sp
import matplotlib as mp

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
print(qv)

lambda_transport = [300,300] #Deffine the number of lambdas from closest to x = 0 (or r = 0) to furthest.
alfa_convex = [30,200] #Deffine the number of convexion coeff from closest to x = 0 (or r = 0) to furthest.
start_run = True

q_v,lam,c1,c2 = sp.symbols ('q_v lambda C_1 C_2')
# if there is no alfa between boundaries set it as 0.
class Plane_block:
    def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex,qv):
        self.Geometry = geometric_par
        self.Temp = geometric_temp
        self.Lambda = lambda_transport
        self.Alfa = alfa_convex
        self.qv = qv
        #Definition of transfer equation
        x = sp.symbols ('x')
        T = -q_v/(2*lam)*x**2 + c1*x + c2
        q_x = q_v*x-lam*c1

    def Solve_for_unique_body (self):
        tw1,tw2,x = sp.symbols('T_{w1} T_{w2} x')

        c1_a = (tw2 - tw1) / (self.Geometry[2]) + (self.qv * self.Geometry[2]) / (2 * self.Lambda[0])
        c1_b = (tw2 - tw1) / (self.Geometry[2]) + (self.qv * self.Geometry[2]) / (2 * self.Lambda[1])

        q_ca = self.Alfa[0]*(self.Temp[0]-tw1)
        q_a = - self.Lambda[0]*c1_a

        q_cb = self.Alfa[1] * (tw2 - self.Temp[1])
        q_b = self.qv * self.Geometry[2] - self.Lambda[1] * c1_b

        tw1_sol = sp.solve((q_ca-q_a),tw1)

        q_b = q_b.subs(tw1,tw1_sol[0])

        tw2_sol = sp.solve((q_cb-q_b),tw2)

        sp.pprint(tw2_sol)

        tw1_sol = tw1_sol[0].subs(tw2,tw2_sol[0])

        return tw1_sol-273.15,tw2_sol[0]-273.15

class Sphere:
    def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex):
        self.Geometry = geometric_par
        self.Temp = geometric_temp
        self.Lambda = lambda_transport
        self.Alfa = alfa_convex

        r = sp.symbols ('r')
        T = -q_v/(6*lam)*r**2 + (c1/r) + c2
        q_r = (1/3)*q_v*r+lam*c1/(r**2)




class Cylinder:
    def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex):
        self.Geometry = geometric_par
        self.Temp = geometric_temp
        self.Lambda = lambda_transport
        self.Alfa = alfa_convex

        r = sp.symbols ('r')
        T = -q_v/(4*lam)*r**2 + c1*sp.log(r) + c2
        q_r = (1/2)*q_v*r-lam*c1/(r)


match type:
    case 'P':
        Body = Plane_block (geometric_par,geometric_temp,lambda_transport,alfa_convex,qv)
        [T1_deg,T2_deg] = Body.Solve_for_unique_body()
        print(T1_deg)
    case 'C':
        Body = Cylinder(geometric_par, geometric_temp, lambda_transport, alfa_convex)
    case 'S':
        Body = Sphere(geometric_par, geometric_temp, lambda_transport, alfa_convex)
    case _:
        print("Please, input a valid type value")
        start_run = False


