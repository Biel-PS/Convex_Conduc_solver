import numpy as np
import simpy as sp
import matplotlib as mp
#geometric_par --> is a numerical array witch contains the geometric data of any case
# if spehere  (S)--> [polar [deg],azimutal [deg],radius [m]]
# if cylinder (C)--> [0, height [m], radius [m]]
# if plane block (P) --> [z length [m], y length [m], x length [m]]
# for the plane block, q is 1D in x direction

geometric_par = [50E-3,10E-3,5E-3]

type = 'P' #S for sphere, C for cylinder, P for plane block (type must be char!)

#geometric_temp --> is a numerical array witch contains the Temperature data of all any case
# All cases --> [innerT [k], OuterT [K]]
# NOTE: innerT is considered to be the yz plane in x = 0 in case of P, center axis in C and center point on S
# OuterT is set at R surfaces in S and C, and the plane zy in x = x lenght in case of type = P.

geometric_temp = [20,200]



lambda_transport = [0,0] #Deffine the number of lambdas from closest to x = 0 (or r = 0) to furthest.
alfa_convex = [0,0] #Deffine the number of convexion coeff from closest to x = 0 (or r = 0) to furthest.
start_run = True
# if there is no alfa between boundaries set it as 0.
class Plane_block:
    def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex):
        self.Geometry = geometric_par
        self.Temp = geometric_temp
        self.Lambda = lambda_transport
        self.Alfa = alfa_convex


class Sphere:
    def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex):
        self.Geometry = geometric_par
        self.Temp = geometric_temp
        self.Lambda = lambda_transport
        self.Alfa = alfa_convex

class Cylinder:
    def __init__(self,geometric_par,geometric_temp,lambda_transport,alfa_convex):
        self.Geometry = geometric_par
        self.Temp = geometric_temp
        self.Lambda = lambda_transport
        self.Alfa = alfa_convex


match type:
    case 'P':
        Body = Plane_block (geometric_par,geometric_temp,lambda_transport,alfa_convex)
    case 'C':
        Body = Cylinder(geometric_par, geometric_temp, lambda_transport, alfa_convex)
    case 'S':
        Body = Sphere(geometric_par, geometric_temp, lambda_transport, alfa_convex)
    case _:
        print("Please, input a valid type value")
        start_run = False

