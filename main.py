import numpy as np
import simpy as sp
#geometric_par --> is a numerical array witch contains the geometric data of any case
# if spehere  (S)--> [polar [deg],azimutal [deg],radius [m]]
# if cylinder (C)--> [0, height [m], radius [m]]
# if plane block (P) --> [z length [m], y length [m], x length [m]]
# for the plane block, q is 1D in x direction

geometric_par = [50E-3,10E-3,5E-3]

type = 'P' #S for sphere, C for cylinder, P for plane block (type must be char!)

#externalT --> is a numerical array witch contains the Temperature data of all any case
# All cases --> [innerT [k], OuterT [K]]
# NOTE: innerT is considered to be the yz plane in x = 0 in case of P, center axis in C and center point on S
# OuterT is set at R surfaces in S and C, and the plane zy in x = x lenght in case of type = P.




def Iter (type,geometric_par,q_V,externalT):