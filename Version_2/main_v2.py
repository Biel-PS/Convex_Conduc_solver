import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import plane_wall as pw

one = False

if one == True:
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

    lambda_transport = 300 #Deffine the number of lambdas from closest to x = 0 (or r = 0) to furthest (one for each iteration).
    alfa_convex = [30,200] #Deffine the number of convexion coeff from closest to x = 0 (or r = 0) to furthest.
    start_run = True
    Position = 0


    q_v,lam,C1,C2 = sp.symbols ('q_v lambda C_1 C_2')

    pos_queue = 0  # pos_queue 0 if first, 1, 2,...
    pros_queue = 0  # pos_queue 0 if first, 1, 2,...
    previous_iter_T = 0
    number_of_bodyes_to_study = 5
    Control_vector = np.zeros(number_of_bodyes_to_study)

    T, q_x = pw.Solve_for_unique_body(geometric_par, geometric_temp, lambda_transport, alfa_convex, qv)

    pw.Plot_T_and_q_x(q_x,T,geometric_par,1E-4)
else:
    lambda_vector = [300,100,20]
    alfa_vector = [200,100,50,300]
    geometry_vector = [15E-3,10E-3,5E-3]
    boundary_temp = [273.15+200,273.15+20]
    qv_vector = [1E7,2E7,3E7]
    pw.Solve_4_multiple(lambda_vector,alfa_vector,geometry_vector,boundary_temp,qv_vector)


