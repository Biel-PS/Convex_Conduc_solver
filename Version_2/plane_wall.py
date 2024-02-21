import numpy as np
import matplotlib as plt
import sympy as sp


def Solve_for_unique_body(Geometry,Temp,Lambda,Alfa,qv):


    tw1, tw2, x = sp.symbols('T_{w1} T_{w2} x')
    x = sp.symbols('x')
    q_v, lam, C1, C2 = sp.symbols('q_v lambda C_1 C_2')

    T = -q_v / (2 * lam) * x ** 2 + C1 * x + C2
    q_x = q_v * x - lam * C1

    c1 = (tw2 - tw1) / (Geometry[2]) + (qv * Geometry[2]) / (2 * Lambda)

    q_ca = Alfa[0] * (Temp[0] - tw1)
    q_a = - Lambda * c1

    q_cb = Alfa[1] * (tw2 - Temp[1])
    q_b = qv * Geometry[2] - Lambda * c1

    tw1_sol = sp.solve((q_ca - q_a), tw1)

    q_b = q_b.subs(tw1, tw1_sol[0])

    tw2_sol = sp.solve((q_cb - q_b), tw2)

    # sp.pprint(tw2_sol)

    tw1_sol = tw1_sol[0].subs(tw2, tw2_sol[0])
    # print(c1)
    c1 = c1.subs([(tw2, tw2_sol[0]), (tw1, tw1_sol)])
    c2 = tw1_sol - 273.15
    # print(c1)
    T = T.subs([(q_v, qv), (lam, Lambda), (C1, c1), (C2, c2)])
    q_x = q_x.subs([(q_v, qv), (lam, Lambda), (C1, c1)])
    return T, q_x


def Plot_T_and_q_x(self, step):
    x = sp.symbols('x')
    pos_Temp_max = 0
    first = True
    num_positions = int((Geometry[2]) / step + 1)
    x_ = np.linspace(0, Geometry[2], num_positions)
    Temp = np.zeros(num_positions)
    q_x = np.zeros(num_positions)
    cache = 0;
    for i in range(0, num_positions):
        Temp[i] = T.subs(x, x_[i])
        q_x[i] = q_x.subs(x, x_[i])
        if cache != 0 and Temp[i - 1] > Temp[i] and first == True:
            pos_Temp_max = i - 1
            first = False
        cache = i
def Solve