import numpy as np
import matplotlib.pyplot as plt
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

    print("Tw2: " + str(tw2_sol[0]-273.15))
    # sp.pprint(tw2_sol)

    tw1_sol = tw1_sol[0].subs(tw2, tw2_sol[0])
    # print(c1)
    print("Tw1: "+ str(tw1_sol-273.15))
    c1 = c1.subs([(tw2, tw2_sol[0]), (tw1, tw1_sol)])
    c2 = tw1_sol-273.15
    # print(c1)
    T = T.subs([(q_v, qv), (lam, Lambda), (C1, c1), (C2, c2)])
    q_x = q_x.subs([(q_v, qv), (lam, Lambda), (C1, c1)])
    return T, q_x

def Solve_4_multiple (lamnda_vector,alfa_vector,geometry_vector,boundary_temp,qv_vector):

    Wall_temp =list(sp.symbols(" ".join(f"Tw_{x}" for x in range(0,len(geometry_vector)*2)), real=True))
    Conv_Temp  = list(sp.symbols(" ".join(f"T_{x}" for x in range(0,len(lamnda_vector)+1)), real=True))


    Conv_Temp[0], Conv_Temp[len(Conv_Temp) - 1] = Conv_Temp[0].subs(Conv_Temp[0],boundary_temp[0]), Conv_Temp[len(Conv_Temp) - 1].subs(Conv_Temp[len(Conv_Temp) - 1],boundary_temp[1])

    Eq_1,Eq_2,Eq_3 = [0]*len(geometry_vector),[0]*len(geometry_vector),[0]*(len(geometry_vector)-1)
    for i in range(0,len(geometry_vector)):
        C = ((Wall_temp[2*i] - Wall_temp[2*i-1])/(geometry_vector[i])+(qv_vector[i]*geometry_vector[i])/(2*lamnda_vector[i]))
        Eq_1[i] = alfa_vector[i]*(Conv_Temp[i] - Wall_temp[2*i-1]) + lamnda_vector[i]*C
        Eq_2[i]= qv_vector[i] * geometry_vector[i] - lamnda_vector[i] * C - alfa_vector[i+1]*(Wall_temp[2*i] - Conv_Temp[i+1])
    for i in range(0,len(geometry_vector)-1):
        C_A = ((Wall_temp[2 * i] - Wall_temp[2 * i - 1]) / (geometry_vector[i]) + (qv_vector[i] * geometry_vector[i]) / (2 * lamnda_vector[i]))
        C_B = ((Wall_temp[2 * (i+1)] - Wall_temp[2 * (i+1) - 1]) / (geometry_vector[(i+1)]) + (qv_vector[(i+1)] * geometry_vector[(i+1)]) / (2 * lamnda_vector[(i+1)]))
        Eq_3[i]= lamnda_vector[i] * C_A - lamnda_vector[i+1] * C_B - qv_vector[i]*geometry_vector[i]
    #print(Eq_1,"\n", Eq_2,"\n", Eq_3)
    Values = np.append(Conv_Temp[1:len(Conv_Temp)-1],Wall_temp)
    All_Eq = np.append([Eq_1,Eq_2],Eq_3)
    #print(All_Eq)
    Values_solved = sp.solve(All_Eq,Values,dict = False,set = True)
    Values_solved_values = list(Values_solved[1])[0]
    Values_solved_index = Values_solved [0]
    print(Values_solved_values,"\n", Values_solved_index)
    #Values_solved = np.append(Values_solved, [Conv_Temp[0], Conv_Temp[len(Conv_Temp) - 1]])


def Plot_T_and_q_x(q_x,T,Geometry, step):
    x = sp.symbols('x')
    pos_Temp_max = 0
    first = True
    num_positions = int((Geometry[2]) / step + 1)
    x_ = np.linspace(0, Geometry[2], num_positions)
    Temp = np.zeros(num_positions)
    q_x_ = np.zeros(num_positions)
    cache = 0

    for i in range(0, num_positions):
        Temp[i] = T.subs(x, x_[i])
        q_x_[i] = q_x.subs(x, x_[i])
        if cache != 0 and Temp[i - 1] > Temp[i] and first == True:
            pos_Temp_max = i - 1
            first = False
        cache = i

    print('Tmax is ' + str(round(Temp[pos_Temp_max], 2)) + ' cº and its found in position ' + str(
        x_[pos_Temp_max] * 1E3) + ' mm')
    print(
        'T at x=0 is T_0= ' + str(round(Temp[0], 2)) + ' cº and T at x= ' + str(x_[cache] * 1E3) + 'mm is T_e = ' + str(
            round(Temp[cache], 2)) + 'ºc')
    print(
        'q at x=0 is q_0=' + str(round(q_x_[0], 2)) + ' w/m2 and q at x= ' + str(x_[cache] * 1E3) + 'mm is q_e = ' + str(
            round(q_x_[cache], 2)) + 'w/m2')
    plt.figure()
    plt.title("Temperature in function of material lenght")
    plt.plot(x_ * 1E3, Temp, linewidth=2, markersize=12, color='red')
    plt.xlabel("Material lenght [mm]")
    plt.ylabel("Temperature [cº]")
    plt.grid()
    plt.figure()
    plt.title("Heat transfer in function of material lenght")
    plt.plot(x_ * 1E3, q_x_)
    plt.xlabel("Material lenght [mm]")
    plt.ylabel("Heat [w/m2]")
    plt.grid()
    plt.show()