from Data_input import *
import numpy as np

delta_r = (R_ext-R_int)/N

r_w = np.zeros(N + 1)
r_p = np.zeros(N+2)

for i in range(0,N+1):
    r_w[i] = R_int + i*delta_r

r_p[0] = R_int
r_p[N+1] = R_ext

for i in range (1,N+1):
    r_p[i] = (r_w[i] + r_w[i-1])/2

