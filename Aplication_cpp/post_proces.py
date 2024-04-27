#POST PROCESSING MODULE by Biel P.S.

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
import csv

T = []
x = []
y = []

with open('Temp_map.csv') as temp_map:
    Temp_map = csv.reader(temp_map, delimiter=',')
    for row in Temp_map:
        T.append(row[0:(len(row)-1)])
    T_l = [0]*(len(T[0])*len(T))
    cont = 0
    for i in range(0,len(T)):
        for j in range(0,len(T[0])):
            T_l[cont] = float(T[i][j])
            cont +=1

    print(T_l)

with open('posx_map.csv') as posx:
    posx_vector = csv.reader(posx, delimiter=',')
    for row in posx_vector:
        x.append(row[0:(len(row)-1)])
    x_l = [0] * (len(x[0]) * len(x))
    cont = 0
    for i in range(0,len(x)):
        for j in range(0,len(x[0])):
            x_l[cont] = float(x[i][j])
            cont += 1
    print(x_l)

with open('posy_map.csv') as posy:
    posy_vector = csv.reader(posy, delimiter=',')
    for row in posy_vector:
        y.append(row[0:(len(row)-1)])
    y_l = [0] * (len(y[0]) * len(y))
    cont = 0
    for i in range(0,len(y)):
        for j in range(0,len(y[0])):
            y_l[cont] = float(y[i][j])
            cont+=1
    print(y_l)

def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} ºC" if plt.rcParams["text.usetex"] else f"{s} ºC"

data = pd.DataFrame({'x': x_l, 'y': y_l, 'T': T_l})
data_pivoted = data.pivot_table(index="y", columns="x",values= "T",sort = False)

ax = sns.heatmap(data_pivoted,annot = False, cmap= "viridis")

cs = plt.contour(data_pivoted, colors='k',levels = 10,corner_mask = True)
plt.clabel(cs,cs.levels, inline = True, fmt = fmt, fontsize = 10)
#plt.plot([len(x_l,1)/2,], [len(y),0], 'go-', label='line 1', linewidth=2)
plt.show()
