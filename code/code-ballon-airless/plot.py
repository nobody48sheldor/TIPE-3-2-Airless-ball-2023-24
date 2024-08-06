import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import os

#style.use('dark_background')

os.system("./main 9.5 40")

with open("data/data_1.txt", 'r') as file:
    x = file.read()

with open("data/data_2.txt", 'r') as file:
    z = file.read()

x = x.split("/")
z = z.split("/")

X = []
Z = []

for t in range(len(x)-1):
    X.append(float(x[t]))
    Z.append(float(z[t]))

# print(X)
print()
print()
# print(Z)

t = np.linspace(0,2,40)

g = 9.81
v_x = 7.27742
v_z = 6.10648

def x_f(t):
    return(v_x*t)

def z_f(t):
    if -g*t*t/2 + v_z*t + 1.9 > 0:
        return(-g*t*t/2 + v_z*t + 1.9)
    else:
        return(0)

X_f = []
Z_f = []

for time in range(len(t)):
    X_f.append(x_f(t[time]))
    Z_f.append(z_f(t[time]))


X_ = []
Z_ = []


for i in range(len(Z)):
    if Z[i] != 0:
        Z_.append(Z[i])
        X_.append(X[i])

plt.plot(X_,Z_)
#plt.scatter(X_f, Z_f, color="red", marker='o')
plt.plot([6.7-0.15, 6.7+0.15], [3.1, 3.1], color="red")
plt.show()
