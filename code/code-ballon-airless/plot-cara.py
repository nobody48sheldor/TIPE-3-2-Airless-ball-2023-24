import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style

#style.use('dark_background')

nt = 100000

with open("data/data_lance-franc-v.txt", 'r') as file:
    v_lf = file.read()
with open("data/data_lance-franc-alpha.txt", 'r') as file:
    alpha_lf = file.read()

with open("data/data_3-points-v.txt", 'r') as file:
    v_3p = file.read()
with open("data/data_3-points-alpha.txt", 'r') as file:
    alpha_3p = file.read()

#with open("data/data_center-v.txt", 'r') as file:
#    v_c = file.read()
#with open("data/data_center-alpha.txt", 'r') as file:
#    alpha_c = file.read()



v_lf = v_lf.split("/",nt)
alpha_lf = alpha_lf.split("/",nt)
v_3p = v_3p.split("/",nt)
alpha_3p = alpha_3p.split("/",nt)
#v_c = v_c.split("/",nt)
#alpha_c = alpha_c.split("/",nt)


for i in range(len(v_lf)-1):
    plt.scatter(float(v_lf[i]), float(alpha_lf[i]), color="red")

for j in range(len(v_3p)-1):
    plt.scatter(float(v_3p[j]), float(alpha_3p[j]), color="blue")

#for j in range(len(v_c)-1):
#    plt.scatter(float(v_c[j]), float(alpha_c[j]), color="green")

plt.savefig("cara-airless.png")
plt.show()
