import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import os

#style.use('dark_background')

g = 9.81
res = [(0,0)]


for v in range(4*20, 10*20):
#for v in range(6*20, 14*20):
#for v in range(7*20, 16*20):
    for alpha in range(10*4, 80*4):
        print("./main {0} {1}".format(v*0.05,alpha*0.25))
        os.system("./main {0} {1}".format(v*0.05,alpha*0.25))

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


        for pos in range(len(X)):
            
            # check if the ball is at the right pos in x

            #if ((3.66-0.2) <= X[pos] <= (3.66+0.2)):
            #if ((7.25-0.2) <= X[pos] <= (7.25+0.2)):
            if ((11.9-0.2) <= X[pos] <= (11.9+0.2)):

            # check if the ball is at the right pos in z

                if (3.1 <= Z[pos] <= (3.1 + 0.2)):

                # check if the ball is going down

                    if Z[pos-1] > Z[pos]:
                        if res[-1] != (v*0.05,alpha*0.25):
                            res.append( (v*0.05,alpha*0.25) )
                            print("BUCKET !")


res.pop(0)
print(res)

with open("data/data_lance-franc-v.txt", 'w') as data:
#with open("data/data_3-points-v.txt", 'w') as data:
#with open("data/data_center-v.txt", 'w') as data:
    for i in res:
        data.write(str(i[0])+"/")

with open("data/data_lance-franc-alpha.txt", 'w') as data:
#with open("data/data_3-points-alpha.txt", 'w') as data:
#with open("data/data_center-alpha.txt", 'w') as data:
    for i in res:
        data.write(str(i[1])+"/")


for i in res:
    plt.scatter(i[0],i[1], color='red')
plt.savefig("renders/lance-franc.png")
#plt.savefig("renders/3-points.png")
#plt.savefig("renders/center.png")
