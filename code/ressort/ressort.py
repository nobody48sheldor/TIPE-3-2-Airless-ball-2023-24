import matplotlib.pyplot as plt
import numpy as np

Y = 9.81*np.linspace(0.11,0.11+(8*0.05),8)
X = np.array([0.15, 0.248, 0.349, 0.47, 0.602, 0.72, 0.85, 0.97])

print("X = ", X)
print("Y = ", Y)
ecart = []

for i in range(len(Y)-1):
    ecart.append(X[i+1]-X[i])

print(ecart)
pas = Y[1] - Y[0]
k = 1/(np.mean(ecart)/(pas))

print("k = ", k, " en N/m")

#Y_droite = k*(X - X[0]) + np.mean(Y) - (k*(X[-1] - X[0])/2)
Y_droite = k*(X - np.mean(X)) + np.mean(Y)

ecart_type = np.sqrt(np.mean([x*x/(pas*pas) for x in ecart]) - (np.mean(ecart)/pas)**2)
print(ecart_type )

plt.plot(X,Y)
plt.scatter(X,Y)
plt.xlabel("élongation (m)")
plt.ylabel("poids  (N)")
plt.plot(X,Y_droite)
plt.show()
