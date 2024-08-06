import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
 
img = mpimg.imread('airless.png')
object = [[False for i in range(800)] for j in range(400)]

for i in range(len(img)):
    for j in range(len(img[i])):
        if img[i][j][3] != 0:
            print("here")
            object[i][j] = True

print(object)
