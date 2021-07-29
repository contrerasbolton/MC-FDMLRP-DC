import numpy as np
import matplotlib.pyplot as plt
import random
import sys

x_location = [25, 85, 70, 60, 30, 120, 20, 130, 70, 170, 150]
y_location = [65, 85, 20, 70, 100, 50, 10, 120, 140, 10, 150]
# 1: 20, 2: 20, 3: 10, 4: 10, 5: 20, 6: 20, 7: 30, 8: 30, 9: 40, 10: 40
coefficients = [20, 20, 10, 10, 20, 20, 30, 30, 40, 40]

block = 20
number = int(sys.argv[1])
depot = int(sys.argv[2])
posCoef = int(sys.argv[3])
dimension = block * number
x_location = x_location[:depot]
y_location = y_location[:depot]
coef = coefficients[posCoef]

# generate plot of the instance
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
major_ticks = np.arange(0, dimension + 1, block)
minor_ticks = np.arange(0, dimension + 1, block)

index_i = [ i for i in range(10, dimension, block)]
index_j = [ j for j in range(20, dimension, block)]

x_label = []
y_label = []
for i in index_i[::-1]:
    for j in index_i:
        x_label.append(j)
        y_label.append(i)
Narea = len(x_label)
area = [i for i in range(0, len(x_label))]

x_coordinates = []
y_coordinates = []
for i in index_j[::-1]:
    for j in index_j:
        x_coordinates.append(j)
        y_coordinates.append(i)
Nnode = len(x_coordinates)
node = [i + depot for i in range(0, len(x_coordinates))]

ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)

# And a corresponding grid
ax.grid(which = 'both')

# Or if you want different settings for the grids:
ax.grid(which = 'minor', alpha = 0.2)
ax.grid(which = 'major', alpha = 0.5)


location = [i for i in range(len(x_location))]

x_location = x_location[:depot]
y_location = y_location[:depot]

print(len(area), len(x_label), len(y_label))
N = len(x_coordinates) + len(x_location)
ciudades = [i for i in range(N)]
arcos = [(i,j) for i in ciudades for j in ciudades]
# formato de distancia EUC_2D

x_coord = x_location + x_coordinates
y_coord = y_location + y_coordinates
distancia = {(i, j): round(np.hypot(x_coord[i] - x_coord[j], y_coord[i] - y_coord[j]) / coef, 2) for i, j in arcos}

plt.xlim((0, dimension))
plt.ylim((0, dimension))
# print(x_coordinates)
# print(y_coordinates)
plt.scatter(x_coordinates, y_coordinates, s = 220)
plt.scatter(x_location, y_location, color = 'red', s = 220)
# plt.scatter(x_label, y_label, color = 'white', s = 0)

# enumerate with text to the areas
for i, txt in enumerate(area):
    plt.annotate(txt, (x_label[i], y_label[i] + 4))

# enumerate with text to the nodes
for i, txt in enumerate(node):
    plt.annotate(txt, (x_coordinates[i], y_coordinates[i]), color = 'white', ha = 'center', va = 'center')

# enumerate with text to the locations
for i, txt in enumerate(location):
    plt.annotate(txt, (x_location[i], y_location[i]), color = 'white', ha = 'center', va = 'center')

plt.show()

print(Narea, Nnode)
matrix = []
# Matrix in matrix
aux = []
for area in range(0, Narea):
    if area != 0 and area % number == 0:
        matrix.append(aux)
        aux = []
        print()
    if area / 10 < 1:
        aux.append(area)
        print("0" + str(area), end = " ")
    else:
        aux.append(area)
        print(area, end = " ")
print()
matrix.append(aux)
print("Coordinates:")
# print(matrix)

node_i = len(x_location)
for k in range(0, number - 1):
    for v in range(0, number - 1):
        # print(node_i, " = ", end = " ")
        print(node_i, 4, end = " ")
        for i in range(k, k + 2):
            for j in range(v, v  + 2):
                print(matrix[i][j], end = " ")
        print()
        node_i += 1
print()
# Cost Matrix
for i in range(N):
    for j in range(N):
        if i != j:
            print(distancia[(i, j)], end = " ")
        else:
            print(99999, end = " ")
    print()
