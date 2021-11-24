import numpy as np
import matplotlib.pyplot as plt
import random
import sys

if len(sys.argv) < 7:
    print("Usage: python3 checkInstance.py area ballons ewatchtower watchtower coef seed")
    exit(1)
block = 20
number = int(sys.argv[1]) # area^2
ballons = int(sys.argv[2]) # number of ballons
ewatchtower = int(sys.argv[3]) # number of watchtower
watchtower = int(sys.argv[4]) # number of existing watchtower
posCoef = int(sys.argv[5]) # coef that divide the distance
seed = int(sys.argv[6]) # seed
depot = ballons + ewatchtower + watchtower
dimension = block * number
coef = posCoef #coefficients[posCoef]
random.seed(number)

#coefficients = [random.randint(1, 50) for i in range(len(x_location))]
#x_location = [25, 85, 70, 60, 30, 120, 20, 130, 70, 170, 150]
#y_location = [65, 85, 20, 70, 100, 50, 10, 120, 140, 10, 150]
# 1: 20, 2: 20, 3: 10, 4: 10, 5: 20, 6: 20, 7: 30, 8: 30, 9: 40, 10: 40
#coefficients = [20, 20, 10, 10, 20, 20, 30, 30, 40, 40]
#print(coefficients)
interval = 20
if dimension / 2 >= depot:
    interval = 10

x_location = random.sample(range(10, dimension, interval), depot)
y_location = random.sample(range(10, dimension, interval), depot)


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

#plot the locations: ballons: red, ewatchtower: green, watchtower: orange
plt.scatter(x_coordinates, y_coordinates, s = 220)
plt.scatter(x_location[:ballons], y_location[:ballons], color = 'red', s = 220)
plt.scatter(x_location[ballons:ballons + ewatchtower], y_location[ballons:ballons + ewatchtower], color = 'green', s = 220)
plt.scatter(x_location[ballons + ewatchtower:ballons + ewatchtower + watchtower], y_location[ballons + ewatchtower:ballons + ewatchtower + watchtower], color = 'orange', s = 220)
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

print(Narea, Nnode)
matrix = []
# Matrix in matrix
aux = []
for Area in range(0, Narea):
    if Area != 0 and Area % number == 0:
        matrix.append(aux)
        aux = []
        print()
    if Area / 10 < 1:
        aux.append(Area)
        print("0" + str(Area), end = " ")
    else:
        aux.append(Area)
        print(Area, end = " ")
print()
matrix.append(aux)
print("Coordinates:")
print(matrix)

nodes = []
node_i = len(x_location)
for k in range(0, number - 1):
    for v in range(0, number - 1):
        temp = []
        # print(node_i, " = ", end = " ")
        # print(node_i, 4, end = " ")
        temp.append(node_i)
        temp.append(4)
        for i in range(k, k + 2):
            for j in range(v, v  + 2):
                # print(matrix[i][j], end = " ")
                temp.append(matrix[i][j])
        # print()
        node_i += 1
        nodes.append(temp)
print()
# # Cost Matrix
# for i in range(N):
#     for j in range(N):
#         if i != j:
#             print(distancia[(i, j)], end = " ")
#         else:
#             print(99999, end = " ")
#     print()


print("# Instance")
print(Narea, ballons, ewatchtower, watchtower, len(nodes))
print("2 10 1 5")
print("150 10 0 6")
print()
for d in range(depot):
    j = int(x_location[d] / block)
    i = abs(int(y_location[d] / block)  - number) - 1
    Blist = []
    eWlist = []
    Wlist = []
    if d < ballons: # ballons
        # size = random.randint(3, 4)
        # print("d", d, x_location[d], y_location[d], i, j, size)
        # print(d, size, end = " ")
        flag = False
        if i - 1 > -1 and j - 1 > -1 and i + 1 < number and j + 1 < number:
            flag = True
        # print(d, i, number, flag)
        # exit(0)

        if i - 1 > -1 and j - 1 > -1:
            if flag:
                if random.random() < 0.5:
                    Blist.append(matrix[i][j])
                    Blist.append(matrix[i][j+1])
                    Blist.append(matrix[i+1][j])
                    Blist.append(matrix[i+1][j+1])
                else:
                    Blist.append(matrix[i][j])
                    Blist.append(matrix[i][j+1])
                    Blist.append(matrix[i+1][j])
                    Blist.append(matrix[i+1][j+1])
            else:
                Blist.append(matrix[i-1][j-1])
                Blist.append(matrix[i-1][j])
                Blist.append(matrix[i][j-1])
                Blist.append(matrix[i][j])
        elif i + 1 < number and j + 1 < number:
            Blist.append(matrix[i][j])
            Blist.append(matrix[i][j+1])
            Blist.append(matrix[i+1][j])
            Blist.append(matrix[i+1][j+1])
        elif i + 1 < number and j - 1 < number:
            Blist.append(matrix[i][j-1])
            Blist.append(matrix[i][j])
            Blist.append(matrix[i+1][j-1])
            Blist.append(matrix[i+1][j])


        # if len(Blist) == 4:
        #     size = random.randint(len(Blist) - 1, len(Blist))
        #     remove = random.randint(0, size)
        #     Blist.pop(remove)
        # else:
        #     size = 3
        print(d, len(Blist), ' '.join(str(e) for e in Blist))


    elif d < ballons + ewatchtower:
        size = random.randint(5, 9)
        # print("e", d, x_location[d], y_location[d], i, j, size)
        #print("\n" + str(d), size)
        if i - 1 > -1 and j - 1 > -1:
            eWlist.append(matrix[i-1][j-1])
        if i - 1 > -1:
            eWlist.append(matrix[i-1][j])
        if i - 1 > -1 and j + 1 < number:
            eWlist.append(matrix[i-1][j+1])
        if j - 1 > -1:
            eWlist.append(matrix[i][j-1])
        eWlist.append(matrix[i][j])
        if j + 1 < number:
            eWlist.append(matrix[i][j+1])
        if i + 1 < number:
            if j - 1 > -1:
                eWlist.append(matrix[i+1][j-1])
            eWlist.append(matrix[i+1][j])
            if j + 1 < number:
                eWlist.append(matrix[i+1][j+1])

        size = random.randint(4, len(eWlist))

        print(d, len(eWlist), ' '.join(str(e) for e in eWlist))
    else:
        size = random.randint(5, 9)
        #print("w", d, x_location[d], y_location[d], i, j, size)
        # print("\n" + str(d), size)
        if i - 1 > -1 and j - 1 > -1:
            Wlist.append(matrix[i-1][j-1])
        if i - 1 > -1:
            Wlist.append(matrix[i-1][j])
        if i - 1 > -1 and j + 1 < number:
            Wlist.append(matrix[i-1][j+1])
        if j - 1 > -1:
            Wlist.append(matrix[i][j-1])
        Wlist.append(matrix[i][j])
        if j + 1 < number:
            Wlist.append(matrix[i][j+1])
        if i + 1 < number:
            if j - 1 > -1:
                Wlist.append(matrix[i+1][j-1])
            Wlist.append(matrix[i+1][j])
            if j + 1 < number:
                Wlist.append(matrix[i+1][j+1])
        size = random.randint(4, len(Wlist))
        print(d, len(Wlist), ' '.join(str(e) for e in Wlist))

for values in nodes:
    print(values[0], values[1], values[2], values[3], values[4], values[5])
print()
# Cost Matrix
for i in range(N):
    for j in range(N):
        if i != j:
            print(distancia[(i, j)], end = " ")
        else:
            print(99999, end = " ")
    print()
plt.show()
