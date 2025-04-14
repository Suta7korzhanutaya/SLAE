import matplotlib.pyplot as pyplot
import numpy as np

with open("res_Gauss_Seidel.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res1 = data[:22]
    n1 = data[22:44]
    time1 = data[44:66]

with open("res_Jacobi.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res2 = data[:2712]
    n2 = data[2712:5424]
    time2 = data[5424:8136]

with open("res_Simple_iteration.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res3 = data[:1082]
    n3 = data[1082:2164]
    time3 = data[2164:3246]

fig, ax = pyplot.subplots(figsize=(16, 12))
ax.grid()

ax.plot(n1, res1, zorder=10, linewidth=1.5, c="#0000FF", label='Gauss_Seidel')
ax.plot(n2, res2, zorder=10, linewidth=1.5, c="#00FF00", label='Jacobi')
ax.plot(n3, res3, zorder=10, linewidth=1.5, c="#DC143C", label='Simple_iteration')

ax.set_title('Методы простой итерации: скорость сходимости по итерациям. Матрица 100 х 100')
ax.legend()
ax.set_xlabel(r'Количество итерация $n$', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')



fig, ax = pyplot.subplots(figsize=(16, 12))
ax.grid()
ax.plot(time1, res1, zorder=10, linewidth=1.5, c="#0000FF", label='Gauss_Seidel')
ax.plot(time2, res2, zorder=10, linewidth=1.5, c="#00FF00", label='Jacobi')
ax.plot(time3, res3, zorder=10, linewidth=1.5, c="#DC143C", label='Simple_iteration')

ax.set_title('Методы простой итерации: скорость сходимости по времени. Матрица 100 х 100')
ax.legend()
ax.set_xlabel(r'Время итерации $t$, мкс', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')
pyplot.show()