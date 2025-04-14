import matplotlib.pyplot as pyplot
import numpy as np

with open("Seidel.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res1 = data[:12]
    n1 = data[12:24]
    time1 = data[24:36]

with open("Jacobi.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res2 = data[:35]
    n2 = data[35:70]
    time2 = data[70:105]

with open("Simple.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res3 = data[:132]
    n3 = data[132:264]
    time3 = data[264:396]

with open("Chebyshev.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res4 = data[:128]
    n4 = data[128:256]
    time4 = data[256:384]



fig, ax = pyplot.subplots(figsize=(16, 12))
ax.grid()

ax.plot(n1, res1, zorder=10, linewidth=1.5, c="#0000FF", label='Gauss_Seidel')
ax.plot(n2, res2, zorder=10, linewidth=1.5, c="#00FF00", label='Jacobi')
ax.plot(n3, res3, zorder=10, linewidth=1.5, c="#DC143C", label='Simple_iteration')
ax.plot(n4, res4, zorder=10, linewidth=1.5, c="#8A2BE2", label='Chebyshev')

ax.set_title('Скорость сходимости по итерациям. Постоянная матрица')
ax.legend()
ax.set_xlabel(r'Количество итерация $n$', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')



fig, ax = pyplot.subplots(figsize=(16, 12))
ax.grid()
ax.plot(time1, res1, zorder=10, linewidth=1.5, c="#0000FF", label='Gauss_Seidel')
ax.plot(time2, res2, zorder=10, linewidth=1.5, c="#00FF00", label='Jacobi')
ax.plot(time3, res3, zorder=10, linewidth=1.5, c="#DC143C", label='Simple_iteration')
ax.plot(time4, res4, zorder=10, linewidth=1.5, c="#8A2BE2", label='Chebyshev')

ax.set_title('Скорость сходимости по времени. Постоянная матрица')
ax.legend()
ax.set_xlabel(r'Время итерации $t$, мкс', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')
pyplot.show()
