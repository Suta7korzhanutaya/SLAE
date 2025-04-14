import matplotlib.pyplot as pyplot
import numpy as np

with open("results.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res1 = data[:24]
    n1 = data[24:48]
    time1 = data[48:72]

    res2 = data[72:85]
    n2 = data[85:98]
    time2 = data[98:111]

    res3 = data[111:142]
    n3 = data[142:173]
    time3 = data[173:204]

    res4 = data[204:236]
    n4 = data[236:268]
    time4 = data[268:300]

    res5 = data[300:309]
    n5 = data[309:318]
    time5 = data[318:327]

    res6 = data[327:336]
    n6 = data[336:345]
    time6 = data[345:354]

    res7 = data[354:383]
    n7 = data[383:412]
    time7 = data[412:441]


fig, ax = pyplot.subplots(figsize=(8, 6))
ax.grid()

ax.plot(n1, res1, zorder=10, linewidth=1.5, c="#8B0000", label='Jacobi')
ax.plot(n2, res2, zorder=10, linewidth=1.5, c="#FF8C00", label='Gauss_Seidel')
ax.plot(n3, res3, zorder=10, linewidth=1.5, c="#00FF00", label='Simple_iteration')
ax.plot(n4, res4, zorder=10, linewidth=1.5, c="#006400", label='Chebyshev')
ax.plot(n5, res5, zorder=10, linewidth=1.5, c="#1E90FF", label='Symmetrized Gauss_Seidel')
ax.plot(n6, res6, zorder=10, linewidth=1.5, c="#000080", label='Chebyshev Sym-Gauss_Seidel')
ax.plot(n7, res7, zorder=10, linewidth=1.5, c="#9400D3", label='Gradient descent')


ax.set_title('Скорость сходимости по итерациям. Постоянная матрица')
ax.legend()
ax.set_xlabel(r'Количество итерация $n$', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')



fig, ax = pyplot.subplots(figsize=(8, 6))
ax.grid()
ax.plot(time1, res1, zorder=10, linewidth=1.5, c="#8B0000", label='Jacobi')
ax.plot(time2, res2, zorder=10, linewidth=1.5, c="#FF8C00", label='Gauss_Seidel')
ax.plot(time3, res3, zorder=10, linewidth=1.5, c="#00FF00", label='Simple_iteration')
ax.plot(time4, res4, zorder=10, linewidth=1.5, c="#006400", label='Chebyshev')
ax.plot(time5, res5, zorder=10, linewidth=1.5, c="#1E90FF", label='Symmetrized Gauss_Seidel')
ax.plot(time6, res6, zorder=10, linewidth=1.5, c="#000080", label='Chebyshev Sym-Gauss_Seidel')
ax.plot(time7, res7, zorder=10, linewidth=1.5, c="#9400D3", label='Gradient descent')

ax.set_title('Скорость сходимости по времени. Постоянная матрица')
ax.legend()
ax.set_xlabel(r'Время итерации $t$, мкс', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')
pyplot.show()