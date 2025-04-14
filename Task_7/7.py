import matplotlib.pyplot as pyplot
import numpy as np

with open("results_GM.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    res1 = data[:39]
    n1 = data[39:78]
    time1 = data[78:117]

    res2 = data[117:181]
    n2 = data[181:245]
    time2 = data[245:309]

    res3 = data[309:346]
    n3 = data[346:383]
    time3 = data[383:420]

    res4 = data[420:436]
    n4 = data[436:452]
    time4 = data[452:468]

fig, ax = pyplot.subplots(figsize=(8, 6))
ax.grid()

ax.plot(n1, res1, zorder=10, linewidth=1.5, c="#006400", label='Simple_iteration')
ax.plot(n2, res2, zorder=10, linewidth=1.5, c="#1E90FF", label='Chebyshev')
ax.plot(n3, res3, zorder=10, linewidth=1.5, c="#000080", label='Gradient descent')
ax.plot(n4, res4, zorder=10, linewidth=1.5, c="#9400D3", label='Сonjugate gradient')

ax.set_title('Скорость сходимости по итерациям. Постоянная матрица')
ax.legend()
ax.set_xlabel(r'Количество итерация $n$', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')



fig, ax = pyplot.subplots(figsize=(8, 6))
ax.grid()
ax.plot(time1, res1, zorder=10, linewidth=1.5, c="#006400", label='Simple_iteration')
ax.plot(time2, res2, zorder=10, linewidth=1.5, c="#1E90FF", label='Chebyshev')
ax.plot(time3, res3, zorder=10, linewidth=1.5, c="#000080", label='Gradient descent')
ax.plot(time4, res4, zorder=10, linewidth=1.5, c="#9400D3", label='Сonjugate gradient')

ax.set_title('Скорость сходимости по времени. Постоянная матрица')
ax.legend()
ax.set_xlabel(r'Время итерации $t$, мкс', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Невязка $r$', fontsize=14, fontname='Georgia')
pyplot.show()