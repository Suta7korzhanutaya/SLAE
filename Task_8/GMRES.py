import matplotlib.pyplot as pyplot
import numpy as np

with open("res_GMRES.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',','.')
    data = np.array([float(elem) for elem in data])
    y1 = data[:99]

    x = np.arange(1, 100, 1)
    y2 = np.full(99, 394)
fig, ax = pyplot.subplots(figsize=(8, 6))
ax.grid()

ax.plot(x, y1, zorder=10, linewidth=1.5, c="#006400", label='CG')
ax.plot(x, y2, zorder=10, linewidth=1.5, c="#1E90FF", label='GMRES')


ax.set_title('Общее решение задачи в зависимости от $m$')
ax.legend()
ax.set_xlabel(r'$m$ - количество прогонок GMRES до автостарта', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'$t$, мкс - Время решения', fontsize=14, fontname='Georgia')
pyplot.savefig("Test_2.png")

