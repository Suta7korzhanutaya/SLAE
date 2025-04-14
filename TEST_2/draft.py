import matplotlib.pyplot as pyplot
import numpy as np

with open("result.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',','.')
    data = np.array([float(elem) for elem in data])
    x1 = data[:132]
    y1 = data[132:264]

    x2 = data[264:283]
    y2 = data[283:302]

    x3 = data[302:321]
    y3 = data[321:340]

    x4 = data[340:352]
    y4 = data[352:364]

    x5 = data[364:368]
    y5 = data[368:372]

fig, ax = pyplot.subplots(figsize=(8, 6))
ax.grid()

ax.plot(x1, y1, zorder=10, linewidth=1.5, c="#006400", label='GD')
ax.plot(x2, y2, zorder=10, linewidth=1.5, c="#1E90FF", label='GD $/tau_{opt}$')
ax.plot(x3, y3, zorder=10, linewidth=1.5, c="#000080", label='SGD')
ax.plot(x4, y4, zorder=10, linewidth=1.5, c="#9400D3", label='GD+')
ax.plot(x5, y5, zorder=10, linewidth=1.5, c="#8B0000", label='CG')

ax.set_title('Проекция шагов методов на плоскость, задаваемую собстенными векторами')
ax.legend()
ax.set_xlabel(r'$/lambda_{max}$', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'$/lambda_{min}$', fontsize=14, fontname='Georgia')
pyplot.savefig("Test_2.png")

