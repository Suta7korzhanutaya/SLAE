import matplotlib.pyplot as pyplot
import numpy as np

with open("time.txt") as file:
    data = [string for string in file.read().split('\n')]
    for i in range(len(data)):
        data[i] = data[i].replace(',', '.')
    data = np.array([float(elem) for elem in data])
    xs = data[:210]
    csr_20 = data[210:420]
    den_20 = data[420:630]
    csr_25 = data[630:840]
    den_25 = data[840:1050]
    csr_30 = data[1050:1260]
    den_30 = data[1260:1470]
    csr_50 = data[1470:1680]
    den_50 = data[1680:1890]

fig, ax = pyplot.subplots(figsize=(16, 12))
ax.grid()

ax.plot(xs, csr_20, zorder=10, linewidth=1.5, c="#1E90FF", label='CSR_20%')
ax.plot(xs, den_20, zorder=10, linewidth=1.5, c="#0000FF", label='Dense_20%')

ax.plot(xs, csr_25, zorder=10, linewidth=1.5, c="#32CD32", label='CSR_25%')
ax.plot(xs, den_25, zorder=10, linewidth=1.5, c="#228B22", label='Dense_25%')

ax.plot(xs, csr_30, zorder=10, linewidth=1.5, c="#FF00FF", label='CSR_30%')
ax.plot(xs, den_30, zorder=10, linewidth=1.5, c="#8A2BE2", label='Dense_30%')

ax.plot(xs, csr_50, zorder=10, linewidth=1.5, c="#FF8C00", label='CSR_50%')
ax.plot(xs, den_50, zorder=10, linewidth=1.5, c="#FF4500", label='Dense_50%')

ax.set_title('CSR и Dense матрицы c % ненулевых элементов')
ax.legend()
ax.set_xlabel(r'Количество элементов в векторе $n$', fontsize=14, fontname='Georgia')
ax.set_ylabel(r'Время умножения матрицы на вектор $t$, мкс', fontsize=14, fontname='Georgia')

pyplot.show()