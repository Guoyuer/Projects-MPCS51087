import numpy as np
import matplotlib.pyplot as plt

filenames = ['strong_scaling_1.png', 'strong_scaling_2.png', 'weak_scaling.png']
titles = ['Speedup vs # of threads, N = 3200', 'Speedup vs # of threads, N = 200',
          'Speedup vs # of threads, N = 800 * nT']

for i in range(3):
    data = np.genfromtxt(f'{i + 1}.csv', dtype=int, delimiter=',')
    x = data[:, 0]
    y = data[:, 1] / data[0, 1]
    fig = plt.figure()
    plt.plot(x, y)
    plt.xlabel('# of threads')
    plt.ylabel('Speedup')
    plt.title(titles[i])
    plt.savefig(filenames[i])
