import numpy as np
import matplotlib.pyplot as plt
import sys

file_name = sys.argv[1]
data = np.genfromtxt(file_name, dtype=float, delimiter=',')
plt.figure()
_file_name = file_name.split('.')[0]
plt.suptitle(_file_name)
plt.title("n=1000")
x = data[:, 0]
y = data[:, 1]
y = y[0] / y
plt.xlabel('nprocs')
plt.ylabel('speedup')
plt.plot(x, y)
for i in range(len(x)):
    plt.text(x[i], y[i], f'({x[i]}, {y[i]:.1f})')
plt.savefig(f'{_file_name}.png')
plt.close()
