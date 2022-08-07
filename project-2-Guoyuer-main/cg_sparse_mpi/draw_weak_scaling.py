from locale import normalize
import numpy as np
import matplotlib.pyplot as plt
import sys

file_name = sys.argv[1]
data = np.genfromtxt(file_name, dtype=float, delimiter=',')
plt.figure()
_file_name = file_name.split('.')[0]
plt.suptitle(_file_name)
plt.title("init n=600")
x = data[:, 0]
y = data[:, 1]

# normalize by number of iteration!!!
y = y / data[:, 2]
# print(y)
y = y[0] / y
plt.xlabel('nprocs')
plt.ylabel('speedup')
plt.plot(x, y)
plt.ylim(0, 1)
for i in range(len(x)):
    plt.text(x[i], y[i], f'({x[i]}, {y[i]:.1f})')
plt.savefig(f'{_file_name}.png')
plt.close()
