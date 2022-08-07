import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('result.csv', dtype=float, delimiter=',',skip_header=1)
plt.figure()
x = data[:, 0]
y_pure = data[:, 1]
y_mix = data[:, 2]
print(y_pure)
print(y_mix)
plt.xlabel('number of nodes')
plt.ylabel('runtime(s)')
plt.ylim(0, 45)

plt.plot(x, y_pure,label='y_pure')
plt.plot(x, y_mix, label= 'y_mix')
plt.legend()
plt.savefig('strong_scaling.png')