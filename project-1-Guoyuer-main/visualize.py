import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt('milestone-2/init_matrix.csv', dtype=float, delimiter=',')

plt.imshow(data)
plt.colorbar()
plt.show()
# plt.savefig('2.png')