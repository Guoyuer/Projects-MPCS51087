from operator import delitem
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

arr = np.loadtxt("window.out")
plt.imshow(arr, cmap='gray')
plt.show()
# plt.savefig('gpu_ball.png', dpi = 300)