import numpy as np
import matplotlib.pyplot as plt
import sys


def main():
    file_name = sys.argv[1]
    data = np.genfromtxt(file_name, dtype=float, delimiter=' ')
    plt.figure()
    _file_name = file_name.split('.')[0]
    plt.title(_file_name)
    plt.imshow(data)
    plt.colorbar()
    plt.clim(0, 1)
    plt.savefig(f'{_file_name}.png')
    plt.close()


if __name__ == "__main__":
    main()
