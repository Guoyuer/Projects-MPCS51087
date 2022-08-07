import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def read_file(filename):
    data = np.genfromtxt(filename, delimiter=' ', skip_header=1)
    # data = data[:, :2]
    return data


def plot_1_frame(filename, i):
    L = 20000
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    data = read_file(filename)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    ax.scatter3D(x, y, z, s=1)
    ax.set_xlim3d(-L, L)
    ax.set_ylim3d(-L, L)
    ax.set_zlim3d(-L, L)
    plt.savefig(f'frame{i}.png', dpi = 200)


def main():
    for i in range(0, 401):
        filename = f'production_step_{i+1}.dat'
        plot_1_frame(filename, i)
        print(f"{filename} saved!")


main()
