import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import sys


def main():
    n = int(sys.argv[1])  # number of particles
    # n_iter = sys.argv[2] # number of iterations

    def read_file():
        data = np.genfromtxt('particles.dat', delimiter=' ', skip_header=1)
        data = data[:, :2]
        return data

    data = read_file()
    L = 10
    fig = plt.figure(figsize=(L*1.2, L*1.2))
    ax = plt.axes(xlim=(-L, L), ylim=(-L, L))
    colors = np.random.random((len(data), 3))
    scatter = ax.scatter(data[:, 0], data[:, 1], c=colors)

    def update(i):
        one_frame = data[i*n:(i+1)*n, :]
        # print(one_frame)
        scatter.set_offsets(one_frame)
        scatter.set_color(colors[0:n, :])
        return scatter,

    anim = animation.FuncAnimation(
        fig, update, interval=20, blit=True, repeat=False, frames=len(data)//n)
    # plt.show()
    anim.save('animation.mp4')


main()
