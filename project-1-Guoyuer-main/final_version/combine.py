import sys
import numpy as np

def main():
    world_size = int(sys.argv[1])
    save_per_step = int(sys.argv[2])
    tot_step = int(sys.argv[3])
    with open('matrix_0_0.txt') as f:
        line = f.readline().rstrip()
        n = int(line.split()[0])

    data = np.zeros(shape=(n, n), dtype=float)
    for step in range(0, tot_step, save_per_step):
        print(step)
        # for i in range(world_size):
        #     filename = f"matrix_{i}_{step}.txt"
        #     with open(filename) as f:
        #         line = f.readline().rstrip()
        #         [_, x, y, div] = map(int, line.split())
        #     print(x, y, div)
        #     sub = np.loadtxt(filename, skiprows=1, )
        #     data[x:x + div, y:y + div] = sub
        # output_name = f"matrix_{step}.csv"
        # np.savetxt(output_name, data, delimiter=',')


if __name__ == "__main__":
    main()
