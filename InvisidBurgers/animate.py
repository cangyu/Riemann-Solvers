import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

if __name__ == '__main__':
    while True:
        ok = False
        f = ''
        while not ok:
            f = input("Enter filename:")
            if os.path.exists(f):
                ok = True
            else:
                if f is '#':
                    exit(0)

        fin = open(f, 'r')
        NumOfStep, NumOfPnt = fin.readline().strip().split()
        NumOfStep = int(NumOfStep)
        NumOfPnt = int(NumOfPnt)

        x = np.array([float(c) for c in fin.readline().strip().split()])

        all_data = []
        for n in range(NumOfStep):
            cur_data = [float (c) for c in fin.readline().strip().split()]
            all_data.append(cur_data)
        animation_data = np.array(all_data)

        fin.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_ylabel(r'$u$')
        ax.set_xlabel('X')
        ax.set_title(os.path.splitext(f)[0])

        line, = ax.plot(x, animation_data[0, :])

        def update(data):
            cur_data = data
            #ax.set_ylim(-0.05, 1.05)
            line.set_ydata(cur_data)
            return line

        a = animation.FuncAnimation(fig, update, animation_data)

        plt.tight_layout()
        plt.show()