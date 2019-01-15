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

        x = np.zeros(NumOfPnt)
        for k in range(NumOfPnt):
            x[k] = float(fin.readline().strip())

        animation_data = np.zeros((NumOfStep, NumOfPnt, 4))
        for n in range(NumOfStep):
            fin.readline()
            for k in range(NumOfPnt):
                animation_data[n][k] = fin.readline().strip().split()

        fin.close()

        fig = plt.figure()
        ax1 = fig.add_subplot(411)
        ax1.set_ylabel(r'$\rho$')

        ax2 = fig.add_subplot(412)
        ax2.set_ylabel(r'$u$')

        ax3 = fig.add_subplot(413)
        ax3.set_xlabel('X')
        ax3.set_ylabel(r'$P$')

        ax4 = fig.add_subplot(414)
        ax4.set_xlabel('X')
        ax4.set_ylabel(r'$e$')

        line1, = ax1.plot(x, animation_data[0, :, 0])
        line2, = ax2.plot(x, animation_data[0, :, 1])
        line3, = ax3.plot(x, animation_data[0, :, 2])
        line4, = ax4.plot(x, animation_data[0, :, 3])

        margin = 0.05


        def update1(data):
            cur_data = data[:, 0]
            bot = np.min(cur_data)
            top = np.max(cur_data)
            height = top - bot
            mh = margin * height
            if top > bot:
                ax1.set_ylim(bot - mh, top + mh)
            line1.set_ydata(cur_data)
            return line1


        def update2(data):
            cur_data = data[:, 1]
            bot = np.min(cur_data)
            top = np.max(cur_data)
            height = top - bot
            mh = margin * height
            if top > bot:
                ax2.set_ylim(bot - mh, top + mh)
            line2.set_ydata(cur_data)
            return line2


        def update3(data):
            cur_data = data[:, 2]
            bot = np.min(cur_data)
            top = np.max(cur_data)
            height = top - bot
            mh = margin * height
            if top > bot:
                ax3.set_ylim(bot - mh, top + mh)
            line3.set_ydata(cur_data)
            return line3

        def update4(data):
            cur_data = data[:, 3]
            bot = np.min(cur_data)
            top = np.max(cur_data)
            height = top - bot
            mh = margin * height
            if top > bot:
                ax4.set_ylim(bot - mh, top + mh)
            line4.set_ydata(cur_data)
            return line4


        a = animation.FuncAnimation(fig, update1, animation_data)
        b = animation.FuncAnimation(fig, update2, animation_data)
        c = animation.FuncAnimation(fig, update3, animation_data)
        d = animation.FuncAnimation(fig, update4, animation_data)

        plt.tight_layout()
        plt.show()
