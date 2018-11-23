import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

data_files = ['Exact.txt', 'CIR.txt', 'Lax-Friedrichs.txt', 'Lax-Wendroff.txt', 'Warming-Beam.txt']

def load_data(f):
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

    return NumOfStep, NumOfPnt, x, animation_data

if __name__ == '__main__':
    NumOfStep, NumOfPnt, x, exact_sol = load_data(data_files[0])
    _, _, _, CIR_sol = load_data(data_files[1])
    _, _, _, LF_sol = load_data(data_files[2])
    _, _, _, LW_sol = load_data(data_files[3])
    _, _, _, WB_sol = load_data(data_files[4])

    fig = plt.figure()
    
    ax1 = fig.add_subplot(221)
    ax1.set_ylabel(r'$u$')
    ax1.set_xlabel('X')
    line11, = ax1.plot(x, exact_sol[0, :], '-', label='Exact')
    line12, = ax1.plot(x, CIR_sol[0, :], '-.', label='CIR')
    ax1.legend()

    def update1(k):
        line11.set_data(x, exact_sol[k, :])
        line12.set_data(x, CIR_sol[k, :])
        return line11, line12

    a1 = animation.FuncAnimation(fig, update1, range(NumOfStep))

    ax2 = fig.add_subplot(222)
    ax2.set_ylabel(r'$u$')
    ax2.set_xlabel('X')
    line21, = ax2.plot(x, exact_sol[0, :], '-', label='Exact')
    line22, = ax2.plot(x, LF_sol[0, :], '-.', label='Lax-Friedrichs')
    ax2.legend()

    def update2(k):
        line21.set_data(x, exact_sol[k, :])
        line22.set_data(x, LF_sol[k, :])
        return line21, line22

    a2 = animation.FuncAnimation(fig, update2, range(NumOfStep))

    ax3 = fig.add_subplot(223)
    ax3.set_ylabel(r'$u$')
    ax3.set_xlabel('X')
    line31, = ax3.plot(x, exact_sol[0, :], '-', label='Exact')
    line32, = ax3.plot(x, LW_sol[0, :], '-.', label='Lax-Wendroff')
    ax3.legend()

    def update3(k):
        line31.set_data(x, exact_sol[k, :])
        line32.set_data(x, LW_sol[k, :])
        return line31, line32

    a3 = animation.FuncAnimation(fig, update3, range(NumOfStep))

    ax4 = fig.add_subplot(224)
    ax4.set_ylabel(r'$u$')
    ax4.set_xlabel('X')
    line41, = ax4.plot(x, exact_sol[0, :], '-', label='Exact')
    line42, = ax4.plot(x, WB_sol[0, :], '-.', label='Warming-Beam')
    ax4.legend()
    
    def update4(k):
        line41.set_data(x, exact_sol[k, :])
        line42.set_data(x, WB_sol[k, :])
        return line41, line42

    a4 = animation.FuncAnimation(fig, update4, range(NumOfStep))

    plt.tight_layout()
    plt.show()
