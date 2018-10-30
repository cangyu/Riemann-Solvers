import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

fin = open('Case 1.txt', 'r')
NumOfStep, NumOfPnt = fin.readline().strip().split()
NumOfStep = int(NumOfStep)
NumOfPnt = int(NumOfPnt)

x = np.zeros(NumOfPnt)
for k in range(NumOfPnt):
    x[k] = float(fin.readline().strip())

animation_data = np.zeros((NumOfStep, NumOfPnt, 3))
for n in range(NumOfStep):
    fin.readline()
    for k in range(NumOfPnt):
        animation_data[n][k] = fin.readline().strip().split()

fin.close()

fig = plt.figure()
ax1 = fig.add_subplot(111)
line, = ax1.plot(x, animation_data[0, :, 0])

def update(data):
    line.set_ydata(data[:,0])
    return line

a = animation.FuncAnimation(fig, update, animation_data)
plt.show()