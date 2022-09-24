import numpy as np
from matplotlib import pyplot

pi = np.pi   #mathモジュールのπを利用



x1 = np.linspace(0, 2*pi, 1000)  #0から2πまでの範囲を100分割したnumpy配列
y1 = np.sin(x1) * 10**-13

x2 = np.linspace(pi / 4, 2*pi, 1000)  #0から2πまでの範囲を100分割したnumpy配列
y2 = 2 * np.sin(x2) * 10**-13

x3 = np.linspace(pi / 2, 2*pi, 1000)  #0から2πまでの範囲を100分割したnumpy配列
y3 = 3 * np.sin(x3) * 10**-13

y = np.array([y1, y2, y3])
x = np.array([x1, x2, x3])


pyplot.plot(x[0], x[1], x[2], y[0], y[1], y[2])
pyplot.show()

