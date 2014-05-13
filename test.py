from scipy.fftpack import fft, rfft, fftfreq
import pylab as plt
import scipy.signal
import numpy as np
import scipy as sp
import os
import glob
import math
from array import array

another_list = [] 
w = 2 * math.pi / 30
list_count = []
for i in range(1000):
	another_list.append(np.sin(i*w))
	list_count.append(i/200.0)

x = array("d",another_list)
y = array("d",list_count)
x_array = np.array(x, np.float64)
y_array = np.array(y, np.float64)
f = np.linspace(0.001, 100, 10000)

pgram = sp.signal.lombscargle(y_array, x_array, f)

plt.subplot(2, 1, 1)
plt.plot(list_count, another_list)
plt.subplot(2, 1, 2)
plt.plot(f, pgram)

plt.show()
