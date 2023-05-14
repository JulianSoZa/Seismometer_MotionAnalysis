import numpy as np
import time
from scipy.fftpack import fft, ifft, fftfreq
import matplotlib.pyplot as plt
from funtions import*

n = 400

tSpan = np.linspace(0.0,1, n)
y = np.sin(88*tSpan) + 5*np.sin(176*tSpan) + 2*np.sin(126*tSpan) - np.random.randn(n)
dt = 1/n

yfft = fft(y) / n  # Normalizada
frq = fftfreq(n, dt)  # Recuperamos las frecuencias
signal = ifft(yfft)*n

harmonics = 5
maximums = np.flip(np.sort(abs(yfft.imag)))
ysfft = []
signals = []

for i in range(harmonics):
    ysfft.append(fft(y)/n)
    ysfft[i][np.where(abs(ysfft[i].imag) != maximums[i*2])] = 0
    signals.append(ifft(ysfft[i])*n)

fig,(ax,ax1) = plt.subplots(2,1)

ax.plot(tSpan, y)
ax.grid()

ax1.vlines(frq, 0, abs(yfft.imag))
ax1.set_xlim(0, 50)
ax1.grid()

fig2,(ax2) = plt.subplots(1,1)
for i in range(harmonics):
    ax2.plot(tSpan, signals[i])
ax2.grid()

plt.tight_layout()
plt.show()
