import numpy as np
import time
from scipy.fftpack import fft, ifft, fftfreq
import matplotlib.pyplot as plt
from funtions import*
import scipy.signal as signal

"""N_order = 1
frqCut = 60
frqSamp = 400
wn = frqCut/frqSamp
b, a = signal.butter(N_order, wn, 'low')
print(b)
print(a)"""

n = 4000

tSpan = np.linspace(0.0,1, n)
y = 3*np.sin(4*2*np.pi*tSpan) + 5*np.sin(6*2*np.pi*tSpan) + 2*np.sin(8*2*np.pi*tSpan) #- np.random.randn(n)
dt = 1/n

yfft = fft(y) / n  # Normalizada
frq = fftfreq(n, dt)  # Recuperamos las frecuencias
harmonics = 3 # Numeros de armonicos de la descomposicion

ysfftTz, signals, harfhzTz = fourierAnalysis.signal_decomposition(harmonics, y, n, yfft, frq)

fig,(ax,ax1) = plt.subplots(2,1)

ax.plot(tSpan, y)
ax.grid()

ax1.vlines(frq, 0, abs(yfft.imag))
ax1.set_xlim(0, 50)
ax1.grid()

fig2,(ax2) = plt.subplots(1,1)
for i in range(harmonics):
    ax2.plot(tSpan, signals[i].real, label = str(round(harfhzTz[i], 3))+' Hz')
ax2.legend()
ax2.grid()

plt.tight_layout()
plt.show()
