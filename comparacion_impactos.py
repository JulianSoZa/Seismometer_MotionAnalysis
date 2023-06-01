import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*
import pandas as pd

magneticField = 8E-3
spirals = 1824
length = 3E-3

#lectura = float(input("¿A que lectura desea acceder? "))
#lectura = float(120)
#lectura = float(248)
#lectura = float(231)
#lectura = float(220)

#lectura = float(228)
#lectura = float(223)

#lectura = float(219)

#lectura = float(218)
lectura = float(210)
#lectura = float(204)
#lectura = float(190)
datos = pd.read_csv("../DATOS/datosantes.csv")
   

y = datos["Aceleracion" + str(lectura)].to_numpy()
x = datos["Tiempo" + str(lectura)].to_numpy()

dt = x[int(len(x))-1]/len(x)

useful = (y).astype(float)
useful_time = (x).astype(float)

acceleration = useful

n = len(acceleration)

frq, yfft = fourierAnalysis.fourier_transform(acceleration, n, dt)

order = 2
accelerationTz = filters.z_transform(acceleration, order)

frqTz, yfftTz = fourierAnalysis.fourier_transform(accelerationTz, n, dt)

harmonics = 2

ysfftTz, signals, harfhzTz = fourierAnalysis.signal_decomposition(harmonics, accelerationTz, n, yfftTz, frqTz)

fig1,(ax,ax1) = plt.subplots(2,1)

ax.plot(useful_time, acceleration, label='Aceleración')
ax.plot(useful_time, accelerationTz, label='Aceleración - Filtro')
ax.legend()
ax.scatter(useful_time, acceleration, s = 12, c = 'black')
ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Aceleracion (m/s^2)')
ax.set_xlim(0, 1)
ax.grid()

ax1.vlines(frq, 0, abs(yfft.imag), colors ='tab:blue', label='Sin Filtro')
ax1.vlines(frqTz, 0, abs(yfftTz.imag), colors = 'tab:orange', label='Filtrada')
ax1.legend()
ax1.set_xlim(0, max(frq))
ax1.set_xlabel('Frecuencia (Hz)')
ax1.set_ylabel('F(w)')
ax1.grid()

plt.tight_layout()

fig3,(ax4) = plt.subplots(1,1)
for i in range(harmonics):
    ax4.plot(useful_time, signals[i].real, label = str(round(harfhzTz[i], 3))+' Hz')
ax4.legend()
ax4.set_xlim(0, 1)
ax4.set_xlabel('Tiempo (s)')
ax4.set_ylabel('Aceleracion (m/s^2)')
ax4.grid()

plt.tight_layout()

plt.show()