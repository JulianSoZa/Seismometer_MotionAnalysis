import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*
import pandas as pd

datos = pd.read_csv("../DATOS/datos_7_08.csv",sep=',',decimal=".")

dat = datos.to_numpy()
y = dat[:,1]
x = dat[:,0]

useful = (y).astype(float)
useful_time = (x).astype(float)
  


frq, transformada = fourierAnalysis.fourier_transform(useful, n = len(y), dt=0.00458308)

fig,(ax,ax1) = plt.subplots(2,1)

ax.vlines(frq, 0, abs(transformada.imag))
ax.set_xlim(0, max(frq))
ax.set_xlabel('Frecuencia (Hz)')
ax.set_ylabel('F(w)')
ax.grid()

ax1.plot(useful_time,useful)
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Aceleracion (m/s^2)')
ax1.grid()
plt.show()


