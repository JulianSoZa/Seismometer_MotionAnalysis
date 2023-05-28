import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*
import pandas as pd

datos = pd.read_csv("../DATOS/DATOS_FINALES/168-172.csv",sep=';',decimal=",")

dat = datos.to_numpy()
y = dat[:,3]
x = dat[:,0]

useful = (y).astype(float)
useful_time = (x).astype(float)
frq, transformada = fourierAnalysis.fourier_transform(useful, n = len(y), dt=0.00458308)

lectura = float(input("¿A que lectura desea acceder? "))

##Filtros digitales -----------------------------------------------

order = 2
acceleration = filters.z_transform(useful, order)

fig,(ax,ax1) = plt.subplots(2,1)

ax.vlines(frq, 0, abs(transformada.imag))
ax.set_xlim(0, max(frq))
ax.set_xlabel('Frecuencia (Hz)')
ax.set_ylabel('F(w)')
ax.grid()

ax1.plot(useful_time,useful)
ax1.plot(useful_time,acceleration)
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Aceleracion (m/s^2)')
ax1.grid()
data_analis(lectura)
plt.show()



