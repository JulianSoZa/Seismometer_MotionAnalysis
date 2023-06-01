import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*
import pandas as pd

magneticField = 8E-3
spirals = 1824
length = 3E-3

datos = pd.read_csv("../DATOS/DATOS_FINALES/168-172.csv",sep=';',decimal=",")

dat = datos.to_numpy()
y = dat[:,1]
x = dat[:,0]

useful = (y).astype(float)
useful = useful - sum(useful)/len(useful)
useful_time = (x).astype(float)
dt = useful_time[len(useful_time)-1]/len(useful_time)
n = len(y)
frq, transformada = fourierAnalysis.fourier_transform(useful, n, dt)

lectura = float(input("¿A que lectura desea acceder? "))

time_data_exp, data_exp = data_analis(lectura)

acceleration = np.repeat(0.0,len(data_exp))
velocity, offsetVelocity = kinematics.velocity_calculation(data_exp, magneticField, spirals, length)

acceleration = kinematics.acceleration_calculation(acceleration, velocity, time_data_exp)

fig,(ax,ax1) = plt.subplots(2,1)

"""useful_time = useful_time[15000:16600]
useful_time = useful_time - useful_time[0]
useful = useful[15000:16600]"""

#ax.plot(useful_time,useful)
ax.plot(time_data_exp, acceleration)
ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Aceleracion (m/s^2)')
ax.grid()

ax1.vlines(frq, 0, abs(transformada.imag))
ax1.set_xlim(0, max(frq))
ax1.set_xlabel('Frecuencia (Hz)')
ax1.set_ylabel('F(w)')
ax1.grid()

plt.show()