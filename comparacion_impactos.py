import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*
import pandas as pd

magneticField = 8E-3
spirals = 1824
length = 7E-3

#lectura = float(input("¿A que lectura desea acceder? "))
lectura = float(220)
datos = pd.read_csv("../DATOS/datosSismometro.csv")
   

y = datos["Aceleracion" + str(lectura)].to_numpy()
x = datos["Tiempo" + str(lectura)].to_numpy()

position = np.repeat(0.0, 1600)

dt = x[1599]/len(x)

useful = (y).astype(float)
useful_time = (x).astype(float)

velocity, offsetVelocity = kinematics.velocity_calculation(useful, magneticField, spirals, length)
velocity = velocity*0.09
order = 2
velocity = filters.z_transform(velocity, order)
position = kinematics.position_calculation(position, velocity, dt)
acceleration = np.repeat(0.0,1600)
acceleration = kinematics.acceleration_calculation(acceleration, velocity, useful_time)

n = 1600

frq, yfft = fourierAnalysis.fourier_transform(acceleration, n, dt)

harmonics = 2

ysfft, signals = fourierAnalysis.signal_decomposition(harmonics, acceleration, n, yfft)


fig2,(ax,ax1) = plt.subplots(2,1)

ax.plot(useful_time, acceleration)
ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Aceleracion (m/s^2)')
ax.grid()

ax1.plot(useful_time, velocity)
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Velocidad (m/s)')
ax1.grid()

fig3,(ax4) = plt.subplots(1,1)
for i in range(harmonics):
    ax4.plot(useful_time, signals[i].real)
ax4.set_xlim(0, 3)
ax4.set_xlabel('Tiempo (s)')
ax4.set_ylabel('velocidad')
ax4.grid()

plt.tight_layout()
plt.show()