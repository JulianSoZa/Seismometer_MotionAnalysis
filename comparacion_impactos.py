import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*
import pandas as pd

magneticField = 8E-3
spirals = 1824
length = 7E-3

lectura = float(input("¿A que lectura desea acceder? "))
datos = pd.read_csv("../DATOS/datosSismometro.csv")
   

y = datos["Aceleracion" + str(lectura+1)].to_numpy()
x = datos["Tiempo" + str(lectura+1)].to_numpy()

position = np.repeat(0.0, 1600)

dt = x[1599]/len(x)

useful = (y).astype(float)
useful_time = (x).astype(float)

velocity, offsetVelocity = kinematics.velocity_calculation(useful, magneticField, spirals, length)
position = kinematics.position_calculation(position, velocity, dt)

fig2,(ax,ax1) = plt.subplots(2,1)

ax.plot(useful_time, position)
ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Posición (m)')
ax.grid()

ax1.plot(useful_time, velocity)
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Velocidad (m/s)')
ax1.grid()
plt.show()