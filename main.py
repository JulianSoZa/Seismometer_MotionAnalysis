import numpy as np
import serial
import time
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from funtions import*

COM = 'COM7'
arduinoSerial = serial.Serial(COM, 9600)

n = 400  # Número de muestras

voltage = np.repeat(0.0,n)
position = np.repeat(0.0,n)
velocity = np.repeat(0.0,n)
acceleration = np.repeat(0.0,n)
samples = np.linspace(0.0,n, n)
tSpan = np.repeat(0.0,n)
tFix = np.repeat(0.0,n)

magneticField = 8E-3
spirals = 1515
length = 2*3.14159*7E-3

mass = 12E-3

#Filtro del ruido ---------------------------------------------------
values = 400
noiseMean = noise_filter(values, arduinoSerial)

#Lectura de los datos ------------------------------------------------
voltage, tFix, tSpan, totalTime, timeValues, dt = data_reading(voltage, arduinoSerial, noiseMean, tFix, tSpan, n)

#Cinematica del movimiento ------------------------------------------

velocity, offsetVelocity = kinematics.velocity_calculation(voltage, magneticField, spirals, length)

position = kinematics.position_calculation(position, offsetVelocity)

acceleration = kinematics.acceleration_calculation(acceleration, velocity, timeValues)

##Transformada de Fourier DFT --------------------------------------

frq, Y = fourier_transform(velocity, n, dt)

#Dinamica del movimiento -----------------------------------------------------



#Graficar ----------------------------------------------------------

fig,((ax,ax1)) = plt.subplots(2,1)

ax.plot(timeValues, velocity)
ax.scatter(timeValues, velocity, s = 12, c = 'black')
ax.set_xlim(0, 1)
ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Velocidad (m/s)')
ax.grid()

ax1.plot(timeValues, acceleration)
ax1.set_xlim(0, 1)
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Aceleracion (m/s^2)')
ax1.grid()

fig2,((ax2,ax3)) = plt.subplots(2,1)

ax2.vlines(frq, 0, abs(Y.imag))
ax2.set_xlim(0, max(frq))
ax2.set_xlabel('Frecuencia (Hz)')
ax2.set_ylabel('F(w)')
ax2.grid()

ax3.plot(timeValues, position)
ax3.set_xlim(0, 1)
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylabel('Posición (m)')
ax3.grid()

plt.tight_layout()
plt.show()
