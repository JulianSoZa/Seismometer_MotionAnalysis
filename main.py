import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*

COM = 'COM7'
arduinoSerial = serial.Serial(COM, 19200)

n = 1600  # Número de muestras

voltage = np.repeat(0.0,n)
position = np.repeat(0.0,n)
velocity = np.repeat(0.0,n)
acceleration = np.repeat(0.0,n)

accelerationTz = np.repeat(0.0,n)
positionTz = np.repeat(0.0,n)

samples = np.linspace(0.0,n, n)
tSpan = np.repeat(0.0,n)
tFix = np.repeat(0.0,n)

magneticField = 8E-3
spirals = 1824
length = 7E-3

mass = 12E-3

#Filtro del ruido ---------------------------------------------------
values = 400
noiseMean = noise_filter(values, arduinoSerial)

#Lectura de los datos ------------------------------------------------
voltage, tFix, tSpan, totalTime, timeValues, dt = data_reading(voltage, arduinoSerial, noiseMean, tFix, tSpan, n)

#Cinematica del movimiento ------------------------------------------

velocity, offsetVelocity = kinematics.velocity_calculation(voltage, magneticField, spirals, length)

position = kinematics.position_calculation(position, offsetVelocity, dt)

acceleration = kinematics.acceleration_calculation(acceleration, velocity, timeValues)

##Transformada de Fourier DFT --------------------------------------

frq, yfft = fourierAnalysis.fourier_transform(velocity, n, dt)

##Filtros digitales -----------------------------------------------

order = 2
velocityTz = filters.z_transform(velocity, order)

frqTz, yfftTz = fourierAnalysis.fourier_transform(velocityTz, n, dt)

accelerationTz = kinematics.acceleration_calculation(accelerationTz, velocityTz, timeValues)

offsetVelocityTz = velocityTz - (sum(velocityTz)/len(velocityTz))

positionTz = kinematics.position_calculation(positionTz, offsetVelocityTz, dt)

harmonics = 3

ysfftTz, signals, harfhzTz = fourierAnalysis.signal_decomposition(harmonics, velocityTz, n, yfftTz, frqTz)


#Dinamica del movimiento -----------------------------------------------------


# Comparación con el acelerómetro

accelerometer_comparison(timeValues, acceleration)

#Graficar ----------------------------------------------------------

fig,(ax,ax1) = plt.subplots(2,1)

ax.plot(timeValues, velocity, label='Velocidad')
ax.plot(timeValues, velocityTz, label='Velocidad - Flltro')
ax.legend()
ax.scatter(timeValues, velocity, s = 12, c = 'black')
ax.set_xlim(0, 3)
ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Velocidad (m/s)')
ax.grid()

ax1.vlines(frq, 0, abs(yfft.imag), colors ='tab:blue', label='Sin Filtro')
ax1.vlines(frqTz, 0, abs(yfftTz.imag), colors = 'tab:orange', label='Filtrada')
ax1.legend()
ax1.set_xlim(0, max(frq))
ax1.set_xlabel('Frecuencia (Hz)')
ax1.set_ylabel('F(w)')
ax1.grid()
plt.tight_layout()

fig2,(ax2,ax3) = plt.subplots(2,1)

#ax1.plot(timeValues, acceleration)
ax2.plot(timeValues, accelerationTz)
ax2.set_xlim(0, 3)
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Aceleracion (m/s^2)')
ax2.grid()

#ax3.plot(timeValues, position)
ax3.plot(timeValues, positionTz)
ax3.set_xlim(0, 3)
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylabel('Posición (m)')
ax3.grid()
plt.tight_layout()

fig3,(ax4) = plt.subplots(1,1)
for i in range(harmonics):
    ax4.plot(timeValues, signals[i].real, label = str(round(harfhzTz[i], 3))+' Hz')
ax4.legend()
ax4.set_xlim(0, 3)
ax4.set_xlabel('Tiempo (s)')
ax4.set_ylabel('velocidad')
ax4.grid()

plt.tight_layout()
plt.show()
