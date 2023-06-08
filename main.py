import numpy as np
import serial
import matplotlib.pyplot as plt
from funtions import*

# --------------- Comunicacion con el arduino ------------
COM = 'COM10' # Puerto del arduino
arduinoSerial = serial.Serial(COM, 19200) # Se establece la comunicacion. (Puerto, Baudios)
##

# --------------- Condiciones / Variables iniciales ----------
n = 1600  # Numero de muestras

samples = np.linspace(0.0,n, n) # Vector de muestreo (Cantidad de muestras)

magneticField = 8E-3 # Magnitud del campo magnetico (T)
spirals = 1824 # Numero de espiras
length = 3E-3 # Diametro del iman

mass = 12E-3 # masa del sistema
##

# ------------- Correccion del offSet ----------------
values = 400 # Muestras 
noiseMean = noise_offset(values, arduinoSerial) # media de los datos
##

# ------------- Lectura de los datos -----------------
voltage, tFix, tSpan, totalTime, timeValues, dt = data_reading(arduinoSerial, noiseMean, n)
# "voltage": Voltaje
# "tFix": Vector de tiempo sin traslacion
# "tSpan": Vector de tiempo con inicio en 0
# "totalTime": tiempo total de lectura
# "timeValues": vector de timepo uniformemente espaciado)
# "dt": diferencial de tiempo
##

# ------------- Cinematica del movimiento -------------- 
velocity, offsetVelocity = kinematics.velocity_calculation(voltage, magneticField, spirals, length)
# "velocity": Velocidad
# "offsetVelocity": correcion del valor constante, (se usa en la integracion)  

position = kinematics.position_calculation(offsetVelocity, dt, n) #Posicion

acceleration = kinematics.acceleration_calculation(velocity, timeValues, n) #Aceleracion
##

# -------------- Analisis de la señal ------------------
# Transformada de Fourier DFT
frq, yfft, fHz = fourierAnalysis.fourier_transform(velocity, n, dt)
print(f'La frecuencia de mayor amplitud sin filtado es: {round(fHz, 2)} Hz.')
# "frq": dominio de las frecuencias
# "yfft": DFT
# "fHz" : Frecuencia de mayor amplitud en el espectro

# ------------- Filtros digitales ----------------------
#Filtro butterworth
order = 2 # orden del filtro
frqCut = 40 #frecuencia de corte
frqSamp = 400 #frecuencia de muestreo
#frqSamp = 1/dt

velocityTz = filters.butterworth(velocity, order, frqCut, frqSamp)

accelerationTz = kinematics.acceleration_calculation(velocityTz, timeValues, n)

offsetVelocityTz = velocityTz - (sum(velocityTz)/len(velocityTz))

positionTz = kinematics.position_calculation(offsetVelocityTz, dt, n)

frqTz, yfftTz, fHzTz = fourierAnalysis.fourier_transform(accelerationTz, n, dt)
print(f'La frecuencia de mayor amplitud con filtado es: {round(fHzTz, 2)} Hz.')
##

# ------------- Descomposicion en Series de Fourier ---------------
harmonics = 3 # Numeros de armonicos de la descomposicion

ysfftTz, signals, harfhzTz = fourierAnalysis.signal_decomposition(harmonics, accelerationTz, n, yfftTz, frqTz)
# "ysfftTz": Vector con la DFT para cada armonico
# "signals": Vector con los valores de cada armonico
# "harfhzTz": Vector con la frecuencia de cada armonico
## 

# ------------- Almacenamiento de los datos -------------------
accelerometer_comparison.data_storage(timeValues, acceleration)
##

#-------------- Graficos de los resultados ------------

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
ax4.set_ylabel('Aceleracion (m/s^2)')
ax4.grid()

plt.tight_layout()
plt.show()
