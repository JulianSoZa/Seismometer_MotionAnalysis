import numpy as np
import serial
import time
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from scipy.interpolate import make_interp_spline

COM = 'COM7'
arduinoSerial = serial.Serial(COM, 9600)

n = 400  # Número de intervalos

voltage = np.repeat(0.0,n)
position = np.repeat(0.0,n)
velocity = np.repeat(0.0,n)
acceleration = np.repeat(0.0,n)
samples = np.linspace(0.0,n, n)
tSpan = np.repeat(0.0,n)
tFix = np.repeat(0.0,n)

magneticField = 8E-3
espiras = 1515
length = 2*3.14159*7E-3

mass = 12E-3

#Filtro del ruido ---------------------------------------------------
values = 400
rui = 0

for i in range(values):
    try:
        sensorValue = (float(arduinoSerial.readline().decode('utf-8')))
    except:
        sensorValue = 0
    rui = rui + sensorValue

meanRuido = rui/values

#Lectura de los datos ------------------------------------------------
print('Comienza')
for i in range(len(voltage)):
    try:
        voltage[i] = (float(arduinoSerial.readline().decode('utf-8'))) - meanRuido
        tFix[i] = time.time()
        tSpan[i] = time.time() - tFix[0]
    except:
        print('Entra')
        voltage[i] = voltage[i-1]
        tFix[i] = time.time()
        tSpan[i] = time.time() - tFix[0]
print('Termina')
tiempoT = tFix[n-1]-tFix[0]
print('Tiempo: ', tiempoT)

#Cinematica del movimiento ------------------------------------------
t = np.linspace(0.0, tiempoT, n)  # Intervalo de tiempo en segundos
dt = tiempoT/n  # Espaciado, 16 puntos por período
velocity = (voltage)/(magneticField*espiras*length*1000)
constP = sum(velocity)/len(velocity)
velocityN = velocity - constP

position[0] = velocityN[0]
for i in range(len(velocityN)-1):
    position[i+1] = position[i] + velocityN[i+1]

acceleration[0] = 0
for i in range(len(velocity)-1):
     acceleration[i+1] = (velocity[i+1] - velocity[i]) / (t[i+1] - t[i])

##Transformada de Fourier DFT --------------------------------------

y = velocity

Y = fft(y) / n  # Normalizada
frq = fftfreq(n, dt)  # Recuperamos las frecuencias

fHz = frq[np.where(abs(Y.imag) == max(abs(Y.imag)))][0]

print('La frecuencia de mayor amplitud es: ', fHz)

#Dinamica del movimiento



#Graficar ----------------------------------------------------------

fig,((ax,ax1)) = plt.subplots(2,1)

ax.plot(t, velocity)
ax.scatter(t, velocity, s = 12, c = 'black')
ax.set_xlim(0, 1)
ax.set_xlabel('Tiempo (s)')
ax.set_ylabel('Velocidad (m/s)')
ax.grid()

ax1.plot(t, acceleration)
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

ax3.plot(t, position)
ax3.set_xlim(0, 1)
ax3.set_xlabel('Tiempo (s)')
ax3.set_ylabel('Posición (m)')
ax3.grid()

plt.tight_layout()
plt.show()