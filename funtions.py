import numpy as np
import time
from scipy.fftpack import fft, fftfreq

def noise_filter(values, arduinoSerial): #Filtro del ruido ---------------------------------------------------
    noise = 0
    for i in range(values):
        try:
            sensorValue = (float(arduinoSerial.readline().decode('utf-8')))
        except:
            sensorValue = 0
        noise = noise + sensorValue
    noiseMean = noise/values

    return noiseMean

def data_reading(voltage, arduinoSerial, noiseMean, tFix, tSpan, n): #Lectura de los datos ------------------------------------------------
    print('Comienza')
    for i in range(len(voltage)):
        try:
            voltage[i] = (float(arduinoSerial.readline().decode('utf-8'))) - noiseMean
            tFix[i] = time.time()
            tSpan[i] = time.time() - tFix[0]
        except:
            print('Entra')
            voltage[i] = voltage[i-1]
            tFix[i] = time.time()
            tSpan[i] = time.time() - tFix[0]
    print('Termina')
    totalTime = tFix[n-1]-tFix[0]
    print('Tiempo: ', totalTime)
    timeValues = np.linspace(0.0, totalTime, n)  # Intervalo de tiempo en segundos
    dt = totalTime/n  # Espaciado, 16 puntos por período

    return voltage, tFix, tSpan, totalTime, timeValues, dt

class kinematics: #Cinematica del movimiento ------------------------------------------
    def velocity_calculation(voltage, magneticField, spirals, length):
        velocity = (voltage)/(magneticField*spirals*length*1000)
        offsetVelocity = velocity - (sum(velocity)/len(velocity))
        return velocity, offsetVelocity
    
    def position_calculation(position, offsetVelocity):
        position[0] = offsetVelocity[0]
        for i in range(len(offsetVelocity)-1):
            position[i+1] = position[i] + offsetVelocity[i+1]          
        return position
    
    def acceleration_calculation(acceleration, velocity, timeValues):
        acceleration[0] = 0
        for i in range(len(velocity)-1):
            acceleration[i+1] = (velocity[i+1] - velocity[i]) / (timeValues[i+1] - timeValues[i])
        return acceleration
    
def fourier_transform(velocity, n, dt):
    Y = fft(velocity) / n  # Normalizada
    frq = fftfreq(n, dt)  # Recuperamos las frecuencias
    fHz = frq[np.where(abs(Y.imag) == max(abs(Y.imag)))][0]
    print('La frecuencia de mayor amplitud es: ', fHz)
    return frq, Y