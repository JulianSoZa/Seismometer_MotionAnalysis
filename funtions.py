import numpy as np
import time
import pandas as pd
from scipy.fftpack import fft, ifft, fftfreq
import scipy.signal as signal

def constants_characterization(A1, A2, t1, t2, mass):
    delta = np.log(A1/A2)
    xi = 1/(np.sqrt(1+(2*np.pi/delta)**2))
    omega_n = delta/(xi*(t2-t1))
    k = mass*(omega_n)**2
    c = 2*xi*np.sqrt(k*mass)

    print('La constante de elasticidad es: ', k)
    print('La constante de amortiguamiento es: ', c)
    print('La frecuencia angular natural es', omega_n)
    print('El factor de amortiguamiento es', xi)

def noise_offset(values, arduinoSerial): #Filtro del ruido ---------------------------------------------------
    noise = 0
    print('Eliminando el offset de la señal')
    for i in range(values):
        try:
            sensorValue = (float(arduinoSerial.readline().decode('utf-8')))
        except:
            sensorValue = 0
        noise = noise + sensorValue
    noiseMean = noise/values
    return noiseMean

def data_reading(voltage, arduinoSerial, noiseMean, tFix, tSpan, n): #Lectura de los datos ------------------------------------------------
    print('Comienza la medicion')
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
    print(f'Tiempo de medicion: {round(totalTime, 4)} seg.')
    timeValues = np.linspace(0.0, totalTime, n)
    dt = totalTime/n

    return voltage, tFix, tSpan, totalTime, timeValues, dt

class kinematics: #Cinematica del movimiento ------------------------------------------
    def velocity_calculation(voltage, magneticField, spirals, length):
        velocity = (voltage)/(magneticField*spirals*length*1000)
        offsetVelocity = velocity - (sum(velocity)/len(velocity))
        return velocity, offsetVelocity
    
    def position_calculation(position, offsetVelocity, dt):
        offsetVelocity = offsetVelocity*dt
        position[0] = offsetVelocity[0]
        for i in range(len(offsetVelocity)-1):
            position[i+1] = position[i] + offsetVelocity[i+1]          
        return position
    
    def acceleration_calculation(acceleration, velocity, timeValues):
        acceleration[0] = 0
        for i in range(len(velocity)-1):
            acceleration[i+1] = (velocity[i+1] - velocity[i]) / (timeValues[i+1] - timeValues[i])
        return acceleration

class fourierAnalysis:
    def fourier_transform(velocity, n, dt):
        yfft = fft(velocity)/n
        frq = fftfreq(n, dt)
        fHz = frq[np.where(abs(yfft.imag) == max(abs(yfft.imag)))][0]
        return frq, yfft, fHz

    def signal_decomposition(harmonics, velocity, n, yfft, frq):
        maximums = np.flip(np.sort(abs(yfft.imag)))
        ysfft = []
        harfhz = []
        signals = []
        for i in range(harmonics):
            ysfft.append(fft(velocity)/n)
            harfhz.append(frq[np.where(abs(ysfft[i].imag) == maximums[i*2])][0])
            ysfft[i][np.where(abs(ysfft[i].imag) != maximums[i*2])] = 0
            signals.append(ifft(ysfft[i])*n)
        return ysfft, signals, harfhz
    
class filters:
    def butterworth(velocity, order, frqCut, frqSamp): 
        if (order == 1): # Filtro de primer orden
            wn = frqCut/frqSamp
            b, a = signal.butter(order, wn, 'low')
            yn = np.repeat(0.0, len(velocity))
            ynn, xn, xnn = (0, 0, 0)
            for i in range(len(velocity)):
                yn[i] = b[0]*xn + b[1]*xnn - a[1]*ynn 
                ynn = yn[i]
                xnn = xn
                xn = velocity[i]
        elif (order == 2): # Filtro de segundo orden
            wn = frqCut/frqSamp
            b, a = signal.butter(order, wn, 'low')
            yn = np.repeat(0.0, len(velocity))
            xn, xnn, xnnn, ynn, ynnn = (0, 0, 0, 0, 0)
            for i in range(len(velocity)):
                yn[i] = b[0]*xn + b[1]*xnn + b[2]*xnnn - a[1]*ynn - a[2]*ynnn
                xnnn = xnn
                xnn = xn
                xn = velocity[i]
                ynnn = ynn
                ynn = yn[i]
        return yn

class accelerometer_comparison:
    def data_storage(timeValues, acceleration):
        try:
            df = pd.read_csv('Datos/datosSismometro.csv')
            i = str(int(len(df.axes[1])/2+1))
            df['Tiempo'+ i] = timeValues
            df['Aceleracion'+ i] = acceleration
            df.to_csv('Datos/datosSismometro.csv', index=False)
        except:
            datos= {'Tiempo1':timeValues,'Aceleracion1':acceleration}
            df = pd.DataFrame(datos)
            df.to_csv('Datos/datosSismometro.csv', index=False)

    def data_extraction(reading):
        data = pd.read_csv("Datos/datosSismometro.csv")
        acceleration = data["Aceleracion" + str(reading+1)].to_numpy().astype(float)
        timeValues = data["Tiempo" + str(reading+1)].to_numpy().astype(float)
        return timeValues, acceleration

    def accelerometer_reading():
        try:
            data = pd.read_csv("Datos/Muestras/1-5.csv",sep=';',decimal=",").to_numpy()
            acceleration = data[:,1].astype(float)
            acceleration = acceleration - sum(acceleration)/len(acceleration)
            timeValues = data[:,0].astype(float)
        except:
            print("""Error en accelerometer_reading(): 
                  * Por favor, asegurese que existen datos por leer *""")
            timeValues, acceleration = (0, 0)
        return timeValues, acceleration