import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq

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
        print('La frecuencia de mayor amplitud es: ', fHz)
        return frq, yfft

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

class filters:
    def z_transform(velocity, order): 
        if (order == 1): # Filtro de primer orden
            yn = np.repeat(0.0, len(velocity))
            ynn, xn = (0, 0)
            for i in range(len(velocity)):
                yn[i] = 0.4*ynn + 0.6*xn # 40 hz
                #yn = 0.07079*ynn + 0.9292*xn
                ynn = yn[i]
                xn = velocity[i]
        elif (order == 2): # Filtro de segundo orden
            yn = np.repeat(0.0, len(velocity))
            xn, xnn, xnnn, ynn, ynnn = (0, 0, 0, 0, 0)
            for i in range(len(velocity)):
                yn[i] = 0.020083365564211*xn + 0.040166731128422*xnn + 0.020083365564211*xnnn + 1.561018075800718*ynn - 0.641351538057563*ynnn # relacion 0.1
                #yn[i] = 0.097631072937818*xn + 0.195262145875635*xnn + 0.097631072937818*xnnn + 0.942809041582063*ynn - 0.333333333333333*ynnn # relacion 0.25
                #yn[i] = 0.067455273889072*xn + 0.134910547778144*xnn + 0.067455273889072*xnnn + 1.142980502539901*ynn - 0.412801598096189*ynnn # relacion 0.2
                xnnn = xnn
                xnn = xn
                xn = velocity[i]
                ynnn = ynn
                ynn = yn[i]
        return yn


def accelerometer_comparison(timeValues, acceleration):

    try:
        df = pd.read_csv('../DATOS/datosSismometro.csv')
        i = str(len(df.axes[1])/2+1)
        df['Tiempo'+ i] = timeValues
        df['Aceleracion'+ i] = acceleration
        df.to_csv('../DATOS/datosSismometro.csv', index=False)
    except:
        datos= {'Tiempo1':timeValues,'Aceleracion1':acceleration}
        df = pd.DataFrame(datos)
        df.to_csv('../DATOS/datosSismometro.csv', index=False)

def data_analis(lectura):
    datos = pd.read_csv("../DATOS/datosSismometro.csv")
    y = datos["Aceleracion" + str(lectura+1)].to_numpy()
    x = datos["Tiempo" + str(lectura+1)].to_numpy()
    dt = x[int(len(x))-1]/len(x)
    useful = (y).astype(float)
    useful_time = (x).astype(float)
    return useful_time, useful