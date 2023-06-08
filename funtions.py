import numpy as np
import time
import pandas as pd
from scipy.fftpack import fft, ifft, fftfreq
import scipy.signal as sg

def constants_characterization(A1, A2, t1, t2, mass):
    #Decrecimiento logaritmico (Solo para sistemas subamortiguados) ------
    delta = np.log(A1/A2) #Valor del decrecimiento logaritmico
    xi = 1/(np.sqrt(1+(2*np.pi/delta)**2)) #Factor de amortiguamiento
    omega_n = delta/(xi*(t2-t1)) #Frecuencia natural
    k = mass*(omega_n)**2 #constante de elasticidad
    c = 2*xi*np.sqrt(k*mass) #constante de amortiguamiento

    print('La constante de elasticidad es: ', k)
    print('La constante de amortiguamiento es: ', c)
    print('La frecuencia angular natural es', omega_n)
    print('El factor de amortiguamiento es', xi)

def noise_offset(values, arduinoSerial):
    #Eliminacion del offset de la señal ----------------------------------
    noise = 0 #Suma de los valores de la señal
    print('Eliminando el offset de la señal')
    for i in range(values):
        try:
            sensorValue = (float(arduinoSerial.readline().decode('utf-8'))) #Lectura de la señal
        except:
            sensorValue = 0
        noise = noise + sensorValue
    noiseMean = noise/values #Valor medio de la señal
    return noiseMean

def data_reading(arduinoSerial, noiseMean, n):
    #Lectura de los datos ------------------------------------------------
    voltage = np.repeat(0.0,n) # Voltaje
    tSpan = np.repeat(0.0,n) # Vector de tiempo con inicio en 0
    tFix = np.repeat(0.0,n) # Vector de tiempo sin traslacion
    print('Comienza la medicion')
    for i in range(len(voltage)):
        try:
            voltage[i] = (float(arduinoSerial.readline().decode('utf-8'))) - noiseMean #Lectura de los datos
            tFix[i] = time.time()
            tSpan[i] = time.time() - tFix[0]
        except:
            print('Entra')
            voltage[i] = voltage[i-1]
            tFix[i] = time.time()
            tSpan[i] = time.time() - tFix[0]
    print('Termina')
    totalTime = tFix[n-1]-tFix[0] #tiempo total de lectura
    print(f'Tiempo de medicion: {round(totalTime, 4)} seg.')
    timeValues = np.linspace(0.0, totalTime, n) #vector de timepo uniformemente espaciado)
    dt = totalTime/n #diferencial de tiempo
    return voltage, tFix, tSpan, totalTime, timeValues, dt

class kinematics: 
    #Cinematica del movimiento --------------------------------------------
    def velocity_calculation(voltage, magneticField, spirals, length):
        # Calculo de la velocidad y eliminacion del offset (uso para la integracion) -----
        velocity = (voltage)/(magneticField*spirals*length*1000) #Velocidad
        offsetVelocity = velocity - (sum(velocity)/len(velocity)) #correcion del valor constante, (se usa en la integracion)
        return velocity, offsetVelocity
    
    def position_calculation(offsetVelocity, dt, n):
        # Calculo de la posicion ---------------------------------------------------------
        position = np.repeat(0.0,n) #Posicion 
        offsetVelocity = offsetVelocity*dt #Discretizacion del area
        position[0] = offsetVelocity[0]
        # Integracion
        for i in range(len(offsetVelocity)-1):
            position[i+1] = position[i] + offsetVelocity[i+1]          
        return position
    
    def acceleration_calculation(velocity, timeValues, n):
        # Calculo de la aceleracion -------------------------------------------------------
        acceleration = np.repeat(0.0,n) # Aceleracion
        acceleration[0] = 0
        # Derivacion
        for i in range(len(velocity)-1):
            acceleration[i+1] = (velocity[i+1] - velocity[i]) / (timeValues[i+1] - timeValues[i])
        return acceleration

class fourierAnalysis:
    #Analisis de la señal ---------------------------------------
    def fourier_transform(signal, n, dt):
        # Transformada de Fourier DFT -----------------------------------------------------
        yfft = fft(signal)/n #DFT
        frq = fftfreq(n, dt) #Dominio de las frecuencias
        fHz = frq[np.where(abs(yfft.imag) == max(abs(yfft.imag)))][0] #Frecuencia de mayor amplitud en el espectro
        return frq, yfft, fHz

    def signal_decomposition(harmonics, signal, n, yfft, frq):
        # Descomposicion en Series de Fourier ---------------------------------------------
        # "Signal": señal a descomponer
        maximums = np.flip(np.sort(abs(yfft.imag))) # Ordena las frecuencias de mayor a menor
        ysfft = [] #Vector con la DFT para cada armonico
        harfhz = [] #Vector con la frecuencia de cada armonico
        signals = [] #Vector con los valores de cada armonico
        for i in range(harmonics):
            ysfft.append(fft(signal)/n) #DFT para toda la señal
            harfhz.append(frq[np.where(abs(ysfft[i].imag) == maximums[i*2])][0])
            ysfft[i][np.where(abs(ysfft[i].imag) != maximums[i*2])] = 0
            signals.append(ifft(ysfft[i])*n)
        return ysfft, signals, harfhz
    
class filters:
    # Filtros digitales -----------------------------------------------
    def butterworth(signal, order, frqCut, frqSamp):
        # Filtro digital butterworth -------------------------------------------------------
        if (order == 1): # Filtro de primer orden
            wn = frqCut/frqSamp #Frecuencia critica (frecuencia de corte / frecuencia de muestreo)
            b, a = sg.butter(order, wn, 'low') #Coeficientes del filtro (b:numerador, a: denominador)
            yn = np.repeat(0.0, len(signal))
            ynn, xn, xnn = (0, 0, 0)
            for i in range(len(signal)):
                yn[i] = b[0]*xn + b[1]*xnn - a[1]*ynn 
                ynn = yn[i]
                xnn = xn
                xn = signal[i]
        elif (order == 2): # Filtro de segundo orden
            wn = frqCut/frqSamp
            b, a = sg.butter(order, wn, 'low')
            yn = np.repeat(0.0, len(signal))
            xn, xnn, xnnn, ynn, ynnn = (0, 0, 0, 0, 0)
            for i in range(len(signal)):
                yn[i] = b[0]*xn + b[1]*xnn + b[2]*xnnn - a[1]*ynn - a[2]*ynnn
                xnnn = xnn
                xnn = xn
                xn = signal[i]
                ynnn = ynn
                ynn = yn[i]
        return yn

class accelerometer_comparison:
    # Comparacion con los datos de un acelerometro ----------------------
    def data_storage(timeValues, acceleration):
        # Almacenamineto de los datos del sismometro -------------------------------------
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
        # Lectura de los datos almacenados del sismometro --------------------------------
        data = pd.read_csv("Datos/datosSismometro.csv")
        acceleration = data["Aceleracion" + str(reading+1)].to_numpy().astype(float)
        timeValues = data["Tiempo" + str(reading+1)].to_numpy().astype(float)
        return timeValues, acceleration

    def accelerometer_reading():
        # Lectura de los datos del acelerometro -----------------------------------------
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