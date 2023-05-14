"""t = np.linspace(0.0, tiempoT, n)  # Intervalo de tiempo en segundos
t2 = np.linspace(0.0, tiempoT + 9*dt, 10*n)
y = np.append(voltage * np.blackman(n), np.zeros(9 * n))

Y = fft(y) / 10*n  # Normalizada
frq = fftfreq(10*n, dt)  # Recuperamos las frecuencias"""

values = 400
rui = 0

for i in range(values):
    print(i)

print(rui)