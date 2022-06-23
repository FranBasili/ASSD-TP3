# Modulacion Delta Sigma - Codificador
# entrada x(t), salida: y[n]
# propuesta:edelros@espol.edu.ec
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as waves

# INGRESO DE DATOS
muestreo, sonido = waves.read('muestra.wav')
fs=44100; fos=fs*16

# rango de observación en segundos
inicia = 0
termina = 0.002

# PROCEDIMIENTO
# Codificar Sigma-Delta
deltaY = 0.1*np.max(sonido) 
#deltaT = 1/muestreo 
deltaT = 1/fos

# Extrae solo una porcion del sonido
donde = int(inicia/deltaT)
# tiempo muestreo de la señal analógica
t = np.arange(inicia,termina,deltaT) 
k = len(t) 

muestra = np.copy(sonido[donde:donde+k])

# Señal Digital
xdigital = np.zeros(k, dtype=float) 
ysalida = np.zeros(k, dtype=int) 

for i in range(1,k):
    diferencia = muestra[i]-xdigital[i-1]
    signo= 1 if diferencia>0 else -1
    xdigital[i] = xdigital[i-1]+signo*deltaY
    ysalida[i] = signo
    
#parametros=np.array([deltaT,deltaY,k])

# SALIDA
#print('parámetros:[deltaT, deltaY, k]')
#print(parametros)
print('datos:')
print(ysalida)
#np.savetxt('deltasigma_parametros.txt',parametros)
#np.savetxt('deltasigma_datos.txt',ysalida,fmt='%i')
#print('... archivos.txt guardados ...')

verdesde=0
verhasta=90

# Graficar
plt.figure(1)       # define la grafica
plt.suptitle('Codificador Delta-Sigma')

plt.subplot(211)    # grafica de 2x1 y subgrafica 1
plt.ylabel('x(t), x[n]')
plt.plot(t[verdesde:verhasta],muestra[verdesde:verhasta], 'g')
plt.step(t[verdesde:verhasta],xdigital[verdesde:verhasta], where='post',color='m') # Puntos x[n]

plt.subplot(212)    # grafica de 2x1 y subgrafica 2
plt.ylabel('y[n]')
#plt.plot(ysalida, 'bo')     # Puntos y[n]
plt.axis((verdesde,verhasta,-1.1,1.1))
puntos=np.arange(verdesde,verhasta,1)     #posicion eje x para escalon
plt.step(puntos[verdesde:verhasta],ysalida[verdesde:verhasta], where='post')

plt.show()