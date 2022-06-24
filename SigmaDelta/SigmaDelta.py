# http://blog.espol.edu.ec/estg1003/sigma-delta-codificador-con-python/

from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wav
import scipy.signal as ss
import scipy.fft as fft

def graph(t, muestra, xi, ydig, ydecim, q):
    # n0=0; deltaN=int(len(t)/1000); nf=n0+deltaN
    # plt.figure(1)
    # plt.suptitle('Conversor $\Sigma \Delta$')

    # # Comp muestra con xi
    # plt.subplot(211)    
    # plt.ylabel('x(t), x[n]')
    # plt.plot(t[n0:nf],muestra[n0:nf], 'g')
    # plt.step(t[n0:nf],xi[n0:nf], where='post',color='m')

    # # Grafica de y[n]
    # plt.subplot(212)
    # plt.ylabel('y[n]')
    # plt.axis((n0,nf,-0.1,1.1))
    # puntos=np.arange(n0,nf,1)
    # puntosdec=np.arange(n0,nf,q); puntosdec=np.append(puntosdec,nf)
    # plt.step(puntos, ydig[n0:nf], where='post', label='y[n]', color='b')
    # plt.step(puntosdec[n0:nf], ydecim[n0:int(nf/q)+2], where='post', label='decimado', color='r')
    # plt.legend()
    # plt.tight_layout()
    # plt.show()

    return 0

def FAA(fb, muestra, t, verbose=False):
    b, a = ss.butter(N=1, Wn=np.pi*fb, btype='lowpass', analog=True)
    _, yout, _ = ss.lsim((b, a), U=muestra, T=t)

    if verbose:
        plt.plot(t, yout)
        plt.plot(t, muestra)
        plt.show()

    return yout

def conversor(k, filtered, deltaY, t, noise=False, fs=44100):
    yanalog=np.zeros(k, dtype=float)
    ydig=np.zeros(k, dtype=int)
    xi=np.zeros(k, dtype=float)
    diferencia=np.zeros(k, dtype=float)

    for i in range(0,k):
        diferencia[i] = filtered[i]-yanalog[i-1] if i>0 else filtered[i]    #Sumador
        xi[i] = xi[i-1]+diferencia[i] if i>0 else diferencia[i]             #Integrador                                               #Ruido   
        ydig[i] = 1 if xi[i]>0 else 0                                       #Cuantizador
        yanalog[i]= deltaY if ydig[i]==1 else -deltaY                       #Realimentación
    return ydig, xi, diferencia

def downsampler(q, ydig):
    #ss.decimate() solo decima hasta orden 13, por lo que divido en rep decimaciones 
    # de orden 8 y una de orden res

    rep=int(np.log(q)/np.log(8)); res=int(q/8**rep)

    if res>1:
        decimated=ss.decimate(ydig,res)
    else:
        decimated=ydig

    for i in range(rep):
        decimated=ss.decimate(decimated,8)
    return decimated

def noiseShapping(k):
    yanalog=np.zeros(k)
    ydig=np.zeros(k)
    xi=np.zeros(k)
    diferencia=np.zeros(k)
    
    ruido=np.random.random(size=k)-0.5
    # print(ruido)

    for i in range(k):
        diferencia[i] = -yanalog[i-1] if i>0 else 0    #Sumador
        xi[i] = xi[i-1]+diferencia[i] if i>0 else diferencia[i]             #Integrador  
        ydig[i] = xi[i]+ruido[i]                                            #Cuantizador
        yanalog[i]= ydig[i]                                                 #Realimentación
    return ydig, ruido   

fb=22000; fs=44100; q=4; fos=fs*q

archivo = 'SigmaDelta\muestra.wav'
fdata, data = wav.read(archivo)
deltaY = 1; deltaT = 1/fos 

tfinal=(len(data)/fdata)/q
t = np.arange(0,tfinal,deltaT)

#FAA
filtered = FAA(fb, data/np.max(data), t, verbose=False)
k=len(filtered)

#Conversor
ydig, xi, diff= conversor(k, filtered, deltaY, t, noise=False, fs=fs)

#DownSampler
ysalida=np.zeros(int(k/q))
ysalida = downsampler(q=q, ydig=ydig)

n0=0; deltaN=int(len(t)/500); nf=n0+deltaN
puntos=(np.arange(n0+1,nf+1,1))/fos

#Grafico entrada-salida
plt.plot(t[n0:nf],filtered[n0:nf], 'g', label= "x(t)")
plt.step(puntos, (ydig[n0:nf]-0.5)*2, where='mid', color='k', alpha=0.2, label= "y[n]")
plt.legend(loc='lower right')
plt.show()

# Grafico I/O+integrador y sumador 
plt.plot(t[n0:nf],filtered[n0:nf], 'green', label= "x(t)")
plt.step(t[n0:nf],xi[n0:nf], 'orange', where='pre', label= "xi(n)")
plt.step(t[n0:nf],diff[n0:nf], 'blue', where='pre', label= "Sum(n)")
plt.step(puntos, (ydig[n0:nf]-0.5)*4, where='pre', color='k', alpha=0.2, label= "y[n]")
plt.legend()
plt.show()

# Grafico Output+decimador
puntos=(np.arange(n0,nf,1))/fos
plt.plot(t[n0:nf],filtered[n0:nf], 'green', label= "x(t)")
plt.step(puntos, (ydig[n0:nf]-0.5)*1.7, where='mid', color='k', alpha=0.2, label= "y[n]")
plt.step(puntos*fos/fs, (ysalida[n0:nf]-0.5)*1.7, where='mid', color='blue', alpha=0.5, label= "Decimada")
plt.xlim(0,0.001)
plt.legend(loc='upper right')
plt.show()

# FFT
noisydata, ruido= noiseShapping(k)                              # Señal con ruido
noiseout_fft = 2.0 / k * np.abs(fft.fft(noisydata)[1:k // 2])   # Espectro salida
noise_fft = 2.0 / k * np.abs(fft.fft(ruido)[1:k // 2])          # Espectro ruido
freqs = fft.fftfreq(k, 1 / fos)[1:k // 2]
plt.plot(freqs, noiseout_fft, 'g', label= "Noise(f)")           # Plot espectro salida
plt.legend()
plt.show()

plt.title("NTF")
plt.plot(freqs,noiseout_fft/noise_fft, label="NTF")             # Plot transferencia
plt.legend()
plt.show() 