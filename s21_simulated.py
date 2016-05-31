import numpy as np 
import matplotlib.pyplot as plt 
import h5py
import math
from scipy.fftpack import fft,ifft
import time
'''Inicializacion de valores'''
npos_r = 8
posx_r = np.linspace(0,1.5,npos_r)
Et = 0.0005


'''Parametros del sistema'''
dr = 0.01 # resolucion de la grilla
c0 = 3. * (10**8)
fre1 = np.linspace(15.,16.,1601) * (10**9)
dr = c0/(2*10**9)

'''Generar grilla'''
y	=	np.arange(4,12+dr,dr)
x 	= 	np.arange(0,1.5+dr,dr)
#y   =   np.linspace(4,8,1601)
#x   =   np.linspace(0,1.5,1601)
maxtam = c0/(2*(fre1[1] - fre1[0]))
tam = np.arange(0,maxtam+dr,dr)
grid1 = np.where(tam>=4)
grid2 = np.where(tam<=12+dr)
grid1 = grid1[0] ; grid2 = grid2[0]
print grid1[0]
print grid2[-1]
indexg = np.arange(grid1[0],grid2[-1]+1)
#print indexg.size
'frecuencia'
#fre1 = fre1[indexg]
'''Escena'''
p0 = np.zeros([y.size,x.size], dtype = complex)
a = np.where(y >= 8.+dr*0)
#print 'importante ',a
a = a[0]

b = np.where(x >= 1.+dr*0)
b = b[0]

p0[a[0],b[0]] = 10**3 
print 'p0 ',np.flipud(p0)
Ft = 1
Fr = 1
'grilla'
R = np.zeros([npos_r,y.size,x.size],dtype = complex)
start = time.time()
for k in range(npos_r):
	for i in range(y.size):
		for j in range(x.size):
			R[k,i,j] = math.hypot(x[j] - posx_r[k] ,y[i] - 0)
end = time.time()
print end - start
print 'importante'
#print R 
Er = np.zeros([npos_r,y.size,x.size])
for k in range(npos_r):
	for i in range(y.size):
		for j in range(x.size):
			Er[k,i,j] = (Et * Ft * Fr * p0[i,j] * np.exp(4j*np.pi*fre1[i]*R[k,i,j] /c0) / R[k,i,j]**2)
#print 4j*np.pi*15*10*(10**9)/c0
print Er.size
Er = sum(Er,0)/npos_r
Er = np.flipud(Er)
plt.imshow(np.absolute(Er),aspect = 'auto',extent=[0, 1.5, 4, 12])
plt.colorbar()
#plt.plot(20*(np.log(np.absolute(Er))))
plt.show()
