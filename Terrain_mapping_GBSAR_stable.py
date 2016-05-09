from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import time
import h5py

'''
RUTAS PARA LA LECTURA DE ARCHIVOS
'''
#FECHA_HORA = '26.04.2016_13.49.48'
FECHA_HORA = ''

#RUTA = '/home/cielosar/Desktop/Datos_hdf5/experimentos_%s/'
RUTA = './'
NOMBRE_ARCHIVO = 'datos_toma_%d.hdf5'
#RUTA_RESULTADOS = "/home/cielosar/Desktop/Resultados_paper/Resultados2"
NOMBRE_RESULTADOS = "/res_08_toma_%d.hdf5"

'''
ESPECIFICACIONES TECNICAS DEL MOTOR:
-LONG_VUELTA : la longitud de desplazamiento de una vuelta del motor a pasos, 66 mm
-PASOS_VUELTA : el numero de pasos por vuelta del motor a pasos
-FACTOR_LONG_PASOS : factor de conversion de una distancia a un numero de pasos
'''
LONG_VUELTA = 66
PASOS_VUELTA = 20000
FACTOR_LONG_PASOS = 1000 * PASOS_VUELTA / LONG_VUELTA

for n_toma in range(1,2):

	print RUTA + FECHA_HORA + NOMBRE_ARCHIVO %(n_toma)
	
	f = h5py.File(RUTA + FECHA_HORA + NOMBRE_ARCHIVO %(n_toma), 'r')
	dset = f['sar_dataset']

	xai = dset.attrs['xi']
	#xaf = dset.attrs['xf']
	dax = dset.attrs['dx']
	nx = dset.attrs['npos']
	xaf = xai + dax * (nx - 1)

	fre_min = dset.attrs['fi'] * 1E9
	fre_max = dset.attrs['ff'] * 1E9
	nfre = dset.attrs['nfre']
	df = (fre_max - fre_min) / (nfre - 1)

	print "Toma numero %d" %n_toma
	print fre_min
	print fre_max
	print nfre
	print xai
	print xaf
	print nx
	print dax

	'''
	Dimensiones de la grilla de la reconstruccion, xai y xaf representan
	las posiciones inicial y final de la antena respectivamente
	'''

	xi =  xai - 1.5
	xf =  xaf + 0.5
	yi =  10.0
	yf =  15.0

	'''
	Velocidad de la luz, inicialmente el R0 lo hemos considerado 0 para
	simplificar un poco los calculos
	'''

	c0  =  299792458.0
	R0 	=  0.0

	'''
	Tamano de los pixeles que buscamos reconstruir:
	-dx y dy representan el las dimensiones en Cross-range y range respectivamente.
	 Ambas de 0.01 m, muchisimo mas finas que las resoluciones teoricas.
	'''

	dx    =  0.01
	dy    =  0.01
	nposx =  int(np.ceil((xf - xi) / dx) + 1)
	nposy =  int(np.ceil((yf - yi) / dy) + 1)
	xf    =  xi + dx * (nposx - 1)
	yf    =  yi + dy * (nposy - 1)
	
	
	#fact  =  10
	#nr 	  =  nfre * fact
	##nr =  2 ** np.ceil(np.log2(nfre))
	
	dr0  = c0 / (2*df*nfre)
	fact = dr0 / min(dx,dy)
	if fact<1:
		fact = 1
	
	#nr =  nfre * fact
	nr = 2 ** int(np.ceil(np.log2(nfre*fact)))
	n  = np.arange(nr)
	B  = df*nr
	dr = c0 / (2*B)

	rn   = dr * n
	R    = dr * nr
	npos = nposx * nposy

	'''
	Lectura de datos
	'''

	s21 = np.empty([nx, nfre], dtype = np.complex64)
	s21 = dset[...]
	f.close()

	'''
	Inicializacion de datos:
	-xa: representa las posiciones de la antena
	-xn: representa las coordenadas en Cross-range de la grilla
	-yn: representa las coordenadas en Range de la grilla
	-Imagen: la imagen compleja
	-Rnk: el vector para almacenar las distancias de cada posicion de la antena
	      a cada posicion de la grilla
	'''

	xa = xai + dax*np.arange(nx)
	xn = xi + dx*np.arange(nposx)
	yn = yi + dy*np.arange(nposy)
	Imagen = np.zeros([nposx * nposy], dtype = complex)
	Rnk = np.zeros([nposx * nposy], dtype = float)

	'''
	-Primero se calculan las distancias Rnk.

	-Luego de haberle dado a Fn la forma de una ifft, se calculan los valores para
	cada posicion de las antenas en el riel. Aqui se hace un zero padding en la
	funcion ifft, para obtener resoluciones mas finas de las que han sido muestreadas.

	-Finalmente, se va formando la imagen, primero como un vector, interpolando
	los valores de Fn obtenidos para que se ajusten a las distancias Rnk de la grilla
	'''

	for k in range(0,nx):
		for y in range(0,nposy):
			for x in range(0, nposx):
				Rnk[x + y * nposx] = np.sqrt((xn[x]-xa[k])**2+(yn[y])**2)
		
		Fn0 = ifft(s21[k,:], n = nr)
		Fn = np.interp(Rnk - R0, rn, np.real(Fn0), period=R) + 1j * np.interp(Rnk - R0, rn, np.imag(Fn0), period=R)
		#fn_re = interpolate.interp1d(rn, np.real(Fn0),kind='cubic')
		#fn_im = interpolate.interp1d(rn, np.imag(Fn0),kind='cubic')
		#Fn = fn_re(Rnk - R0) + 1j * fn_im(Rnk - R0)
		Fn = np.exp(4j * np.pi * (fre_min/c0) * (Rnk - R0)) * Fn 
		
		Imagen += Rnk**2 * Fn

	'''
	Se hace un reshape de la Imagen para que corresponda a la forma de la grilla
	y se normaliza
	'''

	Imagen /= (nfre * (nx))
	Imagen = np.reshape(Imagen, (nposy, nposx))
	Imagen = np.flipud(Imagen)
	Imagen = np.absolute(Imagen)
	#Imagen *= (1.0/Imagen.max())
	#Imagen = 20 * np.log(Imagen)

	fig = plt.figure(1)
	im = plt.imshow(Imagen, cmap = 'jet', aspect = 'auto', extent = [xi,xf,yi,yf], interpolation = 'none')
	cbar = plt.colorbar(im, orientation = 'vertical')
	plt.ylabel('Range (m)', fontsize = 14)
	plt.xlabel('Cross-range (m)', fontsize = 14)
	plt.show()

	#plt.savefig(RUTA %FECHA_HORA + 'Imagen_%d.png' %n_toma)
'''
	Imagen /= (nfre * (nx))
	Imagen = np.reshape(Imagen, (nposy, nposx))
	Imagen = np.flipud(Imagen)

	f = h5py.File(RUTA_RESULTADOS + NOMBRE_RESULTADOS %(n_toma), 'w')
	dset = f.create_dataset("Imagen_compleja", (nposy, nposx), dtype = np.complex64)
	dset.attrs['dx'] = dx
	dset.attrs['dy'] = dy
	dset.attrs['date'] = FECHA_HORA
	dset.attrs['ntoma'] = n_toma
	dset[...] = Imagen
	f.close()

Imagen = np.absolute(Imagen)

fig = plt.figure(1)
im = plt.imshow(Imagen, cmap = 'jet', aspect = 'auto', extent = [xi,xf,yi,yf], interpolation = 'none')
fig.suptitle('Holograma', fontsize = 16)
cbar = plt.colorbar(im, orientation = 'vertical')
cbar.set_label('dB', fontsize = 12)
plt.ylabel('Rango (m)', fontsize = 14)
plt.xlabel('Eje transversal al rango (m)', fontsize = 14)
toc = time.clock()
print "Se demoro %.2f segundos"%(toc - tic)
plt.show()
'''
