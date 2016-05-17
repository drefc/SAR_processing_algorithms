from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
from scipy import interpolate
#from mpl_toolkits.mplot3d import Axes3D
import time
import h5py

FECHA_HORA = '16.05.2016_13.55.28'

RUTA = '/home/cielosar/Desktop/Datos_hdf5/experimentos_%s/'
NOMBRE_ARCHIVO = 'datos_toma_%d.hdf5'
RUTA_RESULTADOS = "/home/cielosar/Desktop/Resultados_paper/Resultados2"
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
    print RUTA % FECHA_HORA + NOMBRE_ARCHIVO %(n_toma)
    f = h5py.File(RUTA % FECHA_HORA + NOMBRE_ARCHIVO %(n_toma), 'r')
    dset = f['sar_dataset']
    fre_min = dset.attrs['fi'] * 1E9
    fre_max = dset.attrs['ff'] * 1E9
    nx = dset.attrs['npos']
    c0 = 299792458.0
    B = fre_max - fre_min
    dr = c0 / (2 * B)
    nfre = dset.attrs['nfre']
    df = (fre_max - fre_min) / (nfre - 1)

    distance = np.arange(start= 0.0, stop = nfre * dr, step = dr)

    print "Toma numero %d" %n_toma
    print fre_min
    print fre_max
    print nfre

    s21 = np.empty([nx, nfre], dtype = np.complex64)
    s21 = dset[100,:]
    s21 = ifft(s21)
    f.close()
    s21 = np.absolute(s21)
    s21 *= (1.0 / s21.max())
    s21 = 20 * np.log10(s21)

    print 'Plotting...'
    print 'Max value %f at %f' %(s21.max(), distance[s21.argmax()])

    fig = plt.figure(1)
#	im = plt.imshow(Imagen, cmap = 'jet', aspect = 'auto', extent = [xi,xf,yi,yf], interpolation = 'none')
#	im = plt.imshow(ImagendB, cmap = 'jet', aspect = 'auto', extent = [xi,xf,yi,yf], interpolation = 'none')
    im = plt.plot(distance, s21)
    plt.xlabel('Range (m)', fontsize = 14)
    plt.ylabel('Relative Radar Reflectivity (dB)', fontsize = 14)

    plt.show()
