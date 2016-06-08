from __future__ import division
import os
import numpy as np
import h5py
from anritsuvna import *
from slp_eth import *

RUTA_GUARDADO = '/home/cielosar/Desktop/Datos_hdf5'
NOMBRE_FOLDER = '/experimentos_%s_%s'
NOMBRE_ARCHIVO = '/datos_toma_%d.hdf5'

LONG_VUELTA = 66
PASOS_VUELTA = 20000
FACTOR_LONG_PASOS = 1000 * PASOS_VUELTA / LONG_VUELTA

xi = 0.0
xf = 1.4
npos = 2
fi = 15.5
ff = 16.5
nfre = 4001
ptx = 'HIGH'
ifbw = 1000
ntomas = 1

xf = int(xf * FACTOR_LONG_PASOS)
xi = int(xi * FACTOR_LONG_PASOS)
#dx = int((xf - xi) / (npos + 1))
dx = int((xf - xi) / (npos - 1))
#xf = xi + dx * (npos - 1)
xf = xi + dx * (npos - 1)

vna = VnaClient()
vna.connect()
vna.send_ifbw(ifbw)
vna.send_number_points(nfre)
vna.send_freq_start(freq_start = fi)
vna.send_freq_stop(freq_stop = ff)
vna.send_power(ptx)
vna.send_select_instrument()
vna.send_cfg()

'''
slp = SLPClient()
slp.connect()
slp.send_to_start()
slp.close()
'''

FECHA   = (time.strftime("%d.%m.%Y"))
HORA    = (time.strftime("%H.%M.%S"))

if not os.path.exists(RUTA_GUARDADO + NOMBRE_FOLDER %(FECHA, HORA)):
    os.makedirs(RUTA_GUARDADO + NOMBRE_FOLDER %(FECHA, HORA))

for i in range(0,ntomas):
    f = h5py.File(RUTA_GUARDADO + NOMBRE_FOLDER %(FECHA, HORA) +
                       NOMBRE_ARCHIVO %(i+1), 'w')
    dset = f.create_dataset("sar_dataset", (npos, nfre), dtype = np.complex64)
    dset.attrs['xi'] = xi / FACTOR_LONG_PASOS
    dset.attrs['xf'] = xf / FACTOR_LONG_PASOS
    dset.attrs['dx'] = dx / FACTOR_LONG_PASOS
    dset.attrs['npos'] = npos
    dset.attrs['fi'] = fi
    dset.attrs['ff'] = ff
    dset.attrs['nfre'] = nfre
    dset.attrs['ptx'] = ptx
    dset.attrs['ifbw'] = ifbw

    for j in range(0, npos):
        if j == 0:
            dset[j,:] = vna.send_sweep()
        #slp.connect()
        #slp.send_move(dx, 'L')
        #slp.close()
        dset[j,:] = vna.send_sweep()
    f.close()
    #slp.connect()
    #slp.send_to_start()
    #slp.close()

vna.close()
