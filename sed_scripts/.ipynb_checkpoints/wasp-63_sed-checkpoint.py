"""
@verison: 1

@author: David Wilson

@date 20230126

Makes the WASP-63 SED
"""

import instruments
import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from astropy.table import Table, vstack
from astropy.io import ascii
import astropy.units as u
import make_mm_sed as sed
# import prepare_cos
# import prepare_stis
# import prepare_model
# import prepare_xmm
# import prepare_euv
# import prepare_chandra
from craftroom import resample
from scipy.interpolate import interp1d
import make_sed_files
from shutil import copyfile
import make_fits
import instruments

path = '/home/david/work/meats/SEDs/draft_hlsp/wasp-63/'  #path where the files are

star = 'WASP-63'
version = 2

airglow =  [1207, 1222, 1300, 1310, 1353, 1356]
stis_gratings = ['G140M','E140M','G140L', 'G230L', 'G230LB', 'G430L']

def make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var'):
    
    # sed_table = Table.read('{}/hlsp_muscles_hst_stis_wasp-63_g230l_v1_component-spec.ecsv'.format(path)) #it needs a sed or everything breaks, fix later
    sed_table = []
    instrument_list = []
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, path, airglow[0:2], instrument_list, airglow[2:], norm=False, remove_negs=True,to_1A=to_1A, trims={'G230L':[1980, 3140]}) #doing remove negatives all the time as there's only one negative flux point
    # print('adding phx')
    
    
    sed_table, instrument_list = sed.add_phoenix_and_g430l(sed_table, path, instrument_list, remove_negs=remove_negs, to_1A=to_1A, trims={'G430L':[3050, 5650]})
    
    
    sed_table, instrument_list = sed.add_euv(sed_table, path, instrument_list, [0, min(sed_table['WAVELENGTH'])], 'sol',to_1A=to_1A)
    
    sed_table.sort(['WAVELENGTH'])
    
    sed_table = sed.add_bolometric_flux(sed_table, path)

    sed_table.meta['WAVEMIN'] = min(sed_table['WAVELENGTH'])
    sed_table.meta['WAVEMAX'] = max(sed_table['WAVELENGTH'])
    sed_table.meta['FLUXMIN'] = min(sed_table['FLUX'])
    sed_table.meta['FLUXMAX'] = max(sed_table['FLUX'])

    
    plt.figure('{}_adapt={}_const={}'.format(star, remove_negs, to_1A))
    plt.step(sed_table['WAVELENGTH'], sed_table['FLUX'], where='mid')
    plt.step(sed_table['WAVELENGTH'], sed_table['ERROR'], where='mid', alpha=0.5)
    plt.yscale('log')
    plt.xscale('log')
    
    make_fits.make_mm_fits(path, sed_table, instrument_list, version,sed_type=sed_type)
    

    

# plt.figure()    
    
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var')
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=True, sed_type='const')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=False, sed_type='adapt-var')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=True, sed_type='adapt-const')

plt.show()
