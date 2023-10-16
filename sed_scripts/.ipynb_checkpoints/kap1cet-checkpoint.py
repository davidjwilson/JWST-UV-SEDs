"""
@verison: 1

@author: David Wilson

@date 20231016

Makes the Kap1Cet. Adding Starcat spectra implemented for the first time.
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
from scipy.interpolate import interp1d
import make_sed_files
from shutil import copyfile
import make_fits
import instruments

path = '/home/david/work/meats/SEDs/draft_hlsp/kap1cet/'  #path where the files are

star = 'kap1Cet'
version = 1
airglow =  [1214, 1217, 1300, 1310, 1353, 1356]

# trims = {'G430L':[3155, 5690], 'G230L':[1710, 3155], 'G140L':[1161, 1710]}

def make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var'):
    
    sed_table = []
    instrument_list = []
    # sed_table, instrument_list = sed.add_stis_and_lya(sed_table, path, airglow[0:2], instrument_list, airglow[2:], norm=False, remove_negs=remove_negs,to_1A=to_1A, trims=trims, lya_max = True, optical=True, Ebv=0.04)
    sed_table, instrument_list = sed.add_starcat(sed_table, path, instrument_list, remove_negs=remove_negs, to_1A=to_1A)
    
    
    # sed_table, instrument_list = sed.add_phoenix(sed_table, path, instrument_list, to_1A=to_1A, ranges = [5690, 1e10])
    
    
    
    # sed_table, instrument_list, gap = sed.add_xray_spectrum(sed_table, path, instrument_list, 'cxo', add_apec = True, find_gap=True, to_1A=to_1A)
    
    
    
    # sed_table, instrument_list = sed.add_euv(sed_table, path, instrument_list, gap, 'dem',to_1A=to_1A)
    
    
    
    sed_table.sort(['WAVELENGTH'])
    
    #deredden
#     sed_table = sed.deredden(sed_table, 0.04, [1161, 5690])
    
    sed_table = sed.add_bolometric_flux(sed_table, path)

    
    
    sed_table.meta['WAVEMIN'] = min(sed_table['WAVELENGTH'])
    sed_table.meta['WAVEMAX'] = max(sed_table['WAVELENGTH'])
    sed_table.meta['FLUXMIN'] = min(sed_table['FLUX'])
    sed_table.meta['FLUXMAX'] = max(sed_table['FLUX'])

    
    plt.figure()
    # if to_1A:
        # print(np.unique(np.diff(sed_table['WAVELENGTH'])))
    plt.step(sed_table['WAVELENGTH'], sed_table['FLUX'], where='mid')
    plt.step(sed_table['WAVELENGTH'], sed_table['ERROR'], where='mid', alpha=0.5)
    plt.yscale('log')
    plt.xscale('log')
    
    
    
#     print(sed_table.meta)
    
    make_fits.make_mm_fits(path, sed_table, instrument_list, version,sed_type=sed_type)
    
    
    # plt.show()
    

    
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var')
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=True, sed_type='const')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=False, sed_type='adapt-var')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=True, sed_type='adapt-const')

plt.show()