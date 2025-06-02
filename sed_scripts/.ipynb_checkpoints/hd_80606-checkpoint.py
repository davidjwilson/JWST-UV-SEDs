"""
@verison: 1

@author: David Wilson

@date 20250602

Makes the HD 80606 sed.
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

path = '/home/david/work/meats/SEDs/draft_hlsp/hd80606/'  #path where the files are

star = 'HD80606'
version = 1
airglow =  [1214, 1217, 1300, 1310, 1353, 1356]

trims = {'G230L':[1715, 3160], 'G140L':[1160, 1715], 'G430L':[3160, 5690]}

def make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var'):
    
    sed_table = []
    instrument_list = []
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, path, airglow[0:2], instrument_list, airglow[2:], norm=True, remove_negs=remove_negs,to_1A=to_1A, trims=trims, lya_max = False, optical=True, Ebv=0.045)
    
    
    
    sed_table, instrument_list = sed.add_phoenix(sed_table, path, instrument_list, to_1A=to_1A)#, ranges = [5690, 1e10])
    
    
    
    sed_table, instrument_list, gap = sed.add_xray_spectrum(sed_table, path, instrument_list, 'xmm', add_apec = True, find_gap=True, to_1A=to_1A, trims=trims, remove_negs=remove_negs, just_apec=True)
    
    
    
    sed_table, instrument_list = sed.add_euv(sed_table, path, instrument_list, gap, 'dem',to_1A=to_1A)
    
    
    
    sed_table.sort(['WAVELENGTH'])
    

    
    sed_table = sed.add_bolometric_flux(sed_table, path)

    
    
    sed_table.meta['WAVEMIN'] = min(sed_table['WAVELENGTH'])
    sed_table.meta['WAVEMAX'] = max(sed_table['WAVELENGTH'])
    sed_table.meta['FLUXMIN'] = min(sed_table['FLUX'])
    sed_table.meta['FLUXMAX'] = max(sed_table['FLUX'])

    
    plt.figure(sed_type)
    # if to_1A:
        # print(np.unique(np.diff(sed_table['WAVELENGTH'])))
    plt.step(sed_table['WAVELENGTH'], sed_table['FLUX'], where='mid')
    plt.step(sed_table['WAVELENGTH'], sed_table['ERROR'], where='mid', alpha=0.5)
    plt.yscale('log')
    plt.xscale('log')
    
    
    
#     print(sed_table.meta)
    
    make_fits.make_mm_fits(path, sed_table, instrument_list, version,sed_type=sed_type)

    sed.update_meta(star.lower(), version, newpath='../fixed_hlsp/', oldpath = '../draft_hlsp/')
    
    
    # plt.show()
    

    
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var')
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=True, sed_type='const')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=False, sed_type='adapt-var')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=True, sed_type='adapt-const')

plt.show()
