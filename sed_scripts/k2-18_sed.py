"""
@verison: 1

@author: David Wilson

@date: 20240212

Makes the K2-18 SED
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
# from craftroom import resample
from scipy.interpolate import interp1d
import make_sed_files
from shutil import copyfile
import make_fits
import instruments



path = '/home/david/work/meats/SEDs/draft_hlsp/k2-18/'  #path where the files are

star = 'K2-18'
version = 1
airglow =  [1207, 1222, 1300, 1310, 1353, 1356]

trims = {'G430L':[3151, 5690], 'G230L':[1700, 3150 ]}

def make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var'):
    
    sed_table = []
    instrument_list = []
    
    
    
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, path, airglow[0:2], instrument_list, airglow[2:], norm=False, remove_negs=remove_negs,to_1A=to_1A, trims=trims, optical=True)
    
    
    proxy_path = '/media/david/2tb_ext_hd/hddata/mega_muscles/muscles_hlsp/gj832/hlsp_muscles_multi_multi_gj832_broadband_v22_var-res-sed.fits'

    d_prox = 1000/201.3252
    d_star = 38.025
    scale = (d_prox/d_star)**2
    sed_table, instrument_list = sed.add_proxy(sed_table, proxy_path, instrument_list, scale, ranges = [1, 1214, 1217, 1699], remove_negs=remove_negs,to_1A=to_1A)
    
    
    sed_table, instrument_list = sed.add_phoenix(sed_table, path, instrument_list, to_1A=to_1A)
    
    # sed_table, instrument_list, gap = sed.add_xray_spectrum(sed_table, path, instrument_list, 'xmm', add_apec = True, find_gap=True, to_1A=to_1A)
    
    #cutting the apec model down, fix later (fixed, wrong dem)
#     apec_mask = (sed_table['WAVELENGTH'] < 50) | (sed_table['WAVELENGTH'] > 200)
#     gap[0] = 40
#     sed_table = sed_table[apec_mask]
    
    
    
    # sed_table, instrument_list = sed.add_euv(sed_table, path, instrument_list, gap, 'dem',to_1A=to_1A)
    
    
    
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
    
#     plt.xlim(1160, 1230)
    plt.xscale('log')
    
    
    
#     print(sed_table.meta)
    
    make_fits.make_mm_fits(path, sed_table, instrument_list, version,sed_type=sed_type)
    
    
    # plt.show()
    

    
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var')
make_sed(path, star, version, norm=False, remove_negs=False, to_1A=True, sed_type='const')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=False, sed_type='adapt-var')
make_sed(path, star, version, norm=False, remove_negs=True, to_1A=True, sed_type='adapt-const')

plt.show()