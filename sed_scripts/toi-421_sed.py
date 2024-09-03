"""
@verison: 1

@author: David Wilson

@date: 20240903

Makes the TOI-421 SED
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



path = '/home/david/work/meats/SEDs/draft_hlsp/toi-421/'  #path where the files are

star = 'TOI-421'
version = 1
airglow =  [1207, 1222, 1300, 1310, 1353, 1356]

trims = {'G430L':[3161, 5690], 'G230L':[1711, 3160 ], 'G140L':[1160, 1710]}

def make_sed(path, star, version, norm=False, remove_negs=False, to_1A=False, sed_type='var'):
    
    sed_table = []
    instrument_list = []
    
    # sed_table, instrument_list = sed.add_cos(path, airglow,  remove_negs=remove_negs,to_1A=to_1A, trims={'G130M':[1080, 1363]})
    
    
    sed_table, instrument_list = sed.add_stis_and_lya(sed_table, path, [1212, 1220], instrument_list, airglow[2:], norm=False, remove_negs=remove_negs,to_1A=to_1A, trims=trims, optical=True, lya_max=False,  Ebv=0.06)
    
    proxy_path = '/home/david/work/meats/SEDs/draft_hlsp/tau_ceti/hlsp_muscles_multi_multi_tau_ceti_broadband_v1_{}-res-sed.fits'.format(sed_type)

    d_prox = 3.65
    d_star = 74.96 # updated distance 
    scale = (d_prox/d_star)**2
    sed_table, instrument_list = sed.add_proxy(sed_table, proxy_path, instrument_list, scale, ranges = [1, 1159], remove_negs=False,to_1A=False)
    
    sed_table, instrument_list = sed.add_phoenix(sed_table, path, instrument_list, to_1A=to_1A)
    
  
    
    sed_table.sort(['WAVELENGTH'])
    
    
    
    sed_table = sed.add_bolometric_flux(sed_table, path)

    
    
    sed_table.meta['WAVEMIN'] = min(sed_table['WAVELENGTH'])
    sed_table.meta['WAVEMAX'] = max(sed_table['WAVELENGTH'])
    sed_table.meta['FLUXMIN'] = min(sed_table['FLUX'])
    sed_table.meta['FLUXMAX'] = max(sed_table['FLUX'])

    # if remove_negs: #fixing some weird errors that show up in the adapt_var
    #     mask = (sed_table['WAVELENGTH'] >= 1175) & (sed_table['WAVELENGTH'] <= 1195) | (sed_table['WAVELENGTH'] >= 1210) & (sed_table['WAVELENGTH'] <= 1226)
    #     args = np.where(mask==True)
    #     new_err = np.median(sed_table['ERROR'][(sed_table['WAVELENGTH'] >= 1200) & (sed_table['WAVELENGTH'] <= 1209)])
    #     sed_table['ERROR'][args] = new_err
        
    
    
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