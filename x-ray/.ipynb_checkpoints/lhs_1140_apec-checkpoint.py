import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from  xspec import *

ntries = 10000

"""
Making an APEC model for LHS 1140 from the result in https://ui.adsabs.harvard.edu/abs/2023AJ....165..200S/abstract

"""
kt1_in = 0.15
kt1e_in = 0.03
kt2_in = 1.14
kt2e_in = 0.37
nh_in = 4.6e-2 #e20 for some reason
fx_in = 0.5 #e-14
fx_e_in = 0.075

x = dict(fx = fx_in, fx_e=fx_e_in, kt1=kt1_in, kt1_e=kt1e_in, kt2=kt2_in, kt2_e=kt2e_in, nh = nh_in, abd=1.0, Star='LHS_1140')



n = 0
fluxes = []
while n < ntries:
    print('Iteration:', n)
    fx = np.random.normal(x['fx'], x['fx_e'])
    if fx < 0.0:
        fx = x['fx']
    kt1 = np.random.normal(x['kt1'], x['kt1_e'])
    if kt1 < 0.008:
        kt1 = 0.008
    kt2 = np.random.normal(x['kt2'], x['kt2_e'])
    if kt2 < 0.008:
        kt2 = 0.008
    abd, nh = x['abd'], x['nh'] * 0.01
    # fx, kt1, kt2, abd, nh = x['fx'], x['kt1'], x['kt2'], x['abd'], x['nh']*0.01

    mod = Model('(apec+apec)*phabs', setPars={1:kt1, 2:abd,
                                                  5:kt2, 6:abd,
                                                  9:nh})
    Plot.xAxis = "angstrom"
    Plot.perHz = False
    Plot.area=True
    AllModels.setEnergies(".3 10. 1000")
    flux = AllModels.calcFlux(".3 10")
    fluxnum = mod.flux[0]
    norm = (fx*1e-14)/fluxnum
    mod.setPars({4:norm})
    mod.setPars({8:norm})

    AllModels.setEnergies("0.1 2.5 2400")
    Plot("model")
    xVals = Plot.x()
    yVals = Plot.model()
    wx = xVals*u.AA
    fx  = (yVals * (u.photon/u.s/u.cm**2/u.AA)).to(u.erg/u.s/u.cm**2/u.AA, equivalencies=u.spectral_density(wx))
  
    
    fluxes.append(fx[::-1])

    n += 1


wavelength = wx[::-1]
fluxes = np.array(fluxes)
fluxmean, fluxstd = np.mean(fluxes, axis=0), np.std(fluxes, axis=0)
savdat = Table((wavelength*u.AA, fluxmean*u.erg/u.s/u.cm**2/u.AA , fluxstd*u.erg/u.s/u.cm**2/u.AA), 
           names=['WAVELENGTH', 'FLUX', 'ERROR'])
ascii.write(savdat, '../models/{}_apec_errs.ecsv'.format(x['Star']), format='ecsv', overwrite=True)

fig, ax = plt.subplots()

ax.step(wavelength, fluxmean, where='mid')
ax.step(wavelength, fluxstd, where='mid', alpha=0.5)
ax.set_yscale('log')

print('DONE')


plt.show()
