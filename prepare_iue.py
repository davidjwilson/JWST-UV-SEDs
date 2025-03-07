import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u

from scipy.interpolate import interpolate
from astropy.units import cds
import astropy.stats as stats
from scipy import stats
from specutils import Spectrum1D
from specutils.manipulation import FluxConservingResampler
from astropy.nddata import StdDevUncertainty


from astropy import time, coordinates as coord
import astropy.constants as const
from datetime import datetime


cds.enable()

"""
@author: David Wilson

version 1 20250306

Finds all IUE files, groups and coadds with senistivity weighted coadd.

Initial version will just do low res.

"""

def wavelength_edges(w):
    """
    Calulates w0 and w1
    """
    diff = np.diff(w)
    diff0 = np.concatenate((np.array([diff[0]]), diff)) #adds an extravalue to make len(diff) = len(w)
    diff1 = np.concatenate((diff, np.array([diff[-1]]))) #adds an extravalue to make len(diff) = len(w)
    w0 = w - diff0/2.
    w1 = w + diff1/2.
    return w0, w1
    
def coadd_iue(spectra, dq_cut=0, rebin=2):
    """
    Coadds IUE spectra using throughput weighting and dq masking. Bins by factor two by default. Makes the other arrays for a MUSCLES HLSP
    """
    if len(spectra) > 0:
        print('there are no spectra here')

            

        fluxes = []
        dqs = []
        waves = []
        counts = []
        nets = []
        exptimes = []
        expstarts = []
        expends = []
        for spec in spectra:
            hdul = fits.open(spec)
            hdr = hdul[0].header
            data = hdul[1].data
            hdul.close()
            fluxes.append(data[0]['FLUX'])
            dqs.append(-1*data[0]['QUALITY'])
            counts.append(data[0]['NET']+data[0]['BACKGROUND']) 
            nets.append(data[0]['NET'])
            w = np.array([data[0]['WAVELENGTH']+ data[0]['DELTAW']*i for i in range(data[0]['NPOINTS'])])
            waves.append(w)
            expstart = hdr['LJD-OBS']- 2400000.5 #convert JD to MJD
            exptime = hdr['LEXPTIME']
            expend = expstart + ((exptime*u.s).to(u.d)).value
            expstarts.append(np.full(len(w)), expstart)
            exptimes.append(np.full(len(w)), exptime)
            exends.append(np.full(len(w)), expend)

            
    
        #mask and coadd data, weighted by throughput
        mask = np.array(dqs) > dq_cut
        if len(spectra) == 1:
            wavelength = w[0][mask]
            coadd_flux = fluxes[0][mask]
            coadd_error = data[0]['SIGMA'][mask]
            exptime = exptimes[0][mask]
            start = expstarts[0][mask]
            end = expends[0][mask]

        else:
        
            fluxes = np.ma.array(fluxes, mask = mask)
            nets = np.ma.array(nets, mask = mask)
            coadd_flux = np.ma.average(fluxes,weights=nets/fluxes, axis=0).data
        
            #sum all counts and calculate S/N
            sumcounts = np.ma.sum(np.ma.array(counts, mask=mask), axis=0)#.data 
            sumnets = np.ma.sum(np.ma.array(nets, mask=mask), axis=0)#.data 
            coadd_error = (sumcounts**0.5 / sumnets) * coadd_flux #/ np.sum(exptimes)
            n = 0
            while np.min(coadd_error) < 0.0 and n < 100: #removing negative error values and replacing them with adjacent points. N count stops infinate loop
                new_error = np.copy(coadd_error)
                for i in range(len(new_error)):
                    if new_error[i] < 0.0:
                        new_error[i] = np.mean((new_error[i-1], new_error[i-1]))
                coadd_error = new_error
                n +=1
                                
        
            wavelength = stats.mode(waves).mode # wavelength arrays are near identical so just use the mode in each case. Experiment with uniform grid?
            exptime = np.ma.sum(np.ma.array(exptimes, mask = mask), axis=0)
            start = np.ma.min(np.ma.arrya(exptstarts, mask = mask), axis=0)
            end = np.ma.max(np.ma.arrya(exptends, mask = mask), axis=0)
        
        
            # rebinning 
        if rebin >0:
            w0, w1 = wavelength[0], wavelength[-1]
            dw = np.median(np.diff(wavelength)) *rebin_pix
            new_wavelength = np.arange(w0, w1, dw)
            fluxcon = FluxConservingResampler(extrapolation_treatment='zero_fill')
            input_spec = Spectrum1D(spectral_axis=wavelength*u.AA, 
                                    flux=coadd_flux*u.erg/u.s/u.cm**2/u.AA ,
                                    uncertainty= StdDevUncertainty(coadd_error))
            new_spec_fluxcon = fluxcon(input_spec, new_wavelength*u.AA)
            new_wavelength = (new_spec_fluxcon.spectral_axis.value)
            coadd_flux = (new_spec_fluxcon.flux.value)
            coadd_error = (1/(new_spec_fluxcon.uncertainty.array**0.5))

            start = interpolate.interpld(wavelength, start, kind='nearest',bounds_error=False, fill_value=0.)(new_wavelength)
            end = interpolate.interpld(wavelength, end, kind='nearest',bounds_error=False, fill_value=0.)(new_wavelength)
            exptime = interpolate.interpld(wavelength, exptime, kind='nearest',bounds_error=False, fill_value=0.)(new_wavelength)

            wavelength = new_wavelength
        
        dq_new = np.zeros(len(wavelength)) # need a dq array for compatability with MUSCLES
        w0, w1 = wavelength_edges(wavelength)
           
        new_data = {'WAVELENGTH':wavelength*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':coadd_flux*u.erg/u.s/u.cm**2/u.AA,
                'ERROR':coadd_error*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq_new,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD}
        return new_data



def sort_mxlos(mxlos):
    """
    Sorts mxlows into short and long- wavelength groups
    """
    swls = []
    lwls = []
    for spec in mxlos:
        camera = fits.getheader(mxlo, 0)['CAMERA']
        if camera == 'SWP':
            swls.append(spec)
        elif camera in ['LWP', 'LWR]':
            lwls.append(spec)
    return dict(SWLO=swls, LWLO=lwls)

def make_metadata(mxlos, new_data, hlsp, normfac, star):
    """
    Makes the metadata for the fits file
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    dates = []
    for x in mxlos:      
        hdr = fits.getheader(x,0)
        dates.append(hdr['LDATEOBS'])
        if star == '':
            starname = hdr['LTARGET']
        else:
            starname = star
    # meats_name = 'MUSCLES Extension for Atmospheric Transmisson Spectroscopy'
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','APERTURE','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN','EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['',hdr['CAMERA'],hdr['DISPERSN'],'',starname,hdr['LTARGRA'],hdr['LTARGRDEC'],hdr['PGM-ID'],hlsp['PROGRAM'],hlsp['PROGRAMSHORT'],hlsp['HLSPLEAD'],
                 'n/a','n/a',min(dates),min(new_data['EXPSTART'].value),max(new_data['EXPEND'].value),max(new_data['EXPTIME'].value),'SUM', 
                min(exptimes), max(exptimes), np.median(exptimes),normfac,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = hdr[name]
        else:
            metadata[name] = filler
    return metadata
    
def save_to_ecsv(data, metadata, save_path, version):
    """
    save the new model to an ecsv file
    """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    file_name = make_component_filename(metadata, version)
    savedat = Table(data, meta=metadata)
    savedat.write(save_path+file_name+'.ecsv', overwrite=True, format='ascii.ecsv')
    print('Spectrum saved as '+file_name+'.ecsv')

    
def save_to_fits(data, metadata, hlsp, dataset_hdu, savepath, version):
    """
    Saves to a MUSCLES-standard fits file
    """
    if os.path.exists(savepath) == False:
        os.mkdir(savepath)
    file_name = make_component_filename(metadata, version, hlsp)
    hdr = fits.Header(metadata)
    primary_hdu = fits.PrimaryHDU(header=hdr)
    hdu = fits.table_to_hdu(Table(data))
    descriptions =['midpoint of the wavelength bin', 'left/blue edge of the wavelength bin','right/red edge of the wavelength bin','average flux over the bin',
                'error on the flux','cumulative exposure time for the bin','data quality flags (HST data only)','modified julian date of start of first exposure', 
                'modified julian date of end of last exposure']
    hdu.header.insert(8, ('EXTNAME','SPECTRUM'))
    hdu.header.insert(9, ('EXTNO',2))
    [hdu.header.insert(i[0]+10, ('TDESC{}s'.format(i[0]), i[1])) for i in enumerate(descriptions)]
    hdul = fits.HDUList([primary_hdu, hdu, dataset_hdu])
    hdul.writeto(savepath+file_name+'.fits', overwrite=True)
    print('Spectrum saved as '+file_name+'.fits')
    

def plot_spectrum(data, metadata):
    fig, ax = plt.subplots()
    ax.step(data['WAVELENGTH'], data['FLUX'], where='mid', label='FLUX')
    ax.step(data['WAVELENGTH'], data['ERROR'], where='mid', alpha=0.5, label='ERROR')
    ax.set_xlabel('Wavelength (\AA)', size=20)
    ax.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)', size=20)
    ax.set_title('{}_{}'.format(metadata['TARGNAME'],metadata['GRATING']))
    ax.legend()
    fig.tight_layout()
    plt.show()


    
def make_iue_spectrum(mxlopath, version, hlsp, savepath = '', plot=False, save_ecsv=False, save_fits=False, return_data=False, return_gratings = False, normfac=1.0, star = '', nclip=5):
    """
    main function
    """
    mxlos = glob.glob('{}/*mxlo.gz'.format(mxlopath))
    hlsp = Table.read(hlsp)[0]
    if len(mxlos) > 0:
        grouped_mxlos = sort_mxlos(mxlos)
        for band in grouped_mixlos.keys():
            if len(grouped_mxlos[band]) > 0:
                print('Working with {} spectra'.format(band))
                data = coadd_iue(grouped_mixlos[band])
                metadata = make_metadata(grouped_mixlos[band], data, hlsp, normfac, star)
            if plot:
                plot_spectrum(data, metadata)
            if save_ecsv:
                save_to_ecsv(data, metadata, savepath, version)
            if save_fits:
                data_set_hdu = make_dataset_extension(x1ds)
                save_to_fits(data, metadata, hlsp, data_set_hdu, savepath, version)
    if return_data:
        return data
    if return_gratings:
        return grouped_mixlos.keys()    
                
        
        
        











