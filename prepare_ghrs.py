#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
from scipy.interpolate import interp1d

import astropy.units as u
import astropy.constants as const
from datetime import date
from astropy.units import cds
cds.enable()

__author__ = 'David Wilson'
__version__ = 0.1
__date__ = '20250310'

"""
Turns GRHS fits files into HLSP files


Takes the standard GHRS files from MAST, where wavelength, flux and error are in separate files, and combines them into on file. If multiple sub exposures are present then they are coadded into one spectrum. 

The calibrated files used are:
c0f = wavelength
c1f = flux
c2f = error
cqf = data quality flag 

As I don't need to coadd for Eps Ind I'm just putting in the error weighted coadd as a place holder for now, can add more sophisticated stuff in the future.

"""

def coadd_flux(f_array, e_array, scale_correct=True):
    """
    Returns a variance-weighted coadd with standard error of the weighted mean (variance weights, scale corrected).
    f_array and e_arrays are collections of flux and error arrays, which should have the same lenth and wavelength scale
    """
    weights = 1 / (e_array**2)
    flux = np.average(f_array, axis =0, weights = weights)
    var = 1 / np.sum(weights, axis=0)
    rcs = np.sum((((flux - f_array)**2) * weights), axis=0) / (len(f_array)-1) #reduced chi-squared
    if scale_correct:
        error = (var * rcs)**0.5
    else:
        error = var**0.5
    return flux,error

def make_plot(rootname, data):
    fig, ax = plt.subplots(num=rootname)
    ax.step(data['WAVELENGTH'], data['FLUX'], where='mid', label='FLUX')
    ax.step(data['WAVELENGTH'], data['ERROR'], where='mid', alpha=0.5, label='ERROR')
    ax.set_xlabel('Wavelength (\AA)')
    ax.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
    ax.set_title(rootname)
    ax.legend()
    fig.tight_layout()
    plt.show()

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
    

def make_ghrs_file(inpath, outpath='', plot=False, save_file=False, return_spectrum=True, find_exposure_times=True):
    """
    This is a basic funtion to turn a set of ghrs files into one fits file. Keeping it here for convienence.
    inpath = where the files from mast are stored
    outpath = where you would like the output files to go.
    """
    
    c0f_list = glob.glob('{}*c0f.fits'.format(inpath))
    roots = []
    for c0f in c0f_list:
        roots.append(fits.getheader(c0f, 0)['ROOTNAME'])
    print('inpath contains c0f (wavelength) files for the following datasets:', roots)
    
    for c0f in c0f_list:
        hdr = fits.getheader(c0f, 0)
        rootname = hdr['ROOTNAME']
        print('working with dataset {}'.format(rootname))
        wavelength_arrays = fits.getdata(c0f, 0)
        c1f = '{}/{}_c1f.fits'.format(inpath, rootname.lower())
        try: 
            flux_arrays = fits.getdata(c1f, 0) #test if this works
        except: 
                print('WARNING: no c1f (flux) file found for dataset {}, skipping'.format(rootname))
                continue
        c2f = '{}/{}_c2f.fits'.format(inpath, rootname.lower())
        try:
            error_arrays = fits.getdata(c2f, 0)
        except: 
            print('WARNING: no c2f (error) file found for dataset {}, skipping'.format(rootname))
            continue
        cqf = '{}/{}_cqf.fits'.format(inpath, rootname.lower())
        try:
            dq_arrays = fits.getdata(cqf, 0)
        except: 
            print('WARNING: no cqf (data quality) file found for dataset {}, skipping'.format(rootname))
            continue
    
        if len(np.shape(wavelength_arrays)) == 1: #check if there are multiple subexposures that need to be combined
            wavelength, flux, error, dq = wavelength_arrays, flux_arrays, error_arrays, dq_arrays
        else:
            print('this dataset has subexposures')
            wavelength = wavelength_arrays[0] #each wavelength array in a subexposure is the same as far as I can tell
            flux, error = coadd_flux(flux_arrays, error_arrays)
            dq = [(np.sum(np.unique(dq_arrays[:,i]))) for i in range(len(dq_arrays[0]))]
        
        filename = '{}_{}_{}_{}.fits'.format(hdr['TARGNAME'].lower(), hdr['INSTRUME'].lower(), hdr['GRATING'].lower(), hdr['ROOTNAME'].lower())
      
        data = Table((wavelength*u.AA, flux*u.erg/u.s/u.cm**2/u.AA, error*u.erg/u.s/u.cm**2/u.AA, dq), names = ['WAVELENGTH', 'FLUX', 'ERROR', 'DQ'])
        
        if plot:
            make_plot(rootname, data)

        if save_file:  #output a basic fits file 
            data_ext = fits.table_to_hdu(data)
            data_ext.header.set('EXTNAME', 'SPECTRUM')
            hdr_new = hdr
            hdr_new.set('FILENAME', filename)
            hdr_new.set('FILETYPE', 'SPECTRUM')
            hdr_new.set('FITSDATE', date.today().strftime('%Y-%m-%d'))
            primary_hdu = fits.PrimaryHDU(header=hdr_new)
            
            hdul = fits.HDUList([primary_hdu, data_ext])
            hdul.writeto('{}/{}'.format(outpath, filename), overwrite=True)
            print('spectrum saved as {}'.format(filename))
            
    if plot:
        plt.show()

    if return_spectrum:
        return data
        
def sort_ghrs(inpath):
    """
    Takes a path to where some ghrs files are stored, sorts them by grating and checks that they have wavelength, flux, error and dq arrays  
    """
    grating_collection = {}
    c0f_list = glob.glob('{}*c0f.fits'.format(inpath))
    # print(c0f_list)
    # roots = []
    for c0f in c0f_list:
        hdr = fits.getheader(c0f, 0)
        rootname = hdr['ROOTNAME']
        cfs = ['{}{}_c1f.fits'.format(inpath, rootname.lower()), '{}{}_c2f.fits'.format(inpath, rootname.lower()), '{}{}_cqf.fits'.format(inpath, rootname.lower())]
        if all([os.path.exists(cf) for cf in cfs]):
            # print('yes')
            grating = hdr['GRATING']
            if grating not in grating_collection.keys():
                grating_collection[grating] = [rootname]
            else:
                grating_collection[grating].append(rootname)

    return grating_collection


def combine_ghrs(inpath, rootnames):
    """
    Coadds a collection of ghrs spectra from the same grating. Input is a list of rootnames that should have already been sorted by grating, and checked to see if the avaialble
    wavelength, flux, error and dq files exist.

    """
    if len(rootnames) > 0:
        spectrum = dict(WAVELENGTH=[], WAVELENGTH0=[], WAVELENGTH1=[], FLUX=[], ERROR=[], DQ=[], EXPTIME=[], EXPSTART=[], EXPEND=[])       
        for rootname in rootnames:
            hdr = fits.getheader('{}/{}_c0f.fits'.format(inpath, rootname.lower()), 0)
            wavelength_arrays = fits.getdata('{}/{}_c0f.fits'.format(inpath, rootname.lower()), 0)
            # print(wavelength_arrays)
            flux_arrays = fits.getdata('{}/{}_c1f.fits'.format(inpath, rootname.lower()), 0)
            error_arrays = fits.getdata('{}/{}_c2f.fits'.format(inpath, rootname.lower()), 0)
            dq_arrays = fits.getdata('{}/{}_cqf.fits'.format(inpath, rootname.lower()), 0)
            # [spectrum['WAVELENGTH'].append(w) for w in wavelength_arrays]
            # [spectrum['FLUX'].append(f) for f in flux_arrays]
            # [spectrum['ERROR'].append(e) for e in error_arrays]
            # [spectrum['DQ'].append(dq) for dq in dq_arrays]
            # spectrum['EXPTIME'].append(np.full(len(wavelength_arrays[0]), hdr['EXPTIME']))
            # spectrum['EXPSTART'].append(np.full(len(wavelength_arrays[0]), hdr['EXPSTART']))
            # spectrum['EXPEND'].append(np.full(len([wavelength_arrays][0]), hdr['EXPEND']))
        if len(np.shape(wavelength_arrays)) == 1: #check if there are multiple subexposures that need to be combined
            spectrum['WAVELENGTH'] =  wavelength_arrays
            spectrum['FLUX'] = flux_arrays
            spectrum['ERROR'] = error_arrays
            spectrum['DQ'] = dq_arrays
            spectrum['EXPTIME']= np.full(len(spectrum['WAVELENGTH']), hdr['EXPTIME'])
            spectrum['EXPSTART'] = np.full(len(spectrum['WAVELENGTH']), hdr['EXPSTART'])
            spectrum['EXPEND']= np.full(len(spectrum['WAVELENGTH']), hdr['EXPEND'])
            w0, w1 = wavelength_edges(spectrum['WAVELENGTH'])
            spectrum['WAVELENGTH0'], spectrum['WAVELENGTH1'] = w0, w1
        else:
            print('You have not written this function yet')

        #add units 
        spectrum['WAVELENGTH'] = spectrum['WAVELENGTH']*u.AA
        spectrum['WAVELENGTH0'] = spectrum['WAVELENGTH0']*u.AA
        spectrum['WAVELENGTH1'] = spectrum['WAVELENGTH1']*u.AA
        spectrum['FLUX'] = spectrum['FLUX']*u.erg/u.s/u.cm**2/u.AA
        spectrum['ERROR'] = spectrum['ERROR']*u.erg/u.s/u.cm**2/u.AA
        spectrum['EXPTIME'] = spectrum['EXPTIME']*u.s
        spectrum['EXPSTART'] = spectrum['EXPSTART']*cds.MJD
        spectrum['EXPEND'] = spectrum['EXPEND']*cds.MJD

        return spectrum

def make_metadata(inpath, rootnames, spectrum, hlsp, normfac, star):
    """
    Makes the metadata for the ecsv files- eventually will be converted into the fits file
    """
    wavelength, flux = spectrum['WAVELENGTH'].value, spectrum['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    exptimes = []
    start_times = []
    end_times = []
    dates = []
    for rootname in rootnames:      
        hdr = fits.getheader('{}/{}_c0f.fits'.format(inpath, rootname.lower()), 0)
        exptimes.append(hdr['EXPTIME'])
        start_times.append(hdr['EXPSTART'])
        end_times.append(hdr['EXPEND'])
        dates.append(hdr['DATE-OBS'])
        if star == '':
            starname = hdr['TARGNAME']
        else:
            starname = star
    # meats_name = 'MUSCLES Extension for Atmospheric Transmisson Spectroscopy'
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','APERTURE','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN','EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['HST','','','',starname,'','','',hlsp['PROGRAM'],hlsp['PROGRAMSHORT'],hlsp['HLSPLEAD'],
                 'n/a','n/a',min(dates),min(start_times),max(end_times),sum(exptimes),'SUM', 
                min(exptimes), max(exptimes), np.median(exptimes),normfac,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = hdr[name]
        else:
            metadata[name] = filler
    return metadata
       
def make_dataset_extension(rootnames):
    """
    Makes a fits extension containg a list of rootnames and dataset IDs used in the spectrum.
    """
    description_text = 'This extension contains a list of HST rootnames (9 character string in HST files downloaded from MAST) and dataset IDs of the exposures used to create this spectrum file. The dataset IDs can be used to directly locate the observations through the MAST HST data archive search interface. Multiple identifiers indicate the spectra were coadded. In some cases the spectra have been rextracted and will be different from those availbale from MAST.' 
    
    dataset_table = Table([rootnames,rootnames], names=[('ROOTNAME'),('DATASET_ID')])
    hdu = fits.table_to_hdu(dataset_table)
    hdu.header.insert(8, ('EXTNAME','SRCSPECS'))
    hdu.header.insert(9, ('EXTNO',3))
    hdu.header['COMMENT'] = description_text
    return hdu
              
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

def make_component_filename(metadata, version, hlsp):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    filename = 'hlsp_{}_{}_{}_{}_{}_v{}_component-spec'.format(hlsp['HLSPNAME'].lower(),metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

def make_ghrs_spectrum(inpath, version, hlsp, savepath = '', plot=False, save_ecsv=False, save_fits=False, return_data=False, return_gratings = False, normfac=1.0, star = ''):
    """
    main function
    """
    grating_collection = sort_ghrs(inpath)
    print(grating_collection)
    hlsp = Table.read(hlsp)[0]
    if len(grating_collection) > 0:
        for grating in grating_collection.keys():
            spectrum = combine_ghrs(inpath, grating_collection[grating])
            metadata = make_metadata(inpath, grating_collection[grating], spectrum, hlsp, normfac, star)
            if plot:
                make_plot(grating, spectrum)
            if save_ecsv:
                save_to_ecsv(spectrum, metadata, savepath, version)
            if save_fits:
                data_set_hdu = make_dataset_extension(grating_collection[grating])
                save_to_fits(spectrum, metadata, hlsp, data_set_hdu, savepath, version)
    if return_data:
        return spectrum
    if return_gratings:
        return gratings_collection.keys()




