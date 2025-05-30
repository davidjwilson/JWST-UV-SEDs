import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.units import cds
cds.enable()

"""

@verison: 4

@author: David Wilson

@date 20250422

Turns AB's Chandra files into HSLP escv and fits files. MEATS version

V4 added bit to take a grating spectrum in a text file.
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
    
def nan_clean(array):
    """
    replace nans in arrays with zeros
    """
    for i in range(len(array)):
        if np.isnan(array[i]) == True:
            array[i] = 0.0
    return array

def apec_to_ecsv(model_path, sed_meta, save_path):
    """
    save the apec model to an ecsv file
    """
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
       
    target = sed_meta['TARGNAME']
    wavelength, bin_width, flux = np.loadtxt(model_path, skiprows=3, unpack=True)
    flux = [fi*1.99e-8/wi for fi, wi in zip(flux,wavelength)]
    savedat = Table([wavelength, flux], names=['WAVELENGTH', 'FLUX'])
    name = target+'apec.txt' 
    ascii.write(savedat, save_path+name, overwrite=True)
        
def build_chandra_data(spectrum_path, hdr0):
    w, bins, counts, counts_err, model  = np.loadtxt(spectrum_path, unpack=True, skiprows=4)
    f = [fi*1.99e-8/wi for fi, wi in zip(counts,w)]
    e = (counts_err/counts)*f
    args = np.argsort(w)
    w, e, bins = w[args], e[args], bins[args] #files are often backwards
    f = np.array(f)[args]
    w0, w1 = w - bins, w+bins    
    start = np.full(len(w), (Time(hdr0['DATE-OBS']).mjd))
    end = np.full(len(w), (Time(hdr0['DATE-END']).mjd))
    exptime = np.full(len(w), ((Time(hdr0['DATE-END']).mjd - (Time(hdr0['DATE-OBS']).mjd))*u.d.to(u.s)))
    dq = np.zeros(len(w), dtype=int)
    f, e = nan_clean(f), nan_clean(e)
    new_data = {'WAVELENGTH':w*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f*u.erg/u.s/u.cm**2/u.AA,
                'ERROR':e*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD}
    return new_data

def build_from_txt(spectrum_path, hdr0):
    """
    Taking a spectrum from an ecsv that has w, f, e in it, this is for the 70 Oph B grating spectrum at first.
    """
    data = Table.read(spectrum_path, format='ascii.basic')
    w, f, e = data['WAVELENGTH'], data['FLUX'], data['ERROR']
    # w, bins, counts, counts_err, model  = np.loadtxt(spectrum_path, unpack=True, skiprows=4)
    # f = [fi*1.99e-8/wi for fi, wi in zip(counts,w)]
    # e = (counts_err/counts)*f
    # args = np.argsort(w)
    # w, e, bins = w[args], e[args], bins[args] #files are often backwards
    # f = np.array(f)[args]
    # w0, w1 = w - bins, w+bins
    w0, w1 = wavelength_edges(w)
    start = np.full(len(w), (Time(hdr0['DATE-OBS']).mjd))
    end = np.full(len(w), (Time(hdr0['DATE-END']).mjd))
    exptime = np.full(len(w), ((Time(hdr0['DATE-END']).mjd - (Time(hdr0['DATE-OBS']).mjd))*u.d.to(u.s)))
    dq = np.zeros(len(w), dtype=int)
    f, e = nan_clean(f), nan_clean(e)
    new_data = {'WAVELENGTH':w*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f*u.erg/u.s/u.cm**2/u.AA,
                'ERROR':e*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD}
    return new_data
                  
def build_chandra_metadata(hdr1, new_data, hlsp, star):
    """
    Makes the metadata for the chandra data table. Version 2 rewriting to not requre the SED metadata so can make x-ray files separately 
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    start, end, exptime = np.min(new_data['EXPSTART'][new_data['EXPSTART']>0].value), np.max(new_data['EXPEND'].value), np.max(new_data['EXPTIME'].value)
    if star == '':
        star = hdr1['OBJECT']
    meta_names =['TELESCOP','INSTRUME','GRATING','DETECTOR','FILTER',
                 'TARGNAME','RA_TARG','DEC_TARG','PROPOSID',
                 'HLSPNAME','HLSPACRN','HLSPLEAD',
                 'PR_INV_L','PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN',
                 'EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['CXO','1','NONE',hdr1['DETNAM'],'NA',
                 star,'1','1',hlsp['PROPOSID'][0],
                 hlsp['PROGRAM'][0],hlsp['PROGRAMSHORT'][0],hlsp['HLSPLEAD'][0],
                 hlsp['PR_INV_L'][0],hlsp['PR_INV_F'][0],'1', start, end, exptime, 'MEAN', exptime, 
                 exptime, exptime, 1.0, min(wavelength), max(wavelength), 'ang', 'vac', specres, waveres,np.min(flux[np.isnan(flux)==False]),
                 np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']  
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == 'sed':
            metadata[name] = sed_meta[name]
     #   elif filler == '0':
      #      metadata[name] = hdr0[name]
        elif filler == '1':
            metadata[name] = hdr1[name]
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

    
def save_to_fits(data, metadata, dataset_hdu, savepath, version):
    """
    Saves to a MUSCLES-standard fits file
    """
    if os.path.exists(savepath) == False:
        os.mkdir(savepath)
    file_name = make_component_filename(metadata, version)
    hdr = fits.Header(metadata)
    primary_hdu = fits.PrimaryHDU(header=hdr)
    hdu = fits.table_to_hdu(Table(data))
    descriptions =['midpoint of the wavelength bin', 'left/blue edge of the wavelength bin','right/red edge of the wavelength bin','average flux over the bin',
                'error on the flux','cumulative exposure time for the bin','data quality flags (HST data only)','modified julian date of start of first exposure', 
                'modified julian date of end of last exposure']
    hdu.header.insert(8, ('EXTNAME','SPECTRUM'))
    hdu.header.insert(9, ('EXTNO',2))
    [hdu.header.insert(i[0]+10, ('TDESC%s' %(i[0]), i[1])) for i in enumerate(descriptions)]
    hdul = fits.HDUList([primary_hdu, hdu, dataset_hdu])
    hdul.writeto(savepath+file_name+'.fits', overwrite=True)
    print('Spectrum saved as '+file_name+'.fits')    
    
    
def make_component_filename(metadata, version):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    filename = 'hlsp_muscles_%s_%s_%s_%s_v%s_component-spec' %(metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

def make_dataset_extension(hdr):
    """
    Makes a fits extension containg a list of rootnames and dataset IDs used in the spectrum.
    """
    description_text = 'This extension contains a list of observation IDs (DATASET_ID used for consistency with HST data) that can be used to locate the data in the CXO archives. CXO data all come from only a single observation (unlike the HST observations), but this extension is retained in place of a keyword for consistency with the HST files. '
    
    rootnames = ['',]
    datasets = [hdr['OBS_ID']]
    dataset_table = Table([rootnames,datasets], names=[('ROOTNAME'),('DATASET_ID')])
    hdu = fits.table_to_hdu(dataset_table)
    hdu.header.insert(8, ('EXTNAME','SRCSPECS'))
    hdu.header.insert(9, ('EXTNO',3))
    hdu.header['COMMENT'] = description_text
    return hdu
    
    
def make_chandra_spectra(chandra_path, evt_path, savepath, version, hlsp, apec_repo='', make_apec=True, save_ecsv=False, save_fits=False, grating =False, star=''):
    hdr = fits.getheader(evt_path, 1)
    if not grating:
        data_file = glob.glob(chandra_path+'*spectrum*')[0]
        data = build_chandra_data(data_file, hdr)
    else:
        data = build_from_txt(chandra_path, hdr)
    metadata = build_chandra_metadata(hdr, data, hlsp, star)
    if make_apec:
        model_file = glob.glob(chandra_path+'*model*')[0]                     
        apec_to_ecsv(model_file, metadata, apec_repo)
    if save_ecsv:
        save_to_ecsv(data, metadata, savepath, version)
    if save_fits:
        data_set_hdu = make_dataset_extension(hdr)
        save_to_fits(data, metadata, data_set_hdu, savepath, version)
 
    
    
