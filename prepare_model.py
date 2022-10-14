import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import griddata, interp1d


"""
@author: David Wilson

@version: 24 

@date :20221014

Turns models into standard MUSCLES file. Added capability for errors
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

def get_model_data(model_path):
    """
    Makes the model data array, assuming an input .txt file with WAVELENGTH and FLUX
    """
    data = Table.read(model_path, format = 'ascii')
    w, f = data['WAVELENGTH'], data['FLUX']
    w0, w1 = wavelength_edges(w)
    new_data = {'WAVELENGTH':w*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f*u.erg/u.s/u.cm**2/u.AA}
    if 'ERROR' in data.dtype.names:
        new_data['ERROR'] = data['ERROR']*u.erg/u.s/u.cm**2/u.AA 
    return new_data
    
def make_model_metadata(new_data, normfac, sed_metadata, model_name):
    """
    Makes the metadata for the lya file -sed metadata is from the total sed so far (should just be COS?)
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['Model',model_name,'NA','','','','','','','','','',normfac,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = sed_metadata[name]
        else:
            metadata[name] = filler
    return metadata

def make_component_filename(metadata, version):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    filename = 'hlsp_muscles_%s_%s_%s_%s_v%s_component-spec' %(metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

def save_to_ecsv(data, metadata, save_path, version):
    """
    save the new model to an ecsv file
    """
    print(save_path)
    if os.path.exists(save_path) == False:
        os.mkdir(save_path)
    file_name = make_component_filename(metadata, version)
    savedat = Table(data, meta=metadata)
    savedat.write(save_path+file_name+'.ecsv', overwrite=True, format='ascii.ecsv')
    print('Spectrum saved as '+file_name+'.ecsv')

    
def model_save_to_fits(data, metadata, savepath, version):
    """
    Saves to a MUSCLES-standard fits file for models
    """
    if os.path.exists(savepath) == False:
        os.mkdir(savepath)
    file_name = make_component_filename(metadata, version)
    hdr = fits.Header(metadata)
    primary_hdu = fits.PrimaryHDU(header=hdr)
    hdu = fits.table_to_hdu(Table(data))
    descriptions =['midpoint of the wavelength bin', 'left/blue edge of the wavelength bin','right/red edge of the wavelength bin','average flux over the bin']
    hdu.header.insert(8, ('EXTNAME','SPECTRUM'))
    hdu.header.insert(9, ('EXTNO',2))
    [hdu.header.insert(i[0]+10, ('TDESC%s' %(i[0]), i[1])) for i in enumerate(descriptions)]
    hdul = fits.HDUList([primary_hdu, hdu])
    hdul.writeto(savepath+file_name+'.fits', overwrite=True)
    print('Spectrum saved as '+file_name+'.fits')

def make_model_spectrum(model_path, version, sed_data ,savepath = '', save_ecsv=False, save_fits=False, normfac=1.0, model_name='LYA-RECONSTRUCTION'):
    """
    Main function.
    """
    data = get_model_data(model_path)
    metadata = make_model_metadata(data,  normfac, sed_data.meta, model_name)
    if save_ecsv:
        save_to_ecsv(data, metadata, savepath, version)
    if save_fits:
        model_save_to_fits(data, metadata, savepath, version)
        
        
