import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import scipy.interpolate as interpolate
# import stisblazefix


"""
Adjusts a STIS echelle spectrum using stis blaze fix and splices each order by uncertainty-weighted coadding. Also returns the ratio between flux measurements in the same wavelength bin.

20221012 now works on files with more than one data extension, adding them to the file name.

20230227 stisblazefix doesn't play well with the new stenv, turning off for now

"""



def make_x1f(x1dpath, outpath, rootname):
    """
    Turn off usetex and run stisblazefix. Turn usetex back on again it is used.
    """
    tex = plt.rcParams["text.usetex"]
    if tex:
        plt.rcParams.update({"text.usetex": False})
    stisblazefix.fluxfix([x1dpath], '{}{}_blazefix.pdf'.format(outpath, rootname))
    if tex:
        plt.rcParams.update({"text.usetex": True})

def spectra_adder(f_array, e_array, scale_correct=True):
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

def echelle_coadd_dq(wavelength, flux, err, dq, nclip =5, find_ratio =True, dq_adjust=False, dq_cut =0):
    """
    combines echelle orders into one spectrum, stiching them together at the overlap 
    """
    #slice dodgy ends off orders (usually 5-10 for stis el40m)
    wavelength = wavelength[:, nclip:-(nclip+1)]
    flux = flux[:, nclip:-(nclip+1)]
    err = err[:, nclip:-(nclip+1)]
    dq = dq[:, nclip:-(nclip+1)]
    
    #new arrays to put the output in
    w_full = np.array([], dtype=float)
    f_full = np.array([], dtype=float)
    e_full = np.array([], dtype=float)
    dq_full = np.array([], dtype=int)
    if find_ratio:
        r_full = np.array([], dtype=float) #ratio between orders

    shape = np.shape(flux)
    order = 0
    while order < (shape[0]):
        
        #first add the part that does not overlap ajacent orders to the final spectrum
        if order == 0: #first and last orders do not overlap at both ends
            overmask = (wavelength[order] > wavelength[order + 1][-1])
        elif order == shape[0]-1:
            overmask = (wavelength[order] < wavelength[order - 1][1])
        else:
            overmask = (wavelength[order] > wavelength[order + 1][-1]) & (wavelength[order] < wavelength[order - 1][1])
        w_full = np.concatenate((w_full, wavelength[order][overmask]))
        f_full = np.concatenate((f_full, flux[order][overmask]))
        e_full = np.concatenate((e_full, err[order][overmask]))
        dq_full = np.concatenate((dq_full, dq[order][overmask]))
        if find_ratio:
            r_full = np.concatenate((r_full, np.full(len(err[order][overmask]), -1)))
  
        if order != shape[0]-1:
            
            #interpolate each order onto the one beneath it, with larger wavelength bins. Code adapted from stisblazefix
            f = interpolate.interp1d(wavelength[order + 1], flux[order + 1], fill_value='extrapolate')
            g = interpolate.interp1d(wavelength[order + 1], err[order + 1], fill_value='extrapolate')
            dqi = interpolate.interp1d(wavelength[order + 1], dq[order + 1], kind='nearest',bounds_error=False, fill_value=0)
            overlap = np.where(wavelength[order] <= wavelength[order + 1][-1])
            f0 = flux[order][overlap]
            f1 = f(wavelength[order][overlap])
            g0 = err[order][overlap]
            g1 = g(wavelength[order][overlap])
            dq0 = dq[order][overlap]
            dq1 = dqi(wavelength[order][overlap])
       
             
            #combine flux and error at overlap and add to final spectrum
            w_av = wavelength[order][overlap]
            if dq_adjust: #removes values with high dq: #THIS DOESN'T REALLY WORK YET
                if dq_cut == 0:
                    dq_cut = 1 #allows zero to be the default
                for i in range(len(wavelength[order][overlap])):
                    if dq0[i] >= dq_cut:
                        g0 *= 100 #make error very large so it doesn't contribute to the coadd
                    if dq1[i] >= dq_cut:
                        g1 *= 100 
                        
            
            
            f_av, e_av = spectra_adder(np.array([f0,f1]),np.array([g0,g1]))
            dq_av = [(np.sum(np.unique(np.array([dq0, dq1])[:,i]))) for i in range(len(dq0))]
            
            
            w_full = np.concatenate((w_full, w_av))
            f_full = np.concatenate((f_full, f_av))
            e_full = np.concatenate((e_full, e_av))
            dq_full = np.concatenate((dq_full, dq_av))
            
            if find_ratio:
                r_full = np.concatenate((r_full, f0/f1))
        order += 1
    
    #stis orders are saved in reverse order, so combined spectra are sorted by the wavelength array
    arr1inds = w_full.argsort()
    sorted_w = w_full[arr1inds]
    sorted_f = f_full[arr1inds]
    sorted_e = e_full[arr1inds]
    sorted_dq = dq_full[arr1inds]
    if find_ratio:
        sorted_r = r_full[arr1inds]
 
    if find_ratio:
        return sorted_w, sorted_f, sorted_e, sorted_dq, sorted_r
    else:
        return sorted_w, sorted_f, sorted_e, sorted_dq
    
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

    
def plot_spectrum(data,rootname):
    plt.figure(rootname)
    mask = data['WAVELENGTH'] > 1160 #stis sensitivity falls off < 1160AA
    plt.step(data['WAVELENGTH'][mask], data['FLUX'][mask], where='mid')
    plt.step(data['WAVELENGTH'][mask], data['ERROR'][mask], where='mid')
    plt.xlabel('Wavelength (\AA)', size=20)
    plt.ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)', size=20)
    plt.tight_layout()
    plt.show()

    
def splice(filepath='data/', outpath='output/', nclip=5, save_fits=True, save_dat=True, plot=True, save_ecsv=False, blazefix = False):
    """
    Applies stisblazefix and splices all echelle x1d files in filepath. Results are saved in outpath.
    """
    if filepath[-1] != '/':
        filepath += '/'
    if outpath[-1] != '/':
        outpath += '/'
    x1ds = glob.glob('{}*x1d.fits'.format(filepath))
    if len(x1ds) == 0:
        print('No x1ds in folder (are you sure you unzipped them?)')
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    for x in x1ds:
        hdr = fits.getheader(x)
        if hdr['INSTRUME'] == 'STIS' and hdr['OPT_ELEM'][0] == 'E': #check it's an echelle:
            rootname = hdr['ROOTNAME']
            print(rootname)
            if blazefix:
                import stisblazefix
                make_x1f(x, outpath, rootname)
            nextend = fits.getheader('{}{}_x1f.fits'.format(filepath, rootname), 0)['NEXTEND']
            for i in range(nextend):
                data = fits.getdata('{}{}_x1f.fits'.format(filepath, rootname), i+1)
                w, f, e, dq, r = echelle_coadd_dq(data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ'], nclip=nclip)
                newdata = Table([w*u.AA, f*u.erg/u.s/u.cm**2/u.AA, e*u.erg/u.s/u.cm**2/u.AA, dq, r], names=['WAVELENGTH', 'FLUX', 'ERROR', 'DQ', 'ORDER_RATIO'])
                if nextend == 1:
                    extname = ''
                else:
                    extname = '_ext{}_'.format(i+1) #if more than one data extension, add the extension number to the file name
                if plot:
                    plot_spectrum(newdata, rootname)
                if save_ecsv:
                    ascii.write(newdata, '{}{}_{}spliced.ecsv'.format(outpath, rootname, extname), format='ecsv', overwrite=True)
                if save_dat:
                    savdat = Table([w,f,e], names=['#WAVELENGTH', 'FLUX', 'ERROR'])
                    ascii.write(savdat, '{}{}_{}spliced.dat'.format(outpath, rootname, extname), format='basic', overwrite=True)
                if save_fits:
                    primary_hdu = fits.PrimaryHDU(header=hdr)
                    hdu = fits.table_to_hdu(newdata)
                    hdu.name='SPECTRUM'
                    hdu.header = fits.getheader(x, i+1)
                    hdul = fits.HDUList([primary_hdu, hdu])
                    hdul.writeto('{}{}_{}spliced.fits'.format(outpath, rootname, extname), overwrite=True)
    print('Done')
            
if __name__ == "__main__": 
    splice()
    