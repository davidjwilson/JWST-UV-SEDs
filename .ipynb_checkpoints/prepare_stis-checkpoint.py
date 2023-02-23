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



cds.enable()

"""
@author: David Wilson

version 5 20230203


Finds all STIS x1d files, groups them by grating, coadds them and saves to file with required metadata. V3 fixes not knowning the exptimes for e140m

V5 updated to make it more adapatable to different programs and new error corrections

"""
def coadd_flux(f_array, e_array, scale_correct=True):
    """
    Returns a variance-weighted coadd with standard error of the weighted mean (variance weights, scale corrected).
    f_array and e_arrays are collections of flux and error arrays, which should have the same lenth and wavelength scale
    """
    if len(f_array) > 1:
        weights = 1 / (e_array**2)
        flux = np.average(f_array, axis =0, weights = weights)
        var = 1 / np.sum(weights, axis=0)
        rcs = np.sum((((flux - f_array)**2) * weights), axis=0) / (len(f_array)-1) #reduced chi-squared
        if scale_correct:
            error = (var * rcs)**0.5
        else:
            error = var**2
    else:
        flux, error = f_array[0,:], e_array[0,:]
    return flux,error

def no_zero_errors_old(flux, error):
    """
    Corrects instances where negative flux measurements have very small errors. Superceeded by make_person_errors
    """
    e_new = error
    for i in range(len(error)):
        if flux[i] < 0.0 and error[i] < 0.1*abs(flux[i]):
            e_new[i] = abs(flux[i])
    return e_new

def no_zero_errors(error):
    """
    Replaces any instances of error == 0 with the nearest value. Stops problems in coadding
    """
    while len(error[error == 0.0]) > 0: #clean out zeros
        for i in range(len(error)-1):
            if error[i] ==0:
                error[i] = error[i+1] #and replace them with the nearest value 
    return error

def make_person_errors(data, hdr):
    """
    Recalculates the error array using a pearson confidence interval - this corrects for the pipeline producing errors with zero
    """
    sensitivity = data['FLUX'] / data['NET'] #sensitivity curve of the spectrum 
    counts = data['GROSS'] * hdr['EXPTIME'] #total counts obtained in each wavelength bin during the exposure
    ci = stats.poisson_conf_interval(counts, interval='pearson') #pearson confidence interval
    ci = np.nan_to_num(ci, nan=0.0)
    count_errors = np.mean(abs(ci -counts), axis = 0) #average the upper and lower errorbars
    new_error = count_errors * sensitivity/hdr['EXPTIME'] #convert error to flux units    
    return new_error

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
                        
            
            
            f_av, e_av = coadd_flux(np.array([f0,f1]),np.array([g0,g1]))
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

def build_wavelength(x1ds):
    """
    builds a wavelength array covering all wavelength ranges in x1d (different cenwaves)
    """
    starts = []
    ends = []
    diffs = []
    for x in x1ds:
        nextend = fits.getheader(x, 0)['NEXTEND']
        for i in range(nextend):
            data = fits.getdata(x, i+1)
            for dt in data:
                w = dt['WAVELENGTH']
                starts.append(min(w))
                ends.append(max(w))
                diffs.append(np.max(np.diff(w)))
    w_new = np.arange(min(starts),max(ends), max(diffs))
    return w_new

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
    array = np.nan_to_num(array, nan=0.0)
    return array
 
def get_ayres_e140m(x1ds): 
    """ 
    reads T Ayres combined e140m files and bulds the data arrays from them. Hacky for now, will replace with my own routines when the new e140m calibrations are available.
    """
    target = fits.getheader(x1ds[0])['TARGNAME']
    savpath = '/home/david/work/muscles/SEDs/common/ayres_e140m/{}_E140M_coadd.sav'.format(target)
    data = readsav(savpath)
    w_new, f_new, e_new, dq_new, exptime = data['wave'], data['flux'], data['photerr'], data['epsilon'], data['texpt']
    return w_new, f_new, e_new, dq_new, exptime
    
def combine_x1ds(x1ds, correct_error=True, nclip=5):
    """
    coadds a collection of x1d fluxes and adds columns for exposure time detials. Input is a list of paths to x1d files with the same grating. Also works for sx1 files

    """
    # if len(x1ds) > 0:
      #   if fits.getheader(x1ds[0])['OPT_ELEM'] == 'E140M':
      # #  print('yes')
      #       w_new, f_new, e_new, dq_new, exptime = get_ayres_e140m(x1ds)
      #       start = []
      #       end = []
      #       for x in x1ds:
      #           hdr = fits.getheader(x,0)
      #           start.append(hdr['TEXPSTRT'])
      #           end.append(hdr['TEXPEND'])
      #       print('start', start)
      #       print('end', end)
      #       start, end = np.min(np.array(start)), np.max(np.array(end))
      #       start, end = np.full(len(w_new), start), np.full(len(w_new), end)

            
                
#         start, end = ,np.full(len(w_new), 0), np.full(len(w_new), 0) 
        # else:
    f_new = []
    e_new = []
    dq_new = []
    exptime = []
    start = []
    end = []
    w_new = build_wavelength(x1ds)
    for x in x1ds:
        hdr = fits.getheader(x,0)
        nextend = hdr['NEXTEND']
        for i in range(nextend):
            hdr1 = fits.getheader(x,i+1)
            data = fits.getdata(x, i+1)
            if hdr['OPT_ELEM'][0] == 'E':
                print(hdr['OPT_ELEM'], 'yes')
                wi, fi, ei, dqi = echelle_coadd_dq(data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ'], nclip=nclip, find_ratio=False)
            else:
                # print(hdr['OPT_ELEM'], 'No')
                data = data[0]
                wi, fi, ei, dqi = data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ']
                if correct_error and hdr['OPT_ELEM'] in ['G140M', 'G140L']:    
                    ei = make_person_errors(data, hdr1)
            fi = interpolate.interp1d(wi, fi, bounds_error=False, fill_value=0.)(w_new)
            ei = interpolate.interp1d(wi, ei, bounds_error=False, fill_value=0.)(w_new)
            dqi =  interpolate.interp1d(wi, dqi, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
            expi = np.full(len(wi), hdr1['EXPTIME'])
            expi = interpolate.interp1d(wi, expi, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
            starti = np.full(len(wi), hdr1['EXPSTART'])
            starti = interpolate.interp1d(wi, starti, kind='nearest',bounds_error=False, fill_value=0.)(w_new)
            endi = np.full(len(wi), hdr['TEXPEND'])
            endi = interpolate.interp1d(wi, endi, kind='nearest',bounds_error=False, fill_value=0.)(w_new)


            f_new.append(fi)
            e_new.append(ei)
            dq_new.append(dqi)
            exptime.append(expi)
            start.append(starti)
            end.append(endi)


    f_new, e_new = coadd_flux(np.array(f_new), np.array(e_new))
    dq_new = np.array(dq_new, dtype=int)
    dq_new = [(np.sum(np.unique(dq_new[:,i]))) for i in range(len(dq_new[0]))]
    exptime = np.sum(np.array(exptime), axis=0)
    start = np.min(np.ma.masked_array(start, mask=[np.array(start) == 0.]), axis=0)
    end = np.max(np.array(end), axis=0)

    # else: #in the case where there's only one available spectrum
    #     data_extension = 1
    #     # if x1ds[0][-8:-5] == 'sx1': #modified 1 off for t1 spectrum, must improve later
    #         # data_extension = 0
    #       #  
    #     #else:
    #     data = fits.getdata(x1ds[0],data_extension)[0]   
    #     hdr = fits.getheader(x1ds[0],0)
    #     w_new, f_new, e_new, dq_new = data['WAVELENGTH'], data['FLUX'], data['ERROR'], data['DQ']
    #     exptime, start, end = np.full(len(data['WAVELENGTH']), hdr['TEXPTIME']), np.full(len(data['WAVELENGTH']), hdr['TEXPSTRT']), np.full(len(data['WAVELENGTH']), hdr['TEXPEND'])
    #     if correct_error and hdr['OPT_ELEM'] in ['G140M', 'G140L']:    
    #                 enew = make_person_errors(data, hdr1)
   
    f_new, e_new = nan_clean(f_new), nan_clean(e_new)
    w0, w1 = wavelength_edges(w_new)
    new_data = {'WAVELENGTH':w_new*u.AA,'WAVELENGTH0':w0*u.AA,'WAVELENGTH1':w1*u.AA,'FLUX':f_new*u.erg/u.s/u.cm**2/u.AA,
                'ERROR':e_new*u.erg/u.s/u.cm**2/u.AA,'EXPTIME':exptime*u.s,'DQ':dq_new,'EXPSTART':start*cds.MJD,'EXPEND':end*cds.MJD}
   

    return new_data


 
def make_metadata(x1ds, new_data, hlsp, normfac, star):
    """
    Makes the metadata for the ecsv files- eventually will be converted into the fits file
    """
    wavelength, flux = new_data['WAVELENGTH'].value, new_data['FLUX'].value
    mid = int(len(wavelength) / 2)
    specres = wavelength[mid]
    waveres = wavelength[mid+1] - wavelength[mid]
    exptimes = []
    start_times = []
    end_times = []
    dates = []
    for x in x1ds:      
        hdr = fits.getheader(x,0)
        exptimes.append(hdr['TEXPTIME'])
        start_times.append(hdr['TEXPSTRT'])
        end_times.append(hdr['TEXPEND'])
        dates.append(hdr['TDATEOBS'])
        if star == '':
            starname = hdr['TARGNAME']
        else:
            starname = star
    # meats_name = 'MUSCLES Extension for Atmospheric Transmisson Spectroscopy'
    meta_names =  ['TELESCOP', 'INSTRUME','GRATING','APERTURE','TARGNAME','RA_TARG','DEC_TARG','PROPOSID','HLSPNAME','HLSPACRN','HLSPLEAD','PR_INV_L',
                   'PR_INV_F','DATE-OBS','EXPSTART','EXPEND','EXPTIME','EXPDEFN','EXPMIN','EXPMAX','EXPMED','NORMFAC','WAVEMIN','WAVEMAX','WAVEUNIT','AIRORVAC','SPECRES','WAVERES','FLUXMIN',
                  'FLUXMAX','FLUXUNIT']
    meta_fill = ['','',hdr['OPT_ELEM'],'',starname,'','','',hlsp['PROGRAM'],hlsp['PROGRAMSHORT'],hlsp['HLSPLEAD'],
                 '','',min(dates),min(start_times),max(end_times),sum(exptimes),'SUM', 
                min(exptimes), max(exptimes), np.median(exptimes),normfac,wavelength[0], wavelength[-1],'ang','vac',specres,waveres,np.min(flux[np.isnan(flux)==False]), np.max(flux[np.isnan(flux)==False]),'erg/s/cm2/ang']
    metadata = {}
    for name, filler in zip(meta_names, meta_fill):
        if filler == '':
            metadata[name] = hdr[name]
        else:
            metadata[name] = filler
    return metadata
    
def setup_list(x1ds):
    """
    Takes all x1ds in selection and sorts them by instrument setup
    """
    gratings = []
    x1ds_by_setup = []
    for x in x1ds:
        hdr = fits.getheader(x,0)
        grating = hdr['OPT_ELEM']
      
        gratings.append(grating)
    setups = np.unique(gratings, axis=0)
    for s in setups:
        collection = []
        for i in range(len(x1ds)):
            if gratings[i] == s:
                collection.append(x1ds[i])
        if len(collection) > 0:
            x1ds_by_setup.append(collection)
    return setups, x1ds_by_setup

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

def stis_clean(x1ds):
    """
    checks that all x1d files are stis spectra
    """
    stis_x1ds =[]
    for x in x1ds:
        if fits.getheader(x,0)['INSTRUME'] == 'STIS':
            stis_x1ds.append(x)
    return stis_x1ds

def make_component_filename(metadata, version, hlsp):
    """
    Construct the filename for a component spectrum ,eg "hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits"hlsp_muscles_hst_stis_gj176_g230l_v22_component-spec.fits
    """
    # filename = 'hlsp_muscles_%s_%s_%s_%s_v%s_component-spec' %(metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    filename = 'hlsp_{}_{}_{}_{}_{}_v{}_component-spec'.format(hlsp['HLSPNAME'].lower(),metadata['TELESCOP'].lower(), metadata['INSTRUME'].lower(), metadata['TARGNAME'].lower(), metadata['GRATING'].lower(), version) 
    return filename

def make_dataset_extension(x1ds):
    """
    Makes a fits extension containg a list of rootnames and dataset IDs used in the spectrum.
    """
    description_text = 'This extension contains a list of HST rootnames (9 character string in HST files downloaded from MAST) and dataset IDs of the exposures used to create this spectrum file. The dataset IDs can be used to directly locate the observations through the MAST HST data archive search interface. Multiple identifiers indicate the spectra were coadded. In some cases the spectra have been rextracted and will be different from those availbale from MAST.' 
    
    rootnames = []
    datasets = []
    for x in x1ds:
        hdr = fits.getheader(x)
        rootnames.append(hdr['ROOTNAME'],)
        datasets.append(hdr['ASN_ID'])
    dataset_table = Table([rootnames,datasets], names=[('ROOTNAME'),('DATASET_ID')])
    hdu = fits.table_to_hdu(dataset_table)
    hdu.header.insert(8, ('EXTNAME','SRCSPECS'))
    hdu.header.insert(9, ('EXTNO',3))
    hdu.header['COMMENT'] = description_text
    return hdu

    
def make_stis_spectrum(x1dpath, version, hlsp, savepath = '', plot=False, save_ecsv=False, save_fits=False, return_data=False, return_gratings = False, normfac=1.0, star = '', nclip=5):
    """
    main function
    """
    all_x1ds = np.hstack((glob.glob('{}*_x1d.fits'.format(x1dpath)), glob.glob('{}*sx1.fits'.format(x1dpath))))
    stis_x1ds = stis_clean(all_x1ds)
    hlsp = Table.read(hlsp)[0]
    if len(stis_x1ds) > 0:
        gratings, x1ds_by_grating = setup_list(stis_x1ds)
        # print(gratings, x1ds_by_grating)
        for x1ds in x1ds_by_grating:
            # print(x1ds)
            data = combine_x1ds(x1ds, nclip)
           # data = [wavelength*u.AA, flux*u.erg/u.s/u.cm**2/u.AA, error*u.erg/u.s/u.cm**2/u.AA, dq]
            metadata = make_metadata(x1ds, data, hlsp, normfac, star)
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
        return gratings


# def test():
#     """
#     testing with GJ 699
#     """
#     x1dpath = '/home/david/work/muscles/MegaMUSCLES/GJ_699/HST/STIS/'
#     version = 2
#     savepath = 'test_files/'
#     make_stis_spectum(x1dpath, version, savepath = savepath, plot=True, save_ecsv=True, save_fits=True)
    
# def e140m_redo():
#     stars = ['GJ15A', 'GJ729']
#     for star in stars:
#         x1dpath = '/media/david/1tb_storage1/emergency_data/mega_muscles/e140m_fits/{}/'.format(star)
#         version = 2
#         savepath = 'e140m_test/'
#         make_stis_spectum(x1dpath, version, savepath = savepath, plot=True, save_ecsv=True, save_fits=True)
    
# e140m_redo()

    
#test()