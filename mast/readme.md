README for Measurements of the Ultraviolet Spectral Characteristics of Low-mass Exoplanetary Systems (MUSCLES)
MAST webpage: https://archive.stsci.edu/hlsp/muscles/
DOI: https://doi.org/10.17909/T9DG6F

# Principal Investigators
Kevin France, Cynthia Froning, Allison Youngblood 

# HLSP authors
Patrick Behr, Girish Duvvuri, Kevin France, Cynthia Froning, Parke Loyd, David Wilson, Allison Youngblood

The MUSCLES HLSP consists of several surveys providing panchromatic spectral energy distributions (SEDs) of an increasing number of stars. The methods and data sources used to construct the SEDs in each survey are mostly similar, but small improvements have been made over time. This readme provides information general to all surveys first, followed by details specific to each survey.

The four main surveys providing SEDs to this HLSP and their version numbers are:
1. MUSCLES, v22: 7 M and 4 K dwarfs (Loyd et al 2016).
2. MUSCLES Extension (MEATS), v24, 11 JWST GTO targets (Behr et al. 2023).
3. Mega-MUSCLES, v25: 12 M dwarfs (Wilson et al. 2025).
4. Mega-MEATS, v26: 24 JWST Cycle 1 and HWO targets (Wilson et al. submitted).

# Complete Target List

| Star | Survey |
|------|--------|
|	GJ 1214	|	MUSCLES	|
|	GJ 876	|	MUSCLES	|
|	GJ 436	|	MUSCLES	|
|	GJ 581	|	MUSCLES	|
|	GJ 667C	|	MUSCLES	|
|	GJ 176	|	MUSCLES	|
|	GJ 832	|	MUSCLES	|
|	HD 85512	|	MUSCLES	|
|	HD 40307	|	MUSCLES	|
|	HD 97658	|	MUSCLES	|
|	eps Eri	|	MUSCLES	|
|	Proximia Cen	|	MUSCLES	|
|	WASP-17	|	MEATS	|
|	HD 149026	|	MEATS	|
|	WASP-127	|	MEATS	|
|	WASP-77 A	|	MEATS	|
|	LTT 9779	|	MEATS	|
|	HAT-P-26	|	MEATS	|
|	HAT-P-12	|	MEATS	|
|	WASP-43	|	MEATS	|
|	L 678-39	|	MEATS	|
|	L 98-59	|	MEATS	|
|	LP 791-18	|	MEATS	|
|	GJ 676 A	|	Mega-MUSCLES	|
|	GJ 15 A	|	Mega-MUSCLES	|
|	GJ 649	|	Mega-MUSCLES	|
|	GJ 674	|	Mega-MUSCLES	|
|	GJ 729	|	Mega-MUSCLES	|
|	GJ 163	|	Mega-MUSCLES	|
|	GJ 1132	|	Mega-MUSCLES	|
|	L 980-5	|	Mega-MUSCLES	|
|	GJ 849	|	Mega-MUSCLES	|
|	GJ 699	|	Mega-MUSCLES	|
|	LHS 2686	|	Mega-MUSCLES	|
|	TRAPPIST-1	|	Mega-MUSCLES	|
|	WASP-63	|	Mega-MEATS	|
|	tau Ceti	|	Mega-MEATS	|
|	WASP-121	|	Mega-MEATS	|
|	HIP 67522	|	Mega-MEATS	|
|	kappa Ceti	|	Mega-MEATS	|
|	Kepler-51	|	Mega-MEATS	|
|	GJ 367	|	Mega-MEATS	|
|	K2-18	|	Mega-MEATS	|
|	TOI-776	|	Mega-MEATS	|
|	GJ 341	|	Mega-MEATS	|
|	TOI-421	|	Mega-MEATS	|
|	TOI-134	|	Mega-MEATS	|
|	WASP-166	|	Mega-MEATS	|
|	NGTS-10	|	Mega-MEATS	|
|	TOI-260	|	Mega-MEATS	|
|	TOI-836	|	Mega-MEATS	|
|	HD 15337	|	Mega-MEATS	|
|	K2-141	|	Mega-MEATS	|
|	TOI-178	|	Mega-MEATS	|
|	HD 80606	|	Mega-MEATS	|
|	HATS-72	|	Mega-MEATS	|
|	LHS 475	|	Mega-MEATS	|
|	eps Indi	|	Mega-MEATS	|
|	70 Oph B 	|	Mega-MEATS	|

# Data Product Description (All surveys)
MUSCLES data products for a given star consist of an SED in four formats, along with the subspectra that the SED is constructed from. The data sources for the SEDs are:

- **X-rays**: *Chandra*/*XMM-Newton*/*Swift* and APEC models (Smith+ 2001; ApJ 556:L91)
- **EUV**: Empirical scaling relation based on Lya flux (Linsky+ 2014; ApJ 780:61) and/or Diferential Emission Measure models (Duvvuri+ 2021; ApJ 913:40)
- **Lya**: Reconstructed spectrum from model fits to Lyα line profiles accounting for the stellar flux and the interstellar medium (Youngblood+ 2016; ApJ 913:40)
- **FUV - blue visible**: *HST* COS and STIS
- **Visible - IR**: synthetic photospheric spectra from PHOENIX atmosphere models (Husser+ 2013; A&A 553:A6) and (Allard 2016)

The exact wavelength covered by each subspectrum vary by SED. 

## Data file naming convention for SEDs:

hlsp_muscles_<tel>_<inst>_<target>_<fgrating>_<ver>_<prodType>.fits

where:

    <tel/model> : telescope used or model, "cxo", "hst", "model", or "xmm",
    <inst/model> : instrument used such as "acis", "cos", "fos" or "stis",  or model such as "apec", "euv-scaling", "phx", "lya-reconstruction", or "multi" ,   
    <target> : target name such as "trappist-1" or "l-980-5",
    <grating> : grating used such as "g130m" or "na" for not applicable for models,
    <ver>: version number, "v22", "v23", "v24," "v25", and
    <prodType> : product type such as "const-res-sed"  or "component-spec".

Primary data file types:
- _const-res-sed.fits 	Panchromatic SED binned to a constant 1 Å resolution. One file per star.
- _var-res-sed.fits 	Panchromatic SED that retains the native instrument resolutions where possible. One file per star.
- _adapt-const-res-sed.fits 	Panchromatic SED binned to a constant 1 Å resolution, downsampled in low signal-to-noise regions to avoid negative fluxes. One file per star.
- _adapt-var-res-sed.fits 	Panchromatic SED that retains the native instrument resolutions where possible, downsampled in low signal-to-noise regions to avoid negative fluxes. One file per star.
 
Ancillary data file products. The exact products vary by star:
- _\[telescope\]_\[instrument\]_\[grating\]*component-spec.fits Observed, coadded spectrum. One file per observed grating
- _apec*component-spec.fits APEC model X-ray spectrum.
- _euv-scaling*component-spec.fits 	Extreme UV spectrum based on empirical scaling relation using Lyα flux. 
- _dem*component-spec.fits 	Extreme UV spectrum based on Differential Emission Measure Model. 
- _lya-reconstruction*component-spec.fits 	Reconstructed spectrum from model fits to Lyα line profiles accounting for the stellar flux and the interstellar medium.
- _phx*component-spec.fits 	Synthetic spectrum in the visible and infrared from PHOENIX models. 

All the spectra have primary extension with a header that gives information on the observation(s) that make up the spectrum. The spectrum itself is in the first extension, named 'SPECTRUM'. For all header values, 'MULTI' indicates multiple values for the field, which will be followed by a numbered list of the various values (e.g. a MULTI value for INSTRUME will be followed by INSTRU00, INSTRU01, ... to specify the various instruments that contributed data to the spectrum).


### Primary Header Keywords
Those keywords that are not self-explanatory are described below. Some keywords are omitted for data products to which they do not apply or are not well-defined (e.g. EXPTIME in the panchromatic SEDs, DATEOBS for model spectra).


- TARGNAME: name of the target star
- RA_TARG : right ascension coordinate of the target (decimal degrees)
- DEC_TARG: declination coordinate of the target (decimal degrees)
- PROPOSID: *HST* proposal number for MUSCLES (always 13650)
- HLSPNAME: name of this HLSP (always 'Measurements of the Ultraviolet Spectral Characteristics of Low-mass Exoplanet Host Stars')
- HLSPACRN: acronym for the HLSP (always 'MUSCLES')
- HLSPLEAD: lead preparer of the data product (always R. O. Parke Loyd)
- PR_INV_L: principle investigator last name (always France)
- PR_INV_F: principle investigator first name (always Kevin)
- DATE-OBS: date of the start of the first observation in YYYY-MM-DDTHH:MM:SS.SSS
- EXPSTART: modified julian date of the start of the first observation
- EXPEND  : modified julian date of the end of the last observation
- EXPTIME : characteristic exposure time in s, determined using the method specified in EXPDEFN
- EXPDEFN : method used to determine EXPTIME value
- EXPMAX  : maximum exposure time for a spectral element in s
- EXPMIN  : minimum nonzero exposure time for a spectral element in s
- EXPMED  : median exposure time in s
- NORMFAC : normalization factor applied to the spectrum before merging into the composite spectrum and/or before use for the Lya reconstruction (see Loyd+ 2106; ApJ in press)
- WAVEMIN : minimum wavelength in the spectrum in Angstroms
- WAVEMAX : maximum wavelength in the spectrum in Angstroms
- AIRORVAC: wavelengths specified in their vacuum or air values (always 'vac')
- SPECRES : wavelength at which WAVERES us specified in Angstroms
- WAVERES : wavelength resolution at SPECRES in Angstroms
- FLUXMIN : minimum flux in the spectrum in erg s-1 cm-2 Angstroms-1
- FLUXMAX : maximum flux in the spectrum in erg s-1 cm-2 Angstroms-1
- BOLOFLUX: bolometric flux measurement (integral of flux >5 Angstroms) for the star, including a blackbody fit for flux >5.5 microns that contributes a few percent to the integral in most cases. Because of the included tail, this is NOT the same value as one will find by integrating the full extend of the SED.


### Spectrum Header Keywords
The spectrum extension of the FITS files contains just the default required keywords for the FITS standard, with the addition of a TDESC keyword associated with each data column (along with the standard TTYPE, TFORM, and TUNIT keywords) that describes the data in that column.


### Spectrum Data Columns
The spectrum extension of the FITS files is a data table that contains these columns (according to the type of spectrum):


All spectra:


- WAVELENGTH : midpoint of the wavelength bin in Angstroms
- WAVELENGTH0: left (blue) edge of the wavelength bin in Angstroms
- WAVELENGTH1: right (red) edge of the wavelength bin in Angstroms
- FLUX       : average flux density in the wavelength bin in erg s-1 cm-2 Angstroms-1


All *except* model spectra:


- DQ         : data quality flags propagated from the original observation (defined by the instrument, NOT defined by MUSCLES -- see docs for that instrument for DQ flag definitions)
- ERROR      : error on the flux density in erg s-1 cm-2 Angstroms-1. Mega-MUSCLES and Mega-MEATS also include estimated uncertainties in model spectra. 
- EXPEND     : modified julian date of the end of the last exposure contributing data to the bin
- EXPSTART   : modified julian date of the start of the first exposure contributing data to the bin
- EXPTIME    : total exposure time of observations contributing to that bin, averaged according to bin widths when two or more elements have been rebinned to a coarser resolution


SEDs:


- BOLOFLUX   : flux density normalized by the bolometric integral of the flux given in the primary header as BOLOFLUX, in Angstroms-1
- BOLOERR    : error on BOLOFLUX
- INSTRUMENT : binary flag(s) identifying the source spectrum or spectra for the data in the bin. These values are defined in the second FITS extension, INSTLGNd. When rebinning, some bins acquire data from multiple sources. In these cases, the binary flags are combined in a bitwise OR sense (e.g. 010 (2) and 100 (4) becomes 110 (6)).
- NORMFAC    : normalization factor applied to the source spectrum when merging into the composite SED
- LNZ_GAM, LNZ_NORM : (MUSCLES only) the width and scale factor the best fit to the Lorentzian wings of the Lya profile as measured in the COS data (see the notes for Version 2.1 earlier in this readme). The equation of the fit is then flux_density = lnz_norm * (lnz_gam/2/((w - w0) + (lnz_gam/2)**2) where w is wavelength and w is the center of the reconstructed Lya profile.


### Component Spectrum SRCSPECS Extension
The second extension of the component spectrum source files, named SRCSPECS, gives identifiers that can be used to locate the observatory-level data used to create the spectrum. For *HST* observations, both DATASET_ID values and ROOTNAME are given, since sometimes multiple files with different ROOTNAMEs are associated with a single observation.


### SED INSTLGND Extension
The SED FITS files have a second extension named INSTLGND that provides the information needed to associate component spectra with the binary keys given in the INSTRUMENT column of the SPECTRUM extension. The included data columns are:


- BITVALUE   : value corresponding to the INSTRUMENT column of the SPECTRUM extension
- GRATING    : grating used for observation
- HLSP_FILE  : name of the corresponding HLSP file
- INSTRUMENT : instrument used for observation
- TELESCOPE  : telescope used for observation

Negative Flux Bins
------------------
Data in regions of low flux contain some bins with negative flux values, as expected due to statistical noise in the signal and subtracted background levels. The physical impossibility of negative flux poses a problem for some applications of these spectra. *We caution the user against simply setting these bins to zero flux. This introduces an upward bias in the flux over the broadband FUV continuum of several percent to over 10% in the worst case.*


We have created "adaptive resolution" SEDs as a solution that relies only on the data (as opposed to further models or interpolations) for applications where negative-flux bins are undesirable. For these spectra, we iteratively averaged negative-flux bins with their neighbors until no negative-flux bins remain. This results in spectra that are effectively downsampled in low-flux regions to improve S/N to a point where the flux over the larger wavelength range can be measured. As with the unaltered SEDs, we provide both `var-res` (variable resolution) and `const-res` (constant 1 Angstrom resolution) versions of these spectra. The `var-res` spectra retain the (often very) coarse resolution in the low flux regions of the SED. In contrast, the 1 Angstrom binning of the `const-res` SEDs implies that areas of the spectrum that have been downsampled to a resolution coarser than 1 Angstrom will then have been upsampled. We caution the user that these regions will appear to have a fidelity higher than their true fidelity.  

Users can also implement their own solutions to the problem of negative-flux bins. An (incomplete) list of solutions that users could implement include replacing the low-flux regions of the SEDs with:


- PHOENIX model output (included in this HLSP). This could serve as a lower limit (i.e. photosphere-only) flux estimate in these regions.
- Data from one of the brighter targets (appropriately renormalized).
- A custom model of stellar upper-atmosphere emission, such as that of Fontenla+ 2016 in the case of GJ 832.

# Details for specific programs

MUSCLES (Version 2.2, v22)
------------------------------------
See Loyd et al. 2016 for full details

- All: The PHOENIX spectra obtained from Husser+ 2013 grid truncate at 5.5 microns. To compute a more accurate bolometric flux, we fit a blackbody spectrum to the data near this edge when computing the full integral of the spectrum. This value is included in the SED header and is used when computing the `boloflux` (flux normalized by the bolometric integral) values.
- All: *HST* STIS data showed systematically lower absolute fluxes and were normalized upward to match corresponding *HST* COS data where statistically justifiable. Note that the Eps Eri *HST* STIS E230M and E230H data are suspect, but no *HST* COS data were available to normalize against.
- All (esp. GJ 832, GJ 876): Flares were observed in some of the data. Data from the flaring times was *not* discarded. Particularly large flares that occurred during the GJ 832 and GJ 876 *HST* COS G130M observations.
- All but Eps Eri: A miscalibration in the wavelength solution of the *HST* COS G230L data by as much as 4 Angstroms (~600 km s-1) appears to be present near the short wavelength edge (~1800 Angstroms).
- Eps Eri, GJ 1214, GJ 876, GJ 436, GJ 581: Wavelengths of the *HST* STIS E230H and G230L data appear miscalibrated at the Mg II 2796,2803 lines by 1-2 Angstroms (165-330 km s-1).
- GJ 1214, GJ 667C: STIS G430L data blueward of ~3850 Angstroms was culled in preference for PHOENIX model because of low S/N. (In all other cases, data was always given preference over model.)
- GJ 667C: slit aperture was significantly off-center in *HST* STIS observations. Normalization factors to match STIS data to COS data are consequently very high.
- GJ 1214, GJ 832, GJ 581, GJ 436: *HST* pipeline extraction of the STIS G140M data failed to locate the spectrum on the detector. The spectrum was "manually" extracted from the two dimensional data (x2d). This was also necessary for the STIS G230L data for GJ 1214.
- GJ 1214: *Chandra* X-ray data provided only upper limits. We used a model fit to archival *XMM-Newton* data instead.
- HD 97658: No X-ray data were collected. We used the HD 85512 data scaled according to the ratio of bolometric flux between the two stars because of their similar levels of Fe XII emission relative to the bolometric flux.
- GJ 581, GJ 876: Two markedly different levels of X-ray activity were observed. We only use data from the more quiescent levels.
- GJ 436: Evidence of a faint second source is present in the *HST* STIS G230L data and may contribute flux to all observations.

In response to popular demand after the announcement of the discovery of a habitable-zone planet orbiting Prox Cen, we added Prox Cen to the MUSCLES sample and created an SED from archival data using methodology consistent with the other MUSCLES SEDs. We highly encourage users of the Prox Cen spectra to read [these notes](https://archive.stsci.edu/missions/hlsp/muscles/gj551/hlsp_muscles_multi_multi_gj551_broadband_v22_reduction-notes.pdf) on the reduction of the Prox Cen spectra, as they address an important issue regarding the stellar effective temperature and bolometric flux.  

Also available for download is a synthetic spectrum of the MUSCLES M2 dwarf GJ 832 from Fontenla et al. 2016. The synthetic spectrum was created from a one-dimensional model of the stellar atmosphere that incorporates non-LTE radiative transfer techniques and many molecular lines. The synthetic spectrum covers 1 Å - 1 mm at a resolving power R = 100,000 (denoted by "r1e5" in the file names). A single resolution element equals two wavelength samples. For consistency with the MUSCLES data products, the synthetic spectrum's wavelength is in Angstroms and the flux is reported as would be observed from Earth (erg/cm2/s/Å) under the assumption that the stellar radius = 0.499 R_solar (Houdebine 2010) and distance = 4.95 pc (van Leeuwen 2007). Note that Houdebine et al. 2016 determine GJ 832's stellar radius to be 0.458 ± 0.039 R_solar. The atomic and molecular transitions appear at their laboratory rest wavelengths (vacuum). Semi-empirical models of the other 10 MUSCLES stars as well as Proxima Centauri and TRAPPIST-1 are coming soon. If you utilize the GJ 832 synthetic spectrum, please cite Fontenla et al. 2016, ApJ, 830, 154. 

Data file naming convention for semi-empirical models for GJ 832:

hlsp_muscles_model_<waveband>_gj832-r1e5_na_v1_synth-spec.fits

where:    
<waveband>:  the waveband of the synthetic spectrum, "xray", "euv", "fuv", "mir", "fir, or "all"

MUSCLES Extension (Version 2.4, v24)
------------------------------------
See Behr et al. 2023 for full details.

  - The MUSCLES Extension SEDs are similar to the MUSCLES SEDs and follow the same file naming conventions. The readme for the MUSCLES SEDs applies in all cases except where detailed below.\n Major differences between MUSCLES and MUSCLES Extension are as follows:


  - FUV and X-ray proxy stars: Due to low signal to noise in FUV and X-ray observations, the MUSCLES Extension SEDs use scaled proxy stars of similar spectral type and age for these regions. The regions where proxy stars are used is typically 5 < λ < 100 Å and 1170 < λ < ~1900 Å. The end of the proxy spectrum varies from target-to-target depending on the S/N of the STIS G230L spectra.

  - Flux Scalings: No adjustments have been made to the STIS flux calibration.
  PHOENIX models: MUSCLES Extension uses PHOENIX photospheric model spectra from the Lyon BT-Settl CIFIST 2011_2015 grid (Allard 2016). The models extend to 5.5 microns.



Mega-MUSCLES (Version 2.5, v25)
-------------------------------------------------------------
See Wilson et al. 2025 for full details.

The Mega-MUSCLES SEDs are similar to all previous SEDs and follow the same file naming conventions. Major differences between MUSCLES and this release of Mega-MUSCLES are as follows:

  - Flux uncertainties are provided for the entire SED, instead of only the observed spectra. For the Lyman alpha profiles and DEMs we use the 1-sigma uncertainty arrays provided by their respective modeling routines. For the X-ray APEC models uncertainties are estimated using a Monte Carlo method based on the uncertainties on their fitted parameters. For the PHX models the uncertainties are based on the uncertainty in Teff, as that was found to dominate the flux uncertainty over other parameters. 

  - Archival data is used more extensively than in previous releases, including both X-ray and archival HST data. All archival data has been carefully checked for data quality and consistency with new data. All datasets used are listed in the relevant dataset extensions.

  - Differential Emission Measure (DEM) models are provided for all SEDs.

  - Flux scalings: No adjustments have been made to the STIS flux calibration as: 1. The COS spectra were generally of too poor a signal to provide the same baseline measurements used in MUSCLES, and 2. The optical G430L STIS spectra were found to be in good agreement with photometry and other, reliably flux calibrated spectra.

  - PHOENIX models: Mega-MUSCLES uses PHOENIX photospheric model spectra from the Lyon BT-Settl CIFIST 2011\_2015 grid (Allard 2016), which extend out to 9995000.0 Angstroms. “const-res” SEDs, which are binned to 1 Angstrom resolution, only extend out to 130000 Angstroms to avoid the file size becoming very large.


Mega-MEATS (Version 2.6, v26)
-------------------------------------------------------------
See Wilson et al. submitted for full details.

- All the changes described for Mega-MUSCLES also apply here. 

  - Proxy spectra are extensively used to replace low signal and unobserved regions of the spectrum. The proxy spectra used are listed in the dataset extensions of the SED. New SEDs for several proxy stars (kappa 1 Ceti, tau Ceti, eps Indi and 70 Oph B) are included as separate HLSPs. The 70 Oph B data products do not include a full SED at time of release as insufficient data is available.
 


# References
The following papers provide details on various aspects of the MUSCLES treasury dataset and its creation. When using the MUSCLES spectra (v22), please cite all or a selection of these works as appropriate to your application, along with the version number of the spectra used.

- [France+ 2016; ApJ 820:89F](https://ui.adsabs.harvard.edu/#abs/2016ApJ...820...89F)
- [Youngblood+ 2016; ApJ 824:101](https://ui.adsabs.harvard.edu/#abs/2016ApJ...824..101Y)
- [Loyd+ 2016; ApJ 824:102](https://ui.adsabs.harvard.edu/#abs/2016ApJ...824..102L)

When using Mega-MUSCLES targets (v23), please cite.
- [Wilson+ 2021; ApJ 911:18W](https://ui.adsabs.harvard.edu/abs/2021ApJ...911...18W/abstract)

When using MUSCLES Extension targets (v24), please cite.
- [Behr+ 2016; AJ 166:35](https://ui.adsabs.harvard.edu/abs/2023AJ....166...35B/abstract)

When using MUSCLES Extension targets (v25), please cite.
- [Wilson et al. 2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...978...85W/abstract)

When using the semi-empirical models, please cite.

- [Fontenla+ 2016; ApJ 830:154F](https://ui.adsabs.harvard.edu/abs/2016ApJ...830..154F/abstract)
 
# Version History

Version 0.0 (v00) 2015/11/13
- Initial internal release

Version 0.1 (v01) 2015/01/04
- corrected a mistake in the normalization of GJ 176 STIS G230L data (~10% change)
- Lya profile for GJ 1214 revised to match a new Lya - Mg II empirical relationship (see Youngblood+ 2016 for details)
- improved artifact correction in the X-ray spectra
- added metadata to headers of CXO component spectra

Version 0.2 (v02) 2015/01/07
- corrected a gap at ~50-100 Angstroms between the APEC model spectrum and EUV estimates accidentally introduced in version 0.1

Version 1.0 (v10) 2015/01/28
- **official public release**
- fixed error in the comment in header of HD 97658 X-ray component spectrum that said the spectrum was scaled from HD 85512 according to the two stars' distances when really it was their bolometric flux
- adaptive resolution spectra added 2015/04


Version 2.0 (v20) 2016/10/28
- Revised noise estimation for photometry. This did not directly affect the normalization values computed, but indirectly did so by potentially affecting which measurements were identified as outliers.
- Improved the filtering of photometry gathered by the Vizier photometry viewer. The photometry table this service returns includes duplicate entries, and the previous pipeline version was overly-conservative in rejecting potential duplicates. The new version does a better job of identifying and retaining unique measurements. This and the above change resulted in changes to the normalization of the PHOENIX and G430L spectra (and thus of the bolometric flux) of <1% in all cases except GJ667C, for which the normalization decreased by 10.1%.
- Corrected error in the propagation of exposure times through the coaddition of Echelle spectrum orders.
- Revised the coaddition of Echelle spectra from coadding all orders from all exposures using 1/error^2 weighting to coadding orders from the same exposure with 1/error^2 wighting then coadding separate exposures using exposure time weighting.
- Moved to double precision for all FITS columns.
- Corrected for the accidental use of an earlier draft of the EUV reconstruction in the SEDs, changing EUV levels in the broadband EUV bins by 0.1-3% (except for GJ1214, where the change was 13-20%).
- Revised optimal splicing algorithm to account for serious data quality flags when evaluating which data to choose within spectral overlaps when assembling stitched SEDs.
- Added Proxima Centauri to the library of spectra.


Version 2.1 (v21) 2017/01/25
- Changed the method for splicing the reconstructed Lya flux into the data. The previous splice range was overly extensive (10 Å), eliminating the FUV continuum and broad extended wings of the Lya line (see Fig 7 of Youngblood+ 2016). The new method retains the COS data up to ±400 km/s of the Lya rest wavelength. A fit to the Lorentzian wings of the Lya profile, detected in the COS data but below the noise level of the STIS data used for the Lya reconstruction, is used to fill the gap between ±400 km/s and the region where the reconstructed Lya profile begins to dominate. The Lya reconstructions do not themselves include Lorentzian wings because the STIS data could not constrain the amplitude of these wings. This change results in increases of order 1% to the integrated Lya flux, on par with expectations based on the Youngblood+ 2016 work. The fit parameters are saved in header keywords (described later in this readme).
- Fixed bug wherein the 'BOLOFLUX' keyword was not being saved in the SED headers.


Version 2.2 (v22) 2018/01/11
- Tweaked the "adaptive binning" spectra: (1) Major emission lines are now excluded from adaptive rebinning and (2) bins that have grown in width to larger than 1.0 AA are resampled to 1.0 AA bins.


Version 2.3 (v23) 2022/03/17
- First release of SEDs for the Mega-MUSCLES targets.


Version 2.4 (v24) 2023/09/06
- First release of SEDs for the MUSCLES Extension targets.


Version 2.5 (v25) 2025/06/XX
- First release of SEDs for the Mega-MUSCLES Extension targets.