import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import glob

stars = os.listdir('draft_hlsp')

for star in stars:
    seds = glob.glob('draft_hlsp/{}/*multi*'.format(star))
    if len(seds) > 0:
        print(star)
        names=fits.getdata(seds[0], 1).names


        for sed in seds:

            plt.figure(figsize=(14, 50))
            filename = os.path.split(sed)[1]
        #     if '{}_check.pdf'.format(filename[:-5]) not in done:
            data = fits.getdata(sed, 1)
            for i, name in enumerate(names[1:]):
                plt.subplot(12, 1, i+1)
                if i == 0:
                     plt.title(filename.replace('_', ' '))
        #              plt.title(filename)
                plt.step(data['WAVELENGTH'], data[name], where='mid')
                plt.xscale('log')
                if name not in ['EXPTIME', 'EXPSTART', 'EXPEND']:
                    plt.yscale('log')
                if name == 'FLUX':
                    plt.ylim(min(data[name]))
        #         if i == len(names)-2:
                plt.xlabel('WAVELENGTH (\AA)')
                plt.ylabel(name)
                plt.xlim(5, 1e7)

            plt.tight_layout()
            plt.savefig('plots/diagnosis_plots/{}_check.pdf'.format(filename[:-5]))
    #         plt.show()
            plt.close()