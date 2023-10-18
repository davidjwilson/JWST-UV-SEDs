"""
Makes an old looking SED
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

from functools import partial
from random import gauss, randrange
from PIL import Image, ImageFilter
plt.rcParams["font.family"] = "sans-serif"



def add_noise_to_image(im, perc=20, blur=0.5):
    gaussian = partial(gauss, 0.50, 0.02)
    width, height = im.size
    for _ in range(width*height * perc // 100):
        noise = int(gaussian() * 255)
        x, y = randrange(width), randrange(height)
        r, g, b, a = im.getpixel((x, y))
        im.putpixel((x, y),
                    (min(r+noise, 255), min(g+noise, 255), min(b+noise, 255)))

    im = im.filter(ImageFilter.GaussianBlur(blur))
    return im

star = 'kap1cet'

plt.style.use('/home/david/work/misc/old-style.mplstyle')
data = fits.getdata('draft_hlsp/{}/hlsp_muscles_multi_multi_{}_broadband_v1_adapt-const-res-sed.fits'.format(star, star), 1)
w, f, e= data['WAVELENGTH'], data['FLUX'], data['ERROR']


fig, ax = plt.subplots(figsize=(10, 5))
ax.step(w, f, where='mid', c='k')
ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlabel(' '.join('WAVELENGTH')+' [ $\mathrm{\AA}$ ]')
ax.set_ylabel('F L U X [ E R G S$^{-1}$ C M$^{-2}$ $\mathrm{\AA}^{-1}$ ]')
ax.set_title(' '.join('{}'.format(star.upper().replace('_',''))))

labely = max(f)*5
efac = 0.3
labelfac = 1.5

ax.errorbar(np.array([7,70]),np.array([labely, labely]), yerr= [[efac*labely,efac*labely],[0,0]], c ='k')
ax.annotate('XMM', (40, labelfac*labely),  ha='center')

ax.errorbar(np.array([72,1120]),np.array([labely, labely]), yerr= [[efac*labely,efac*labely],[0,0]], c ='k')
ax.annotate('DEM', (190, labelfac*labely))

ax.errorbar(np.array([1160, 2840]),np.array([labely, labely]), yerr= [[efac*labely,efac*labely],[0,0]], c ='k')
ax.annotate('HST/LYA', (1800, labelfac*labely), ha='center')

ax.errorbar(np.array([2950,  105000]),np.array([labely, labely]), yerr= [[efac*labely,efac*labely],[0,0]], c ='k')
ax.annotate('PHOENIX', (15000, labelfac*labely), ha='center')






ax.set_facecolor('#E8DCB8')

ax.set_ylim(min(f)*0.1, max(f)*100)

fig.tight_layout()


fig.savefig('plots/old_{}.png'.format(star), dpi=100, facecolor='#E8DCB8')






im = Image.open('plots/old_{}.png'.format(star))
im = add_noise_to_image(im, 5)
im.save('plots/old_{}_noise.png'.format(star))
im.show()