
import numpy as np
import matplotlib.cm as cm

import scipy.stats as stats


import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.limits
import flares_utility.plt
import flares_utility.analyse as analyse
import flares_utility.stats

x_limits = [7.5, 9.99]

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

print(flares.tags)
print(flares.zeds)

# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Pure_Stellar', 'dataset': 'FUV', 'name': 'PSFUV', 'log10': True})


D = {}
s = {}
s = {}

for tag, z in zip(flares.tags, flares.zeds):

    # --- get quantities (and weights and deltas)
    D[z] = flares.get_datasets(tag, quantities)
    D[z]['AFUV'] = 2.5*(D[z]['log10PSFUV']-D[z]['log10FUV'])

    s[z] = D[z]['log10Mstar_30']>x_limits[0]




x = 'log10Mstar_30'


y = 'log10FUV'

limits = flares_utility.limits.limits


limits[y] = [27.25, 29.5]
limits[x] = x_limits


# cmap = cmr.get_sub_cmap('cmr.sapphire', 0.4, 1.0)
cmap = cmr.guppy_r

fig, axes = flares_utility.plt.linear_redshift(D, flares.zeds, x, y, s, limits = limits, scatter_colour_quantity = 'AFUV', scatter_cmap = cmap, bins = 15, rows = 2, lowz=True, add_weighted_range = False)

for z, ax in zip(flares.zeds, axes):

    ax.set_yticks(np.arange(27.5, 29.5, 0.5))

    for m, ls in zip([30, 29, 28, 27],[':','-.','--','-']):

        log10L = np.log10(phot.m_to_lum(m, z))
        print(log10L)
        ax.axhline(log10L, lw=1, c='k', alpha=0.2, ls=ls)
        ax.text(9.75, log10L+0.05, rf'$\rm {m}$', fontsize = 7, c='k', alpha=0.2)


    wm1 = flares_utility.stats.weighted_quantile(D[z]['log10Mstar_30'][s[z]], [0.5], D[z]['weight'][s[z]])
    wm2 = flares_utility.stats.weighted_quantile(D[z]['log10FUV'][s[z]], [0.5], D[z]['weight'][s[z]])

    f = lambda x: (x-wm1) + wm2
    ax.plot([7.5, 10], f(np.array([7.5,10])), alpha=0.15, c='k', lw=3)

    # f = lambda x: (x-np.median(D[z]['log10Mstar_30'][s[z]])) + np.median(D[z]['log10FUV'][s[z]])
    # ax.plot([7.5, 10], f(np.array([7.5,10])), alpha=0.15, c='k', lw=3)






fig.savefig(f'figs/L.pdf')
