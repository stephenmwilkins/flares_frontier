
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


import flare.obs.literature.beta as beta_observations

x_limits = [27.5, 29.5]

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

print(flares.tags)
print(flares.zeds)

# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar', 'log10': True})
quantities.append({'path': 'Galaxy/SFR_aperture/30', 'dataset': f'50Myr', 'name': 'SFR', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Binggeli', 'dataset': f'DustModelI', 'name': f'BB', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Pure_Stellar', 'dataset': 'FUV', 'name': 'PSFUV', 'log10': True})



D = {}
s = {}
s = {}

for tag, z in zip(flares.tags, flares.zeds):

    # --- get quantities (and weights and deltas)
    D[z] = flares.get_datasets(tag, quantities)
    D[z]['log10sSFR'] = D[z]['log10SFR']-D[z]['log10Mstar']+9.0

    s[z] = D[z]['log10FUV']>x_limits[0]



limits = flares_utility.limits.limits
limits['log10FUV'] = x_limits
limits['log10BB'] = [-0.3, 0.5]

fig, axes = flares_utility.plt.linear_redshift(D, flares.zeds, 'log10FUV', 'log10BB', s, limits = limits, scatter_colour_quantity = 'log10sSFR', scatter_cmap = cmr.get_sub_cmap('cmr.pride', 0.15, 0.85), bins = 10, rows = 1)






fig.savefig(f'figs/balmer.pdf')
