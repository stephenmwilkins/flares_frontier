
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
quantities.append({'path': f'Galaxy/SFR_aperture/30', 'dataset': f'50Myr', 'name': f'SFR_50', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Pure_Stellar', 'dataset': 'FUV', 'name': 'IntrinsicFUV', 'log10': True})
quantities.append({'path': f'Galaxy/StellarAges', 'dataset': 'MassWeightedStellarAge', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/Metallicity', 'dataset': 'MassWeightedGasZ', 'name': None, 'log10': True})



D = {}
s = {}
s = {}

for tag, z in zip(flares.tags, flares.zeds):

    # --- get quantities (and weights and deltas)
    D[z] = flares.get_datasets(tag, quantities)
    D[z]['log10sSFR'] = D[z]['log10SFR_50'] - D[z]['log10Mstar_30'] + 9.0
    D[z]['AFUV'] = 2.5*(D[z]['log10IntrinsicFUV'] - D[z]['log10FUV'])

    s[z] = D[z]['log10Mstar_30']>x_limits[0]

    print(z, np.median(D[z]['AFUV']))


x = 'log10Mstar_30'


properties = ['log10sSFR','log10MassWeightedStellarAge', 'log10MassWeightedGasZ','AFUV']

limits = flares_utility.limits.limits

limits['log10sSFR'] = [0.51, 1.59]
limits['log10MassWeightedStellarAge'] = [1.21, 2.19]
limits['log10MassWeightedGasZ'] = [-3.99, -1.51]
limits['log10FUV'] = [27.5, 29.5]

print(limits)

limits[x] = x_limits


cmap = cmr.get_sub_cmap('cmr.sapphire', 0.4, 1.0)

fig, axes = flares_utility.plt.linear_redshift_mcol(D, flares.zeds, x, properties, s, limits = limits, scatter_colour_quantity = 'log10FUV', scatter_cmap = cmap, add_linear_fit = False, bins = 10, add_zevo = True, fontsize = 8)



fig.savefig(f'figs/physical_properties.pdf')
