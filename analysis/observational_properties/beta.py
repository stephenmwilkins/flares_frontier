
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

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'DustModelI', 'name': f'beta', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Pure_Stellar', 'dataset': 'FUV', 'name': 'PSFUV', 'log10': True})



D = {}
s = {}
s = {}

for tag, z in zip(flares.tags, flares.zeds):

    # --- get quantities (and weights and deltas)
    D[z] = flares.get_datasets(tag, quantities)
    D[z]['AFUV'] = 2.5*(D[z]['log10PSFUV']-D[z]['log10FUV'])

    s[z] = D[z]['log10FUV']>x_limits[0]



limits = flares_utility.limits.limits
limits['log10FUV'] = x_limits
limits['log10Mstar_30'] = [7.5, 10.]

fig, axes = flares_utility.plt.linear_redshift(D, flares.zeds, 'log10FUV', 'beta', s, limits = limits, scatter_colour_quantity = 'AFUV', scatter_cmap = cmr.guppy_r, bins = 20, rows = 2)



# --- add observational comparisons

for ax, z in zip(axes.flatten(), flares.zeds):
    print(z)
    if z in beta_observations.observed.keys():
        print(z)
        for obs in beta_observations.observed[z]:
            print(obs)
            if obs.dt == 'individual':
                print('here')
                ax.errorbar(obs.log10Luv, obs.beta, xerr = obs.log10Luv_err, yerr = obs.beta_err, c='k',fmt='o',elinewidth=1, ms=3, label = rf'$\rm {obs.label}$')

    ax.legend(fontsize=7, handletextpad = 0.0, loc = 'upper left')




fig.savefig(f'figs/beta.pdf')
