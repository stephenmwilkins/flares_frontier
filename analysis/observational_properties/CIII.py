
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


import flare.obs.literature.CIII as CIII_observations

x_limits = [28.0, 29.6]

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/CIII1907', 'dataset': f'EW', 'name': f'CIII1907_EW', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/CIII1909', 'dataset': f'EW', 'name': f'CIII1909_EW', 'log10': False})

D = {}
s = {}
s = {}

for tag, z in zip(flares.tags, flares.zeds):

    # --- get quantities (and weights and deltas)
    D[z] = flares.get_datasets(tag, quantities)
    D[z]['CIIIEW'] = D[z]['CIII1907_EW']+D[z]['CIII1909_EW']

    s[z] = D[z]['log10FUV']>x_limits[0]



limits = flares_utility.limits.limits
limits['log10FUV'] = x_limits
limits['CIIIEW'] = [0., 40.]




# fig, axes = flares_utility.plt.linear_redshift(D, flares.zeds, 'log10FUV', 'CIII_EW', s, limits = limits, scatter_colour_quantity = 'CIII_EW', scatter_cmap = cmr.guppy_r, bins = 10, add_zevo = True)

# fig, ax = flares_utility.plt.zevo(D, flares.zeds, 'log10FUV', 'CIIIEW', s, limits = limits, bins = 10)

hist_bins = np.arange(0, 30, 2)

fig, ax, hax = flares_utility.plt.zevo_whist(D, flares.zeds, 'log10FUV', 'CIIIEW', s, limits = limits, bins = 10, hist_bins = hist_bins)



o = CIII_observations.Llerena2021_B

print(o.MFUV)

ax.errorbar(o.log10LFUV, o.CIII_EW, xerr=o.log10LFUV_binw/2, yerr=o.CIII_EW_err, c='k', label = rf'$\rm {o.label}$', fmt='o',elinewidth=1, ms=5, capsize=2)

o = CIII_observations.Topping2021

ax.errorbar(o.log10LFUV, o.CIII_EW, yerr=o.CIII_EW_err, uplims = o.CIII_EW_uplims, c='0.5', label = rf'$\rm {o.label}$',fmt='o',elinewidth=1, ms=3, capsize=2)

o = CIII_observations.Stark2017

ax.errorbar(o.log10LFUV, o.CIII_EW, yerr=o.CIII_EW_err, uplims = o.CIII_EW_uplims, c='0.5', label = rf'$\rm {o.label}$',fmt='d',elinewidth=1, ms=3, capsize=2)


o = CIII_observations.Stark2015

ax.errorbar(o.log10LFUV, o.CIII_EW, yerr=o.CIII_EW_err, uplims = o.CIII_EW_uplims, c='0.5', label = rf'$\rm {o.label}$',fmt='s',elinewidth=1, ms=3, capsize=2)

ax.legend(fontsize=7, labelspacing = 0.4)




ax.set_xticks(np.arange(28., 29.5, 0.5))


# hax.hist(Maseda2017.CIII_EW, bins=hist_bins, orientation='horizontal', color='k', alpha=0.2, density=True, zorder=0, label = r'$\rm Maseda+17$')
# hax.legend(fontsize = 6)
hax.axis('off')

# --- add observational comparisons

# for ax, z in zip(axes, flares.zeds):
#     print(z)
#     if z in beta_observations.observed.keys():
#         print(z)
#         for obs in beta_observations.observed[z]:
#             if obs.dt == 'io':
#                 ax.errorbar(obs.log10Luv, obs.beta, xerr = obs.log10Luv_err, yerr = obs.beta_err, c='k',fmt='o',elinewidth=1, ms=3, label = rf'$\rm {obs.label}$')
#
#     ax.legend(fontsize=7, handletextpad = 0.0, loc = 'upper left')
#



fig.savefig(f'figs/CIII.pdf')
