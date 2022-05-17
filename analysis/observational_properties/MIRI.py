
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines

import scipy.stats as stats
from scipy.stats import binned_statistic

import cmasher as cmr

import h5py

from flare.utils import bin_centres
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


for t in ['Pure_Stellar','Intrinsic', 'DustModelI']:
    for filter in ['NIRCAM/F444W', 'MIRI/F770W']:
        inst, f = filter.split('/')
        quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Flux/{t}/JWST/{inst}', 'dataset': f, 'name': f'{f}_{t}', 'log10': True})




left  = 0.15
width = 0.7
height = 0.8
bottom = 0.15

fig = plt.figure(figsize = (3.5, 3.0))

ax = fig.add_axes((left, bottom, width, height))


bin_edges = np.arange(27.5, 30., 0.25)

for z,c in zip(flares.zeds, cmr.take_cmap_colors('cmr.gem_r', len(flares.zeds))):

    # --- get quantities (and weights and deltas)
    D = flares.get_datasets(flares.tag_from_zed[z], quantities)

    s = D['log10FUV']>28.0

    for t, ls in zip(['Pure_Stellar', 'DustModelI'],[':','-']):

        Y = D[f'log10F770W_{t}']-D[f'log10F444W_{t}']

        Y, bin_edges, _ = binned_statistic(D['log10FUV'], Y, statistic='median', bins = bin_edges)

        ax.plot(bin_centres(bin_edges), Y, c=c, lw=1, ls=ls)




# handles = [mlines.Line2D([], [], color='0.5', ls=ls, lw=1, label=label) for ls, label in zip([':','--','-'],[r'$\rm stellar$',r'$\rm stellar\ +\ nebular$', r'$\rm stellar\ +\ nebular\ +\ dust$'])]

handles = [mlines.Line2D([], [], color='0.5', ls=ls, lw=1, label=label) for ls, label in zip([':','-'],[r'$\rm stellar$', r'$\rm stellar\ +\ nebular\ +\ dust$'])]

legend1 = plt.legend(loc ='upper right', handles=handles, fontsize = 7, labelspacing = 0.1)

handles = [mlines.Line2D([], [], color=c, ls='-', lw=1, label=rf'$\rm z={z:.0f}$') for z, c in zip(flares.zeds, cmr.take_cmap_colors('cmr.gem_r', len(flares.zeds)))]

legend2 = plt.legend(loc ='upper left', handles=handles, fontsize = 7, labelspacing = 0.1)

ax.add_artist(legend1)
ax.add_artist(legend2)


# ax.xaxis.set_major_locator(plt.MaxNLocator(6))
ax.set_xticks(np.arange(28, 31, 1.0))

ylim = np.array([-0.1, 0.45])
ax.set_ylim(ylim)

ax.set_xlabel(r'$\rm log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(r'$\rm log_{10}(f_{F770W}/f_{F444W})$')

ax2 = ax.twinx()
ax2.set_ylim(2.5*ylim)
ax2.set_ylabel(r'$\rm F444W-F770W$')


fig.savefig(f'figs/miri.pdf')
