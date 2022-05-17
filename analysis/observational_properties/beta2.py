
import numpy as np

import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines

# import scipy.stats as stats
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

x_limits = [27.5, 29.5]

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)



# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'DustModelI', 'name': f'DustModelI_beta', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'Pure_Stellar', 'name': f'Pure_Stellar_beta', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'Intrinsic', 'name': f'Intrinsic_beta', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})


z = 10


fig, ax = fplt.simple()


bin_edges = np.arange(27.5, 30., 0.25)

for z,c in zip(flares.zeds, cmr.take_cmap_colors('cmr.gem_r', len(flares.zeds))):

    # --- get quantities (and weights and deltas)
    D = flares.get_datasets(flares.tag_from_zed[z], quantities)

    s = D['log10FUV']>28.5

    for t, ls in zip(['Pure_Stellar','Intrinsic','DustModelI'],[':','--','-']):

        Y, bin_edges, _ = binned_statistic(D['log10FUV'], D[f'{t}_beta'], statistic='median', bins = bin_edges)

        ax.plot(bin_centres(bin_edges), Y, c=c, lw=1, ls=ls)

        print(f'{t} {np.median(D[f"{t}_beta"][s]):.2f} {np.median(D[f"{t}_beta"][~s]):.2f}')


handles = [mlines.Line2D([], [], color='0.5', ls=ls, lw=1, label=label) for ls, label in zip([':','--','-'],[r'$\rm stellar$',r'$\rm stellar\ +\ nebular$', r'$\rm stellar\ +\ nebular\ +\ dust$'])]

handles += [mlines.Line2D([], [], color=c, ls='-', lw=1, label=rf'$\rm z={z:.0f}$') for z, c in zip(flares.zeds, cmr.take_cmap_colors('cmr.gem_r', len(flares.zeds)))]

ax.legend(loc ='upper left', handles=handles, fontsize = 7, labelspacing = 0.1)



# W16 = []
# W16.append({'id': 'GN-z10-1', 'M_FUV': -21.6, 'beta': [-2.1, 0.3]})
# W16.append({'id': 'GN-z10-2', 'M_FUV': -20.7, 'beta': [-2.5, 0.7]})
# W16.append({'id': 'GN-z10-3', 'M_FUV': -20.6, 'beta': [-2.3, 0.5]})
# W16.append({'id': 'GS-z10-1', 'M_FUV': -20.6, 'beta': [-1.9, 0.5]})
# W16.append({'id': 'MACS1149-JD1', 'M_FUV': -19.4, 'beta': [-2.1, 0.3]})
#
# X = [np.log10(M_to_lum(o['M_FUV'])) for o in W16]
# Y = [o['beta'][0] for o in W16]
# Yerr = [o['beta'][1] for o in W16]
#
# ax.errorbar(X, Y, yerr=Yerr, c='k',fmt='o',elinewidth=1, ms=3, label = r'$\rm W16$')
#
# # ax.errorbar(np.median(X), -2.1, yerr=0.3, c='r',fmt='o',elinewidth=1, ms=5, label = r'$\rm W16\ stack$')
#



# ax.xaxis.set_major_locator(plt.MaxNLocator(6))
ax.set_xticks(np.arange(28, 31, 1.0))
ax.set_ylim([-2.9, -1.51])

ax.set_xlabel(r'$\rm log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(r'$\rm\beta$')
fig.savefig(f'figs/beta2.pdf')
