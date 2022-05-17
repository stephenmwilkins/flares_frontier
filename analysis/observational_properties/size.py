
# import sys
# import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import cmasher as cmr
import h5py
from scipy.stats import binned_statistic

import flare.plt as fplt
from flare.utils import bin_centres
import flares_utility.analyse as analyse

fig, ax = fplt.simple()


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares-sizes-results.hdf5'
f = h5py.File(filename)


f.visit(print)

redshifts = [10,11,12]

bin_edges = np.arange(27.5, 30., 0.25)





for t, ls in zip(['Observed', 'Intrinsic'], ['-','--']):
    for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', 6)):

        tag = flares.tag_from_zed[z]


        clim = f[f'{tag}/FAKE.TH.FUV/{t}/Luminosity']>f[f'{tag}/FAKE.TH.FUV/{t}'].attrs['Complete_Luminosity']

        log10LFUV = np.log10(f[f'{tag}/FAKE.TH.FUV/{t}/Luminosity'][clim])

        if t == 'Observed': log10r = np.log10(f[f'{tag}/FAKE.TH.FUV/{t}/HLR_Pixel_0.5'][clim])
        if t == 'Intrinsic': log10r = np.log10(f[f'{tag}/FAKE.TH.FUV/{t}/HLR_0.5'][clim])

        # print(np.median(log10LFUV))
        # print(np.median(log10r))

        Y, _ , _ = binned_statistic(log10LFUV, log10r, statistic='median', bins = bin_edges)

        ax.plot(bin_centres(bin_edges), Y, c=c, lw=1, ls=ls)



handles = [mlines.Line2D([], [], color='0.5', ls=ls, lw=1, label=rf'$\rm {t}$') for t, ls in zip(['Observed', 'Intrinsic'], ['-','--'])]

legend1 = plt.legend(loc ='upper right', handles=handles, fontsize = 7, labelspacing = 0.1)

handles = [mlines.Line2D([], [], color=c, ls='-', lw=1, label=rf'$\rm z={z:.0f}$') for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(flares.zeds)))]

legend2 = plt.legend(loc ='upper left', handles=handles, fontsize = 7, labelspacing = 0.1)

ax.add_artist(legend1)
ax.add_artist(legend2)



ax.set_xticks(np.arange(28, 31, 1.0))
ax.set_ylim([-1.25, 0.24])

ax.set_xlabel(r'$\rm log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(r'$\rm log_{10}(r_{0.5}/pkpc)$')
fig.savefig(f'figs/size.pdf')
