
import sys
import numpy as np
import pickle

import matplotlib.pyplot as plt
from matplotlib import cm
import cmasher as cmr



import scipy.stats

import flare

from flare.LF import models, evo, completeness
from flare.LF.literature import UV

from flare.photom import m_to_flux
import flare.plt as fplt
import flare.stats
from flare.surveys import jwst


z_r = [8.0, 14.4]



fig, ax = fplt.simple()


# ax.fill_between([8,9],[0]*2, [10000]*2, color='k', alpha=0.05)



line_styles = [':','-.','--','-','-','-']


redshift_limits = [7, 15]
log10L_limits = [27.5, 31.]
dz = 0.1
dlog10L = 0.1

m = getattr(UV, 'FLARES_binned')()
bin_edges, bin_centres, volume, N_tot = m.N(redshift_limits = redshift_limits, log10L_limits = log10L_limits, dz = dz, dlog10L = dlog10L, return_volumes = True)





f = 'Webb.NIRCam.F277W'


surveys = [jwst.PEARLS, jwst.CEERS, jwst.COSMOS_Web, jwst.NGDEEP, jwst.PANORAMIC, jwst.PRIMER, jwst.JADES, jwst.Cy1][::-1]

lws = [2.0]+[1.0]*7

colors = cmr.take_cmap_colors('cmr.lavender', len(surveys))

lss = ['-','--','-.',':']*5

for i, (Survey, color, lw, ls) in enumerate(zip(surveys, colors, lws, lss)):

    print(Survey.name)

    N = np.zeros((*N_tot.shape, len(Survey.fields_)))

    for k, field in enumerate(Survey.fields_):

        flux_limit = field['depths'][f]*(10/5)  # <--- assumes we can detect anything at >10\sigma
        area = field['area'] # arcmin2

        # --- simple completeness
        completeness = flare.LF.completeness.completeness_cut(bin_centres, flux_limit, cosmo = flare.default_cosmo())

        # --- get expected number of galaxies
        N_ = np.multiply(N_tot, completeness) * area
        N[:,:,k] = N_ # apply completeness correction

    n = np.sum(np.sum(N, axis=2), axis=0)[::-1]

    cn = np.cumsum(n)
    s = Survey.name.replace(' ',r'\ ')
    ax.plot(bin_centres['z'][::-1], cn, c=color, alpha = 0.7, lw=lw, ls = ls, label = rf"$\rm {s}$")

    for z in [10.,12.,14]:
        print(z, np.interp(z, bin_centres['z'], cn[::-1]))





#
ax.legend(loc= 'lower left', fontsize=8, labelspacing = 0.05)
#
#



ax.set_ylim([0.11, 5000])
ax.set_xlim(z_r)
ax.set_yscale('log')
ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm N(>z,z<15)$')

# from matplotlib.ticker import ScalarFormatter
# ax.yaxis.set_major_formatter(ScalarFormatter())

fig.savefig(f'figs/CN.pdf')
