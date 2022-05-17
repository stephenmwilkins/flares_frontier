
import sys
import numpy as np
import pickle

import matplotlib.pyplot as plt
from matplotlib import cm
import cmasher as cmr



import scipy.stats

import flare

from flare.LF import models, evo, literature, completeness

import flare.photom
import flare.plt as fplt
import flare.stats

import flags


z_r = [8.0, 14.4]



fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


# ax.fill_between([8,9],[0]*2, [10000]*2, color='k', alpha=0.05)



line_styles = [':','-.','--','-','-','-']


redshift_limits = [7, 15]
log10L_limits = [27.5, 31.]
dz = 0.1
dlog10L = 0.1

m = getattr(literature, 'FLARES_binned')()
bin_edges, bin_centres, volume, N_tot = m.N(redshift_limits = redshift_limits, log10L_limits = log10L_limits, dz = dz, dlog10L = dlog10L, return_volumes = True)








surveys = ['PEARLS', 'CEERS', 'COSMOS-Web', 'NGDEEP', 'PANORAMIC',  'PRIMER',  'JADES', 'Webb All Cycle 1'][::-1]

lws = [2.0]+[1.0]*7

colors = cmr.take_cmap_colors('cmr.lavender', len(surveys))


for i, (survey, color, lw) in enumerate(zip(surveys, colors, lws)):

    print(survey)

    Survey = flags.Surveys[survey]

    N = np.zeros((*N_tot.shape, len(Survey)))

    for k,sub in enumerate(Survey):

        flux_limit = flare.photom.m_to_flux(sub['depths_abmag']['Webb.NIRCam.F150W'])*(10/5)  # <--- assumes we can detect anything at >10\sigma
        area = sub['area'] # arcmin2

        # --- simple completeness
        completeness = flare.LF.completeness.completeness_cut(bin_centres, flux_limit, cosmo = flare.default_cosmo())

        # --- get expected number of galaxies
        N_ = np.multiply(N_tot, completeness) * area
        N[:,:,k] = N_ # apply completeness correction

    n = np.sum(np.sum(N, axis=2), axis=0)[::-1]

    cn = np.cumsum(n)
    s = survey.replace(' ',r'\ ')
    ax.plot(bin_centres['z'][::-1], cn, c=color, alpha = 0.7, lw=lw, ls = '-', label = rf"$\rm {s}$")






#
ax.legend(loc= 'lower left', fontsize=8, labelspacing = 0.05)
#
#



ax.set_ylim([0.11, 3000])
ax.set_xlim(z_r)
ax.set_yscale('log')
ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm N(>z,z<15)$')

# from matplotlib.ticker import ScalarFormatter
# ax.yaxis.set_major_formatter(ScalarFormatter())

fig.savefig(f'figs/surveys.pdf')
