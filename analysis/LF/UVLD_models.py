


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.lines as mlines
import cmasher as cmr

import sys
import os
from termcolor import colored

from flare.photom import  M_to_lum, lum_to_M
import flare.plt as fplt
from flare.LF import models, evo, plots
from flare.LF.literature import UV, SFR



print(colored('test', 'red'))


fig, ax = fplt.simple() # UV luminosity density


M_limits = [-18, -19, -20]
colors = cmr.take_cmap_colors('cmr.flamingo', len(M_limits), cmap_range=(0.15, 0.85))

offsets = np.linspace(-0.1, 0.1, len(M_limits))




for M_limit, color, offset in zip(M_limits, colors, offsets):

    L_limit = M_to_lum(M_limit)
    mod_handles = []

    for model, ls in zip(['FLARES_binned','Bluetides','Ma2019','Mason2015'], ['-','-.','--',':']):

        m = getattr(UV, model)()

        mod_handles.append(mlines.Line2D([], [], color='0.5', lw=1, linestyle=ls, label=rf'$\rm {m.name} $'))

        rhos = []


        for z in m.redshifts:
            rho = m.density(z, L_limit)
            print(z, rho, np.log10(rho))
            rhos.append(np.log10(rho))

        ax.plot(m.redshifts, rhos, ls = ls, lw=1, c=color)




# Create a legend for the first line.
mod_legend = plt.legend(handles=mod_handles, loc='lower left', fontsize=8)

# Add the legend manually to the current Axes.
plt.gca().add_artist(mod_legend)



# ax.legend(fontsize = 8, title = r'$\rm \rho_{FUV}\ (L_{FUV}>10^{28}\ erg\ s^{-1}\ Hz^{-1}$')

ax.set_xlim([6.5, 15])
ax.set_ylim([23, 26.1])

ax.legend(fontsize = 8)
ax.set_xlabel(r'$\rm z $')

ax.set_ylabel(r'$\rm log_{10}(\rho_{FUV}/erg\ s^{-1}\ Hz^{-1}\ Mpc^{-3}) $')

fig.savefig('figs/UVLD_models.pdf')
