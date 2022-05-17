


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.lines as mlines
import cmasher as cmr

import sys
import os


from flare.photom import  M_to_lum, lum_to_M
import flare.plt as fplt
from flare.LF import models, evo, plots
from flare.LF.literature import UV, SFR






fig, ax = fplt.simple() # UV luminosity density


M_limits = [-18, -19, -20]
colors = cmr.take_cmap_colors('cmr.flamingo', len(M_limits), cmap_range=(0.15, 0.85))

offsets = np.linspace(-0.1, 0.1, len(M_limits))




for M_limit, color, offset in zip(M_limits, colors, offsets):


    # for model, ls in zip(['FLARES_binned','Bluetides','Ma2018'], ['-','-.','--']):

    L_limit = M_to_lum(M_limit)

    m = getattr(UV, 'FLARES_binned')()

    rhos = []


    for z in m.redshifts:
        rho = m.density(z, L_limit)
        print(z, rho, np.log10(rho))
        rhos.append(np.log10(rho))

    ax.plot(m.redshifts, rhos, label = rf'$\rm M_{{FUV}}<{M_limit}$', ls = '-', lw=1, c=color)



    # --- observations

    handles = []

    for model, ms in zip(['Bouwens2021','Mcleod16'], ['o', 'd']):

        m = getattr(UV, model)()

        handles.append(mlines.Line2D([], [], color='0.5', marker=ms, linestyle='None', markersize=3, label=rf'$\rm {m.name} $'))

        for z in m.redshifts:

            rho = m.density(z, L_limit)
            ax.plot([z+offset]*2, m.density_range(z, L_limit), c=color, lw=1)

            # --- sample the observed uncertainties

            ax.scatter(z+offset, np.log10(np.array(rho)), s = 10, c=[color], marker=ms)



# Create a legend for the first line.
obs_legend = plt.legend(handles=handles, loc='lower left', fontsize=8)

# Add the legend manually to the current Axes.
plt.gca().add_artist(obs_legend)




# ax.legend(fontsize = 8, title = r'$\rm \rho_{FUV}\ (L_{FUV}>10^{28}\ erg\ s^{-1}\ Hz^{-1}$')

ax.set_xlim([6.5, 15])
ax.set_ylim([23, 26.1])

ax.legend(fontsize = 8)
ax.set_xlabel(r'$\rm z $')

ax.set_ylabel(r'$\rm log_{10}(\rho_{FUV}/erg\ s^{-1}\ Hz^{-1}\ Mpc^{-3}) $')

fig.savefig('figs/UVLD.pdf')
