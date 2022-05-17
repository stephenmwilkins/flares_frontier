

# Create a model SED


import numpy as np
import matplotlib.pyplot as plt

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import cmasher as cmr
import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh
import flare.filters

import flare.plt as fplt
from interrogator.sed.core import rebin



# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

# SPS = models.SPS('P2/ModSalpeter_100')
SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')


# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)

sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -3., 'log10M*': 8.})

# plt.imshow(sfzh)
# plt.show()


print('star formation rate: {0}'.format(sfr))

SED = SPS.get_Lnu(sfzh, {'fesc': 0.0}, dust = False)
# SED = SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': 0.0}, dust = ('simple', {'slope':-1}))


left  = 0.15
height = 0.65
bottom = 0.15
width = 0.8
hheight = 0.15

fig = plt.figure(figsize = (3.5, 3.5))

ax = fig.add_axes((left, hheight+bottom, width, height))
filt_ax = fig.add_axes([left, bottom, width, hheight])

cosmo = flare.default_cosmo()


for z in [10, 15]:

    SED.total.get_fnu(cosmo, z) # calculate observer frame wavelength
    SED.stellar.get_fnu(cosmo, z) # calculate observer frame wavelength

    l, f = SED.total.lamz, SED.total.fnu
    l, f = rebin(SED.total.lamz, SED.total.fnu, 10)

    ax.plot(l/10000., np.log10(f), zorder = 1, lw=1, c='0.8') # plot SED
    ax.text(0.1216*(1+z)*0.98, 0.6, rf'$\rm z={z}$', c='0.6', fontsize=10, ha='center', rotation=90) # --- plot filter transmission curve


filters = flare.filters.NIRCam_W + flare.filters.MIRI
F = flare.filters.add_filters(filters) # --- NOTE: need to give it the redshifted


colors = cmr.take_cmap_colors('cmr.neon', len(F['filters']))

for f, c, os in zip(F['filters'], colors, [0, 1]*10):

    pivwv = F[f].pivwv()/10000

    if pivwv > 1 and pivwv<17:

        filt_ax.plot(F[f].lam/1E4, F[f].T, lw = 1, zorder = 1, c=c) # --- plot filter transmission curve

        label = f.split('.')[-1]

        filt_ax.text(pivwv, 0.55+os*0.1, label, c=c, fontsize=4.5, ha='center') # --- plot filter transmission curve




filters = ['Hubble.WFC3.f160w']
F = flare.filters.add_filters(filters) # --- NOTE: need to give it the redshifted


for f in F['filters']:

    filt_ax.plot(F[f].lam/1E4, F[f].T, lw = 2, zorder = 1, c=k, alpha=0.2, zorder=0) # --- plot filter transmission curve

    label = f.split('.')[-1]

    filt_ax.text(pivwv, 0.1, label, c='0.5', fontsize=4.5, ha='center') # --- plot filter transmission curve





# hax.set_xlim(xlim)
# hax.set_xticks([])
# # hax.set_yticks([])
# hax.set_ylabel(r'$\rm N_{\rm sim}$')

xlim = [0.9, 15.]

ax.set_xticks([])
ax.set_xlim(xlim)
ax.set_ylim([0.51, 1.49])

filt_ax.set_xlim(xlim)
filt_ax.set_ylim([0.,0.75])

# ax.legend(fontsize=8)

filt_ax.set_xlabel(r'$\rm \lambda/\mu m$')
ax.set_ylabel(r'$\rm \log_{10}(f_{\nu}/nJy\cdot 10^{8}\ M_{\odot})$')

fig.savefig(f'figs/sed_linear.pdf')

fig.clf()




# SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
# for f in filters: plt.scatter(F[f].pivwv(), np.log10(SED.total.Fnu[f]), edgecolor = 'k', zorder = 2, label = f)
#
# plt.xlim([5000.,50000.])
#
# mx = np.max(np.log10(SED.total.fnu))
# plt.ylim([mx-4., mx+0.3])
# plt.show()
