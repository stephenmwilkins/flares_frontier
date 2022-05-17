

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


filters = flare.filters.NIRCam_W + flare.filters.MIRI


colors = cmr.take_cmap_colors('cmr.neon', len(filters))

for z in [10, 15]:

    SED.total.get_fnu(cosmo, z) # calculate observer frame wavelength
    SED.stellar.get_fnu(cosmo, z) # calculate observer frame wavelength



    l, f = SED.total.lamz, SED.total.fnu
    l, f = rebin(SED.total.lamz, SED.total.fnu, 10)

    ax.plot(np.log10(l)-4, np.log10(f), zorder = 1, lw=1, c='0.8') # plot SED
    ax.text(np.log10(0.1216*(1+z))-0.03, 0.6, rf'$\rm z={z}$', c='0.6', fontsize=10, ha='center', rotation=90) # --- plot filter transmission curve

    # ax.plot([0.0, np.log10(0.35*(1+z))], [np.interp(1500., SED.total.lam, np.log10(SED.total.fnu))]*2, c='k', lw=3, alpha=0.2)


    F = flare.filters.add_filters(filters, new_lam = SED.total.lamz)
    SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)


    # add beta
    beta = SED.total.return_beta()
    l = np.array([0.125, 0.3])
    fnu = np.interp(2000., SED.total.lam, SED.total.fnu)*(l/0.2)**(beta+2.0)
    ax.plot(np.log10(l)+np.log10(1+z), np.log10(fnu),c='b',alpha=0.3)

    # add BB
    l = np.array([0.35, 0.45])
    fnu = np.interp(l*1E4, SED.total.lam, SED.total.fnu)
    print(fnu)
    ax.plot(np.log10(l)+np.log10(1+z), np.log10(fnu),c='r',alpha=0.3)

    for f, c, os in zip(F['filters'], colors, [0, 1]*10):

        pivwv = F[f].pivwv()/10000

        if z==15:
            ax.scatter(np.log10(pivwv), np.log10(SED.total.Fnu[f]), s=10, c=[c], zorder = 2)
        else:
            ax.scatter(np.log10(pivwv), np.log10(SED.total.Fnu[f]), edgecolor='0.3', s=13, lw=1, c=[c], zorder = 2)




F = flare.filters.add_filters(filters) # --- NOTE: need to give it the redshifted

for f, c, os in zip(F['filters'], colors, [0, 1]*10):

    pivwv = F[f].pivwv()/10000

    if pivwv > 1 and pivwv<17:

        filt_ax.plot(np.log10(F[f].lam)-4, F[f].T, lw = 1, zorder = 1, c=c) # --- plot filter transmission curve

        label = f.split('.')[-1]

        filt_ax.text(np.log10(F[f].pivwv())-4, 0.55+os*0.1, label, c=c, fontsize=4.5, ha='center') # --- plot filter transmission curve



filters = ['Hubble.WFC3.f125w', 'Hubble.WFC3.f160w', 'Spitzer.IRAC.ch1', 'Spitzer.IRAC.ch2']
F = flare.filters.add_filters(filters) # --- NOTE: need to give it the redshifted


for f in F['filters']:

    filt_ax.plot(np.log10(F[f].lam)-4, F[f].T, lw = 2, c='k', alpha=0.2, zorder=0) # --- plot filter transmission curve

    label = f.split('.')[-1]

    filt_ax.text(np.log10(F[f].pivwv())-4, 0.1, label, c='0.5', fontsize=4.5, ha='center') # --- plot filter transmission curve




# hax.set_xlim(xlim)
# hax.set_xticks([])
# # hax.set_yticks([])
# hax.set_ylabel(r'$\rm N_{\rm sim}$')

xlim = [0.0, 1.25]

ax.set_xticks([])
ax.set_xlim(xlim)
ax.set_ylim([0.51, 1.49])

filt_ax.set_xlim(xlim)
filt_ax.set_ylim([0.,0.75])

# ax.legend(fontsize=8)

filt_ax.set_xlabel(r'$\rm \log_{10}(\lambda/\mu m) $')
ax.set_ylabel(r'$\rm \log_{10}(f_{\nu}/nJy\cdot 10^{8}\ M_{\odot})$')

fig.savefig(f'figs/sed.pdf')

fig.clf()




# SED.total.get_Fnu(F) # generates Fnu (broad band fluxes)
# for f in filters: plt.scatter(F[f].pivwv(), np.log10(SED.total.Fnu[f]), edgecolor = 'k', zorder = 2, label = f)
#
# plt.xlim([5000.,50000.])
#
# mx = np.max(np.log10(SED.total.fnu))
# plt.ylim([mx-4., mx+0.3])
# plt.show()
