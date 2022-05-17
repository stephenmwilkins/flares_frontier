
import sys
import numpy as np
import pickle

import matplotlib.pyplot as plt
from matplotlib import cm

import cmasher as cmr

from astropy.io import ascii

import flare

from flare.LF import models, evo, completeness
from flare.LF.literature import UV

from flare.photom import m_to_flux
import flare.plt as fplt
import flare.stats
from flare.surveys import jwst



survey = jwst.Cy1
f = 'Webb.NIRCam.F277W'

redshifts = [10,11,12,13,14,15][::-1]
# redshifts = [15]
X_limits = [27.9, 29.6]
Y_limits = [ -8.5,-1.01]


N = len(redshifts)
left = 0.075
top = 0.9
bottom = 0.2
right = 0.975
panel_width = (right-left)/N
panel_height = top-bottom
fig, axes = plt.subplots(1, N, figsize = (7,7/(panel_height/panel_width)), sharey = True, sharex = True)
plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)



for z, c, ax in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts)), axes):

    ax.text(0.5, 1.02, rf'$\rm z={z}$', horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize = 7)

    print(z)

    z_r = [z-0.5,z+0.5]

    # --- plot binned version

    m = UV.FLARES_binned() # binned

    ax.step(m.log10L, np.log10(m.phi[z]), where = 'mid', label = 'original binned', lw=2, c=c, alpha = 0.2) # raw data


    # --- plot Schechter version with redshift range

    binw = 0.01
    bin_edges = np.arange(27, 30.5, binw)
    bin_centres = models.bin_centres(bin_edges)

    m = UV.FLARES() # Schechter

    phi = m.phi_binned(z, bin_edges)/binw
    phi_low = m.phi_binned(z_r[0], bin_edges)/binw
    phi_high = m.phi_binned(z_r[1], bin_edges)/binw

    # ax.fill_between(bin_centres, np.log10(phi_low), np.log10(phi_high), alpha = 0.2, color='0.7')
    # ax.plot(bin_centres, np.log10(phi))
    ax.plot(bin_centres, np.log10(phi), lw=1, alpha=0.5, c=c)



    # --------------------------------------------------------
    # --- make LFE model

    m = UV.FLARES()

    log10L_limits = [27.05, 31.]
    dz = 0.1
    dlog10L = 0.1

    bin_edges, bin_centres, volume, N_tot = m.N(redshift_limits = z_r, log10L_limits = log10L_limits, dz = dz, dlog10L = dlog10L, return_volumes = True)

    N_fields = len(survey.fields_)

    N = np.zeros((*N_tot.shape, N_fields))
    N_CV = np.zeros((*N_tot.shape, N_fields))
    V = np.zeros((*N_tot.shape, N_fields))

    total_area = 0.0

    for k,field in enumerate(survey.fields_):

        flux_limit = field['depths'][f]*(10/5)  # <--- assumes we can detect anything at >10\sigma
        area = field['area'] # arcmin2
        area_sd = area / 3600.  # Area in square degrees
        area_sr = (np.pi / 180.) ** 2 * area_sd  # Area in steradian
        total_area += area

        # --- simple completeness
        completeness = flare.LF.completeness.completeness_cut(bin_centres, flux_limit, cosmo = flare.default_cosmo())

        # --- get expected number of galaxies
        N_ = np.multiply(N_tot, completeness) * area
        N[:,:, k] = N_ # apply completeness correction
        N_CV[:,:, k] = N_

        # --- get volume
        V_ = np.multiply(np.reshape(np.repeat(volume, 40), (40,10)), completeness) * area_sr
        V[:,:, k] = V_



    N = np.sum(N, axis=2)
    N_CV = np.sum(N_CV, axis=2)
    V = np.sum(V, axis=2)

    # --- select redshifts of interest
    #


    n = np.sum(N, axis=1)
    v = np.sum(V, axis=1)
    n_cv = np.sum(N_CV, axis=1)



    b_c = bin_centres['log10L']


    # --- re-bin
    n = n[0::2] + n[1::2]
    n_cv = n_cv[0::2] + n_cv[1::2]
    v = np.mean(np.array([v[0::2], v[1::2]]), axis=0)
    b_c = (bin_centres['log10L'][0::2] + bin_centres['log10L'][1::2])/2
    dlog10L = 0.2

    # -- rebin again
    if np.sum(n)<50:
        # --- re-bin
        n = n[0::2] + n[1::2]
        n_cv = n_cv[0::2] + n_cv[1::2]
        v = np.mean(np.array([v[0::2], v[1::2]]), axis=0)
        b_c = (b_c[0::2] + b_c[1::2])/2
        dlog10L = 0.4


    n_up = np.array([flare.stats.poisson_confidence_interval(n_gal)[1] for n_gal in n_cv])
    n_down = np.array([flare.stats.poisson_confidence_interval(n_gal)[0] for n_gal in n_cv])

    phi = n_cv/(dlog10L*v)
    phi_upper = n_up/(dlog10L*v)
    phi_lower = n_down/(dlog10L*v)


    c = 'k'
    s = n>0.5
    ax.scatter(b_c[s], np.log10(phi[s]), c=c, s=10, alpha = 1)
    for b,u,l in zip(b_c[s], phi_upper[s], phi_lower[s]):
        ax.plot([b]*2, [np.log10(u),np.log10(l)], c=c, lw=1, alpha = 1)

    # ax.fill_between(b_c, np.log10(phi_upper),np.log10(phi_lower), color=c, alpha = 0.5)

    ax.set_xlim(X_limits)
    ax.set_ylim(Y_limits)

    # ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    # ax.set_xticks(np.arange(28, 31, 1.0))


axes[0].set_ylabel(r'$\rm\log_{10}[\phi/Mpc^{-3}\ dex^{-1}]$', fontsize = 9)
fig.text(left+(right-left)/2, 0.04, r'$\rm \log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$', ha='center', fontsize = 9)

fig.savefig(f'figs/LF_forecast.pdf')
