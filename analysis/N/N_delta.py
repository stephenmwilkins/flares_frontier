
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot

import flares_utility.analyse as analyse










filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

# redshifts = [11,12,13,14,15]
# redshifts = [11,13,15]
redshifts = flares.zeds

print(redshifts)


V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3

xlim = [-0.3, 0.35]
ylim = [-0.2, 2.99]


# --- all redshifts on one figure


left  = 0.15
height = 0.5
bottom = 0.1
width = 0.7
hheight = 0.35


fig = plt.figure(figsize = (3.5, 5.5))

ax = fig.add_axes((left, bottom, width, height))
hax = fig.add_axes([left, bottom+height, width, hheight])
hax2 = hax.twinx()


deltas = flares.deltas
log10deltas = np.log10(deltas+1)

print(np.min(deltas), np.max(deltas))
print(np.min(log10deltas), np.max(log10deltas))

# hax.hist(log10deltas, bins = np.arange(-0.3, 0.35, 0.025), color = '0.5', histtype=u'step')
hax.hist(log10deltas, bins = np.arange(-0.3, 0.35, 0.05), color = 'k', alpha=0.3, log=False)



for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

    tag = flares.tag_from_zed[z]

    X = flares.load_dataset(tag, f'Galaxy', 'Mstar')


    N = []
    n = []
    for sim, delta in zip(flares.sims, flares.deltas):
        x = np.log10(np.array(X[sim])) +10.
        N.append(np.sum(x>7.5))
        n += [delta]*np.sum(x>7.5)

    N = np.array(N)
    N_avg = np.sum(N*flares.weights)/np.sum(flares.weights)



    # ax.axhline(np.log10(N_avg), c=c, alpha = 0.3)
    # ax.scatter(0.0, np.log10(N_avg), c=c, alpha = 1.0, s=20, marker = '*')

    ax.scatter(np.log10(1+flares.deltas), np.log10(N), s=z-7, c=[c], label = rf'$\rm z={z:.0f}$')


    # h, bin_edges = np.histogram(np.log10(np.array(n)+1), bins = np.arange(-0.3, 0.35, 0.05))
    # bin_centres = (bin_edges[1:]+bin_edges[:-1])/2.
    # hax2.scatter(bin_centres, np.log10(h), c=[c], alpha = 0.2, zorder = 0)

    hax2.hist(np.log10(np.array(n)+1), bins = np.arange(-0.3, 0.35, 0.05), color = c, histtype=u'step', log=False)

    # hax.hist(n, bins = np.arange(-0.3, 0.35, 0.025), color = c, histtype=u'step', log=True)



hax.set_xlim(xlim)
hax.set_xticks([])
# hax.set_yticks([])
hax.set_ylabel(r'$\rm N_{\rm sim}$')
hax2.set_ylabel(r'$\rm N(M_{\star}>10^{7.5}\ M_{\odot})$')

ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.legend(fontsize=8)

ax.set_xlabel(r'$\rm\log_{10}[1+\delta_{14}(z=4.7)]$')
ax.set_ylabel(r'$\rm \log_{10}[N(M_{\star}>10^{7.5}\ M_{\odot})]$')

fig.savefig(f'figs/N_delta.pdf')

fig.clf()



# --- all redshifts on one figure but normalised by average density


# fig, ax = fplt.simple()
#
# for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):
#
#     tag = flares.tag_from_zed[z]
#
#     X = flares.load_dataset(tag, f'Galaxy', 'Mstar')
#
#
#     N = []
#     for sim in flares.sims:
#         x = np.log10(np.array(X[sim])) +10.
#         N.append(np.sum(x>8.))
#
#     N = np.array(N)
#     N_avg = np.sum(N*flares.weights)/np.sum(flares.weights)
#
#     ax.scatter(np.log10(1+flares.deltas), np.log10(N/N_avg), s=5, c=[c], label = rf'$\rm z={z}$')
#
#
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
#
# ax.legend(fontsize=8)
#
# ax.set_xlabel(r'$\rm\log_{10}(1+\delta)$')
# ax.set_ylabel(r'$\rm \log_{10}[N(M_{\star}>10^{8}\ M_{\odot})/\bar{N}]$')
#
# fig.savefig(f'figs/N_delta_mean.pdf')
#
# fig.clf()



# --- individual redshifts

# for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):
#
#     fig, ax = fplt.simple()
#
#     tag = flares.tag_from_zed[z]
#
#     X = flares.load_dataset(tag, f'Galaxy', 'Mstar')
#
#
#     N = []
#     for sim in flares.sims:
#         x = np.log10(np.array(X[sim])) +10.
#         N.append(np.sum(x>8.))
#
#     N = np.array(N)
#     N_avg = np.sum(N*flares.weights)/np.sum(flares.weights)
#
#     ax.axhline(np.log10(N_avg), c=c, alpha = 0.3)
#     ax.scatter(np.log10(1+flares.deltas), np.log10(N), s=5, c=[c], label = rf'$\rm z={z}$')
#
#
#     ax.set_ylim([-0.75, 1.7])
#
#     ax.legend(fontsize=8)
#
#     ax.set_xlabel(r'$\rm\log_{10}(1+\delta)$')
#     ax.set_ylabel(r'$\rm \log_{10}[N(M_{\star}>10^{8}\ M_{\odot})]$')
#
#     fig.savefig(f'figs/N_delta_{z}.pdf')
#
#     fig.clf()
