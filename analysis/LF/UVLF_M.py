
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot

import flares_utility.analyse as analyse




redshifts = [11,12,13,14,15]
redshifts = [13]

X_limits = [-24.5,-18.]






filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

print(flares.tags)






fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))





# --- FLARES


# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3

binw = 0.25
bin_edges = np.arange(*X_limits, binw)
bin_centres = bin_edges[:-1]+binw/2

for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):


    tag = flares.tag_from_zed[z]

    for phot_type, ls in zip(['DustModelI', 'Intrinsic'], ['-',':']):

        ## ---- get data
        phi = np.zeros(len(bin_centres))
        N = np.zeros(len(bin_centres))

        X = flares.load_dataset(tag, f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{phot_type}/', 'FUV')




        for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

            x = phot.lum_to_M(np.array(X[sim]))
            # x = x[x>0.0]

            # print(np.median(x))

            N_temp, _ = np.histogram(x, bins = bin_edges)

            N += N_temp

            phi_temp = (N_temp / V) / binw

            phi += phi_temp * w


        # print(phot_type, np.log10(phi))

        if ls == '-':
            ax.plot(bin_centres, np.log10(phi), ls = ls, c=c, label = rf'$\rm z={z}$', lw=2, alpha=0.5)
        else:
            ax.plot(bin_centres, np.log10(phi), ls = ls, c=c, lw=1)




# --- observations

M = [-23.5, [0.5, 0.5]]
phi = [3.7E-8, [3E-8, 9.4E-8]]

c= 'k'
ax.scatter(M[0], np.log10(phi[0]), lw=1, c=c)
ax.plot([M[0]-M[1][0], M[0]+M[1][1]], [np.log10(phi[0])]*2, lw=1, c=c)
ax.plot([M[0]]*2, [np.log10(phi[0]-phi[1][0]), np.log10(phi[0]+phi[1][1])], lw=1, c=c)




ax.legend(labelspacing=0.05,fontsize=8)
ax.set_xlim(X_limits)
# ax.set_ylim([-12., -1.51])


ax.set_xlabel(r'$\rm M_{FUV}$')
ax.set_ylabel(r'$\rm\log_{10}[\phi/Mpc^{-3}\ mag^{-1}]$')


fig.savefig(f'figs/UVLF_M.pdf')


fig.clf()
