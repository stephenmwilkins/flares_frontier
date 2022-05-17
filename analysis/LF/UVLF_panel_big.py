
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.analyse as analyse

from flare.LF.literature import UV
from flare.LF.models import bin_centres as _bin_centres

redshifts = [10,11,12,13,14,15][::-1]

binw = 0.1
X_limits = [27.75, 30.5]
Y_limits = [ -7.9,-2.01]
# binw = 0.5
# X_limits = [27.9, 32]
# # X_limits = [27.9, 30.5]



filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

print(flares.tags)




N = len(redshifts)
left = 0.1
top = 0.95
bottom = 0.1
right = 0.9
rows = 2
panel_width = (right-left)/int(N/2)
panel_height = top-bottom
fig, axes = plt.subplots(2, int(N/2), figsize = (7,rows*7/(panel_height/panel_width)), sharey = True, sharex = True)
plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.1)
axes = axes.flatten()



# --- FLARES


# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3


bin_edges = np.arange(*X_limits, binw)
bin_centres = bin_edges[:-1]+binw/2


imarker = 0
markers = ['o','d','h','p','s','D']


line_styles = ['-','--','-.',':']
used_ls = {}

for z, c, ax in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts)), axes):


    tag = flares.tag_from_zed[z]

    for phot_type, ls in zip(['DustModelI', 'Pure_Stellar'], ['-','--']):

        ## ---- get data
        phi = np.zeros(len(bin_centres))
        N = np.zeros(len(bin_centres))

        X = flares.load_dataset(tag, f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{phot_type}/', 'FUV')


        for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

            x = np.log10(np.array(X[sim]))
            x = x[x>0.0]

            N_temp, _ = np.histogram(x, bins = bin_edges)

            N += N_temp

            phi_temp = (N_temp / V) / binw

            phi += phi_temp * w


        print(phot_type, np.log10(phi))

        s = N>2



        if z==10 and phot_type=='DustModelI':
            for ax_ in axes:
                ax_.plot(bin_centres[s], np.log10(phi[s]), ls = ls, c='k', alpha=0.2, zorder=0, lw=2)



        ax.plot(bin_centres[s], np.log10(phi[s]), ls = ls, c=c, lw=1, zorder=1)



    ax.text(0.5, 1.02, rf'$\rm z={z}$', horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize = 7)

    # --- add observations


    if z in UV.observed_binned.keys():

        for lf in UV.observed_binned[z]: # doesn't preserve marker style



            if lf.M_binw[z][0]<0.2:
                ax.plot(lf.log10L[z], lf.log10phi[z], color='0.4', lw=2, label = rf'$\rm {lf.name}$')
                ax.fill_between(lf.log10L[z], lf.log10phi[z]+lf.log10phi_err_low[z], lf.log10phi[z]+lf.log10phi_err_up[z], color='k',alpha=0.1)
            else:
                yerr = [lf.log10phi_err_low[z], lf.log10phi_err_up[z]]
                xerr = lf.log10L_binw[z]/2.
                ax.errorbar(lf.log10L[z], lf.log10phi[z], xerr=xerr, yerr = yerr, label = rf'$\rm {lf.name}$', color = '0.4',fmt='o',elinewidth=1, lw=1, marker = markers[imarker], ms=5)
                imarker += 1


    if z in UV.theoretical_parameterised.keys():

        for lf in UV.theoretical_parameterised[z]: # doesn't preserve marker style

            # --- this stops repeating the line style
            nm = lf().nameref
            first_time = True
            if nm in used_ls.keys():
                ls = used_ls[nm]
                first_time = False
            else:
                ls = line_styles[0]
                used_ls[nm] = ls
                line_styles.pop(0)

            print(line_styles)
            print(used_ls)

            _binw = 0.01
            _bin_edges = np.arange(27, 30.5, _binw)
            _phi = lf().phi_binned(z, _bin_edges)/_binw

            if first_time:
                ax.plot(_bin_centres(_bin_edges), np.log10(_phi), ls=ls, lw=1, alpha=0.5, c='0.4', label = rf'$\rm {nm}$')
            else:
                ax.plot(_bin_centres(_bin_edges), np.log10(_phi), ls=ls, lw=1, alpha=0.5, c='0.4')



    ax.legend(loc ='lower left', fontsize = 7, handletextpad = 0.4)

    ax.set_xlim(X_limits)
    ax.set_ylim(Y_limits)
    ax.set_xticks(np.arange(28, 31, 1.0))







handles = [mlines.Line2D([], [], color='0.5', ls=ls, lw=1, label=label) for ls, label in zip(['-','--'],[r'$\rm observed$',r'$\rm intrinsic$'])]
handles += [mlines.Line2D([], [], alpha=0.2, color='k', ls='-', lw=2, label=r'$\rm z=10$')]

axes[0].legend(loc ='upper right', handles=handles, fontsize = 8, labelspacing = 0.2)




# axes[0].set_ylabel(r'$\rm\log_{10}[\phi/Mpc^{-3}\ dex^{-1}]$', fontsize = 9)

fig.text(0.04, bottom+(top-bottom)/2, r'$\rm\log_{10}[\phi/Mpc^{-3}\ dex^{-1}]$', rotation = 90, va='center', fontsize = 9)

fig.text(left+(right-left)/2, 0.04, r'$\rm \log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$', ha='center', fontsize = 9)


fig.savefig(f'figs/UVLF_panel.pdf')


fig.clf()
