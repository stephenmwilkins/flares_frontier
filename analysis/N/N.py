


import numpy as np

import h5py

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines
import cmasher as cmr
import flare.photom as phot

import flares_utility.analyse as analyse

import flare.plt as fplt


fig, ax = fplt.simple()


# FLARES
FLARES_redshifts = np.arange(5,16,1)
N = {}
N[15.0]={ 7.0 : 475, 8.0 : 6, 9.0 : 0, 10.0 : 0, 11.0 : 0}
N[14.0]={ 7.0 : 1528, 8.0 : 30, 9.0 : 2, 10.0 : 0, 11.0 : 0}
N[13.0]={ 7.0 : 3893, 8.0 : 83, 9.0 : 2, 10.0 : 0, 11.0 : 0}
N[12.0]={ 7.0 : 8806, 8.0 : 190, 9.0 : 7, 10.0 : 0, 11.0 : 0}
N[11.0]={ 7.0 : 17195, 8.0 : 412, 9.0 : 30, 10.0 : 0, 11.0 : 0}
N[10.0]={ 7.0 : 30945, 8.0 : 892, 9.0 : 73, 10.0 : 3, 11.0 : 0}
N[9.0]={ 7.0 : 51539, 8.0 : 2033, 9.0 : 200, 10.0 : 8, 11.0 : 0}
N[8.0]={ 7.0 : 82414, 8.0 : 4453, 9.0 : 441, 10.0 : 25, 11.0 : 0}
N[7.0]={ 7.0 : 124474, 8.0 : 9359, 9.0 : 1012, 10.0 : 114, 11.0 : 0}
N[6.0]={ 7.0 : 177887, 8.0 : 18526, 9.0 : 2237, 10.0 : 355, 11.0 : 1}
N[5.0]={ 7.0 : 237437, 8.0 : 33415, 9.0 : 4991, 10.0 : 1033, 11.0 : 16}
N[4.77]={ 7.0 : 251053, 8.0 : 37731, 9.0 : 5956, 10.0 : 1279, 11.0 : 21}


EAGLE = h5py.File('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/EAGLE_REF_sp_info.hdf5', 'r')
snap = {5:'008_z005p037',6: '006_z005p971',7:'005_z007p050',8:'004_z008p075',9:'003_z008p988',10:'002_z009p993'}


EAGLE_redshifts = np.arange(5,11,1)

limits = np.arange(8,12,1)
c = 'k'

for limit, c in zip(limits, cmr.take_cmap_colors('cmr.torch', len(limits), cmap_range = [0.15, 0.85])):

    N_ = [np.sum(np.log10(EAGLE[snap[z]]['Galaxy/Mstar_30'])>(limit-10)) for z in EAGLE_redshifts]

    ax.plot(EAGLE_redshifts, np.log10(N_), ls = '--', c=c, alpha = 0.6)

    N_ = [N[z][limit] for z in FLARES_redshifts]

    ax.plot(FLARES_redshifts, np.log10(N_), ls = '-', c=c, alpha = 1.0, lw=1, label = rf'$\rm \log_{{10}}(M_{{\star}}/M_{{\odot}})>{limit}$')





handles = [mlines.Line2D([], [], color='0.5', ls=ls, lw=1, label=label) for ls, label in zip(['-','--'],[r'$\rm FLARES$',r'$\rm EAGLE$'])]

legend2 = plt.legend(loc ='upper left', handles=handles, fontsize = 8, labelspacing = 0.2)

plt.gca().add_artist(legend2)

ax.legend(loc ='upper right', fontsize = 8)


ax.set_xlim([5,15])
ax.set_ylim([-0.2, 5.5])
ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm \log_{10}(N)$')

fig.savefig(f'figs/N.pdf')

fig.clf()
