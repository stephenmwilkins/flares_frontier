
import numpy as np
import scipy.stats as stats

import h5py

import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.limits
import flares_utility.table as table
import flares_utility.analyse as analyse


import flare.obs.literature.beta as beta_observations


bins = np.arange(27.7, 29.4, 0.2)

filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'DustModelI', 'name': f'beta', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': 'LFUV', 'log10': True})


D = {}

for tag, z in zip(flares.tags, flares.zeds):
    # --- get quantities (and weights and deltas)
    D[z] = flares.get_datasets(tag, quantities)

table.redshift(D, flares.zeds, 'log10LFUV', 'beta', filename = '../../flares_frontier_data/beta.dat', bins = bins)


# z = 10
#
# # --- get quantities (and weights and deltas)
# D = flares.get_datasets(flares.tag_from_zed[z], quantities)
#
# table.simple(D, 'log10FUV', 'beta', bins = bins)
