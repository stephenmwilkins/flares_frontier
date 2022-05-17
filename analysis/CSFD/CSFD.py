
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot

import flares_utility.analyse as analyse













filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

# redshifts = [11,12,13,14,15]
# redshifts = [11,13,15]
redshifts = flares.zeds

print(redshifts)


# --- FLARES


# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3
totV = V*40
print(flares.radius)


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy/SFR', 'dataset': '10Myr', 'name': 'SFR', 'log10': True}) # SFR averaged over 100 Myr
# quantities.append({'path': 'Galaxy', 'dataset': 'Mstar', 'name': None, 'log10': True}) # Stellar mass

# --- get quantities (and weights and deltas)


for z in redshifts:
    D = flares.get_datasets(flares.tag_from_zed[z], quantities)
    # print(z, np.min(D['SFR']),np.median(D['SFR']), np.max(D['SFR']))
    s = D['SFR']>10.

    sfrd = np.sum(D['SFR'][s]*D['weight'][s])/(totV)
    # sfrd = np.sum(D['SFR'][s])/(totV)

    print(z, sfrd, np.log10(sfrd))
