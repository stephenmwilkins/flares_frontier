
import numpy as np

import h5py

import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.analyse as analyse


from astropy.io import ascii



redshifts = [10,11,12,13,14,15][::-1]

binw = 0.2
X_limits = [27.7, 30.]


filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3


bin_edges = np.arange(*X_limits, binw)
bin_centres = bin_edges[:-1]+binw/2


table = {'log10LFUV': bin_centres}
formats = {'log10LFUV': '%.1f'}

for z in redshifts:

    tag = flares.tag_from_zed[z]

    for phot_type in ['DustModelI']: # 'Pure_Stellar'

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


    table['log10phi_z'+str(z)] = np.log10(np.array(phi))
    table['log10phi_z'+str(z)][table['log10phi_z'+str(z)]==-np.inf] = None
    formats['log10phi_z'+str(z)] = '%.2f'


ascii.write(table, '../../flares_frontier_data/UVLF.dat', formats = formats, overwrite = True)
ascii.write(table, Writer=ascii.Latex, formats = formats)
