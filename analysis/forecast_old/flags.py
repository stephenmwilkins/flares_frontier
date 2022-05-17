
import os
import numpy as np
import pickle
import flare.photom




c = 'dodgerblue'



filters = [f'Webb.NIRCam.{f}' for f in ['F115W','F150W','F200W','F277W','F356W','F444W']]


from pathlib import Path

base_path = Path(__file__).parent




# --- other surveys

Surveys = {}

Surveys['JADES'] = []
Surveys['JADES'].append({'observatory': 'Webb', 'survey': 'JADES', 'sub': 'deep', 'area': 46., 'depths_abmag': {'Webb.NIRCam.F150W': 30.7, 'Webb.NIRCam.F200W': 30.7, 'Webb.NIRCam.F356W': 30.2}, 'public': False})
Surveys['JADES'].append({'observatory': 'Webb', 'survey': 'JADES', 'sub': 'medium', 'area': 144., 'depths_abmag': {'Webb.NIRCam.F150W': 29.7, 'Webb.NIRCam.F200W': 29.8, 'Webb.NIRCam.F356W': 29.4}, 'public': False})

# --- needs updating
Surveys['PANORAMIC'] = [{'observatory': 'Webb', 'survey': 'CEERS', 'sub':'', 'area': 150*9.1, 'depths_abmag': {'Webb.NIRCam.F115W': 28.1, 'Webb.NIRCam.F150W':28.29, 'Webb.NIRCam.F200W': 28.46, 'Webb.NIRCam.F356W': 28.03, 'Webb.NIRCam.F356W': 28.09, 'Webb.NIRCam.F444W': 27.81}, 'public': True}]

Surveys['CEERS'] = [{'observatory': 'Webb', 'survey': 'CEERS', 'sub':'', 'area': 10*9.1, 'depths_abmag': {'Webb.NIRCam.F150W': 28.9, 'Webb.NIRCam.F200W': 28.95, 'Webb.NIRCam.F356W': 28.95}, 'public': True}]

# --- needs updating
Surveys['PEARLS'] = []
Surveys['PEARLS'].append({'observatory': 'Webb', 'survey': 'CEERS', 'sub':'', 'area': 7*9.1, 'depths_abmag': {'Webb.NIRCam.F150W': 28.9, 'Webb.NIRCam.F200W': 28.95, 'Webb.NIRCam.F356W': 28.95}, 'public': True})


Surveys['COSMOS-Web'] = []
Surveys['COSMOS-Web'] += [{'observatory': 'Webb', 'survey': 'COSMOS-Web', 'sub':'', 'area': 96, 'depths_abmag': {'Hubble.ACS.f814w': 27.2, 'Webb.NIRCam.F115W': 26.74, 'Webb.NIRCam.F150W': 26.99, 'Webb.NIRCam.F277W': 27.43, 'Webb.NIRCam.F444W': 27.10}, 'public': True}]
Surveys['COSMOS-Web'] += [{'observatory': 'Webb', 'survey': 'COSMOS-Web', 'sub':'', 'area': 1071, 'depths_abmag': {'Hubble.ACS.f814w': 27.2, 'Webb.NIRCam.F115W': 27.13, 'Webb.NIRCam.F150W': 27.38, 'Webb.NIRCam.F277W': 27.82, 'Webb.NIRCam.F444W': 27.49}, 'public': True}]
Surveys['COSMOS-Web'] += [{'observatory': 'Webb', 'survey': 'OSMOS-Web', 'sub':'', 'area': 41, 'depths_abmag': {'Hubble.ACS.f814w': 27.2, 'Webb.NIRCam.F115W': 27.36, 'Webb.NIRCam.F150W': 27.61, 'Webb.NIRCam.F277W': 28.05, 'Webb.NIRCam.F444W': 27.72}, 'public': True}]
Surveys['COSMOS-Web'] += [{'observatory': 'Webb', 'survey': 'OSMOS-Web', 'sub':'', 'area': 753, 'depths_abmag': {'Hubble.ACS.f814w': 27.2, 'Webb.NIRCam.F115W': 27.52, 'Webb.NIRCam.F150W': 27.77, 'Webb.NIRCam.F277W': 28.21 }, 'public': True}]


Surveys['PRIMER'] = [{'observatory': 'Webb', 'survey': 'PRIMER', 'sub':'COSMOS-Shallow', 'area': 144.2, 'depths_abmag': {'Webb.NIRCam.F090W': 28.33 ,'Webb.NIRCam.F115W': 28.61, 'Webb.NIRCam.F150W': 28.81, 'Webb.NIRCam.F200W': 28.89, 'Webb.NIRCam.F277W': 28.85, 'Webb.NIRCam.F356W': 28.79, 'Webb.NIRCam.F410M': 28.04, 'Webb.NIRCam.F444W': 28.32}, 'public': True}, {'observatory': 'Webb', 'survey': 'PRIMER', 'sub':'COSMOS-Medium', 'area': 108.3, 'depths_abmag': {'Webb.NIRCam.F090W': 28.57,'Webb.NIRCam.F115W': 28.84, 'Webb.NIRCam.F150W': 29.03, 'Webb.NIRCam.F200W': 29.11, 'Webb.NIRCam.F277W': 29.08, 'Webb.NIRCam.F356W': 29.02, 'Webb.NIRCam.F410M': 28.27, 'Webb.NIRCam.F444W': 28.54}, 'public': True}, {'observatory': 'Webb', 'survey': 'PRIMER', 'sub':'COSMOS-Deep', 'area': 33.4, 'depths_abmag': {'Webb.NIRCam.F090W': 28.96,'Webb.NIRCam.F115W': 29.23, 'Webb.NIRCam.F150W': 29.43, 'Webb.NIRCam.F200W': 29.51, 'Webb.NIRCam.F277W': 29.47, 'Webb.NIRCam.F356W': 29.42, 'Webb.NIRCam.F410M': 28.65, 'Webb.NIRCam.F444W': 28.93}, 'public': True}, {'observatory': 'Webb', 'survey': 'PRIMER', 'sub':'UDS-Shallow', 'area': 234.02, 'depths_abmag': {'Webb.NIRCam.F090W': 27.92,'Webb.NIRCam.F115W': 28.21, 'Webb.NIRCam.F150W': 28.42, 'Webb.NIRCam.F200W': 28.48, 'Webb.NIRCam.F277W': 28.38, 'Webb.NIRCam.F356W': 28.37, 'Webb.NIRCam.F410M': 27.57, 'Webb.NIRCam.F444W': 27.79}, 'public': True}, {'observatory': 'Webb', 'survey': 'PRIMER', 'sub':'UDS-Medium', 'area': 175.17, 'depths_abmag': {'Webb.NIRCam.F090W': 28.32,'Webb.NIRCam.F115W': 28.62, 'Webb.NIRCam.F150W': 28.82, 'Webb.NIRCam.F200W': 28.69, 'Webb.NIRCam.F277W': 28.85, 'Webb.NIRCam.F356W': 28.77, 'Webb.NIRCam.F410M': 27.97, 'Webb.NIRCam.F444W': 28.21}, 'public': True}]


Surveys['NGDEEP'] = [{'observatory': 'Webb', 'survey': 'NGDEEP', 'sub':'', 'area': 2*9.1, 'depths_abmag': {'Webb.NIRCam.F115W': 30.90, 'Webb.NIRCam.F150W': 30.62, 'Webb.NIRCam.F200W': 30.62, 'Webb.NIRCam.F277W': 30.72, 'Webb.NIRCam.F356W': 30.70, 'Webb.NIRCam.F444W': 30.56}, 'public': True}]











Surveys['FLAGS'] = Surveys['NGDEEP'] + Surveys['COSMOS-Web']+ Surveys['PANORAMIC'] + Surveys['PRIMER'] + Surveys['CEERS']
Surveys['Webb Public Cycle 1'] = Surveys['FLAGS']
Surveys['Webb All Cycle 1'] = Surveys['FLAGS'] + Surveys['JADES']
Surveys['Webb/Cy1'] = Surveys['Webb All Cycle 1']
