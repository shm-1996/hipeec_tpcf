#Data Handling
import numpy as np
import itertools
import argparse 
import PhysicalConstantsCGS as const    
import pickle
from scipy.optimize import curve_fit

#System modules
import time
import os
import sys
import warnings
import shutil
import subprocess

#Astropy Data Handling
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
import scipy

#Visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmasher as cmr
mpl.style.use('classic')
mpl.rc_file('/Users/shm/.matplotlib/matplotlibrc',
    use_default_template=False)


##### Some global quantities #########
arcsec_to_degree = 1./3600.
list_of_galaxies = ['NGC_1614','NGC_3256','NGC_34','NGC_3690','NGC_4194',
        'NGC_6052']
       

catalog_links = \
['https://archive.stsci.edu/hlsps/hipeec/ngc1614/hlsp_hipeec_hst_acs-wfc-wfc3-u\
vis_ngc1614_multiband_v1_pdmod-sbext.dat',
'https://archive.stsci.edu/hlsps/hipeec/ngc3256/hlsp_hipeec_hst_acs-wfc-acs-hr\
c-wfc3-uvis_ngc3256_multiband_v1_pdmod-sbext.dat',
'https://archive.stsci.edu/hlsps/hipeec/ngc34/hlsp_hipeec_hst_acs-wfc-wfc3-uvis\
_ngc34_multiband_v1_pdmod-sbext.dat',
'https://archive.stsci.edu/hlsps/hipeec/ngc3690/hlsp_hipeec_hst_wfc3-uvis_ngc36\
90_multiband_v1_pdmod-sbext.dat',
'https://archive.stsci.edu/hlsps/hipeec/ngc4194/hlsp_hipeec_hst_acs-wfc-wfc3-uv\
is_ngc4194_multiband_v1_pdmod-sbext.dat',
'https://archive.stsci.edu/hlsps/hipeec/ngc6052/hlsp_hipeec_hst_wfc3-uvis_ngc60\
52_multiband_v1_pdmod-sbext.dat']

mosaic_links = \
['https://archive.stsci.edu/hlsps/hipeec/ngc1614/hlsp_hipeec_hst_acs-wfc_ngc161\
4_f435w_v1_drc.fits',
'https://archive.stsci.edu/hlsps/hipeec/ngc3256/hlsp_hipeec_hst_acs-wfc_ngc3256\
_f435w_v1_drc.fits',
'https://archive.stsci.edu/hlsps/hipeec/ngc34/hlsp_hipeec_hst_acs-wfc_ngc34_f43\
5w_v1_drc.fits',
'https://archive.stsci.edu/hlsps/hipeec/ngc3690/hlsp_hipeec_hst_wfc3-uvis_ngc36\
90_f438w_v1_drc.fits',
'https://archive.stsci.edu/hlsps/hipeec/ngc4194/hlsp_hipeec_hst_acs-wfc_ngc4194\
_f814w_v1_drc.fits',
'https://archive.stsci.edu/hlsps/hipeec/ngc6052/hlsp_hipeec_hst_wfc3-uvis_ngc60\
52_f438w_v1_drc.fits']

# Pickle Data Handling
def saveObj(obj, name):
    """
    Save a pickle object.

    INPUTS:
    ----------
    obj      - the name of the data object to be saved
    name     - the name of the pickle object

    """

    os.system("touch " + name + ".pkl")
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        #pickle.dump(obj, f, protocol=2)


def loadObj(name):
    """
    Load a pickle object.

    INPUTS:
    ----------
    name     - the name of the pickle object

    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def load_galaxy(galaxy_name,method='masked',directory=None):
    """
    Load a galaxy pickle object. 
    """
    if(method not in ['uniform','masked','masked_radial']):
        raise ValueError("This method does not exist. Allowed values are "+
            "'Uniform', 'Masked', and 'Masked_Radial'.")
    galaxy_name = galaxy_name.upper()
    method = method.lower()
    if(directory is None):
        directory = os.path.abspath('../Results/Galaxies/{}'.format(galaxy_name))
    
    if(method == 'masked'):
        galaxy_class = loadObj(directory+
                       '/Masked/{}_summary'.format(galaxy_name))

    elif(method == 'masked_radial'):
        galaxy_class = loadObj(directory+
                       '/Masked_Radial/{}_summary'.format(galaxy_name))

    else: 
        galaxy_class = loadObj(directory+
                       '/Uniform/{}_summary'.format(galaxy_name))

    return galaxy_class
