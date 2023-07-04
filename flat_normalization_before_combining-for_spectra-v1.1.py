
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
from astropy.io import fits
import os
import warnings
warnings.filterwarnings('ignore')
from astropy import units as u
from pathlib import Path
import ccdproc as ccdp
from astropy.nddata import CCDData
from astropy.stats import mad_std
from argparse import ArgumentParser


# In[77]:


'''
clearup
'''

parser = ArgumentParser()
parser.add_argument("-input", help="images you want to combine, e.g. 'D_B_202001030-00*flat30s.*fits'", dest = 'input')
parser.add_argument("-method",  help="method: it can be 'mean' or 'median' ", dest="method", )
# parser.add_argument("-output", help="set a name for output combined flat image, e.g. Flat_30S.fits ", dest="output", )
parser.add_argument("-xp", help='nth x-axis pixel for setting a range for calculating the mean or median',                   dest='x_pixels')

args = parser.parse_args()

input_file = args.input
method1 = args.method
x_pixels = int(args.x_pixels)
# output_file = args.output

xp_upper = int(x_pixels + x_pixels * 0.11764705882352941)
xp_lower = int(x_pixels - x_pixels * 0.11764705882352941)


'''list the flat images'''
flat_only = ccdp.ImageFileCollection('.', glob_include = input_file)

for ind,v  in enumerate(flat_only.files):
    
    
    print(ind, v)
    flat = fits.open(v)
    
    primary_flat = flat[0].data
    
#     x_pixels = 1700
    
    profile_flat = primary_flat[0:,x_pixels]
    argmax_flat = np.argmax(profile_flat)
    
    if method1 == 'median':
        print('median:', np.nanmedian(primary_flat[argmax_flat][xp_lower:xp_upper]))
        correct = np.nanmedian(primary_flat[argmax_flat][xp_lower:xp_upper])
        
    elif method1 == 'mean':
        print('mean:', np.nanmean(primary_flat[argmax_flat][xp_lower:xp_upper]))
        correct = np.nanmean(primary_flat[argmax_flat][xp_lower:xp_upper])
        
    flat[0].data = flat[0].data / correct
    
    print('save Nor_'+v)
    flat.writeto('Nor_'+v)

