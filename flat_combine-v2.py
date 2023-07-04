
# coding: utf-8

# In[2]:


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


# In[1]:


# cd /sciproc/disk5/carey1/LOT/SLT/slt20201030/bias-dark


# In[1]:


# def inv_median(a):
#     return 1 / np.median(a)


# In[7]:


parser = ArgumentParser()
parser.add_argument("-input", help="images you want to combine, e.g. 'Flat_30S*.fits'", dest = 'input')
parser.add_argument("-method",  help="method: it can be 'average' or 'median' ", dest="method", )
parser.add_argument("-output", help="set a name for output combined flat image, e.g. Flat_30S.fits ", dest="output", )
args = parser.parse_args()

input_file = args.input
method1 = args.method
output_file = args.output


'''list the flat images'''
flat_only = ccdp.ImageFileCollection('.', glob_include = input_file)

'''Combine flat images'''
# combined_flat = ccdp.combine(flat_only.files, scale=inv_median,                 method= method1,                 sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,                 sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,                 mem_limit=350e6,unit='adu')
combined_flat = ccdp.combine(flat_only.files,     method= method1,   sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,                 sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,                 mem_limit=350e6,unit='adu')

'''add header combined'''
combined_flat.meta['combined'] = True   

'''save the combined dark image'''
combined_flat.write(output_file)

print(output_file,' saved')


# In[11]:


# combined_bias.write.help()

