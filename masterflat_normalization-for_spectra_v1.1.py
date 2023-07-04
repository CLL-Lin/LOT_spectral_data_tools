
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
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


# In[2]:


def Gauss(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + b



# In[13]:


# cd /sciproc/disk5/carey1/LOT/LISA/LOT20201031/hclin


# In[14]:


# cd flat_20201030/


# In[15]:


# ls


# In[12]:


# def inv_median(a):
#     return 1 / np.median(a)


# In[ ]:


parser = ArgumentParser()
parser.add_argument("-input", help="Flat images you want to normalize, e.g. 'Flat_30S*fits' ", dest = 'input')
# parser.add_argument("-method",  help="method: it can be 'average' or 'median' ", dest="method", )
# parser.add_argument("-output", help="set a name for output combined flat image, e.g. Flat_30S.fits ", dest="output", )
args = parser.parse_args()

input_file = args.input
# method1 = args.method
# output_file = args.output


'''list the flat images'''
flat_only = ccdp.ImageFileCollection('.', glob_include = input_file)


for ind,v in enumerate(flat_only.files):
    
    print(ind, flat_only.files[ind])
    flat = fits.open(flat_only.files[ind])
    primary_flat = flat[0].data
    primary_flat_for_nor = flat[0].data
    
    flat_median = []

    '''find the median'''
    print('find the median')
    for i in range(len(primary_flat[0])):

        x_pixels = i

        profile_flat = primary_flat[0:,x_pixels]
        argmax_flat = np.argmax(profile_flat)

        y_axis_pixels_flat = np.arange(0, len(profile_flat),1)

        #Gaussian fit
        x = y_axis_pixels_flat[np.abs(argmax_flat-1000) : argmax_flat+1000]
        # y = profile[argmax-1000 : argmax+1000]
        y = profile_flat[np.abs(argmax_flat-1000) : argmax_flat+1000]
        y_min  = np.min(y)
        mean = sum(x*y)/sum(y)
        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
        popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma, y_min])


        ind1 = int(popt[1] - 1.25 *  popt[2])
        ind2 = int(popt[1] + 1.25 *  popt[2])
    #     print(np.median(profile_flat[profile_flat>=3]))

    #     f_median = np.median(profile_flat[profile_flat > np.std(profile_flat) * 2.3])
        f_median = np.median(profile_flat[ind1: ind2+1])


    #     primary_flat_for_nor[0:, x_pixels] = primary_flat_for_nor[0:, x_pixels] / f_median
        flat_median.append(f_median)
    
    
    '''Fit median trend along with x-axis'''
    # print('Fit median trend along with x-axis')
    # x_pixels_array = np.arange(0, len(flat_median))
    # fit_f = np.polyfit(x_pixels_array, flat_median, 25)
    # fit_fn = np.poly1d(fit_f)
    
    '''normalize the flat'''
    print('normalize the flat')
    for a in range(len(flat_median)):
    
        x_pixels = a
        primary_flat_for_nor[0:, x_pixels] = primary_flat_for_nor[0:, x_pixels] / flat_median[x_pixels]
        
    
    '''save the normalized flat file'''
    print('save file: Nor_'+flat_only.files[ind])
    flat.writeto('Nor_'+flat_only.files[ind])

