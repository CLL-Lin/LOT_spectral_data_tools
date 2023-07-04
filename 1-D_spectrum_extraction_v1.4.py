#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
from scipy.optimize import curve_fit
import pandas as pd
import warnings
from argparse import ArgumentParser
import os
warnings.filterwarnings('ignore')


# In[1]:


# cd /sciproc/disk5/carey1/LOT/SLT/slt20201030/bias-dark


# In[2]:


def Gauss(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2)) + b


# In[7]:


parser = ArgumentParser()
parser.add_argument("-input", help="images you want to extract 1-D spectra, e.g. 'D_B_EPIC210651981*.fits' or 'D_B_EPIC210651981_01.fits'", dest = 'input')
parser.add_argument("-flat",  help="the slit master flat image for flat correciton ", dest="flat" )
parser.add_argument("-ow", help="True = overwrite, False = Don't overwrite", dest="overwrite", default = False)
# parser.add_argument("-output", help="set a name for output .csv file contains the data of 1-D extracted spectrum, e.g. Flat_30S.fits ", dest="output", )
args = parser.parse_args()

input_file = args.input
input_flat = args.flat
overwrite = args.overwrite
print(input_file)


if os.path.exists(input_file[:-5]+'-1D_extracted_spectrum.csv') == True:
    
    if overwrite == 'True':


        '''readin the fits file'''
        hdus = pyfits.open(input_file)
        flat = pyfits.open(input_flat)

        primary = hdus[0].data  
        primary_flat = flat[0].data

        '''extract spectrum started at 700th pixel (for LISA only)'''
        x_axis_pixels = np.arange(700,len(primary[0]),1)
        
        
        '''fit optimal centers of apertures'''
        amax = []
        for i,v in enumerate(x_axis_pixels):

            x_pixels = v

            profile = primary[0:,x_pixels]
            profile_flat = primary_flat[0:,x_pixels]
        #     profile_flat = primary_flat.sum(axis=1)

            argmax = np.argmax(profile)
        #     print(argmax)
            amax.append(argmax)
            
            
        fit_apc = np.polyfit(x_axis_pixels, amax, 1)
        fit_ap = np.poly1d(fit_apc)
        

        '''extract 1-D spectrum'''
        extra_flux_1 = [] 
        lower_aperture_for_lamp = [] # aperture information for lamp wavelenght calibration
        upper_aperture_for_lamp = []

        for i,v in enumerate(x_axis_pixels):

            x_pixels = v

            profile = primary[0:,x_pixels]
            profile_flat = primary_flat[0:,x_pixels]
#             profile_flat = primary_flat.sum(axis=1)
#             argmax = np.argmax(profile)
            argmax = int (fit_ap(x_pixels))
            y_axis_pixels = np.arange(0, len(profile),1)

            '''Flat-corrected intensity profile along y-axis'''
            profile_flat_corr = profile/profile_flat

            try:
                '''
                Gaussian fit for flux profile
                '''
                x = y_axis_pixels[argmax-100 : argmax+100]
                # y = profile[argmax-50 : argmax+50]
                y = profile_flat_corr[argmax-100 : argmax+100]
                y_min = mini = np.min(y)
                mean = sum(x*y)/sum(y)
                sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
        #         popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma, y_min], maxfev= int(1e6))
                popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma, y_min])

                sigma_fr = np.abs(popt[2])
                mean_fr = popt[1]
                # print(popt)
                
                '''goodness of fit test'''
                chi_squ = np.sum((y - Gauss(x, *popt))**2                 / Gauss(x, *popt) )

                dof = len(Gauss(x, *popt)) - len(popt) #dof

                red_chisq = chi_squ/dof

                '''
                extract flux in aperture
                '''

                lower_aperture = mean_fr - 3 * sigma_fr
                upper_aperture = mean_fr + 3 * sigma_fr

                lower_aperture_for_lamp.append(lower_aperture)
                upper_aperture_for_lamp.append(upper_aperture)

                lower_background1 = mean_fr - 10 * sigma_fr
                lower_background2 = (mean_fr - 10 * sigma_fr) + 3 * sigma_fr

                upper_background1 = (mean_fr + 10 * sigma_fr) - 3 * sigma_fr
                upper_background2 =  mean_fr + 10 * sigma_fr



                profile_flat_corr = pd.Series(profile_flat_corr)
                profile = pd.Series(profile)
                y_axis_pixels = pd.Series(y_axis_pixels)

                lower_background = profile_flat_corr[y_axis_pixels >= lower_background1][y_axis_pixels <= lower_background2].values
                upper_background = profile_flat_corr[y_axis_pixels >= upper_background1][y_axis_pixels <= upper_background2].values
        #         lower_background = profile[y_axis_pixels >= lower_background1][y_axis_pixels <= lower_background2].values
        #         upper_background = profile[y_axis_pixels >= upper_background1][y_axis_pixels <= upper_background2].values

                background = (np.nanmedian(lower_background)+ np.nanmedian(upper_background)) / 2

                signal_in_ap = profile_flat_corr [y_axis_pixels>=lower_aperture][y_axis_pixels<=upper_aperture] - background
        #         signal_in_ap = profile [y_axis_pixels>=lower_aperture][y_axis_pixels<=upper_aperture] - background

                signal_in_ap_Sum = signal_in_ap.sum()
                print(v, signal_in_ap_Sum)

                if signal_in_ap_Sum <= 0 or red_chisq >= 100:
        #             extra_flux.append(np.nan)
                    extra_flux_1.append(np.nan)
                else:
        #             extra_flux.append(signal_in_ap_Sum)
                    extra_flux_1.append(signal_in_ap_Sum)

            except (ValueError, RuntimeError):
        #         extra_flux.append(np.nan)
                extra_flux_1.append(np.nan)

                lower_aperture_for_lamp.append(np.nan)
                upper_aperture_for_lamp.append(np.nan)


        '''save the extracted result as a .csv file'''
        re1 = pd.DataFrame([], columns=['x-axis-pixel','count/s','lower_Aperture','upper_Aperture'])
        re1['x-axis-pixel'] = pd.Series(x_axis_pixels)
        re1['count/s'] = pd.Series(extra_flux_1)  / hdus[0].header['EXPTIME']
        re1['lower_Aperture'] = pd.Series(lower_aperture_for_lamp)
        re1['upper_Aperture'] = pd.Series(upper_aperture_for_lamp)

        re1.to_csv(input_file[:-5]+'-1D_extracted_spectrum.csv')
        
        print(input_file[:-5]+'-1D_extracted_spectrum.csv',' saved')
        
    elif overwrite == False:
        print(input_file[:-5]+'-1D_extracted_spectrum.csv',' already exists')
        
else:
    
    '''readin the fits file'''
    hdus = pyfits.open(input_file)
    flat = pyfits.open(input_flat)

    primary = hdus[0].data  
    primary_flat = flat[0].data

    '''extract spectrum for 700th pixel (for LISA only)'''
    x_axis_pixels = np.arange(700,len(primary[0]),1)

    '''fit optimal centers of apertures'''
    amax = []
    for i,v in enumerate(x_axis_pixels):

        x_pixels = v

        profile = primary[0:,x_pixels]
        profile_flat = primary_flat[0:,x_pixels]

        argmax = np.argmax(profile)
        amax.append(argmax)


    fit_apc = np.polyfit(x_axis_pixels, amax, 1)
    fit_ap = np.poly1d(fit_apc)
    
    
    '''extract 1-D spectrum'''
    extra_flux_1 = [] 
    lower_aperture_for_lamp = [] # aperture information for lamp wavelenght calibration
    upper_aperture_for_lamp = []

    for i,v in enumerate(x_axis_pixels):

        x_pixels = v

        profile = primary[0:,x_pixels]
        profile_flat = primary_flat[0:,x_pixels]
#         profile_flat = primary_flat.sum(axis=1)
#         argmax = np.argmax(profile)
        argmax = int (fit_ap(x_pixels))
        y_axis_pixels = np.arange(0, len(profile),1)

        '''Flat-corrected intensity profile along y-axis'''
        profile_flat_corr = profile/profile_flat

        try:
            '''
            Gaussian fit for flux profile
            '''
            x = y_axis_pixels[argmax-100 : argmax+100]
            # y = profile[argmax-50 : argmax+50]
            y = profile_flat_corr[argmax-100 : argmax+100]
            y_min = mini = np.min(y)
            mean = sum(x*y)/sum(y)
            sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
    #         popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma, y_min], maxfev= int(1e6))
            popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma, y_min])

            sigma_fr = np.abs(popt[2])
            mean_fr = popt[1]
            # print(popt)
            
            '''goodness of fit test'''
            chi_squ = np.sum((y - Gauss(x, *popt))**2             / Gauss(x, *popt) )

            dof = len(Gauss(x, *popt)) - len(popt) #dof

            red_chisq = chi_squ/dof


            '''
            extract flux in aperture
            '''

            lower_aperture = mean_fr - 3 * sigma_fr
            upper_aperture = mean_fr + 3 * sigma_fr

            lower_aperture_for_lamp.append(lower_aperture)
            upper_aperture_for_lamp.append(upper_aperture)

            lower_background1 = mean_fr - 10 * sigma_fr
            lower_background2 = (mean_fr - 10 * sigma_fr) + 3 * sigma_fr

            upper_background1 = (mean_fr + 10 * sigma_fr) - 3 * sigma_fr
            upper_background2 =  mean_fr + 10 * sigma_fr



            profile_flat_corr = pd.Series(profile_flat_corr)
            profile = pd.Series(profile)
            y_axis_pixels = pd.Series(y_axis_pixels)

            lower_background = profile_flat_corr[y_axis_pixels >= lower_background1][y_axis_pixels <= lower_background2].values
            upper_background = profile_flat_corr[y_axis_pixels >= upper_background1][y_axis_pixels <= upper_background2].values
    #         lower_background = profile[y_axis_pixels >= lower_background1][y_axis_pixels <= lower_background2].values
    #         upper_background = profile[y_axis_pixels >= upper_background1][y_axis_pixels <= upper_background2].values

            background = (np.nanmedian(lower_background)+ np.nanmedian(upper_background)) / 2

            signal_in_ap = profile_flat_corr [y_axis_pixels>=lower_aperture][y_axis_pixels<=upper_aperture] - background
    #         signal_in_ap = profile [y_axis_pixels>=lower_aperture][y_axis_pixels<=upper_aperture] - background

            signal_in_ap_Sum = signal_in_ap.sum()
            print(v, signal_in_ap_Sum)

            if signal_in_ap_Sum <= 0 or red_chisq >= 100:
    #             extra_flux.append(np.nan)
                extra_flux_1.append(np.nan)
            else:
    #             extra_flux.append(signal_in_ap_Sum)
                extra_flux_1.append(signal_in_ap_Sum)

        except (ValueError, RuntimeError):
    #         extra_flux.append(np.nan)
            extra_flux_1.append(np.nan)

            lower_aperture_for_lamp.append(np.nan)
            upper_aperture_for_lamp.append(np.nan)


    '''save the extracted result as a .csv file'''
    re1 = pd.DataFrame([], columns=['x-axis-pixel','count/s','lower_Aperture','upper_Aperture'])
    re1['x-axis-pixel'] = pd.Series(x_axis_pixels)
    re1['count/s'] = pd.Series(extra_flux_1) / hdus[0].header['EXPTIME']
    re1['lower_Aperture'] = pd.Series(lower_aperture_for_lamp)
    re1['upper_Aperture'] = pd.Series(upper_aperture_for_lamp)

    re1.to_csv(input_file[:-5]+'-1D_extracted_spectrum.csv')

    print(input_file[:-5]+'-1D_extracted_spectrum.csv',' saved')
    
    
    
    count_rate = re1['count/s'].values
    x_pixel = re1['x-axis-pixel'].values
#     wavelength = re1['wavelength(\\AA)'].values

    plt.xlabel('X-axis Pixels')
    plt.ylabel('Reltive flux')
    plt.plot(x_pixel, count_rate)
#     plt.plot(wavelength, count_rate)

#     plt.xlabel('$\lambda (\AA)$')
#     plt.ylabel('Reltive flux')
    plt.grid(True)
    plt.savefig(input_file[:-5]+'-1D_extracted_spectrum.png', format='png')


# In[11]:


# combined_bias.write.help()


# In[ ]:




