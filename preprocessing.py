#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:31:17 2019

@author: astro_15
"""

from pathlib import Path
from astropy.io import fits
import pandas as pd

#%%
allfits = list(Path("./AO2019-2/2019-10-26").glob("*.fit"))
print(allfits)

#%%
fits.open(allfits[65])[0].header
#%%
with fits.open(allfits[65], 'update') as f:
    for hdu in f:
            hdu.header['OBJECT'] = 'Neptune'
    

#Neptune.fit file didn't have header 'object' (took image as single when observing), so I added neptune.

#%%
fpaths = dict(bias=[], dark={}, flat=[], comp=[], objt={})

for fpath in allfits:
    hdr = fits.getheader(fpath)
    imagetyp = hdr["IMAGETYP"].lower() #lower() - somunja
    exptime = float(hdr["EXPTIME"])
    obs_object = hdr["OBJECT"]
    if imagetyp.startswith('bias frame') and exptime == 0. and obs_object.lower() == 'bias':
        fpaths['bias'].append(fpath)
    elif imagetyp.startswith('dark frame'): #elif - if not, if
        if exptime == 5 or exptime == 10 or exptime == 20 or exptime == 60:
            try:
                fpaths['dark'][exptime]
            except KeyError:
                fpaths['dark'][exptime] = []
            fpaths['dark'][exptime].append(fpath)
    elif imagetyp.startswith('flat field'):
        fpaths['flat'].append(fpath)
    elif imagetyp.startswith('light frame') and obs_object.lower() == 'comp':
        fpaths['comp'].append(fpath)
    else:
        if obs_object.lower() == 'neptune' or obs_object.lower() =='uranus' or obs_object.lower() =='uranus20' or obs_object.lower() =='hr9087' or obs_object.lower() =='hr718' :
            try:
                fpaths['objt'][obs_object.lower()]
            except KeyError:
                fpaths['objt'][obs_object.lower()] = []
            fpaths['objt'][obs_object.lower()].append(fpath)

allfits.sort()

fpaths

# sorting into bias, dark(exptime), flat, comp, and objt(neptune, uranus, standard stars)
# took Uranus 2 times. uranus in 10 seconds and uranus20 in 20 seconds
# I just preprocessed both uranus and uranus20, but considering values, uranus20 will be used for furthur process.

#%%
from astropy.nddata import CCDData
from ccdproc import Combiner as combiner
from ccdproc import trim_image as trim
import astropy.units as u

b = []

for i in range (0, 9):
    hdul = fits.open(fpaths['bias'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    if i == 0:
        newhdr = hdr

    ccd = CCDData(data, unit=u.adu)
    ccd = trim(ccd, fits_section='[900:1100, 150:1750]')
    # trim images (range chosen by eyes)
    b.append(ccd)
master_bias = combiner(b).median_combine()

print(master_bias)

# combining bias

#%%
try: 
    _ = newhdr["PROCESS"]
except KeyError:
    newhdr["PROCESS"] = ''

newdata = master_bias
newhdr["PROCESS"] += "B"
newhdr.add_history(f"Bias median combined")


newhdul = fits.PrimaryHDU(data = newdata, header = newhdr)
newhdul.writeto(Path("./AO2019-2")/"image_bias.fits", overwrite=True)

# making combined bias fits file
# adding history header


#%%
hdul = fits.open("./AO2019-2/image_bias.fits")
hdr = hdul[0].header
data = hdul[0].data
hdr

# just identifying master_bias
#%%
d = []

for i in range (0, 9):
    hdul = fits.open(fpaths['dark'][5.0][i])
    hdr = hdul[0].header
    data = hdul[0].data
    if i == 0:
        newhdr5 = hdr

    ccd = CCDData(data, unit=u.adu)
    ccd = trim(ccd, fits_section='[900:1100, 150:1750]')
    d.append(ccd)
    
mas_dark_5 = combiner(d).median_combine()

master_dark_5 = mas_dark_5.data - master_bias.data

print(master_dark_5)

#%%
d = []

for i in range (0, 9):
    hdul = fits.open(fpaths['dark'][10.0][i])
    hdr = hdul[0].header
    data = hdul[0].data
    if i == 0:
        newhdr10 = hdr

    ccd = CCDData(data, unit=u.adu)
    ccd = trim(ccd, fits_section='[900:1100, 150:1750]')
    d.append(ccd)
    
mas_dark_10 = combiner(d).median_combine()

master_dark_10 = mas_dark_10.data - master_bias.data

print(master_dark_10)

#%%
d = []

for i in range (0, 9):
    hdul = fits.open(fpaths['dark'][20.0][i])
    hdr = hdul[0].header
    data = hdul[0].data
    if i == 0:
        newhdr20 = hdr

    ccd = CCDData(data, unit=u.adu)
    ccd = trim(ccd, fits_section='[900:1100, 150:1750]')
    d.append(ccd)
    
mas_dark_20 = combiner(d).median_combine()

master_dark_20 = mas_dark_20.data - master_bias.data

print(master_dark_20)

#%%
d = []

for i in range (0, 9):
    hdul = fits.open(fpaths['dark'][60.0][i])
    hdr = hdul[0].header
    data = hdul[0].data
    if i == 0:
        newhdr60 = hdr

    ccd = CCDData(data, unit=u.adu)
    ccd = trim(ccd, fits_section='[900:1100, 150:1750]')
    d.append(ccd)
    
mas_dark_60 = combiner(d).median_combine()

master_dark_60 = mas_dark_60.data - master_bias.data

print(master_dark_60)

# combining dark in each exptime, then subtracting master_bias

#%%
try: 
    _ = newhdr5["PROCESS"]
    _ = newhdr10["PROCESS"]
    _ = newhdr20["PROCESS"]
    _ = newhdr60["PROCESS"]
except KeyError:
    newhdr5["PROCESS"] = ''
    newhdr10["PROCESS"] = ''
    newhdr20["PROCESS"] = ''
    newhdr60["PROCESS"] = ''

newdata5 = master_dark_5
newdata10 = master_dark_10
newdata20 = master_dark_20
newdata60 = master_dark_60

newhdr5["PROCESS"] += "D"
newhdr10["PROCESS"] += "D"
newhdr20["PROCESS"] += "D"
newhdr60["PROCESS"] += "D"

newhdr5.add_history(f"Dark median combined and Bias subtracted")
newhdr10.add_history(f"Dark median combined and Bias subtracted")
newhdr20.add_history(f"Dark median combined and Bias subtracted")
newhdr60.add_history(f"Dark median combined and Bias subtracted")


newhdul5 = fits.PrimaryHDU(data = newdata5, header = newhdr5)
newhdul5.writeto(Path("./AO2019-2")/"image_dark_5.fits", overwrite=True)
newhdul10 = fits.PrimaryHDU(data = newdata10, header = newhdr10)
newhdul10.writeto(Path("./AO2019-2")/"image_dark_10.fits", overwrite=True)
newhdul20 = fits.PrimaryHDU(data = newdata20, header = newhdr20)
newhdul20.writeto(Path("./AO2019-2")/"image_dark_20.fits", overwrite=True)
newhdul60 = fits.PrimaryHDU(data = newdata60, header = newhdr60)
newhdul60.writeto(Path("./AO2019-2")/"image_dark_60.fits", overwrite=True)

# making master_dark_exptime fits file
#%%

hdul = fits.open("./AO2019-2/image_dark_5.fits")
hdr = hdul[0].header
data = hdul[0].data
hdr
#just identifying one image
#%%
fpaths

#%%
#flat = 20 seconds
import numpy as np
s =[]

for i in range (0, 9):
    hdul = fits.open(fpaths['flat'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    ccd = CCDData(data, unit=u.adu)
    ccd = trim(ccd, fits_section='[900:1100, 150:1750]')
    ccd = ccd.data - master_bias.data - master_dark_20.data
    s.append(ccd)
    if i ==0:
        newhdr = hdr

print(s)

flat_mean = np.mean(s)
print(flat_mean)
# flat_mean is a number
flat_norm = s / flat_mean
print(flat_norm)
master_flat = np.median(flat_norm, axis=0)
print(master_flat)
# flat_norm is an array

np.where(master_flat==0)

# making master_flat by each bias and dark subtracted flat, and normalizing and median combining.
#%%
#master_flat[294, 1986] = 1.0
#master_flat[804,4] = 1.0
#master_flat[806,18] = 1.0
#master_flat[809,31] = 1.0
#master_flat[841,0] = 1.0
#master_flat[851,28] = 1.0
#master_flat[857,30] = 1.0
#master_flat[874,22] = 1.0
#master_flat[875,21] = 1.0
#master_flat[891,26] = 1.0

# when I proceeded without trimming image, some pixels showed 0 value, so I made those values as 1,
# but process with trimming did not show 0 values, so I left this cell out.
#%%
try: 
    _ = newhdr["PROCESS"]
except KeyError:
    newhdr["PROCESS"] = ''

newdata = master_flat
newhdr["PROCESS"] += "F"
newhdr.add_history(f"flat subtracted and normalized and median combined")


newhdul = fits.PrimaryHDU(data = newdata, header = newhdr)
newhdul.writeto(Path("./AO2019-2")/"image_flat.fits", overwrite=True)

# making master_flat fits file

#%%
hdul = fits.open("./AO2019-2/image_flat.fits")
hdr = hdul[0].header
hdr

# again checking

#%%
fpaths

#%%

for i in range (0, 9):
    hdul = fits.open(fpaths['objt']['neptune'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    data = CCDData(data, unit=u.adu)
    data = trim(data, fits_section='[900:1100, 150:1750]')
    try: 
        _ = hdr["PROCESS"]
    except KeyError:
        hdr["PROCESS"] = ''
    
    data_reduced = data - master_bias.data
    hdr["PROCESS"] += "B"
    hdr.add_history(f"Bias subtracted by master_bias")
              
    hdu = fits.PrimaryHDU(data=data_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"neptune_0_bxx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"neptune_1_bxx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"neptune_2_bxx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"neptune_3_bxx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"neptune_4_bxx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"neptune_5_bxx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"neptune_6_bxx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"neptune_7_bxx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"neptune_8_bxx.fits", overwrite=True)
        
    data_d_reduced = data_reduced - master_dark_60
    hdr["PROCESS"] += "D"
    hdr.add_history(f"Dark subtracted by master_dark_60")
    
    hdu = fits.PrimaryHDU(data=data_d_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"neptune_0_bdx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"neptune_1_bdx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"neptune_2_bdx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"neptune_3_bdx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"neptune_4_bdx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"neptune_5_bdx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"neptune_6_bdx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"neptune_7_bdx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"neptune_8_bdx.fits", overwrite=True)
        
    data_corrected = data_d_reduced/master_flat
    hdr["PROCESS"] += "F"
    hdr.add_history(f"Flat corrected by master_flat")
    
    hdu = fits.PrimaryHDU(data=data_corrected, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"neptune_0_bdf.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"neptune_1_bdf.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"neptune_2_bdf.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"neptune_3_bdf.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"neptune_4_bdf.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"neptune_5_bdf.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"neptune_6_bdf.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"neptune_7_bdf.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"neptune_8_bdf.fits", overwrite=True)
    
#preprocessing neptune at first.
#in each 9 images, bias subtraction - dark subtraction - flat correction took place,
#and each process was saved as bxx - bdx - bdf file.
#%%
from matplotlib import pyplot as plt

fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharex=False, sharey=False, gridspec_kw=None)

for i, fpath in enumerate([Path("./AO2019-2/2019-10-26/Neptune-0005.fit"),
                           Path("./AO2019-2/neptune_0_bxx.fits"),
                           Path("./AO2019-2/neptune_0_bdx.fits"),
                           Path("./AO2019-2/neptune_0_bdf.fits")]):
    data = fits.getdata(fpath)
    vmin, vmax = np.percentile(data, [10, 90])
    axs[i].imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    axs[i].set(title=fpath)
    print(data)
plt.tight_layout()
plt.show()

# printing results.
# I don't know why, but neptune-0005 file comes at first place
#%%

# now each preprocessing for uranus and standard stars and comprasion lamp.
# dark exptime matches objt exptime

for i in range (0, 9):
    hdul = fits.open(fpaths['objt']['uranus'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    data = CCDData(data, unit=u.adu)
    data = trim(data, fits_section='[900:1100, 150:1750]')
    try: 
        _ = hdr["PROCESS"]
    except KeyError:
        hdr["PROCESS"] = ''
    
    data_reduced = data - master_bias.data
    hdr["PROCESS"] += "B"
    hdr.add_history(f"Bias subtracted by master_bias")
              
    hdu = fits.PrimaryHDU(data=data_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"uranus_0_bxx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"uranus_1_bxx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"uranus_2_bxx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"uranus_3_bxx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"uranus_4_bxx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"uranus_5_bxx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"uranus_6_bxx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"uranus_7_bxx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"uranus_8_bxx.fits", overwrite=True)
        
    data_d_reduced = data_reduced - master_dark_10
    hdr["PROCESS"] += "D"
    hdr.add_history(f"Dark subtracted by master_dark_10")
    
    hdu = fits.PrimaryHDU(data=data_d_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"uranus_0_bdx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"uranus_1_bdx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"uranus_2_bdx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"uranus_3_bdx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"uranus_4_bdx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"uranus_5_bdx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"uranus_6_bdx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"uranus_7_bdx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"uranus_8_bdx.fits", overwrite=True)
        
    data_corrected = data_d_reduced/master_flat
    hdr["PROCESS"] += "F"
    hdr.add_history(f"Flat corrected by master_flat")
    
    hdu = fits.PrimaryHDU(data=data_corrected, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"uranus_0_bdf.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"uranus_1_bdf.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"urnaus_2_bdf.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"uranus_3_bdf.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"uranus_4_bdf.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"uranus_5_bdf.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"uranus_6_bdf.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"uranus_7_bdf.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"uranus_8_bdf.fits", overwrite=True)
        
#%%
fpaths
#%%
        

fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharex=False, sharey=False, gridspec_kw=None)

for i, fpath in enumerate([Path("./AO2019-2/2019-10-26/Uranus-0004.fit"),
                           Path("./AO2019-2/uranus_0_bxx.fits"),
                           Path("./AO2019-2/uranus_0_bdx.fits"),
                           Path("./AO2019-2/uranus_0_bdf.fits")]):
    data = fits.getdata(fpath)
    vmin, vmax = np.percentile(data, [10, 90])
    axs[i].imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    axs[i].set(title=fpath)
    print(data)
plt.tight_layout()
plt.show()

#%%

for i in range (0, 9):
    hdul = fits.open(fpaths['objt']['uranus20'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    data = CCDData(data, unit=u.adu)
    data = trim(data, fits_section='[900:1100, 150:1750]')
    try: 
        _ = hdr["PROCESS"]
    except KeyError:
        hdr["PROCESS"] = ''
    
    data_reduced = data - master_bias.data
    hdr["PROCESS"] += "B"
    hdr.add_history(f"Bias subtracted by master_bias")
              
    hdu = fits.PrimaryHDU(data=data_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"uranus20_0_bxx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"uranus20_1_bxx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"uranus20_2_bxx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"uranus20_3_bxx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"uranus20_4_bxx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"uranus20_5_bxx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"uranus20_6_bxx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"uranus20_7_bxx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"uranus20_8_bxx.fits", overwrite=True)
        
    data_d_reduced = data_reduced - master_dark_20
    hdr["PROCESS"] += "D"
    hdr.add_history(f"Dark subtracted by master_dark_20")
    
    hdu = fits.PrimaryHDU(data=data_d_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"uranus20_0_bdx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"uranus20_1_bdx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"uranus20_2_bdx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"uranus20_3_bdx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"uranus20_4_bdx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"uranus20_5_bdx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"uranus20_6_bdx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"uranus20_7_bdx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"uranus20_8_bdx.fits", overwrite=True)
        
    data_corrected = data_d_reduced/master_flat
    hdr["PROCESS"] += "F"
    hdr.add_history(f"Flat corrected by master_flat")
    
    hdu = fits.PrimaryHDU(data=data_corrected, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"uranus20_0_bdf.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"uranus20_1_bdf.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"urnaus20_2_bdf.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"uranus20_3_bdf.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"uranus20_4_bdf.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"uranus20_5_bdf.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"uranus20_6_bdf.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"uranus20_7_bdf.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"uranus20_8_bdf.fits", overwrite=True)
        
#%%
fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharex=False, sharey=False, gridspec_kw=None)

for i, fpath in enumerate([Path("./AO2019-2/2019-10-26/Uranus20-0003.fit"),
                           Path("./AO2019-2/uranus20_0_bxx.fits"),
                           Path("./AO2019-2/uranus20_0_bdx.fits"),
                           Path("./AO2019-2/uranus20_0_bdf.fits")]):
    data = fits.getdata(fpath)
    vmin, vmax = np.percentile(data, [10, 90])
    axs[i].imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    axs[i].set(title=fpath)
    print(data)
plt.tight_layout()
plt.show()

#%%

for i in range (0, 9):
    hdul = fits.open(fpaths['objt']['hr9087'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    data = CCDData(data, unit=u.adu)
    data = trim(data, fits_section='[900:1100, 150:1750]')
    try: 
        _ = hdr["PROCESS"]
    except KeyError:
        hdr["PROCESS"] = ''
    
    data_reduced = data - master_bias.data
    hdr["PROCESS"] += "B"
    hdr.add_history(f"Bias subtracted by master_bias")
              
    hdu = fits.PrimaryHDU(data=data_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"hr9087_0_bxx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"hr9087_1_bxx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"hr9087_2_bxx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"hr9087_3_bxx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"hr9087_4_bxx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"hr9087_5_bxx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"hr9087_6_bxx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"hr9087_7_bxx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"hr9087_8_bxx.fits", overwrite=True)
        
    data_d_reduced = data_reduced - master_dark_10
    hdr["PROCESS"] += "D"
    hdr.add_history(f"Dark subtracted by master_dark_10")
    
    hdu = fits.PrimaryHDU(data=data_d_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"hr9087_0_bdx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"hr9087_1_bdx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"hr9087_2_bdx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"hr9087_3_bdx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"hr9087_4_bdx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"hr9087_5_bdx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"hr9087_6_bdx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"hr9087_7_bdx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"hr9087_8_bdx.fits", overwrite=True)
        
    data_corrected = data_d_reduced/master_flat
    hdr["PROCESS"] += "F"
    hdr.add_history(f"Flat corrected by master_flat")
    
    hdu = fits.PrimaryHDU(data=data_corrected, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"hr9087_0_bdf.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"hr9087_1_bdf.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"hr9087_2_bdf.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"hr9087_3_bdf.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"hr9087_4_bdf.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"hr9087_5_bdf.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"hr9087_6_bdf.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"hr9087_7_bdf.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"hr9087_8_bdf.fits", overwrite=True)
        

#%%
fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharex=False, sharey=False, gridspec_kw=None)

for i, fpath in enumerate([Path("./AO2019-2/2019-10-26/HR9087-0001.fit"),
                           Path("./AO2019-2/hr9087_0_bxx.fits"),
                           Path("./AO2019-2/hr9087_0_bdx.fits"),
                           Path("./AO2019-2/hr9087_0_bdf.fits")]):
    data = fits.getdata(fpath)
    vmin, vmax = np.percentile(data, [10, 90])
    axs[i].imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    axs[i].set(title=fpath)
    print(data)
plt.tight_layout()
plt.show()


#%%

for i in range (0, 9):
    hdul = fits.open(fpaths['objt']['hr718'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    data = CCDData(data, unit=u.adu)
    data = trim(data, fits_section='[900:1100, 150:1750]')
    try: 
        _ = hdr["PROCESS"]
    except KeyError:
        hdr["PROCESS"] = ''
    
    data_reduced = data - master_bias.data
    hdr["PROCESS"] += "B"
    hdr.add_history(f"Bias subtracted by master_bias")
              
    hdu = fits.PrimaryHDU(data=data_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"hr718_0_bxx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"hr718_1_bxx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"hr718_2_bxx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"hr718_3_bxx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"hr718_4_bxx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"hr718_5_bxx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"hr718_6_bxx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"hr718_7_bxx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"hr718_8_bxx.fits", overwrite=True)
        
    data_d_reduced = data_reduced - master_dark_5
    hdr["PROCESS"] += "D"
    hdr.add_history(f"Dark subtracted by master_dark_5")
    
    hdu = fits.PrimaryHDU(data=data_d_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"hr718_0_bdx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"hr718_1_bdx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"hr718_2_bdx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"hr718_3_bdx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"hr718_4_bdx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"hr718_5_bdx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"hr718_6_bdx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"hr718_7_bdx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"hr718_8_bdx.fits", overwrite=True)
        
    data_corrected = data_d_reduced/master_flat
    hdr["PROCESS"] += "F"
    hdr.add_history(f"Flat corrected by master_flat")
    
    hdu = fits.PrimaryHDU(data=data_corrected, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"hr718_0_bdf.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"hr718_1_bdf.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"hr718_2_bdf.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"hr718_3_bdf.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"hr718_4_bdf.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"hr718_5_bdf.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"hr718_6_bdf.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"hr718_7_bdf.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"hr718_8_bdf.fits", overwrite=True)
        

#%%
fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharex=False, sharey=False, gridspec_kw=None)

for i, fpath in enumerate([Path("./AO2019-2/2019-10-26/HR718-0001.fit"),
                           Path("./AO2019-2/hr718_0_bxx.fits"),
                           Path("./AO2019-2/hr718_0_bdx.fits"),
                           Path("./AO2019-2/hr718_0_bdf.fits")]):
    data = fits.getdata(fpath)
    vmin, vmax = np.percentile(data, [10, 90])
    axs[i].imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    axs[i].set(title=fpath)
    print(data)
plt.tight_layout()
plt.show()

#%%
for i in range (0, 9):
    hdul = fits.open(fpaths['comp'][i])
    hdr = hdul[0].header
    data = hdul[0].data
    data = CCDData(data, unit=u.adu)
    data = trim(data, fits_section='[900:1100, 150:1750]')
    try: 
        _ = hdr["PROCESS"]
    except KeyError:
        hdr["PROCESS"] = ''
    
    data_reduced = data - master_bias.data
    hdr["PROCESS"] += "B"
    hdr.add_history(f"Bias subtracted by master_bias")
              
    hdu = fits.PrimaryHDU(data=data_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"comp_0_bxx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"comp_1_bxx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"comp_2_bxx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"comp_3_bxx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"comp_4_bxx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"comp_5_bxx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"comp_6_bxx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"comp_7_bxx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"comp_8_bxx.fits", overwrite=True)
        
    data_d_reduced = data_reduced - master_dark_10
    hdr["PROCESS"] += "D"
    hdr.add_history(f"Dark subtracted by master_dark_10")
    
    hdu = fits.PrimaryHDU(data=data_d_reduced, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"comp_0_bdx.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"comp_1_bdx.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"comp_2_bdx.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"comp_3_bdx.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"comp_4_bdx.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"comp_5_bdx.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"comp_6_bdx.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"comp_7_bdx.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"comp_8_bdx.fits", overwrite=True)
        
    data_corrected = data_d_reduced/master_flat
    hdr["PROCESS"] += "F"
    hdr.add_history(f"Flat corrected by master_flat")
    
    hdu = fits.PrimaryHDU(data=data_corrected, header=hdr)
    if i ==0:
        hdu.writeto(Path("./AO2019-2")/"comp_0_bdf.fits", overwrite=True)
    if i ==1:
        hdu.writeto(Path("./AO2019-2")/"comp_1_bdf.fits", overwrite=True)
    if i ==2:
        hdu.writeto(Path("./AO2019-2")/"comp_2_bdf.fits", overwrite=True)
    if i ==3:
        hdu.writeto(Path("./AO2019-2")/"comp_3_bdf.fits", overwrite=True)
    if i ==4:
        hdu.writeto(Path("./AO2019-2")/"comp_4_bdf.fits", overwrite=True)
    if i ==5:
        hdu.writeto(Path("./AO2019-2")/"comp_5_bdf.fits", overwrite=True)
    if i ==6:
        hdu.writeto(Path("./AO2019-2")/"comp_6_bdf.fits", overwrite=True)
    if i ==7:
        hdu.writeto(Path("./AO2019-2")/"comp_7_bdf.fits", overwrite=True)
    if i ==8:
        hdu.writeto(Path("./AO2019-2")/"comp_8_bdf.fits", overwrite=True)


#%%
fig, axs = plt.subplots(4, 1, figsize=(40, 40), sharex=False, sharey=False, gridspec_kw=None)

for i, fpath in enumerate([Path("./AO2019-2/2019-10-26/Comp-0008.fit"),
                           Path("./AO2019-2/comp_0_bxx.fits"),
                           Path("./AO2019-2/comp_0_bdx.fits"),
                           Path("./AO2019-2/comp_0_bdf.fits")]):
    data = fits.getdata(fpath)
    vmin, vmax = np.percentile(data, [10, 90])
    axs[i].imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    axs[i].set(title=fpath)
    print(data)
plt.tight_layout()
plt.show()

#%%

#adding images of master bias, dark and flat
fig, axs = plt.subplots(1, 6, figsize=(50, 10), sharex=False, sharey=False, gridspec_kw=None)

for i, fpath in enumerate([Path("./AO2019-2/image_bias.fits"),
                           Path("./AO2019-2/image_dark_5.fits"),
                           Path("./AO2019-2/image_dark_10.fits"),
                           Path("./AO2019-2/image_dark_20.fits"),
                           Path("./AO2019-2/image_dark_60.fits"),
                           Path("./AO2019-2/image_flat.fits")]):
    data = fits.getdata(fpath)
    vmin, vmax = np.percentile(data, [10, 90])
    axs[i].imshow(data, origin='lower', vmin=vmin, vmax=vmax)
    axs[i].set(title=fpath)
    print(data)
plt.tight_layout()
plt.show()
