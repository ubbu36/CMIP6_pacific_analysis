
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:50:34 2020

@author: ullaheede
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
from mpl_toolkits.basemap import Basemap
import cartopy.feature as cfeature
from matplotlib import colorbar, colors
from matplotlib.cm import get_cmap
from matplotlib import colorbar, colors
import glob as glob
import os
import cmocean
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from pylab import *
import matplotlib.gridspec as gridspec

filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'psl_4xCO2*.nc'))
filelist=sorted(filelist,key=str.lower)

filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'psl_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist=mylist-mylist_c
mylist=mylist.sel(year=slice(0,9)).mean('year')
lat=mylist.lat
weights=np.cos(lat* np.pi / 180.)
mylist_to=(mylist*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,300))

mylist_ga=mylist_to.mean('lon',skipna=True)
mylist=mylist-mylist_ga
mylist=xr.Dataset.to_array(mylist)
mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist_c=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist=mylist-mylist_c
    mylist=mylist.sel(year=slice(0,9)).mean('year')
    lat=mylist.lat
    weights=np.cos(lat* np.pi / 180.)
    mylist_to=(mylist*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
    mylist_to=mylist_to.sel(lon=slice(60,300))

    mylist_ga=mylist_to.mean('lon',skipna=True)
    mylist=mylist-mylist_ga
    mylist=xr.Dataset.to_array(mylist)



    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat1=mylist_out.sel(variable='psl',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3=mylist_out.sel(variable='psl',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')


mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist=mylist-mylist_c
mylist=mylist.sel(year=slice(140,149)).mean('year')
lat=mylist.lat
weights=np.cos(lat* np.pi / 180.)
mylist_to=(mylist*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,300))

mylist_ga=mylist_to.mean('lon',skipna=True)
mylist=mylist-mylist_ga
mylist=xr.Dataset.to_array(mylist)
mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist_c=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist=mylist-mylist_c
    mylist=mylist.sel(year=slice(140,149)).mean('year')
    lat=mylist.lat
    weights=np.cos(lat* np.pi / 180.)
    mylist_to=(mylist*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
    mylist_to=mylist_to.sel(lon=slice(60,300))
#    mylist_to=(mylist*weights).sel(lat=slice(-90,90)).sum('lat',skipna=True) / weights.sel(lat=slice(-90,90)).sum('lat',skipna=True)
#    mylist_to=mylist_to.sel(lon=slice(60,300))
    mylist_ga=mylist_to.mean('lon',skipna=True)
    mylist=mylist-mylist_ga
    mylist=xr.Dataset.to_array(mylist)



    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat11=mylist_out.sel(variable='psl',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat33=mylist_out.sel(variable='psl',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')




lon=mylist.lon
lat=mylist.lat
#%%
filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'uas_4xCO2*.nc'))
filelist=sorted(filelist,key=str.lower)

filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'uas_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist=xr.Dataset.to_array(mylist)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist_c=xr.Dataset.to_array(mylist_c)
mylist=mylist.sel(year=slice(0,1)).mean('year')


mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist=xr.Dataset.to_array(mylist)
    mylist_c=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist_c=xr.Dataset.to_array(mylist_c)
    mylist=mylist-mylist_c
    mylist=mylist.sel(year=slice(0,1)).mean('year')




    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat1_uas=mylist_out.sel(variable='uas',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_uas=mylist_out.sel(variable='uas',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')


mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist=xr.Dataset.to_array(mylist)
mylist_c=xr.Dataset.to_array(mylist_c)
mylist=mylist-mylist_c
mylist=mylist.sel(year=slice(140,149)).mean('year')


mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist=xr.Dataset.to_array(mylist)
    mylist_c=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist_c=xr.Dataset.to_array(mylist_c)
    mylist=mylist-mylist_c
    mylist=mylist.sel(year=slice(140,149)).mean('year')




    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat11_uas=mylist_out.sel(variable='uas',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat33_uas=mylist_out.sel(variable='uas',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')


#%%

filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'vas_4xCO2*.nc'))
filelist=sorted(filelist,key=str.lower)

filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'vas_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist=xr.Dataset.to_array(mylist)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist_c=xr.Dataset.to_array(mylist_c)
mylist=mylist-mylist_c
mylist=mylist.sel(year=slice(0,1)).mean('year')


mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist=xr.Dataset.to_array(mylist)
    mylist_c=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist_c=xr.Dataset.to_array(mylist_c)
    mylist=mylist-mylist_c
    mylist=mylist.sel(year=slice(0,1)).mean('year')



    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat1_vas=mylist_out.sel(variable='vas',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_vas=mylist_out.sel(variable='vas',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')


mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist=xr.Dataset.to_array(mylist)
mylist_c=xr.Dataset.to_array(mylist_c)
mylist=mylist-mylist_c
mylist=mylist.sel(year=slice(140,149)).mean('year')


mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist=xr.Dataset.to_array(mylist)
    mylist_c=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist_c=xr.Dataset.to_array(mylist_c)
    mylist=mylist-mylist_c
    mylist=mylist.sel(year=slice(140,149)).mean('year')




    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat11_vas=mylist_out.sel(variable='vas',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat33_vas=mylist_out.sel(variable='vas',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')



lon=mylist.lon
lat=mylist.lat

lon_q=np.linspace(-180,178,num=359)
lat_q=mylist_cat1_vas.lat

#%%
#cmap = cmocean.cm.curl
cmap='BrBG_r'
k1=-150
k2=150
levels = np.arange(k1, k2, 1)

plt.rcParams.update({'font.size': 22})


fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(25, 8),subplot_kw={'projection': ccrs.PlateCarree(180)})

axlist = axarr.flatten()


plt.figtext(0.12, 0.86, 'a)')
plt.figtext(0.12, 0.45, 'c)')
plt.figtext(0.55, 0.85, 'b)')
plt.figtext(0.55, 0.45, 'd)')



cf1=axlist[0].contourf(lon, lat, mylist_cat1,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
q1 = axlist[0].quiver(lon_q[::5],lat_q[::5],mylist_cat1_uas[::5,::5],mylist_cat1_vas[::5,::5],transform=ccrs.PlateCarree(180),scale=23,headwidth=4)

axlist[0].set_title('OT category: SLP pattern, years 0:10', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[0], orientation='vertical', shrink=0.8, pad=0.05, ticks=[-100, -50, 0, 50, 100])
cb1.set_label('$\Delta$SLP (Pa)',fontsize=19)

axlist[0].set_extent([-120, 120, -35, 35], ccrs.PlateCarree(180))
axlist[0].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[0].add_feature(cfeature.COASTLINE, zorder=100, edgecolor='k')
gl = axlist[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-30,-15,0,15,30])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}


cf1=axlist[1].contourf(lon, lat, mylist_cat11,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
q2 = axlist[1].quiver(lon_q[::5],lat_q[::5],mylist_cat11_uas[::5,::5],mylist_cat11_vas[::5,::5],transform=ccrs.PlateCarree(180),scale=48,headwidth=4)

axlist[1].set_title('OT category: SLP pattern, years 140:150', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[1], orientation='vertical', shrink=0.8, pad=0.05, ticks=[-100, -50, 0, 50, 100])
cb1.set_label('$\Delta$SLP (Pa)',fontsize=19)
axlist[1].set_extent([-120, 120, -35, 35], ccrs.PlateCarree(180))
axlist[1].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[1].add_feature(cfeature.COASTLINE, zorder=100, edgecolor='k')
gl = axlist[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-30,-15,0,15,30])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}


cf1=axlist[2].contourf(lon, lat, mylist_cat3,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
q3 = axlist[2].quiver(lon_q[::5],lat_q[::5],mylist_cat3_uas[::5,::5],mylist_cat3_vas[::5,::5],transform=ccrs.PlateCarree(180),scale=23,headwidth=4)

axlist[2].set_title('EP category: SLP pattern, years 0:10', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[2], orientation='vertical', shrink=0.8, pad=0.05, ticks=[-100, -50, 0, 50, 100])
cb1.set_label('$\Delta$SLP (Pa)')
axlist[2].set_extent([-120, 120, -35, 35], ccrs.PlateCarree(180))
axlist[2].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[2].add_feature(cfeature.COASTLINE, zorder=100, edgecolor='k')
gl = axlist[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-30,-15,0,15,30])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

cf1=axlist[3].contourf(lon, lat, mylist_cat33,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
q4 = axlist[3].quiver(lon_q[::5],lat_q[::5],mylist_cat33_uas[::5,::5],mylist_cat33_vas[::5,::5],transform=ccrs.PlateCarree(180),scale=48,headwidth=4)

axlist[3].set_title('EP category: SLP pattern, years 140:150', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[3], orientation='vertical', shrink=0.8, pad=0.05, ticks=[-100, -50, 0, 50, 100])
cb1.set_label('$\Delta$SLP (Pa)')
axlist[3].set_extent([-120, 120, -35,35], ccrs.PlateCarree(180))
axlist[3].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[3].add_feature(cfeature.COASTLINE, zorder=100, edgecolor='k')
gl = axlist[3].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-30,-15,0,15,30])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

qk = plt.quiverkey(q1, 0.39, 0.87, 1, r'$1  m s^{-1}$', labelpos='E',
                   coordinates='figure')

qk = plt.quiverkey(q2, 0.82, 0.87, 1, r'$1  m s^{-1}$', labelpos='E',
                   coordinates='figure')

qk = plt.quiverkey(q3, 0.39, 0.46, 1, r'$1  m s^{-1}$', labelpos='E',
                   coordinates='figure')

qk = plt.quiverkey(q4, 0.82, 0.46, 1, r'$1  m s^{-1}$', labelpos='E',
                   coordinates='figure')