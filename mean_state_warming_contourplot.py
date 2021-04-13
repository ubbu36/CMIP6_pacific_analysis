
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
from matplotlib.ticker import FuncFormatter

model_names=['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','BSS-ESM1','CanESM5','CESM2','CESM2-WACCM',\
             'CNRM-CM6','CNRM-ESM2-1','FGOALS-f3-L','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H',\
             'HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A','MCM-UA-1-0','MIROC-ES2L','MIROC6',\
                 'MRI-ESM2','NESM3','UKESM1-0-LL']



model_names_cat1=['CanESM5','CESM2','CESM2-WACCM','GFDL-CM4','GFDL-ESM4','GISS-E2-1-H','HadGEM3-GC3-LL','MIROC6']

model_names_cat2=['ACCESS-ESM1-5','FGOALS-f3-L','GISS-E2-1-G','MRI-ESM2']

model_names_cat3=['ACCESS-CM2','BCC-CSM2-MR','BCC-ESM1','CNRM-CM6','CNRM-ESM2-1','IPSL-CM6','NESM3','UKESM1-0-LL']

obs=xr.open_dataset('/Users/ullaklintheede/Downloads/ersst.v4.1854-2020.nc')

ds_out = xr.Dataset({'lat': (['lat'], np.arange(-88, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 359, 1)),
                    }
                   )
regridder = xe.Regridder(obs, ds_out, 'bilinear')
mylist_obs_regrid = regridder(obs)



mask_ocean = 1 * np.ones((mylist_obs_regrid.dims['lat'], mylist_obs_regrid.dims['lon'])) * np.isfinite(mylist_obs_regrid.sst.isel(time=0))
test=mask_ocean.where(mask_ocean == 1)
mask_ocean=mask_ocean.isel(lev=0)



mylist = xr.open_dataset('/Volumes/Armor_CMIP6/4xCO2_ts.nc')


mylist_control = xr.open_dataset('/Volumes/Armor_CMIP6/control_timemean_ts_1deg.nc')


lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)


mylist = mylist_control
mylist = mylist*test
weights=np.cos(lat* np.pi / 180.)*test
mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_ga=mylist_to.mean('lon',skipna=True)

#mylist=mylist-mylist_ga
mylist=mylist-273
mylist=xr.Dataset.to_array(mylist)
mylist_cat1=mylist.sel(variable='ts',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
mylist_cat1_q=(mylist.sel(variable='ts',new_dim=[0,13,14,15,32,33,36])).quantile([0.2,0.8], dim="new_dim")
mylist_cat1_q=mylist_cat1_q.rename({'quantile':'ts'})
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3=mylist.sel(variable='ts',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')
mylist_cat3_q=(mylist.sel(variable='ts',new_dim=[7,8,9,37,12,16,18,20,21,31])).quantile([0.2,0.8], dim="new_dim")
mylist_cat3_q=mylist_cat3_q.rename({'quantile':'ts'})
#%%


filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'psl_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist_out=mylist_c

for x in range(1,len(filelist_c)):
    mylist=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_out=xr.Dataset.to_array(mylist_out)

lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)


mylist_cat1_psl=mylist_out.sel(new_dim=[0,13,14,15,32,33,36]).mean('new_dim')

#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_psl=mylist_out.sel(new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')

#%%


filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'uas_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist_c=xr.Dataset.to_array(mylist_c)
mylist_out=mylist_c

for x in range(1,len(filelist_c)):
    mylist=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist=xr.Dataset.to_array(mylist)
 #   mylist=mylist.sel(variable='uas')
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')



lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)


mylist_cat1_uas=mylist_out.sel(new_dim=[0,13,14,15,32,33,36]).mean('new_dim')

#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_uas=mylist_out.sel(new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')

#%%


filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'vas_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist_c=xr.Dataset.to_array(mylist_c)
mylist_out=mylist_c

for x in range(1,len(filelist_c)):
    mylist=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist=xr.Dataset.to_array(mylist)
 #   mylist=mylist.sel(variable='uas')
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')



lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)


mylist_cat1_vas=mylist_out.sel(new_dim=[0,13,14,15,32,33,36]).mean('new_dim')

#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_vas=mylist_out.sel(new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')

#%%
plot0 = mylist_cat1.sel(lev=0)


plot1 = mylist_cat1_psl.sel(variable='psl')
#plot2 = mylist_cat2.isel(year=slice(0,9)).sel(lev=0).mean('year')
#plot3 = mylist_cat2.isel(year=slice(140,149)).sel(lev=0).mean('year')
plot4 = mylist_cat3.sel(lev=0)


plot5 = mylist_cat3_psl.sel(variable='psl')

lon=plot0.lon
lat=plot0.lat

lon1=plot1.lon
lat1=plot1.lat

lon_q=np.linspace(-180,178,num=359)
lat_q=mylist_cat3_vas.lat


cmap = cmocean.cm.balance
cmap1 = 'BrBG_r'
levels = np.arange(15, 35, 1)
levels1 = np.arange(100000, 103000, 50)

k1=15
k2=35

p1=100000
p2=103000

plt.rcParams.update({'font.size': 22})
fig, axarr = plt.subplots(nrows=3, ncols=2, figsize=(25, 15),subplot_kw={'projection': ccrs.PlateCarree(180)})
plt.figtext(0.125, 0.865, 'a)')
plt.figtext(0.55, 0.845, 'b)')
plt.figtext(0.125, 0.6, 'c)')
plt.figtext(0.55, 0.575, 'd)')
plt.figtext(0.125, 0.34, 'e)')
plt.figtext(0.55, 0.31, 'f)')
axlist = axarr.flatten()

cf1=axlist[0].contourf(lon, lat, plot0,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)

axlist[0].set_title('OT category: mean state SST', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[0], orientation='vertical', shrink=0.8, pad=0.02,ticks=[15, 20, 25, 30])
cb1.set_label('SST ($^o$C)')

axlist[0].set_extent([-120, 120, -41, 41], ccrs.PlateCarree(180))
axlist[0].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[0].add_feature(cfeature.LAND, zorder=100, edgecolor='k')
gl = axlist[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([80, 120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

cf1=axlist[1].contourf(lon1, lat1, plot1,levels1, extend="both",
             transform=ccrs.PlateCarree(),vmin=p1,vmax=p2, cmap=cmap1)
q1 = axlist[1].quiver(lon_q[::5],lat_q[::5],mylist_cat1_uas.sel(variable='uas')[::5,::5],mylist_cat1_vas.sel(variable='vas')[::5,::5],transform=ccrs.PlateCarree(180),scale=180,headwidth=4)
comma_fmt = FuncFormatter(lambda x, p: format(int(x), ','))
axlist[1].set_title('OT category: mean state SLP', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[1], orientation='vertical', shrink=0.6, pad=0.02, ticks=[100500,101000,101500,102000,102500],format=comma_fmt)
cb1.set_label('SLP (Pa)')
axlist[1].set_extent([-120, 120, -30, 30], ccrs.PlateCarree(180))
axlist[1].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[1].add_feature(cfeature.LAND, zorder=100, edgecolor='k')

gl = axlist[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60, 120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}
# cf1=axlist[2].contourf(lon, lat, plot2,levels, extend="both",
#              transform=ccrs.PlateCarree(),vmin=-1.6,vmax=1.6, cmap=get_cmap("coolwarm"))

# axlist[2].set_title('Category 2 warming pattern, years 1:10', fontsize=20)
# cb1 = fig.colorbar(cf1, ax=axlist[2], orientation='vertical', shrink=0.74, pad=0.05)

# axlist[2].set_extent([-100, 120, -40, 40], ccrs.PlateCarree(180))
# axlist[2].coastlines()
# #plt.colorbar()
# #plt.clim(-0.7,0.7)
# axlist[2].add_feature(cfeature.LAND, zorder=100, edgecolor='k')

# cf1=axlist[3].contourf(lon, lat, plot3,levels, extend="both",
#              transform=ccrs.PlateCarree(),vmin=-2,vmax=2, cmap=get_cmap("coolwarm"))

# axlist[3].set_title('Category 2 warming pattern, years 140:150', fontsize=20)
# cb1 = fig.colorbar(cf1, ax=axlist[3], orientation='vertical', shrink=0.74, pad=0.05)

# axlist[3].set_extent([-100, 120, -40, 40], ccrs.PlateCarree(180))
# axlist[3].coastlines()
# #plt.colorbar()
# #plt.clim(-0.7,0.7)
# axlist[3].add_feature(cfeature.LAND, zorder=100, edgecolor='k')

cf1=axlist[2].contourf(lon, lat, plot4,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)

axlist[2].set_title('EP category: mean state SST', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[2], orientation='vertical', shrink=0.8, pad=0.02,ticks=[15, 20, 25, 30])
cb1.set_label('SST ($^o$C)')
axlist[2].set_extent([-120, 120, -41, 41], ccrs.PlateCarree(180))
axlist[2].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[2].add_feature(cfeature.LAND, zorder=100, edgecolor='k')

gl = axlist[2].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([80, 120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

cf1=axlist[3].contourf(lon1, lat1, plot5,levels1, extend="both",
             transform=ccrs.PlateCarree(),vmin=p1,vmax=p2, cmap=cmap1)
q2= axlist[3].quiver(lon_q[::5],lat_q[::5],mylist_cat3_uas.sel(variable='uas')[::5,::5],mylist_cat3_vas.sel(variable='vas')[::5,::5],transform=ccrs.PlateCarree(180),scale=180,headwidth=4)

axlist[3].set_title('EP category: mean state SLP', fontsize=20)


cb1 = fig.colorbar(cf1, ax=axlist[3], orientation='vertical', shrink=0.6, pad=0.02, ticks=[100500,101000,101500,102000,102500],format=comma_fmt)

cb1.set_label('SLP (Pa)')


axlist[3].set_extent([-120, 120, -30,30], ccrs.PlateCarree(180))
axlist[3].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[3].add_feature(cfeature.LAND, zorder=100, edgecolor='k')
gl = axlist[3].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

levels = np.arange(-1, 1, 0.1)
levels1 = np.arange(-100, 100, 10)

k1=-1
k2=1

p1=-100
p2=100

plot6 = mylist_cat1.sel(lev=0)-mylist_cat3.sel(lev=0)
plot7 = mylist_cat1_psl.sel(variable='psl')-mylist_cat3_psl.sel(variable='psl')
uas_dif=mylist_cat1_uas-mylist_cat3_uas
vas_dif=mylist_cat1_vas-mylist_cat3_vas

cf1=axlist[4].contourf(lon, lat, plot6,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)

axlist[4].set_title('OT minus EP, $\Delta$SST', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[4], orientation='vertical', shrink=0.8, pad=0.04)
cb1.set_label('$\Delta$SST ($^o$C)')
axlist[4].set_extent([-120, 120, -41, 41], ccrs.PlateCarree(180))
axlist[4].coastlines()
#plt.coorbar()
#plt.clim(-0.7,0.7)
axlist[4].add_feature(cfeature.LAND, zorder=100, edgecolor='k')

gl = axlist[4].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([80, 120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

cf1=axlist[5].contourf(lon1, lat1, plot7,levels1, extend="both",
             transform=ccrs.PlateCarree(),vmin=p1,vmax=p2, cmap=cmap1)

axlist[5].set_title('OT minus EP, $\Delta$SLP', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[5], orientation='vertical', shrink=0.6, pad=0.02,ticks=[-80,-40,0,40,80])
cb1.set_label('$\Delta$SLP (Pa)')
q3 = axlist[5].quiver(lon_q[::5],lat_q[::5],uas_dif.sel(variable='uas')[::5,::5],vas_dif.sel(variable='vas')[::5,::5],transform=ccrs.PlateCarree(180),scale=50,headwidth=4)

axlist[5].set_extent([-120, 120, -30,30], ccrs.PlateCarree(180))
axlist[5].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[5].add_feature(cfeature.LAND, zorder=100, edgecolor='k')
gl = axlist[5].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60, 120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

qk = plt.quiverkey(q1, 0.80, 0.85, 5, r'$5  m s^{-1}$', labelpos='E',
                   coordinates='figure')

qk = plt.quiverkey(q2, 0.80, 0.58, 5, r'$5  m s^{-1}$', labelpos='E',
                   coordinates='figure')

qk = plt.quiverkey(q3, 0.79, 0.315, 1, r'$1  m s^{-1}$', labelpos='E',
                   coordinates='figure')

