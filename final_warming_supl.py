
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
import cmocean
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import glob
import os

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


mylist = mylist-mylist_control
mylist = mylist*test
weights=np.cos(lat* np.pi / 180.)*test
mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
#mylist_to=(mylist*weights).sel(lat=slice(-90,90)).sum('lat',skipna=True) / weights.sel(lat=slice(-90,90)).sum('lat',skipna=True)

mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_ga=mylist_to.mean('lon',skipna=True)

mylist=mylist-mylist_ga
mylist=xr.Dataset.to_array(mylist)
mylist_cat1=mylist.sel(variable='ts').mean('new_dim')

mylist_cat1_qq=mylist.sel(variable='ts',year=slice(140,149)).mean('year').quantile([0.2,0.8], dim="new_dim")
mylist_cat1_qq=mylist_cat1_qq.rename({'quantile':'ts'})

#%%

mylist_1pct=xr.open_dataset('/Volumes/Armor_CMIP6/1pct_ts.nc')


mylist_control = xr.open_dataset('/Volumes/Armor_CMIP6/control_timemean_ts_1deg.nc')


lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)


mylist = mylist_1pct-mylist_control
mylist = mylist*test
weights=np.cos(lat* np.pi / 180.)*test
mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
#mylist_to=(mylist*weights).sel(lat=slice(-90,90)).sum('lat',skipna=True) / weights.sel(lat=slice(-90,90)).sum('lat',skipna=True)

mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_ga=mylist_to.mean('lon',skipna=True)

mylist=mylist-mylist_ga
mylist=xr.Dataset.to_array(mylist)
mylist_cat1_1pct=mylist.sel(variable='ts').mean('new_dim')

mylist_cat1_qq_1pct=mylist.sel(variable='ts',year=slice(140,149)).mean('year').quantile([0.2,0.8], dim="new_dim")
mylist_cat1_qq_1pct=mylist_cat1_qq_1pct.rename({'quantile':'ts'})

#%%

filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_ssp370*.nc'))
filelist=sorted(filelist,key=str.lower)
ts_control=mylist_control['ts']
control_ssp370subset=ts_control.isel(new_dim=[0,1,2,3,4,5,7,9,12,13,14,15,18,20,21,25,26,27,28,29,30,31,32,34,39])
control_ssp370subset=control_ssp370subset.assign_coords(new_dim=range(0,len(filelist)))


mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist=mylist['ts']

mylist=mylist.mean('ens_member')
mylist_c=control_ssp370subset.sel(new_dim=0)
mylist_c=mylist_c.reset_coords(drop=True)
mylist=mylist-mylist_c
mylist=mylist*test
lat=mylist.lat
weights=np.cos(lat* np.pi / 180.)*test

mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_ga=mylist_to.mean('lon',skipna=True)
mylist=mylist-mylist_ga
#mylist=xr.Dataset.to_array(mylist)
mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist=mylist['ts']
    mylist=mylist.mean('ens_member')
    mylist_c=control_ssp370subset.isel(new_dim=x)
    mylist_c=mylist_c.reset_coords(drop=True)

    mylist=mylist-mylist_c
    mylist=mylist*test
    lat=mylist.lat
    weights=np.cos(lat* np.pi / 180.)*test

    mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
    mylist_to=mylist_to.sel(lon=slice(60,300))
    mylist_ga=mylist_to.mean('lon',skipna=True)
    mylist=mylist-mylist_ga

#    mylist=xr.Dataset.to_array(mylist)



    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat1_ssp3=mylist_out.mean('new_dim')


mylist_cat1_qq_ssp3=mylist_out.sel(year=slice(2080,2100)).mean('year').quantile([0.2,0.8], dim="new_dim")
mylist_cat1_qq_ssp3=mylist_cat1_qq_ssp3.rename({'quantile':'ts'})

#%%
filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_ssp585*.nc'))
filelist=sorted(filelist,key=str.lower)
ts_control=mylist_control['ts']
control_ssp585subset=ts_control.isel(new_dim=[0,1,2,4,5,7,9,11,12,13,14,15,17,18,19,20,21,23,24,25,26,27,28,29,30,31,33,34,35,38,39])
control_ssp585subset=control_ssp585subset.assign_coords(new_dim=range(0,len(filelist)))

mylist=xr.open_dataset(filelist[x],decode_cf=False)
mylist=mylist['ts']
mylist=mylist.mean('ens_member')
mylist_c=control_ssp585subset.isel(new_dim=0)
mylist_c=mylist_c.reset_coords(drop=True)
mylist=mylist-mylist_c
mylist=mylist*test
lat=mylist.lat
weights=np.cos(lat* np.pi / 180.)*test

mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_ga=mylist_to.mean('lon',skipna=True)
mylist=mylist-mylist_ga
#mylist=xr.Dataset.to_array(mylist)
mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist=mylist['ts']
    mylist=mylist.mean('ens_member')
    mylist_c=control_ssp585subset.isel(new_dim=x)
    mylist_c=mylist_c.reset_coords(drop=True)
    mylist=mylist-mylist_c
    mylist=mylist*test
    lat=mylist.lat
    weights=np.cos(lat* np.pi / 180.)*test

    mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
    mylist_to=mylist_to.sel(lon=slice(60,300))
    mylist_ga=mylist_to.mean('lon',skipna=True)
    mylist=mylist-mylist_ga
 #   mylist=mylist-mylist_ga
#    mylist=xr.Dataset.to_array(mylist)



    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

mylist_cat1_ssp5=mylist_out.mean('new_dim')


mylist_cat1_qq_ssp5=mylist_out.sel(year=slice(2080,2100)).mean('year').quantile([0.2,0.8], dim="new_dim")
mylist_cat1_qq_ssp5=mylist_cat1_qq_ssp5.rename({'quantile':'ts'})

#%%
plot0 = mylist_cat1.isel(year=slice(140,149)).sel(lev=0).mean('year')
mylist_cat1_pos=plot0.where(plot0>0)+mylist_cat1_qq[0].where(mylist_cat1_qq[0]>0)
mylist_cat1_neg=plot0.where(plot0<0)+mylist_cat1_qq[1].where(mylist_cat1_qq[1]<0)


plot1 = mylist_cat1_1pct.isel(year=slice(140,149)).sel(lev=0).mean('year')
mylist_cat1_pos_pct=plot1.where(plot1>0)+mylist_cat1_qq_1pct[0].where(mylist_cat1_qq_1pct[0]>0)
mylist_cat1_neg_pct=plot1.where(plot1<0)+mylist_cat1_qq_1pct[1].where(mylist_cat1_qq_1pct[1]<0)

#plot2 = mylist_cat2.isel(year=slice(0,9)).sel(lev=0).mean('year')
#plot3 = mylist_cat2.isel(year=slice(140,149)).sel(lev=0).mean('year')
plot3 = mylist_cat1_ssp3.sel(year=slice(2080,2100)).isel(lev=0).mean('year')
mylist_cat1_pos_ssp3=plot1.where(plot3>0)+mylist_cat1_qq_ssp3[0].where(mylist_cat1_qq_ssp3[0]>0)
mylist_cat1_neg_ssp3=plot1.where(plot3<0)+mylist_cat1_qq_ssp3[1].where(mylist_cat1_qq_ssp3[1]<0)

plot4 = mylist_cat1_ssp5.sel(year=slice(2080,2100)).isel(lev=0).mean('year')
mylist_cat1_pos_ssp5=plot1.where(plot4>0)+mylist_cat1_qq_ssp5[0].where(mylist_cat1_qq_ssp5[0]>0)
mylist_cat1_neg_ssp5=plot1.where(plot4<0)+mylist_cat1_qq_ssp5[1].where(mylist_cat1_qq_ssp5[1]<0)



lon=plot0.lon
lat=plot0.lat

cmap = cmocean.cm.balance


plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'hatch.color': 'grey'})    
#%%
import matplotlib.patches as mpatches
k1=-1.6
k2=1.6

levels = np.arange(k1, k2, 0.1)

fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(25, 10),subplot_kw={'projection': ccrs.PlateCarree(180)})
plt.figtext(0.12, 0.86, 'a)')
plt.figtext(0.12, 0.45, 'c)')
plt.figtext(0.55, 0.85, 'b)')
plt.figtext(0.55, 0.45, 'd)')

axlist = axarr.flatten()

cf1=axlist[0].contourf(lon, lat, plot0,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
cs = axlist[0].contourf(lon, lat, mylist_cat1_pos.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])
cs = axlist[0].contourf(lon, lat, mylist_cat1_neg.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])           
  
axlist[0].set_title('4xCO$_2$: SST pattern, years 140:150', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[0], orientation='vertical', shrink=0.8, pad=0.02)
cb1.set_label('$\Delta$SST ($^o$C)')

axlist[0].set_extent([-120, 120, -41, 41], ccrs.PlateCarree(180))
axlist[0].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[0].add_feature(cfeature.LAND, zorder=100, edgecolor='k')

cr=    axlist[0].add_patch(mpatches.Rectangle(xy=[-100, -5], width=70, height=10,
                                    facecolor='none', edgecolor='black',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(180))
                 )    

cr=    axlist[0].add_patch(mpatches.Rectangle(xy=[0, -5], width=100, height=10,
                                    facecolor='none', edgecolor='black',
                                    linewidth=2,
                                    transform=ccrs.PlateCarree(180))
                 )  


gl = axlist[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

cf1=axlist[1].contourf(lon, lat, plot1,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
cs = axlist[1].contourf(lon, lat, mylist_cat1_neg_pct.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])  
cs = axlist[1].contourf(lon, lat, mylist_cat1_pos_pct.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.']) 
axlist[1].set_title('1pct CO$_2$: SST pattern, years 140:150', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[1], orientation='vertical', shrink=0.8, pad=0.02)
cb1.set_label('$\Delta$SST ($^o$C)')
axlist[1].set_extent([-120, 120, -41, 41], ccrs.PlateCarree(180))
axlist[1].coastlines()
#plt.colorbar()
#plt.clim(-0.7,0.7)
axlist[1].add_feature(cfeature.LAND, zorder=100, edgecolor='k')

gl = axlist[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlines = False
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-20,0,20,40])
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

cf1=axlist[2].contourf(lon, lat, plot3,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
cs = axlist[2].contourf(lon, lat, mylist_cat1_neg_ssp3.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])  
cs = axlist[2].contourf(lon, lat, mylist_cat1_pos_ssp3.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])  


axlist[2].set_title('ssp370: SST pattern, years 2080:2100', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[2], orientation='vertical', shrink=0.8, pad=0.02)
cb1.set_label('$\Delta$SST ($^o$C)')
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
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

cf1=axlist[3].contourf(lon, lat, plot4,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
cs = axlist[3].contourf(lon, lat, mylist_cat1_neg_ssp5.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])  
cs = axlist[3].contourf(lon, lat, mylist_cat1_pos_ssp5.isel(lev=0),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.']) 
axlist[3].set_title('ssp585: SST pattern, years 2080:2100', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[3], orientation='vertical', shrink=0.8, pad=0.02)
cb1.set_label('$\Delta$SST ($^o$C)')
axlist[3].set_extent([-120, 120, -41,41], ccrs.PlateCarree(180))
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