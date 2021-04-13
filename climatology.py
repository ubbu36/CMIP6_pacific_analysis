
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
import matplotlib.patches as patches
import xarray.ufuncs as xu
from pylab import *
import matplotlib.gridspec as gridspec

model_names=['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','BSS-ESM1','CanESM5','CESM2','CESM2-WACCM',\
             'CNRM-CM6','CNRM-ESM2-1','FGOALS-f3-L','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H',\
             'HadGEM3-GC31-LL','INM-CM4-8','INM-CM5-0','IPSL-CM6A','MCM-UA-1-0','MIROC-ES2L','MIROC6',\
                 'MRI-ESM2','NESM3','UKESM1-0-LL']



model_names_cat1=['CanESM5','CESM2','CESM2-WACCM','GFDL-CM4','GFDL-ESM4','GISS-E2-1-H','HadGEM3-GC3-LL','MIROC6']

model_names_cat2=['ACCESS-ESM1-5','FGOALS-f3-L','GISS-E2-1-G','MRI-ESM2']

model_names_cat3=['ACCESS-CM2','BCC-CSM2-MR','BCC-ESM1','CNRM-CM6','CNRM-ESM2-1','IPSL-CM6','NESM3','UKESM1-0-LL']

obs=xr.open_dataset('/Users/ullaklintheede/Downloads/ersst.v4.1854-2020.nc')
obs=obs.groupby('time.year').mean('time')
obs=obs.sel(year=slice(1870,2000)).mean('year')

ds_out = xr.Dataset({'lat': (['lat'], np.arange(-88, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 359, 1)),
                    }
                   )
regridder = xe.Regridder(obs, ds_out, 'bilinear')
mylist_obs_regrid = regridder(obs)



mask_ocean = 1 * np.ones((mylist_obs_regrid.dims['lat'], mylist_obs_regrid.dims['lon'])) * np.isfinite(mylist_obs_regrid.sst)
test=mask_ocean.where(mask_ocean == 1)
mask_ocean=mask_ocean.isel(lev=0)
mylist_obs_regrid  = mylist_obs_regrid.sel(lev=0).rename({'sst': 'ts'})
mylist_obs_regrid=mylist_obs_regrid['ts']
#%%
mylist = xr.open_dataset('/Volumes/Armor_CMIP6/4xCO2_ts.nc')


mylist_control = xr.open_dataset('/Volumes/Armor_CMIP6/control_timemean_ts_1deg.nc')


lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)


mylist = mylist_control
#mylist = mylist*test
#weights=np.cos(lat* np.pi / 180.)*test
mylist_to=(mylist*weights).sel(lat=slice(-40,40)).sum('lat',skipna=True) / weights.sel(lat=slice(-40,40)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_ga=mylist_to.mean('lon',skipna=True)

#mylist=mylist-mylist_ga
mylist=mylist-273
mylist=xr.Dataset.to_array(mylist)
#mylist_cat1=mylist.sel(variable='ts',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
mylist_cat1=mylist.sel(variable='ts',new_dim=[0,13,14,15,32,33,36]).mean('new_dim')
mylist_cat1_diff=mylist_cat1-mylist_obs_regrid
mylist_cat1_q=(mylist.sel(variable='ts',new_dim=[0,13,14,15,32,33,36])-mylist_obs_regrid).quantile([0.2,0.8], dim="new_dim")
mylist_cat1_q=mylist_cat1_q.rename({'quantile':'ts'})

profile_cat1=mylist.sel(variable='ts',new_dim=[0,13,14,15,32,33,36],lat=(slice(-5,5))).mean('lat')
profile_cat1=(mylist.sel(variable='ts',new_dim=[0,13,14,15,32,33,36],lat=(slice(-5,5)))*test.sel()).mean('lat')

#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
#mylist_cat3=mylist.sel(variable='ts',new_dim=[1,4,5,6,7,8,9,12,16,26,27,29,37]).mean('new_dim')
mylist_cat3=mylist.sel(variable='ts',new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')
mylist_cat3_diff=mylist_cat3-mylist_obs_regrid
mylist_cat3_q=(mylist.sel(variable='ts',new_dim=[7,8,9,37,12,16,18,20,21,31])-mylist_obs_regrid).quantile([0.2,0.8], dim="new_dim")
mylist_cat3_q=mylist_cat3_q.rename({'quantile':'ts'})

profile_cat3=(mylist.sel(variable='ts',new_dim=[7,8,9,37,12,16,18,20,21,31],lat=(slice(-5,5)))*test.sel()).mean('lat')

observed_profile=mylist_obs_regrid.sel(lat=slice(-5,5)).mean('lat')

profile_cmip=(mylist.sel(variable='ts',lat=(slice(-5,5)))*test.sel()).mean('lat').mean('new_dim')
#%%

obs2=xr.open_dataset('/Users/ullaklintheede/Downloads/slp.mon.mean.nc')
obs2=obs2.groupby('time.year').mean('time')*100
obs2=obs2.sel(year=slice(1950,2000)).mean('year')
obs2=obs2.rename({'slp':'psl'})
obs2=obs2['psl']
#obs2=obs2.reindex(lat=list(reversed(obs.lat)))

ds_out = xr.Dataset({'lat': (['lat'], np.arange(-88, 87, 1.0)),
                     'lon': (['lon'], np.arange(0, 357, 1)),
                    }
                   )
regridder = xe.Regridder(obs2, ds_out, 'bilinear')
mylist_obs_regrid_psl = regridder(obs2)



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
mylist_cat1_psl_diff=mylist_cat1_psl-mylist_obs_regrid_psl
mylist_cat1_q_psl=(mylist_out.sel(variable='psl',new_dim=[0,13,14,15,32,33,36])-mylist_obs_regrid_psl).quantile([0.2,0.8], dim="new_dim")

mylist_cat1_q_psl=mylist_cat1_q_psl.rename({'quantile':'psl'})
#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_psl=mylist_out.sel(new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim')
mylist_cat3_psl_diff=mylist_cat3_psl-mylist_obs_regrid_psl

mylist_cat3_q_psl=(mylist_out.sel(variable='psl',new_dim=[7,8,9,37,12,16,18,20,21,31])-mylist_obs_regrid_psl).quantile([0.2,0.8], dim="new_dim")
mylist_cat3_q_psl=mylist_cat3_q_psl.rename({'quantile':'psl'})
#%%
obs3=xr.open_dataset('/Users/ullaklintheede/Downloads/uwnd.mon.mean.nc')
obs3=obs3.groupby('time.year').mean('time')
obs3=obs3.sel(year=slice(1950,2000)).mean('year')
obs3=obs3.rename({'uwnd':'uas'})
obs3=obs3['uas']
#obs2=obs2.reindex(lat=list(reversed(obs.lat)))

ds_out = xr.Dataset({'lat': (['lat'], np.arange(-90, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 357, 1)),
                    }
                   )
regridder = xe.Regridder(obs2, ds_out, 'bilinear')
mylist_obs_regrid_uas = regridder(obs3)


filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'uas_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist_c=xr.Dataset.to_array(mylist_c)
mylist_out=mylist_c

for x in range(1,len(filelist_c)):
    mylist=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist=mylist['uas']
 #   mylist=mylist.sel(variable='uas')
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')



lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)

mylist_uas=mylist_out

mylist_cat1_uas=mylist_out.sel(new_dim=[0,13,14,15,32,33,36]).mean('new_dim').sel(variable='uas')
mylist_cat1_uas.lon.astype(dtype='float64')

mylist_cat1_uas_diff=mylist_cat1_uas-mylist_obs_regrid_uas
mylist_cat1_uas_diff=mylist_cat1_uas_diff
profile_cat1_uas=(mylist_out.sel(new_dim=[0,13,14,15,32,33,36],lat=(slice(-5,5)))*test).mean('lat',skipna='False')

#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_uas=mylist_out.sel(new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim').sel(variable='uas')

mylist_cat3_uas_diff=mylist_cat3_uas-mylist_obs_regrid_uas
mylist_cat3_uas_diff=mylist_cat3_uas_diff
profile_cat3_uas=(mylist_out.sel(new_dim=[7,8,9,37,12,16,18,20,21,31],lat=(slice(-5,5)))*test).mean('lat',skipna='False')

profile_cat2_uas=mylist_out.sel(new_dim=[1,2,3,4,5,6,10,11,17,19,22,23,24,25,26,27,28,29,30,34,35,37,38,39],lat=(slice(-5,5))).mean('lat')


observed_profile_uas=(mylist_obs_regrid_uas.sel(lat=slice(-5,5))*test).mean('lat',skipna='False')

profile_cmip_uas=(mylist_out.sel(lat=(slice(-5,5)))*test).mean('lat',skipna='False').mean('new_dim')


#%%

obs4=xr.open_dataset('/Users/ullaklintheede/Downloads/vwnd.mon.mean.nc')
obs4=obs4.groupby('time.year').mean('time')
obs4=obs4.sel(year=slice(1950,2000)).mean('year')
obs4=obs4.rename({'vwnd':'vas'})
obs4=obs4['vas']
#obs2=obs2.reindex(lat=list(reversed(obs.lat)))

ds_out = xr.Dataset({'lat': (['lat'], np.arange(-90, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 357, 1)),
                    }
                   )
regridder = xe.Regridder(obs2, ds_out, 'bilinear')
mylist_obs_regrid_vas = regridder(obs4)


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

mylist_vas=mylist_out

lon=mylist.lon
lat=mylist.lat

weights=np.cos(lat* np.pi / 180.)


mylist_cat1_vas=mylist_out.sel(new_dim=[0,13,14,15,32,33,36]).mean('new_dim').sel(variable='vas')
mylist_cat1_vas_diff=mylist_cat1_vas-mylist_obs_regrid_vas

#mylist_cat2=mylist.sel(variable='ts',new_dim=[2,10,11,24,25,27,38,39]).mean('new_dim')
mylist_cat3_vas=mylist_out.sel(new_dim=[7,8,9,37,12,16,18,20,21,31]).mean('new_dim').sel(variable='vas')
mylist_cat3_vas_diff=mylist_cat3_vas-mylist_obs_regrid_vas

#%%
mylist_ws=0.0012*mylist_uas.sel(variable='uas')*xu.sqrt(mylist_uas.sel(variable='uas')**2+mylist_vas.sel(variable='vas')**2)
obs_ws=0.0012*mylist_obs_regrid_uas*xu.sqrt(mylist_obs_regrid_vas**2+mylist_obs_regrid_uas**2)

ws_profile_cat1=mylist_ws.sel(new_dim=[0,13,14,15,32,33,36]).sel(lat=slice(-5,5)).mean('lat')
ws_profile_cat3=mylist_ws.sel(new_dim=[7,8,9,37,12,16,18,20,21,31]).sel(lat=slice(-5,5)).mean('lat')
ws_profile_cmip=mylist_ws.mean('new_dim').sel(lat=slice(-5,5)).mean('lat')
ws_profile_obs=obs_ws.sel(lat=slice(-5,5)).mean('lat')

#%%
plot0 = mylist_cat1_diff
mylist_cat1_pos=plot0.where(plot0>0)+mylist_cat1_q[0].where(mylist_cat1_q[0]>0)
mylist_cat1_neg=plot0.where(plot0<0)+mylist_cat1_q[1].where(mylist_cat1_q[1]<0)

plot1 = mylist_cat1_psl_diff
mylist_cat1_pos_psl=plot1.where(plot1>0)+mylist_cat1_q_psl[0].where(mylist_cat1_q_psl[0]>0)
mylist_cat1_neg_psl=plot1.where(plot1<0)+mylist_cat1_q_psl[1].where(mylist_cat1_q_psl[1]<0)

#plot2 = mylist_cat2.isel(year=slice(0,9)).sel(lev=0).mean('year')
#plot3 = mylist_cat2.isel(year=slice(140,149)).sel(lev=0).mean('year')
plot4 = mylist_cat3_diff
mylist_cat3_pos=plot4.where(plot4>0)+mylist_cat3_q[0].where(mylist_cat3_q[0]>0)
mylist_cat3_neg=plot4.where(plot4<0)+mylist_cat3_q[1].where(mylist_cat3_q[1]<0)

plot5 = mylist_cat3_psl_diff
mylist_cat3_pos_psl=plot5.where(plot5>0)+mylist_cat3_q_psl[0].where(mylist_cat3_q_psl[0]>0)
mylist_cat3_neg_psl=plot5.where(plot5<0)+mylist_cat3_q_psl[1].where(mylist_cat3_q_psl[1]<0)

lon=plot0.lon
lat=plot0.lat

lon1=plot1.lon
lat1=plot1.lat

lon_q=np.linspace(-180,178,num=359)
lat_q=mylist_cat3_vas.lat


cmap = cmocean.cm.balance
cmap1 = 'BrBG_r'
levels = np.arange(-3, 3, 0.1)
levels1 = np.arange(-200, 200, 10)

k1=-3
k2=3

p1=-200
p2=200

plt.rcParams.update({'font.size': 22})
#%% root mean square error

test=xu.square(mylist_cat1_diff.sel(lon=slice(60, 280),lat=slice(-30, 30)))
rmse=xu.sqrt(test)
rmse1=np.around(rmse, decimals=2)

mylist_to=(rmse1*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,280))
mylist_ga_cat1=mylist_to.mean('lon',skipna=True)

test=xu.square(mylist_cat3_diff.sel(lon=slice(60, 280),lat=slice(-30, 30)))
rmse=xu.sqrt(test)
rmse1=np.around(rmse, decimals=2)

mylist_to=(rmse1*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,280))
mylist_ga_cat3=mylist_to.mean('lon',skipna=True)

test=xu.square(mylist_cat1_psl_diff.sel(lon=slice(60, 280),lat=slice(-30, 30)))
rmse=xu.sqrt(test)
rmse1=np.around(rmse, decimals=2)

mylist_to=(rmse1*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,280))
mylist_ga_cat1_psl=mylist_to.mean('lon',skipna=True)

test=xu.square(mylist_cat3_psl_diff.sel(lon=slice(60, 280),lat=slice(-30, 30)))
rmse=xu.sqrt(test)
rmse1=np.around(rmse, decimals=2)

mylist_to=(rmse1*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,280))
mylist_ga_cat3_psl=mylist_to.mean('lon',skipna=True)

#%%
plt.rcParams.update({'hatch.color': 'grey'})    
fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(25, 10),subplot_kw={'projection': ccrs.PlateCarree(180)})
plt.figtext(0.125, 0.865, 'a)')
plt.figtext(0.55, 0.825, 'b)')
plt.figtext(0.125, 0.45, 'c)')
plt.figtext(0.55, 0.4, 'd)')
#plt.figtext(0.125, 0.34, 'e)')
#plt.figtext(0.55, 0.31, 'f)')
axlist = axarr.flatten()

cf1=axlist[0].contourf(mylist_cat1_diff.lon, mylist_cat1_diff.lat, mylist_cat1_diff,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
cs = axlist[0].contourf(mylist_cat3_pos.lon, mylist_cat3_pos.lat, mylist_cat1_pos,alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])
cs = axlist[0].contourf(mylist_cat1_neg.lon, mylist_cat1_neg.lat, mylist_cat1_neg,alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])           
      
axlist[0].set_title('OT category: bias relative to observed SST', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[0], orientation='vertical', shrink=0.8, pad=0.02,ticks=[-3, -1.5, 0, 1.5, 3])
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
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
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

cf1=axlist[2].contourf(mylist_cat3_diff.lon, mylist_cat3_diff.lat, mylist_cat3_diff,levels, extend="both",
             transform=ccrs.PlateCarree(),vmin=k1,vmax=k2, cmap=cmap)
cs = axlist[2].contourf(mylist_cat3_pos.lon, mylist_cat3_pos.lat, mylist_cat3_pos,alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])
cs = axlist[2].contourf(mylist_cat3_pos.lon, mylist_cat3_pos.lat, mylist_cat3_neg,alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])           
      
axlist[2].set_title('EP category: bias relative to observed SST', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[2], orientation='vertical', shrink=0.8, pad=0.02,ticks=[-3, -1.5, 0, 1.5, 3])
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
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

plt.rcParams.update({'hatch.color': 'black'})   

cf1=axlist[1].contourf(mylist_cat1_psl_diff.lon, mylist_cat1_psl_diff.lat, mylist_cat1_psl_diff.sel(variable='psl'),levels1, extend="both",
             transform=ccrs.PlateCarree(),vmin=p1,vmax=p2, cmap=cmap1)
cs = axlist[1].contourf(mylist_cat1_pos_psl.lon, mylist_cat1_pos_psl.lat, mylist_cat1_pos_psl.sel(variable='psl'),alpha=0.01, extend="both",
                        
             transform=ccrs.PlateCarree(),hatches=['.'])
cs = axlist[1].contourf(mylist_cat1_neg_psl.lon, mylist_cat1_neg_psl.lat, mylist_cat1_neg_psl.sel(variable='psl'),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])           

#q1 = axlist[1].quiver(mylist_cat1_uas_diff.lon[::5],mylist_cat1_uas_diff.lat[::5],mylist_cat1_uas_diff[::5,::5],mylist_cat1_vas_diff[::5,::5],transform=ccrs.PlateCarree(180),scale=50,headwidth=4)
comma_fmt = FuncFormatter(lambda x, p: format(int(x), ','))
axlist[1].set_title('OT category: bias relative to observed SLP', fontsize=20)
cb1 = fig.colorbar(cf1, ax=axlist[1], orientation='vertical', shrink=0.6, pad=0.02, ticks=[-200,-100,0,100,200],format=comma_fmt)
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
gl.xlocator = mticker.FixedLocator([60,120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

cf1=axlist[3].contourf(mylist_cat3_psl_diff.lon, mylist_cat3_psl_diff.lat, mylist_cat3_psl_diff.sel(variable='psl'),levels1, extend="both",
             transform=ccrs.PlateCarree(),vmin=p1,vmax=p2, cmap=cmap1)
#q1 = axlist[3].quiver(mylist_cat3_uas_diff.lon[::5],mylist_cat3_uas_diff.lat[::5],mylist_cat3_uas_diff[::5,::5],mylist_cat3_vas_diff[::5,::5],transform=ccrs.PlateCarree(180),scale=50,headwidth=4)
cs = axlist[3].contourf(mylist_cat3_pos_psl.lon, mylist_cat3_pos_psl.lat, mylist_cat3_pos_psl.sel(variable='psl'),alpha=0.01, extend="both",
                        
             transform=ccrs.PlateCarree(),hatches=['.'])
cs = axlist[3].contourf(mylist_cat3_neg_psl.lon, mylist_cat3_neg_psl.lat, mylist_cat3_neg_psl.sel(variable='psl'),alpha=0.01, extend="both",
             transform=ccrs.PlateCarree(),hatches=['.'])  
axlist[3].set_title('EP category: bias relative to observed SLP', fontsize=20)


cb1 = fig.colorbar(cf1, ax=axlist[3], orientation='vertical', shrink=0.6, pad=0.02, ticks=[-200,-100,0,100,200],format=comma_fmt)

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
gl.xlocator = mticker.FixedLocator([60, 120, 180, -120, -60])
gl.ylocator = mticker.FixedLocator([-40,-20,0,20,40])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 19, 'color': 'gray'}
gl.ylabel_style = {'size': 19, 'color': 'gray'}

#%% version wind stress
plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'axes.titlesize': 'medium'})
fig = figure(figsize=(18,10))
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0, 0:1])
ax2 = plt.subplot(gs[0,1:2])
ax3 = plt.subplot(gs[1,0:1])
ax4 = plt.subplot(gs[1,1:2])


fig = gcf()
gs.tight_layout(fig,h_pad=3,w_pad=2.5)
axlist = [ax1,ax2,ax3,ax4]

plt.figtext(0.05, 0.965, 'a)')
plt.figtext(0.55, 0.965, 'b)')
plt.figtext(0.05, 0.455, 'c)')
plt.figtext(0.55, 0.455, 'd)')
#plt.figtext(0.125, 0.34, 'e)')
#plt.figtext(0.55, 0.31, 'f)')
#axlist = axarr.flatten()
    

for x in range(10):
    te_c=ws_profile_cat3.sel(new_dim=x)
    axlist[3].plot(te_c,'grey',alpha=0.5)
    a,=axlist[3].plot(ws_profile_obs,'black',label='observed',linewidth=2.0)
    b,=axlist[3].plot(ws_profile_cat3.mean('new_dim'),'red',label='EP category mean',linewidth=2.0)
    c,=axlist[3].plot(ws_profile_cmip,'green', linestyle='dashed',label='CMIP6 mean',linewidth=2.0)
    axlist[3].set_title('Equatorial zonal wind stress, EP-category')
    axlist[3].set_xlabel('longitude')
    axlist[3].set_ylabel('N m$^{-2}$')
    axlist[3].set_xlim(40, 280)
    axlist[3].set_ylim(-0.08, 0.03)
    l=axlist[3].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)

rect = patches.Rectangle((95, -0.05), 45, 0.06, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[3].add_patch(rect)

for x in range(7):
    te_c=ws_profile_cat1.sel(new_dim=x)
    axlist[1].plot(te_c,'grey',alpha=0.5)
    a,=axlist[1].plot(ws_profile_obs,'black',label='observed',linewidth=2.0)
    b,=axlist[1].plot(ws_profile_cat1.mean('new_dim'),'blue',label='OT category mean',linewidth=2.0)
    c,=axlist[1].plot(ws_profile_cmip,'green', linestyle='dashed',label='CMIP6 mean',linewidth=2.0)
    axlist[1].set_title('Equatorial zonal wind stress, OT-category')
 #   axlist[1].set_xlabel('longitude')
    axlist[1].set_ylabel('N m$^{-2}$')
    axlist[1].set_xlim(40, 280)
    axlist[1].set_ylim(-0.08, 0.03)
    l=axlist[1].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)
    
rect = patches.Rectangle((95, -0.05), 45, 0.06, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[1].add_patch(rect)

for x in range(10):
    te_c=profile_cat3.sel(new_dim=x)
    axlist[2].plot(te_c,'grey',alpha=0.5)
    a,=axlist[2].plot(observed_profile,'black',label='observed',linewidth=2.0)
    b,=axlist[2].plot(profile_cat3.mean('new_dim'),'red',label='EP category mean',linewidth=2.0)
    c,=axlist[2].plot(profile_cmip,'green', linestyle='dashed',label='CMIP6 mean',linewidth=2.0)
    axlist[2].set_title('Equatorial SST, EP-category')
    axlist[2].set_xlabel('longitude')
    axlist[2].set_ylabel('$^{o}$ C')
    axlist[2].set_xlim(40, 280)
    axlist[2].set_ylim(23, 31)
 #   plt.legend(handles=[a,b])
    l=axlist[2].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)
    
rect = patches.Rectangle((95, 23), 44, 7.5, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[2].add_patch(rect)

for x in range(7):
    te_c=profile_cat1.sel(new_dim=x)
    axlist[0].plot(te_c,'grey',alpha=0.5)
    a,=axlist[0].plot(observed_profile,'black',label='observed',linewidth=2.0)
    b,=axlist[0].plot(profile_cat1.mean('new_dim'),'blue',label='OT category mean',linewidth=2.0)
    c,=axlist[0].plot(profile_cmip,'green', linestyle='dashed',label='CMIP6 mean',linewidth=2.0)
    axlist[0].set_title('Equatorial SST, OT-category')
 #   axlist[0].set_xlabel('longitude')
    axlist[0].set_ylabel('$^{o}$ C')
    axlist[0].set_xlim(40, 280)
    axlist[0].set_ylim(23, 31)
 #   plt.legend(handles=[a,b])
    l=axlist[0].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)
rect = patches.Rectangle((95, 23), 44, 7.5, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[0].add_patch(rect)


#%% version wind speed
plt.rcParams.update({'font.size': 20})
fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(20, 10))
plt.figtext(0.125, 0.885, 'a)')
plt.figtext(0.55, 0.885, 'b)')
plt.figtext(0.125, 0.48, 'c)')
plt.figtext(0.55, 0.48, 'd)')
#plt.figtext(0.125, 0.34, 'e)')
#plt.figtext(0.55, 0.31, 'f)')
axlist = axarr.flatten()
    

for x in range(10):
    te_c=profile_cat3_uas.sel(new_dim=x,variable='uas')
    axlist[3].plot(te_c,'red',alpha=0.5)
    a,=axlist[3].plot(observed_profile_uas,'black',label='observed')
    b,=axlist[3].plot(profile_cat3_uas.sel(variable='uas').mean('new_dim'),'red',label='EP category mean')
    c,=axlist[3].plot(profile_cmip_uas.sel(variable='uas'),'green', linestyle='dashed',label='CMIP6 mean',alpha=0.3)
    axlist[3].set_title('Equatorial zonal wind, EP-category')
    axlist[3].set_xlabel('longitude')
    axlist[3].set_ylabel('m s$^{-1}$')
    axlist[3].set_xlim(40, 280)
    plt.ylim(-7.5, 4)
    l=axlist[3].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)

rect = patches.Rectangle((95, -6), 45, 8, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[3].add_patch(rect)

for x in range(7):
    te_c=profile_cat1_uas.sel(new_dim=x,variable='uas')
    axlist[1].plot(te_c,'blue',alpha=0.5)
    a,=axlist[1].plot(observed_profile_uas,'black',label='observed')
    b,=axlist[1].plot(profile_cat1_uas.sel(variable='uas').mean('new_dim'),'blue',label='OT category mean')
    c,=axlist[1].plot(profile_cmip_uas.sel(variable='uas'),'green', linestyle='dashed',label='CMIP6 mean',alpha=0.3)
    axlist[1].set_title('Equatorial zonal wind, OT-category')
 #   axlist[1].set_xlabel('longitude')
    axlist[1].set_ylabel('m s$^{-1}$')
    axlist[1].set_xlim(40, 280)
    plt.ylim(-7.5, 4)
    l=axlist[1].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)
    
rect = patches.Rectangle((95, -6), 45, 8, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[1].add_patch(rect)

for x in range(10):
    te_c=profile_cat3.sel(new_dim=x)
    axlist[2].plot(te_c,'red',alpha=0.5)
    a,=axlist[2].plot(observed_profile,'black',label='observed')
    b,=axlist[2].plot(profile_cat3.mean('new_dim'),'red',label='EP category mean')
    c,=axlist[2].plot(profile_cmip,'green', linestyle='dashed',label='CMIP6 mean',alpha=0.3)
    axlist[2].set_title('Equatorial SST, EP-category')
    axlist[2].set_xlabel('longitude')
    axlist[2].set_ylabel('$^{o}$ C')
    axlist[2].set_xlim(40, 280)
    axlist[2].set_ylim(23, 31)
 #   plt.legend(handles=[a,b])
    l=axlist[2].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)
    
rect = patches.Rectangle((95, 23), 44, 7.5, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[2].add_patch(rect)

for x in range(7):
    te_c=profile_cat1.sel(new_dim=x)
    axlist[0].plot(te_c,'blue',alpha=0.5)
    a,=axlist[0].plot(observed_profile,'black',label='observed')
    b,=axlist[0].plot(profile_cat1.mean('new_dim'),'blue',label='OT category mean')
    c,=axlist[0].plot(profile_cmip,'green', linestyle='dashed',label='CMIP6 mean',alpha=0.3)
    axlist[0].set_title('Equatorial SST, OT-category')
 #   axlist[0].set_xlabel('longitude')
    axlist[0].set_ylabel('$^{o}$ C')
    axlist[0].set_xlim(40, 280)
    axlist[0].set_ylim(23, 31)
 #   plt.legend(handles=[a,b])
    l=axlist[0].legend(handles=[a,b,c], prop={'size': 18})
    l.set_zorder(51)
rect = patches.Rectangle((95, 23), 44, 7.5, linewidth=1, edgecolor='none', facecolor='white',zorder=50)
ax = plt.gca()
axlist[0].add_patch(rect)

