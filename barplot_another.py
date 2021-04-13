
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 17:41:59 2020

@author: ullaheede
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
import glob as glob
import os

from pylab import *
import matplotlib.gridspec as gridspec

e1=180
e2=280

w1=80
w2=150

model_names=['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','BCC-ESM1','CAMS-CSM1-0','CanESM5','CAS-ESM2-0','CESM2','CESM2-FV2','CESM2-WACCM','CESM2-WACCM-FV2',\
             'CIESM','CMCC-CM2-SR5','CNRM-CM6','CNRM-CM6-HR','CNRM-ESM2-1','E3SM','FGOALS-f3-L','FGOALS-g3','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H',\
             'HadGEM3-GC31-LL','HadGEM3-GC3-MM','INM-CM4-8','INM-CM5-0','IPSL-CM6A','KACE-1-0-G','MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM-1-2-HAM','MPI-ESM1-2-LR',\
                 'MRI-ESM2','NESM3','NorCPM1','SAM0-UNICON','TaiESM1','UKESM1-0-LL']
   
model_names1=['ersstv4','ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','BCC-ESM1','CAMS-CSM1-0','CanESM5','CAS-ESM2-0','CESM2','CESM2-FV2','CESM2-WACCM','CESM2-WACCM-FV2',\
             'CIESM','CMCC-CM2-SR5','CNRM-CM6','CNRM-CM6-HR','CNRM-ESM2-1','E3SM','FGOALS-f3-L','FGOALS-g3','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H',\
             'HadGEM3-GC31-LL','HadGEM3-GC3-MM','INM-CM4-8','INM-CM5-0','IPSL-CM6A','KACE-1-0-G','MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM-1-2-HAM','MPI-ESM1-2-LR',\
                 'MRI-ESM2','NESM3','NorCPM1','SAM0-UNICON','TaiESM1','UKESM1-0-LL']    
obs=xr.open_dataset('/Users/ullaklintheede/Downloads/ersst.v4.1854-2020.nc')
ts_obs=obs['sst']

ts_obs_a=ts_obs.groupby('time.year').mean('time',skipna=True)

east_obs=ts_obs_a.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon').mean('lev')
west_obs=ts_obs_a.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon').mean('lev')
grad_obs=west_obs-east_obs
grad_obs=grad_obs-grad_obs.sel(year=slice(1950,1970)).mean('year')

obs1=grad_obs.sel(year=slice(2000,2014)).mean('year')
std_obs=obs1*0
    
    
mylist_control=xr.open_dataset('/Volumes/Armor_CMIP6/control_timemean_ts_1deg.nc')
ts_control=mylist_control['ts']
mylist_4xCO2=xr.open_dataset('/Volumes/Armor_CMIP6/4xCO2_ts.nc')
#mylist_4xCO2=xr.open_dataset('/Volumes/Armor_CMIP6/1ptCO2_ts_anomaly.nc')
ts_4xCO2=mylist_4xCO2['ts']

ts_anom=ts_4xCO2-ts_control
#ts_anom=ts_4xCO2
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')

grad=west-east
grad_pct=grad.assign_coords(new_dim=range(1,len(model_names)+1))

abrupt1 = grad.sel(year=slice(0,25)).mean('year')
abrupt2 = grad.sel(year=slice(100,149)).mean('year')


###################################################################

mylist_control=xr.open_dataset('/Volumes/Armor_CMIP6/control_timemean_ts_1deg.nc')
ts_control=mylist_control['ts']
#mylist_4xCO2=xr.open_dataset('/Volumes/Armor_CMIP6/4xCO2_ts.nc')
mylist_4xCO2=xr.open_dataset('/Volumes/Armor_CMIP6/1pct_ts.nc')
ts_4xCO2=mylist_4xCO2['ts']

ts_anom=ts_4xCO2-ts_control

east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')

grad=west-east
grad=grad.assign_coords(new_dim=range(1,len(model_names)+1))

abrupt1_pct = grad.sel(year=slice(20,80)).mean('year')
abrupt2_pct = grad.sel(year=slice(100,149)).mean('year')

#################################################################

  

filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_historical*.nc'))
filelist=sorted(filelist,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
ts_anom=mylist['ts']
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
grad0=west-east
grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
abrupt1_hist=grad1.mean('ens_member')
std1_hist=grad1.std('ens_member')
grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
abrupt2_hist=grad2.mean('ens_member')
std2_hist=grad2.std('ens_member') 

for x in range(1,len(filelist)):
#for x in range(11):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    ts_anom=mylist['ts']

    east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    grad0=west-east
    grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
    abrupt1_hist=xr.concat([abrupt1_hist,grad1.mean('ens_member')],'new_dim')
    std1_hist=xr.concat([std1_hist,grad1.std('ens_member')],'new_dim')
    grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
    abrupt2_hist=xr.concat([abrupt2_hist,grad2.mean('ens_member')],'new_dim')
    std2_hist=xr.concat([std2_hist,grad1.std('ens_member')],'new_dim')
    
abrupt1_hist=abrupt1_hist.assign_coords(new_dim=range(1,len(model_names)+1))
abrupt2_hist=abrupt2_hist.assign_coords(new_dim=range(1,len(model_names)+1))
std1_hist=std1_hist.assign_coords(new_dim=range(1,len(model_names)+1))
std2_hist=std2_hist.assign_coords(new_dim=range(1,len(model_names)+1))

#################################################################



filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_ssp585*.nc'))
filelist=sorted(filelist,key=str.lower)
control_ssp585subset=ts_control.isel(new_dim=[0,1,2,4,5,7,9,11,12,13,14,15,17,18,19,20,21,23,24,25,26,27,28,29,30,31,33,34,35,38,39])
control_ssp585subset=control_ssp585subset.assign_coords(new_dim=range(0,len(filelist)))

mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist_control=control_ssp585subset.sel(new_dim=0)
ts_anom=mylist['ts']
east=mylist_control.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=mylist_control.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
gradC=west-east
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
grad0=west-east
grad1=grad0.sel(year=slice(2015,2040)).mean('year')-gradC
abrupt1_ssp=grad1.mean('ens_member')
std1_ssp=grad1.std('ens_member')
grad2=grad0.sel(year=slice(2080,2100)).mean('year')-gradC
abrupt2_ssp=grad2.mean('ens_member')
std2_ssp=grad2.std('ens_member') 

for x in range(1,len(filelist)):
#for x in range(11):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist_control=control_ssp585subset.sel(new_dim=x)
    ts_anom=mylist['ts']
    east=mylist_control.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=mylist_control.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    gradC=west-east
    east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    grad0=west-east
    grad1=grad0.sel(year=slice(2015,2040)).mean('year')-gradC
    abrupt1_ssp=xr.concat([abrupt1_ssp,grad1.mean('ens_member')],'new_dim')
    std1_ssp=xr.concat([std1_ssp,grad1.std('ens_member')],'new_dim')
    grad2=grad0.sel(year=slice(2080,2100)).mean('year')-gradC
    abrupt2_ssp=xr.concat([abrupt2_ssp,grad2.mean('ens_member')],'new_dim')
    std2_ssp=xr.concat([std2_ssp,grad1.std('ens_member')],'new_dim')
    

ssp_models=[0,1,2,4,5,7,9,11,12,13,14,15,17,18,19,20,21,23,24,25,26,27,28,29,30,31,33,34,35,38,39]
ssp_models1=[x+1 for x in ssp_models]
abrupt1_ssp=abrupt1_ssp.assign_coords(new_dim=ssp_models1)
abrupt2_ssp=abrupt2_ssp.assign_coords(new_dim=ssp_models1)
std1_ssp=std1_ssp.assign_coords(new_dim=ssp_models1)
std2_ssp=std2_ssp.assign_coords(new_dim=ssp_models1)
#################################################################



filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_ssp370*.nc'))
filelist=sorted(filelist,key=str.lower)
control_ssp370subset=ts_control.isel(new_dim=[0,1,2,3,4,5,7,9,12,13,14,15,18,20,21,25,26,27,28,29,30,31,32,34,39])
control_ssp370subset=control_ssp370subset.assign_coords(new_dim=range(0,len(filelist)))

mylist=xr.open_dataset(filelist[0],decode_cf=False)
ts_anom=mylist['ts']
mylist_control=control_ssp370subset.sel(new_dim=0)
east=mylist_control.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=mylist_control.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
gradC=west-east
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
grad0=west-east
grad1=grad0.sel(year=slice(2015,2040)).mean('year')-gradC
abrupt1_ssp3=grad1.mean('ens_member')
std1_ssp3=grad1.std('ens_member')
grad2=grad0.sel(year=slice(2080,2100)).mean('year')-gradC
abrupt2_ssp3=grad2.mean('ens_member')
std2_ssp3=grad2.std('ens_member') 

for x in range(1,len(filelist)):
#for x in range(11):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    ts_anom=mylist['ts']
    mylist_control=control_ssp370subset.sel(new_dim=x)
    east=mylist_control.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=mylist_control.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    gradC=west-east
    east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    grad0=west-east
    grad1=grad0.sel(year=slice(2015,2040)).mean('year')-gradC
    abrupt1_ssp3=xr.concat([abrupt1_ssp3,grad1.mean('ens_member')],'new_dim')
    std1_ssp3=xr.concat([std1_ssp3,grad1.std('ens_member')],'new_dim')
    grad2=grad0.sel(year=slice(2080,2100)).mean('year')-gradC
    abrupt2_ssp3=xr.concat([abrupt2_ssp3,grad2.mean('ens_member')],'new_dim')
    std2_ssp3=xr.concat([std2_ssp3,grad1.std('ens_member')],'new_dim')
    

ssp3_models=[0,1,2,3,4,5,7,9,12,13,14,15,18,20,21,25,26,27,28,29,30,31,32,34,39]
ssp3_models1=[x+1 for x in ssp3_models]
abrupt1_ssp3=abrupt1_ssp3.assign_coords(new_dim=ssp3_models1)
abrupt2_ssp3=abrupt2_ssp3.assign_coords(new_dim=ssp3_models1)
std1_ssp3=std1_ssp3.assign_coords(new_dim=ssp3_models1)
std2_ssp3=std2_ssp3.assign_coords(new_dim=ssp3_models1)
#################################################################

filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_GHGonly*.nc'))
filelist=sorted(filelist,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
ts_anom=mylist['ts']
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
grad0=west-east
grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
abrupt1_ghg=grad1.mean('ens_member')
std1_ghg=grad1.std('ens_member')
grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
abrupt2_ghg=grad2.mean('ens_member')
std2_ghg=grad2.std('ens_member') 

for x in range(1,len(filelist)):
#for x in range(11):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    ts_anom=mylist['ts']

    east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    grad0=west-east
    grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
    abrupt1_ghg=xr.concat([abrupt1_ghg,grad1.mean('ens_member')],'new_dim')
    std1_ghg=xr.concat([std1_ghg,grad1.std('ens_member')],'new_dim')
    grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
    abrupt2_ghg=xr.concat([abrupt2_ghg,grad2.mean('ens_member')],'new_dim')
    std2_ghg=xr.concat([std2_ghg,grad1.std('ens_member')],'new_dim')

GHGonly_models=[1,2,5,7,13,18,20,21,23,27,31,32,34,39]
GHGonly_models1=[x+1 for x in GHGonly_models]

abrupt1_ghg=abrupt1_ghg.assign_coords(new_dim=GHGonly_models1)
abrupt2_ghg=abrupt2_ghg.assign_coords(new_dim=GHGonly_models1)
std1_ghg=std1_ghg.assign_coords(new_dim=GHGonly_models1)
stdt2_ghg=std2_ghg.assign_coords(new_dim=GHGonly_models1)

#%%
plt.rcParams.update({'hatch.color': '0.1'})  
plt.rcParams.update({'font.size': 40})
x = np.arange(len(model_names))+1  # the label locations
x3=0
width = 0.25  # the width of the bars
fig = figure(figsize=(47,25))
gs = gridspec.GridSpec(2, 1)
ax1 = plt.subplot(gs[0, 0:1])
ax2 = plt.subplot(gs[1,0:1])
#ax3 = plt.subplot(gs[2,0:1])

fig = gcf()
gs.tight_layout(fig,h_pad=12,w_pad=1)
ax = [ax1,ax2]
plt.figtext(0.04, 0.98, 'a)')
plt.figtext(0.04, 0.375, 'b)')
#plt.figtext(0.04, 0.215, 'c)')

ax[0].bar(x - width, abrupt1, width, label='year 0-25, 4xCO$_2$')
ax[0].bar(x, abrupt1_pct, width, label='year 20-80, 1pct CO$_2$', hatch="\\")
ax[0].bar(x + width, abrupt2, width, label='year 100-150, 4xCO$_2$', hatch="/")

ax[0].bar(x3 + width, obs1, width, label='year 2000-2014 minus 1950-1970, observed', hatch="x")
# Add some text for labels, title and custom x-axis tick labels, etc.
ax[0].set_ylabel('$\Delta$ T ($^o$C)')
ax[0].set_title('Pacific zonal SST gradient change, hypothetical CO$_2$ scenarios',fontsize=53)
ax[0].set_xticks(range(len(model_names1)))
ax[0].set_xticklabels(model_names1,rotation='vertical')
ax[0].legend(ncol=2,fontsize=38)



x = abrupt1_hist.new_dim  # the label locations
x1 = abrupt1_ssp.new_dim
x15 = abrupt1_ssp3.new_dim  
x2 =abrupt1_ghg.new_dim
x3=0
width = 0.25  # the width of the bars

# ax[1].bar(x - width, abrupt1_hist, width, yerr=std1_hist, label='year 2000-2014 minus 1950-1970')
# ax[1].bar(x15, abrupt1_ssp3, width, yerr=std1_ssp3, label='year 2015-2035, ssp370')
# ax[1].bar(x15 + width, abrupt2_ssp3, width, yerr=std1_ssp3, label='year 2080-2100, ssp370')
# ax[1].bar(x3, obs1, width, label='year 2000-2014 minus 1950-1970, obs')

# ax[1].set_ylabel('$\Delta$ K')
# ax[1].set_title('SST gradient, historical and future projection')
# ax[1].set_xticks(range(len(model_names1)))
# ax[1].legend()

ax[1].bar(x - width, abrupt1_hist, width, yerr=std1_hist, error_kw=dict(lw=5),label='year 2000-2014 minus 1950-1970')
#ax[1].bar(x2 - width, abrupt1_ghg, width, yerr=std1_ghg, label='year 2000-2014 minus 1950-1970, GHGonly')
ax[1].bar(x1, abrupt1_ssp, width, yerr=std1_ssp, error_kw=dict(lw=5), label='year 2015-2035, ssp585', hatch="\\")
ax[1].bar(x1 + width, abrupt2_ssp, width, yerr=std1_ssp, error_kw=dict(lw=5), label='year 2080-2100, ssp585', hatch="/")
ax[1].bar(x3, obs1, width, label='year 2000-2014 minus 1950-1970, obs', hatch="x")


# Add some text for labels, title and custom x-axis tick labels, etc.
ax[1].set_ylabel('$\Delta$ T ($^o$C)')
ax[1].set_title('Pacific zonal SST gradient change, historical and future projections',fontsize=53)
ax[1].set_xticks(range(len(model_names1)))
ax[1].set_xticklabels(model_names1,rotation='vertical')
ax[1].legend(ncol=2,fontsize=38)

# ax[2].bar(x - width, abrupt1_hist, width, yerr=std1_hist, error_kw=dict(lw=5),label='year 2000-2014 minus 1950-1970')
# #ax[1].bar(x2 - width, abrupt1_ghg, width, yerr=std1_ghg, label='year 2000-2014 minus 1950-1970, GHGonly')
# ax[2].bar(x15, abrupt1_ssp3, width, yerr=std1_ssp3, error_kw=dict(lw=5), label='year 2015-2035, ssp370')
# ax[2].bar(x15 + width, abrupt2_ssp3, width, yerr=std1_ssp3, error_kw=dict(lw=5), label='year 2080-2100, ssp370')
# ax[2].bar(x3, obs1, width, label='year 2000-2014 minus 1950-1970, obs')


# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax[2].set_ylabel('$\Delta$ T ($^o$C)')
# ax[2].set_title('SST gradient, historical and future projection',fontsize=53)
# ax[2].set_xticks(range(len(model_names1)))
# ax[2].set_xticklabels(model_names1,rotation='vertical')
# ax[2].legend(ncol=2,fontsize=38)

