#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 11:57:59 2020

@author: ullaheede
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
import glob
import os
from pylab import *
import matplotlib.gridspec as gridspec
from scipy.stats import kendalltau, pearsonr, spearmanr

def kendall_pval(x,y):
        return kendalltau(x,y)[1]
    
def pearsonr_pval(x,y):
        return pearsonr(x,y)[1]
    
def spearmanr_pval(x,y):
        return spearmanr(x,y)[1]

def kendall_r(x,y):
        return kendalltau(x,y)[0]
    
def pearsonr_r(x,y):
        return pearsonr(x,y)[0]
    
def spearmanr_r(x,y):
        return spearmanr(x,y)[0]

plt.rcParams.update({'font.size': 22})

model_names=['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','BCC-ESM1','CAMS-CSM1-0','CanESM5','CAS-ESM2-0','CESM2','CESM2-FV','CESM2-WACCM','CESM2-WACCM-FV2',\
             'CIESM','CMCC-CM2-SR5','CNRM-CM6','CNRM-CM6-HR','CNRM-ESM2-1','E3SM','FGOALS-f3-L','FGOALS-g3','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H',\
             'HadGEM3-GC31-LL','HadGEM3-GC3-MM','INM-CM4-8','INM-CM5-0','IPSL-CM6A','KACE-1-0-G','MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM-1-2-HAM','MPI-ESM1-2-LR',\
                 'MRI-ESM2','NESM3','NorCPM1','SAM0-UNICORN','TaiESM1','UKESM1-0-LL']
    
mylist_control=xr.open_dataset('/Volumes/Armor_CMIP6/control_timemean_ts_1deg.nc')
ts_control=mylist_control['ts']
mylist_4xCO2=xr.open_dataset('/Volumes/Armor_CMIP6/4xCO2_ts.nc')
mylist_1pct=xr.open_dataset('/Volumes/Armor_CMIP6/1pct_ts.nc')
ts_4xCO2=mylist_4xCO2['ts']
ts_1pct=mylist_1pct['ts']

e1=180
e2=280

w1=80
w2=150

ts_anom=ts_4xCO2-ts_control
ts_pct=ts_1pct-ts_control
#ts_anom=ts_4xCO2
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')

grad=west-east
grad=grad.assign_coords(new_dim=range(0,len(model_names)))

east_pct=ts_pct.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west_pct=ts_pct.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')

grad_pct=west_pct-east_pct
grad_pct=grad_pct.assign_coords(new_dim=range(0,len(model_names)))


abrupt1_4x = grad.sel(year=slice(0,25)).mean('year')
abrupt2_4x = grad.sel(year=slice(100,149)).mean('year')

abrupt1_pct = grad_pct.sel(year=slice(20,80)).mean('year')
abrupt2_pct = grad_pct.sel(year=slice(20,40)).mean('year')

#%%

filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'psl_4xCO2*.nc'))
filelist=sorted(filelist,key=str.lower)

filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'psl_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist=mylist['psl']
mylist_c=xr.open_dataset(filelist_c[0],decode_cf=False)
mylist_c=mylist_c['psl']
mylist=mylist-mylist_c

lat=mylist.lat
weights=np.cos(lat* np.pi / 180.)
# mylist_to=(mylist*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
# mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_to=(mylist*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
mylist_to=mylist_to.sel(lon=slice(60,300))
mylist_ga=mylist_to.mean('lon',skipna=True)
#mylist=mylist-mylist_ga
#mylist=xr.Dataset.to_array(mylist)
mylist_out=mylist


for x in range(1,len(filelist)):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist_c=xr.open_dataset(filelist_c[x],decode_cf=False)
    mylist=mylist-mylist_c
    mylist=mylist['psl']
#    mylist=mylist.sel(year=slice(0,1)).mean('year')
    lat=mylist.lat
    weights=np.cos(lat* np.pi / 180.)

    mylist_to=(mylist*weights).sel(lat=slice(-30,30)).sum('lat',skipna=True) / weights.sel(lat=slice(-30,30)).sum('lat',skipna=True)
    mylist_to=mylist_to.sel(lon=slice(60,300))
    mylist_ga=mylist_to.mean('lon',skipna=True)
 #   mylist=mylist-mylist_ga
#    mylist=xr.Dataset.to_array(mylist)


    
    mylist_out=xr.concat([mylist_out,mylist],'new_dim')

east=mylist_out.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=mylist_out.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')

grad_psl=east-west
grad_psl=grad_psl.assign_coords(new_dim=range(0,len(model_names)))

abrupt1_psl = grad_psl.sel(year=slice(0,25)).mean('year')
abrupt2_psl = grad_psl.sel(year=slice(100,149)).mean('year')
#%%


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

#for x in range(1,len(filelist)):
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
ssp_models1=[x for x in ssp_models]
abrupt1_ssp_ts=abrupt1_ssp.assign_coords(new_dim=ssp_models1)
abrupt2_ssp_ts=abrupt2_ssp.assign_coords(new_dim=ssp_models1)
std1_ssp_ts=std1_ssp.assign_coords(new_dim=ssp_models1)
std2_ssp_ts=std2_ssp.assign_coords(new_dim=ssp_models1)
#%%

filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'psl_ssp585*.nc'))
filelist=sorted(filelist,key=str.lower)

filelist_c = glob.glob(os.path.join('/Volumes/Armor_CMIP6/psl_control_SSP585subset/', 'psl_control*.nc'))
filelist_c=sorted(filelist_c,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
mylist_control=xr.open_dataset(filelist_c[0],decode_cf=False)
ts_anom=mylist['psl']
east=mylist_control.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=mylist_control.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
gradC=east-west
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
grad0=east-west
grad1=grad0.sel(year=slice(2015,2040)).mean('year')-gradC
abrupt1_ssp=grad1.mean('ens_member')
std1_ssp=grad1.std('ens_member')
grad2=grad0.sel(year=slice(2080,2100)).mean('year')-gradC
abrupt2_ssp=grad2.mean('ens_member')
std2_ssp=grad2.std('ens_member') 

#for x in range(1,len(filelist)):
for x in range(1,len(filelist)):
#for x in range(11):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    mylist_control=xr.open_dataset(filelist_c[x],decode_cf=False)
    ts_anom=mylist['psl']
    east=mylist_control.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=mylist_control.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    gradC=east-west
    east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    grad0=east-west
    grad1=grad0.sel(year=slice(2015,2040)).mean('year')-gradC
    abrupt1_ssp=xr.concat([abrupt1_ssp,grad1.mean('ens_member')],'new_dim')
    std1_ssp=xr.concat([std1_ssp,grad1.std('ens_member')],'new_dim')
    grad2=grad0.sel(year=slice(2080,2100)).mean('year')-gradC
    abrupt2_ssp=xr.concat([abrupt2_ssp,grad2.mean('ens_member')],'new_dim')
    std2_ssp=xr.concat([std2_ssp,grad1.std('ens_member')],'new_dim')
    

ssp_models=[0,1,2,4,5,7,9,11,12,13,14,15,17,18,19,20,21,23,24,25,26,27,28,29,30,31,33,34,35,38,39]
ssp_models1=[x for x in ssp_models]
abrupt1_ssp_psl=abrupt1_ssp.assign_coords(new_dim=ssp_models1)
abrupt2_ssp_psl=abrupt2_ssp.assign_coords(new_dim=ssp_models1)
std1_ssp_psl=std1_ssp.assign_coords(new_dim=ssp_models1)
std2_ssp_psl=std2_ssp.assign_coords(new_dim=ssp_models1)

markerlist='o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'
colorlist='grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'

#%%
abrupt1_4x_C=pd.Series(abrupt2_4x)
abrupt2_4x_C_psl=pd.Series(abrupt2_psl)

abrupt2_ssp_C=pd.Series(abrupt2_ssp_ts)
abrupt2_ssp_psl=abrupt2_ssp_psl['psl']
abrupt2_ssp_C_psl=pd.Series(abrupt2_ssp_psl)


corr1=abrupt1_4x_C.corr(abrupt2_4x_C_psl,method=pearsonr_r)
pval1 = abrupt1_4x_C.corr(abrupt2_4x_C_psl,method=pearsonr_pval)

corr2=abrupt2_ssp_C.corr(abrupt2_ssp_C_psl,method=pearsonr_r)
pval2 = abrupt2_ssp_C.corr(abrupt2_ssp_C_psl,method=pearsonr_pval)

#%%
plt.rcParams.update({'font.size': 25})
fig = figure(figsize=(20,9))
gs = gridspec.GridSpec(1, 2)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0,1])


fig = gcf()
gs.tight_layout(fig,h_pad=4,w_pad=3)
axlist = [ax1,ax2]
ax = axlist
#ax[0].axvline(x=0.24, color='grey', linestyle='dashed',alpha=0.5)
#ax[0].axvline(x=-0.24, color='grey', linestyle='dashed',alpha=0.5)
for x in range(0,len(model_names)):
#for x in range(40,41):
    ax[0].scatter(abrupt2_4x.sel(new_dim=x),abrupt2_psl.sel(new_dim=x), marker=markerlist[x],s=400,c=colorlist[x])

ax[0].set_title('abrupt4xCO$_2$,years 100-150',fontsize=29)
ax[0].set_xlabel('SST gradient change, $^o$C')
ax[0].set_ylabel('SLP gradient change, Pa')

for x in ssp_models1:
#for x in range(40,41):
    ax[1].scatter(abrupt2_ssp_ts.sel(new_dim=x),abrupt2_ssp_psl.sel(new_dim=x), marker=markerlist[x],s=400,c=colorlist[x])

ax[1].set_title('ssp585 years 2080-2100',fontsize=29)
ax[1].set_xlabel('SST gradient change, $^o$C')
ax[1].set_ylabel('SLP gradient change, Pa')
#ax[0].set_ylim(-2.5,0.7)

new_dim_values=abrupt1_4x.new_dim.values
plt.figtext(0.07,0.89,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
plt.figtext(0.58,0.89,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))

plt.figtext(0.05,0.96,'a)')
plt.figtext(0.56,0.96,'b)')

