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
grad=grad.assign_coords(new_dim=range(1,len(model_names)+1))

east_pct=ts_pct.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west_pct=ts_pct.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')

grad_pct=west_pct-east_pct
grad_pct=grad_pct.assign_coords(new_dim=range(1,len(model_names)+1))


abrupt1_4x = grad.sel(year=slice(0,25)).mean('year')
abrupt2_4x = grad.sel(year=slice(100,149)).mean('year')

abrupt1_pct = grad_pct.sel(year=slice(20,80)).mean('year')
abrupt2_pct = grad_pct.sel(year=slice(20,40)).mean('year')
#####################HISTORICAL####################################


filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_historical*.nc'))
filelist=sorted(filelist,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
ts_anom=mylist['ts']
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
grad0=west-east
grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
abrupt1_hist=grad1.mean('ens_member')
grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
abrupt2_hist=grad2.mean('ens_member')
    

for x in range(1,len(filelist)):
#for x in range(11):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    ts_anom=mylist['ts']

    east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    grad0=west-east
    grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
    abrupt1_hist=xr.concat([abrupt1_hist,grad1.mean('ens_member')],'new_dim')
    grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
    abrupt2_hist=xr.concat([abrupt2_hist,grad2.mean('ens_member')],'new_dim')

abrupt1_hist=abrupt1_hist.assign_coords(new_dim=range(1,len(model_names)+1))
abrupt2_hist=abrupt2_hist.assign_coords(new_dim=range(1,len(model_names)+1))

#########################HISTORIAL GHGonly################################
    
filelist = glob.glob(os.path.join('/Volumes/Armor_CMIP6/', 'ts_GHGonly*.nc'))
filelist=sorted(filelist,key=str.lower)

mylist=xr.open_dataset(filelist[0],decode_cf=False)
ts_anom=mylist['ts']
east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
grad0=west-east
grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
abrupt1_ghg=grad1.mean('ens_member')
grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
abrupt2_ghg=grad2.mean('ens_member')
    

for x in range(1,len(filelist)):
#for x in range(11):
    mylist=xr.open_dataset(filelist[x],decode_cf=False)
    ts_anom=mylist['ts']

    east=ts_anom.sel(lat=slice(-5,5),lon=slice(e1,e2)).mean('lat').mean('lon')
    west=ts_anom.sel(lat=slice(-5,5),lon=slice(w1,w2)).mean('lat').mean('lon')
    grad0=west-east
    grad1=grad0.sel(year=slice(2000,2014)).mean('year')-grad0.sel(year=slice(1950,1970)).mean('year')
    abrupt1_ghg=xr.concat([abrupt1_ghg,grad1.mean('ens_member')],'new_dim')
    grad2=grad0.sel(year=slice(1990,2014)).mean('year')-grad0.sel(year=slice(1850,1880)).mean('year')
    abrupt2_ghg=xr.concat([abrupt2_ghg,grad2.mean('ens_member')],'new_dim')

GHGonly_models=[1,2,5,7,13,18,20,21,23,27,31,32,34,39]
GHGonly_models1=[x+1 for x in GHGonly_models]

abrupt1_ghg=abrupt1_ghg.assign_coords(new_dim=GHGonly_models1)
abrupt2_ghg=abrupt2_ghg.assign_coords(new_dim=GHGonly_models1)

#########################SSP 585################################
    

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
ssp_models1=[x+1 for x in ssp_models]
abrupt1_ssp=abrupt1_ssp.assign_coords(new_dim=ssp_models1)
abrupt2_ssp=abrupt2_ssp.assign_coords(new_dim=ssp_models1)
std1_ssp=std1_ssp.assign_coords(new_dim=ssp_models1)
std2_ssp=std2_ssp.assign_coords(new_dim=ssp_models1)

#########################SSP 370################################
    
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
#######################################################################

abrupt1_4x_ghgsubset=abrupt1_4x.isel(new_dim=GHGonly_models)
abrupt1_4x_sspsubset=abrupt1_4x.isel(new_dim=ssp_models)
abrupt1_4x_ssp3subset=abrupt1_4x.isel(new_dim=ssp3_models)

######################################################################
abrupt1_4x_C=pd.Series(abrupt1_4x)
abrupt2_4x_C=pd.Series(abrupt2_4x)
abrupt1_pct_C=pd.Series(abrupt1_pct)
abrupt1_pct_C=pd.Series(abrupt1_pct)
abrupt1_hist_C=pd.Series(abrupt1_hist)
abrupt2_hist_C=pd.Series(abrupt2_hist)
abrupt1_ghg_C=pd.Series(abrupt1_ghg)
abrupt2_ghg_C=pd.Series(abrupt2_ghg)
abrupt1_ssp_C=pd.Series(abrupt1_ssp)
abrupt2_ssp_C=pd.Series(abrupt2_ssp)
abrupt1_ssp3_C=pd.Series(abrupt1_ssp3)
abrupt2_ssp3_C=pd.Series(abrupt2_ssp3)

abrupt1_4x_ghgsubset_C=pd.Series(abrupt1_4x_ghgsubset)
abrupt1_4x_sspsubset_C=pd.Series(abrupt1_4x_sspsubset)
abrupt1_4x_ssp3subset_C=pd.Series(abrupt1_4x_ssp3subset)


corr1=abrupt1_4x_C.corr(abrupt2_4x_C,method=pearsonr_r)
pval1 = abrupt1_4x_C.corr(abrupt2_4x_C,method=pearsonr_pval)
corr2=abrupt1_4x_C.corr(abrupt1_pct_C,method=pearsonr_r)
pval2=abrupt1_4x_C.corr(abrupt1_pct_C,method=pearsonr_pval)
corr3=abrupt1_4x_C.corr(abrupt2_hist_C,method=pearsonr_r)
pval3=abrupt1_4x_C.corr(abrupt2_hist_C,method=pearsonr_pval)
corr35=abrupt1_4x_ghgsubset_C.corr(abrupt2_ghg_C,method=pearsonr_r)
pval35=abrupt1_4x_ghgsubset_C.corr(abrupt2_ghg_C,method=pearsonr_pval)
corr4=abrupt1_4x_sspsubset_C.corr(abrupt2_ssp_C,method=pearsonr_r)
pval4=abrupt1_4x_sspsubset_C.corr(abrupt2_ssp_C,method=pearsonr_pval)
corr5=abrupt1_4x_ssp3subset_C.corr(abrupt2_ssp3_C,method=pearsonr_r)
pval5=abrupt1_4x_ssp3subset_C.corr(abrupt2_ssp3_C,method=pearsonr_pval)
#########################################################
#markerlist='o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'
#colorlist='grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
#    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'
    
markerlist='o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_','o','v','^','<','>','1','2','3','4','s','p','P','*','h','+','x','X','D','|','_'
colorlist='grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k','grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k',\
    'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k', 'grey','brown','orange','olive','green','cyan','blue','purple','pink','red','k'
#########################################################3
#%%
plt.rcParams.update({'font.size': 25})
fig = figure(figsize=(28,16))
gs = gridspec.GridSpec(2, 3)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])
ax4 = plt.subplot(gs[1,0])
ax5 = plt.subplot(gs[1,1])
ax6 = plt.subplot(gs[1,2])

fig = gcf()
gs.tight_layout(fig,h_pad=4,w_pad=3)
axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
ax = axlist
ax[0].axvline(x=0.24, color='grey', linestyle='dashed',alpha=0.5)
ax[0].axvline(x=-0.24, color='grey', linestyle='dashed',alpha=0.5)
for x in range(1,len(model_names)+1):
#for x in range(40,41):
    ax[0].scatter(abrupt1_4x.sel(new_dim=x),abrupt2_4x.sel(new_dim=x), marker=markerlist[x-1],s=400,c=colorlist[x-1])

ax[0].set_title('abrupt4xCO$_2$ initial vs long-term',fontsize=29)
ax[0].set_xlabel('abrupt4xCO$_2$, year 0-25')
ax[0].set_ylabel('abrupt4xCO$_2$, year 100-150')
ax[0].set_ylim(-2.5,0.7)

new_dim_values=abrupt1_4x.new_dim.values
#for i, txt in enumerate(new_dim_values):
#    ax[0].annotate(txt, (abrupt1_4x[i], abrupt2_4x[i]))


ax[1].axvline(x=0.24, color='grey', linestyle='dashed',alpha=0.5)
ax[1].axvline(x=-0.24, color='grey', linestyle='dashed',alpha=0.5)
for x in range(1,len(model_names)+1):
    ax[1].scatter(abrupt1_4x.sel(new_dim=x),abrupt1_pct.sel(new_dim=x), marker=markerlist[x-1],s=400,c=colorlist[x-1])
ax[1].set_title('abrupt4xCO$_2$ vs 1pct',fontsize=29)
ax[1].set_xlabel('abrupt4xCO$_2$, year 0-25')
ax[1].set_ylabel('1pct, year 20-80')
ax[1].set_ylim(-1,0.4)
new_dim_values=abrupt1_4x.new_dim.values
#for i, txt in enumerate(new_dim_values):
#    ax[1].annotate(txt, (abrupt1_4x[i], abrupt1_pct[i]))

ax[5].axvline(x=0.24, color='grey', linestyle='dashed',alpha=0.5)
ax[5].axvline(x=-0.24, color='grey', linestyle='dashed',alpha=0.5)
for x in ssp_models1:
    ax[5].scatter(abrupt1_4x.sel(new_dim=x),abrupt2_ssp.sel(new_dim=x), marker=markerlist[x-1],s=400,c=colorlist[x-1])
ax[5].set_title('abrupt4xCO$_2$ vs ssp 585',fontsize=29)
ax[5].set_xlabel('abrupt4xCO$_2$, year 0-25')
ax[5].set_ylabel('ssp, year 2080-2100')
ax[5].set_ylim(-1.25,0.35)
new_dim_values=abrupt1_4x_sspsubset.new_dim.values

ax[4].axvline(x=0.24, color='grey', linestyle='dashed',alpha=0.5)
ax[4].axvline(x=-0.24, color='grey', linestyle='dashed',alpha=0.5)
for x in ssp3_models1:
    ax[4].scatter(abrupt1_4x.sel(new_dim=x),abrupt2_ssp3.sel(new_dim=x), marker=markerlist[x-1],s=400,c=colorlist[x-1])
ax[4].set_title('abrupt4xCO$_2$ vs ssp 370',fontsize=29)
ax[4].set_xlabel('abrupt4xCO$_2$, year 0-25')
ax[4].set_ylabel('ssp, year 2080-2100')
ax[4].set_ylim(-1.5,0.35)
new_dim_values=abrupt1_4x_sspsubset.new_dim.values


ax[3].axvline(x=0.24, color='grey', linestyle='dashed',alpha=0.5)
ax[3].axvline(x=-0.24, color='grey', linestyle='dashed',alpha=0.5)
for x in range(1,len(model_names)+1):
    ax[3].scatter(abrupt1_4x.sel(new_dim=x),abrupt2_hist.sel(new_dim=x), marker=markerlist[x-1],s=400,c=colorlist[x-1])
ax[3].set_title('abrupt4xCO$_2$ vs historical',fontsize=29)
ax[3].set_xlabel('abrupt4xCO$_2$, year 0-25')
ax[3].set_ylabel('historical, year 1990-2014')
ax[3].set_ylim(-0.6,0.7)

ax[2].axvline(x=0.24, color='grey', linestyle='dashed',alpha=0.5)
ax[2].axvline(x=-0.24, color='grey', linestyle='dashed',alpha=0.5)

for x in GHGonly_models1:
    ax[2].scatter(abrupt1_4x.sel(new_dim=x),abrupt2_ghg.sel(new_dim=x), marker=markerlist[x-1],s=400,c=colorlist[x-1])
ax[2].set_title('abrupt4xCO$_2$ vs histGHGonly',fontsize=29)
ax[2].set_xlabel('abrupt4xCO$_2$, year 0-25')
ax[2].set_ylabel('historical, year 1990-2014')
ax[2].set_ylim(-0.6,0.1)



plt.figtext(0.07,0.41,'R='+str("%.2f" % corr3)+',  p='+str("%.2f" % pval3))
plt.figtext(0.40,0.41,'R='+str("%.2f" % corr5)+',  p='+str("%.2f" % pval5))
plt.figtext(0.74,0.41,'R='+str("%.2f" % corr4)+',  p='+str("%.2f" % pval4))

plt.figtext(0.07,0.93,'R='+str("%.2f" % corr1)+',  p='+str("%.2f" % pval1))
plt.figtext(0.40,0.93,'R='+str("%.2f" % corr2)+',  p='+str("%.2f" % pval2))
plt.figtext(0.74,0.93,'R='+str("%.2f" % corr35)+',  p='+str("%.2f" % pval35))

plt.figtext(0.035,0.455,'d)')
plt.figtext(0.38,0.455,'e)')
plt.figtext(0.71,0.455,'f)')

plt.figtext(0.035,0.975,'a)')
plt.figtext(0.38,0.975,'b)')
plt.figtext(0.71,0.975,'c)')
