#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# module import
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
model_names=['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','BCC-ESM1','CAMS-CSM1-0','CanESM5','CAS-ESM2-0','CESM2','CESM2-FV','CESM2-WACCM','CESM2-WACCM-FV2',\
             'CIESM','CMCC-CM2-SR5','CNRM-CM6','CNRM-CM6-HR','CNRM-ESM2-1','E3SM','FGOALS-f3-L','FGOALS-g3','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-H',\
             'HadGEM3-GC31-LL','HadGEM3-GC3-MM','INM-CM4-8','INM-CM5-0','IPSL-CM6A','KACE-1-0-G','MCM-UA-1-0','MIROC-ES2L','MIROC6','MPI-ESM-1-2-HAM','MPI-ESM1-2-LR',\
                 'MRI-ESM2','NESM3','NorCPM1','SAM0-UNICORN','TaiESM1','UKESM1-0-LL']

def regrid_anomaly(forcing,a):
#control 


 #   uas_control= control['U']

#4xCO2

    uas_4xCO2=forcing['ts']
   # uas_4xCO2=forcing['U']
 
    uas_4xCO2_anom=uas_4xCO2#-control_timemean


   #uas_4xCO2_anom_an=uas_4xCO2_anom
    uas_4xCO2_anom_an=uas_4xCO2_anom.groupby('time.year').mean('time')

    ds_out = xr.Dataset({'lat': (['lat'], np.arange(-90, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 359, 1)),
                    }
                   )

    regridder = xe.Regridder(uas_4xCO2_anom_an, ds_out, 'bilinear')

    uas_regrid = regridder(uas_4xCO2_anom_an)
    uas_regrid = uas_regrid.assign_coords(year=list(range(a)))
    
    return uas_regrid



### load and concatenate data ###
control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-CM2_piControl_r1i1p1f1_gn_095001-144912.nc')
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-CM2_1pctCO2_r1i1p1f1_gn_095001-109912.nc')
a=int(forcing.sizes['time']/12)

output=regrid_anomaly(forcing,a)
mylist=output


control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-ESM1-5_piControl_r1i1p1f1_gn_010101-060012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-ESM1-5_piControl_r1i1p1f1_gn_060101-100012.nc')
control=xr.concat([control1,control2],'time')
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-ESM1-5_1pctCO2_r1i1p1f1_gn_010101-025012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-CSM2-MR_piControl_r1i1p1f1_gn_185001-244912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-CSM2-MR_1pctCO2_r1i1p1f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)
uas_BCC_ESM1=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,uas_BCC_ESM1], 'new_dim')

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-ESM1_piControl_r1i1p1f1_gn_185001-230012.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-ESM1_1pctCO2_r1i1p1f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)
uas_BCC_ESM1=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,uas_BCC_ESM1], 'new_dim')

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CAMS-CSM1-0_piControl_r1i1p1f1_gn_290001-314912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CAMS-CSM1-0_piControl_r1i1p1f1_gn_315001-339912.nc')
control=xr.concat([control1,control2],'time')
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CAMS-CSM1-0_1pctCO2_r1i1p1f1_gn_303001-317912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')


control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_piControl_r1i1p1f1_gn_600101-620012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_piControl_r1i1p2f1_gn_560101-580012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_piControl_r1i1p2f1_gn_580101-600012.nc')
control=xr.concat([control1,control2,control3],'time')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_1pctCO2_r2i1p2f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)

uas_CanESM5=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,uas_CanESM5], 'new_dim')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CAS-ESM2-0_1pctCO2_r1i1p1f1_gn_025001-040212.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')


control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_piControl_r1i1p1f1_gn_000101-009912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_piControl_r1i1p1f1_gn_010001-019912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_piControl_r1i1p1f1_gn_020001-029912.nc')
control=xr.concat([control1,control2,control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_1pctCO2_r1i1p1f1_gn_000101-005012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_1pctCO2_r1i1p1f1_gn_005101-010012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_1pctCO2_r1i1p1f1_gn_010101-015012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output],'new_dim')

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_piControl_r1i1p1f1_gn_000101-005012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_piControl_r1i1p1f1_gn_005101-010012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_piControl_r1i1p1f1_gn_010101-015012.nc')
control=xr.concat([control1,control2,control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_1pctCO2_r1i1p1f1_gn_000101-005012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_1pctCO2_r1i1p1f1_gn_005101-010012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_1pctCO2_r1i1p1f1_gn_010101-015012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output],'new_dim')
                  
control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_piControl_r1i1p1f1_gn_000101-009912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_piControl_r1i1p1f1_gn_010001-019912.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_1pctCO2_r1i1p1f1_gn_000101-004912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_1pctCO2_r1i1p1f1_gn_005001-009912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_1pctCO2_r1i1p1f1_gn_010001-015012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-004912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_005001-009912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_010001-014912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_015001-019912.nc')
control=xr.concat([control1,control2,control3,control4],'time')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_1pctCO2_r1i1p1f1_gn_000101-005012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_1pctCO2_r1i1p1f1_gn_005101-010012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_1pctCO2_r1i1p1f1_gn_010101-015012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CIESM_1pctCO2_r1i1p1f1_gr_000101-015012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CMCC-CM2-SR5_1pctCO2_r1i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1_piControl_r1i1p1f2_gr_185001-234912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1_1pctCO2_r1i1p1f2_gr_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1-HR_piControl_r1i1p1f2_gr_185001-214912.nc')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1-HR_1pctCO2_r1i1p1f2_gr_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1-HR_1pctCO2_r1i1p1f2_gr_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-ESM2-1_piControl_r1i1p1f2_gr_185001-234912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-ESM2-1_1pctCO2_r1i1p1f2_gr_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_010101-012512.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_012601-015012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_015101-017512.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_017601-020012.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_020101-022512.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_022601-025012.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6],'time')

forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_1pctCO2_r1i1p1f1_gr_000101-002512.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_1pctCO2_r1i1p1f1_gr_002601-005012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_1pctCO2_r1i1p1f1_gr_005101-007512.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_1pctCO2_r1i1p1f1_gr_007601-010012.nc')
forcing5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_1pctCO2_r1i1p1f1_gr_010101-012512.nc')
forcing6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_1pctCO2_r1i1p1f1_gr_012601-015012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3,forcing4,forcing5,forcing6],'time')

a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')



control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-f3-L_piControl_r1i1p1f1_gr_060001-116012.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-f3-L_1pctCO2_r1i1p1f1_gr_185001-200912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_037001-037912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_038001-038912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_039001-039912.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_040001-040912.nc')
forcing5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_041001-041912.nc')
forcing6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_042001-042912.nc')
forcing7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_043001-043912.nc')
forcing8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_044001-044912.nc')
forcing9=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_045001-045912.nc')
forcing10=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_046001-046912.nc')
forcing11=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_047001-047912.nc')
forcing12=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_048001-048912.nc')
forcing13=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_049001-049912.nc')
forcing14=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_051001-051912.nc')
forcing15=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_1pctCO2_r1i1p1f1_gn_052001-052612.nc')
forcing=xr.concat([forcing1,forcing2, forcing3,forcing4,forcing5,forcing6, forcing7, forcing8,forcing9,forcing10,forcing11,forcing12,forcing13,forcing14,forcing15],'time')

a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')


control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_055101-065012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_055101-065012.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_000101-010012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_010101-015012.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_piControl_r1i1p1f1_gr1_040101-050012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_piControl_r1i1p1f1_gr1_040101-050012.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_1pctCO2_r1i1p1f1_gr1_000101-010012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_1pctCO2_r1i1p1f1_gr1_010101-015012.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_piControl_r1i1p1f1_gn_415001-420012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_piControl_r1i1p1f1_gn_420101-425012.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_piControl_r1i1p1f1_gn_425101-430012.nc')
control=xr.concat([control1,control2, control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_1pctCO2_r102i1p1f1_gn_185001-190012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_1pctCO2_r102i1p1f1_gn_190101-195012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_1pctCO2_r102i1p1f1_gn_195101-200012.nc')

forcing=xr.concat([forcing1,forcing2, forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_piControl_r1i1p1f1_gn_318001-323012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_piControl_r1i1p1f1_gn_323101-328012.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_piControl_r1i1p1f1_gn_328101-333012.nc')
control=xr.concat([control1,control2, control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_1pctCO2_r1i1p1f1_gn_185001-190012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_1pctCO2_r1i1p1f1_gn_190101-195012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_1pctCO2_r1i1p1f1_gn_195101-200012.nc')

forcing=xr.concat([forcing1,forcing2, forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

#forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_1pctCO2_r1i1p1f1_gn_185001-187512.nc')
#forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_1pctCO2_r1i1p1f1_gn_187601-190012.nc')
#forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_1pctCO2_r1i1p1f1_gn_190101-192512.nc')
#forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_1pctCO2_r1i1p1f1_gn_192601-195012.nc')
#forcing5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_1pctCO2_r1i1p1f1_gn_195101-197512.nc')
#forcing6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_1pctCO2_r1i1p1f1_gn_197601-200012.nc')
#forcing=xr.concat([forcing1,forcing2, forcing3,forcing4,forcing5,forcing6],'time')
#a=int(forcing.sizes['time']/12)
#output=regrid_anomaly(forcing,a)
#mylist=xr.concat([mylist,output], 'new_dim')
#del output


control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_225001-234912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_225001-234912.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_1pctCO2_r1i1p1f3_gn_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_1pctCO2_r1i1p1f3_gn_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_185001-186912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_187001-188912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_189001-190912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_191001-192912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_193001-194912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_195001-196912.nc')
control7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_197001-198912.nc')
control8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_201001-202912.nc')
control9=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_203001-204912.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6,control7,control8,control9],'time')

forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_185001-186912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_187001-188912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_189001-190912.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_191001-192912.nc')
forcing5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_193001-194912.nc')
forcing6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_195001-196912.nc')
forcing7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_197001-198912.nc')
forcing8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_1pctCO2_r1i1p1f3_gn_199001-199912.nc')
forcing=xr.concat([forcing1,forcing2,forcing3,forcing4,forcing5,forcing6,forcing7,forcing8],'time')

a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')



#control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_piControl_r1i1p1f1_gr1_185001-194912.nc')
#control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_piControl_r1i1p1f1_gr1_195001-204912.nc')
#control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_piControl_r1i1p1f1_gr1_205001-214912.nc')
#control=xr.concat([control1,control2,control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_1pctCO2_r1i1p1f1_gr1_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_1pctCO2_r1i1p1f1_gr1_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

#control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_piControl_r1i1p1f1_gr1_199601-209512.nc')
#control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_piControl_r1i1p1f1_gr1_209601-219512.nc')
#control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_1pctCO2_r1i1p1f1_gr1_185001-188312.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_1pctCO2_r1i1p1f1_gr1_188401-196212.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_285001-304912.nc')
#control = xr.decode_cf(control)
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_IPSL-CM6A-LR_1pctCO2_r1i1p1f1_gr_185001-199912.nc',use_cftime=True)
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_KACE-1-0-G_1pctCO2_r1i1p1f1_gr_185001-200012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_1pctCO2_r1i1p1f1_gn_000101-010012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_1pctCO2_r1i1p1f1_gn_010101-020012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_1pctCO2_r1i1p1f1_gn_020101-030012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3],'time')
forcing = forcing.rename({'longitude': 'lon', 'latitude': 'lat'})
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC-ES2L_piControl_r1i1p1f2_gn_225001-234912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC-ES2L_1pctCO2_r1i1p1f2_gn_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_piControl_r1i1p1f1_gn_330001-339912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_piControl_r1i1p1f1_gn_340001-349912.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_piControl_r1i1p1f1_gn_390001-399912.nc')
control=xr.concat([control1, control2,control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_1pctCO2_r1i1p1f1_gn_320001-329912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_1pctCO2_r1i1p1f1_gn_330001-334912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_185001-186912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_187001-188912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_189001-190912.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_191001-192912.nc')
forcing5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_193001-194912.nc')
forcing6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_195001-196912.nc')
forcing7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_197001-198912.nc')
forcing8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_1pctCO2_r1i1p1f1_gn_199001-199912.nc')
forcing=xr.concat([forcing1,forcing2,forcing3,forcing4,forcing5,forcing6,forcing7,forcing8],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')


control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_185001-186912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_187001-188912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_189001-190912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_191001-192912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_193001-194912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_195001-196912.nc')
control7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_197001-198912.nc')
control8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_199001-200912.nc')
control9=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_201001-202912.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6,control7,control8,control9],'time')

forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_185001-186912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_187001-188912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_189001-190912.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_191001-192912.nc')
forcing5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_193001-194912.nc')
forcing6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_195001-196912.nc')
forcing7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_197001-198912.nc')
forcing8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_199001-200912.nc')
forcing9=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_201001-201412.nc')
forcing=xr.concat([forcing1,forcing2,forcing3,forcing4,forcing5,forcing6,forcing7,forcing8,forcing9],'time')

a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MRI-ESM2-0_piControl_r1i1p1f1_gn_185001-255012.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MRI-ESM2-0_1pctCO2_r1i1p1f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_piControl_r1i1p1f1_gn_050001-059912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_piControl_r1i1p1f1_gn_060001-069912.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_piControl_r1i1p1f1_gn_070001-079912.nc')
control=xr.concat([control1, control2,control3],'time')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_1pctCO2_r1i1p1f1_gn_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NorCPM1_1pctCO2_r1i1p1f1_gn_000101-016412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output


forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_185001-185912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_186001-186912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_187001-187912.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_188001-188912.nc')
forcing5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_189001-189912.nc')
forcing6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_190001-190912.nc')
forcing7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_191001-191912.nc')
forcing8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_192001-192912.nc')
forcing9=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_193001-193912.nc')
forcing10=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_194001-194912.nc')
forcing11=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_195001-195912.nc')
forcing12=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_196001-196912.nc')
forcing13=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_197001-197912.nc')
forcing14=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_198001-198912.nc')
forcing15=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_1pctCO2_r1i1p1f1_gn_199001-199912.nc')

forcing=xr.concat([forcing1,forcing2, forcing3,forcing4,forcing5,forcing6, forcing7, forcing8,forcing9,forcing10,forcing11,forcing12,forcing13,forcing14,forcing15],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_TaiESM1_piControl_r1i1p1f1_gn_020101-030012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_TaiESM1_piControl_r1i1p1f1_gn_030101-040012.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_TaiESM1_piControl_r1i1p1f1_gn_040101-050012.nc')
control=xr.concat([control1, control2,control3],'time')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_TaiESM1_1pctCO2_r1i1p1f1_gn_000101-015012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_UKESM1-0-LL_piControl_r1i1p1f2_gn_255001-264912.nc')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_UKESM1-0-LL_1pctCO2_r2i1p1f2_gn_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_UKESM1-0-LL_1pctCO2_r2i1p1f2_gn_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'new_dim')
del output

xr.concat([mylist], pd.Index(list(model_names), name='new_dim'))

mylist.to_netcdf('/Volumes/Armor_CMIP6/1pct_ts.nc')
