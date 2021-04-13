#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 18:29:28 2020

@author: ullaheede
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

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


def regrid_anomaly(control):
#control 

    uas_control= control['ts']
 #   uas_control= control['U']

#4xCO2

  #  uas_4xCO2=forcing['ts']
   # uas_4xCO2=forcing['U']

    control_timemean=uas_control.mean("time")
    #uas_4xCO2_anom=uas_4xCO2#-control_timemean


   #uas_4xCO2_anom_an=uas_4xCO2_anom
  #  uas_4xCO2_anom_an=uas_4xCO2_anom.groupby('time.year').mean('time')

    ds_out = xr.Dataset({'lat': (['lat'], np.arange(-88, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 359, 1)),
                    }
                   )

    regridder = xe.Regridder(control_timemean, ds_out, 'bilinear')

    uas_regrid = regridder(control_timemean)
   # uas_regrid = uas_regrid.assign_coords(year=list(range(a)))
    
    return uas_regrid


### load and concatenate data ###
control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-CM2_piControl_r1i1p1f1_gn_095001-144912.nc')
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-CM2_abrupt-4xCO2_r1i1p1f1_gn_095001-109912.nc')
a=int(forcing.sizes['time']/12)

output=regrid_anomaly(control)
mylist=output


control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-ESM1-5_piControl_r1i1p1f1_gn_010101-060012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-ESM1-5_piControl_r1i1p1f1_gn_060101-100012.nc')
control=xr.concat([control1,control2],'time')
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-ESM1-5_abrupt-4xCO2_r1i1p1f1_gn_010101-025012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-CSM2-MR_piControl_r1i1p1f1_gn_185001-244912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-CSM2-MR_abrupt-4xCO2_r1i1p1f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)
uas_BCC_ESM1=regrid_anomaly(control)
mylist=xr.concat([mylist,uas_BCC_ESM1], 'new_dim')

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-ESM1_piControl_r1i1p1f1_gn_185001-230012.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_BCC-ESM1_abrupt-4xCO2_r1i1p1f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)
uas_BCC_ESM1=regrid_anomaly(control)
mylist=xr.concat([mylist,uas_BCC_ESM1], 'new_dim')

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CAMS-CSM1-0_piControl_r1i1p1f1_gn_290001-314912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CAMS-CSM1-0_piControl_r1i1p1f1_gn_315001-339912.nc')
control=xr.concat([control1,control2],'time')
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_ACCESS-ESM1-5_abrupt-4xCO2_r1i1p1f1_gn_010101-025012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_piControl_r1i1p1f1_gn_600101-620012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_piControl_r1i1p2f1_gn_560101-580012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_piControl_r1i1p2f1_gn_580101-600012.nc')
control=xr.concat([control1,control2,control3],'time')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CanESM5_abrupt-4xCO2_r1i1p1f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)

uas_CanESM5=regrid_anomaly(control)
mylist=xr.concat([mylist,uas_CanESM5], 'new_dim')

control= xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CAS-ESM2-0_piControl_r1i1p1f1_gn_000101-054912.nc')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')


control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_piControl_r1i1p1f1_gn_000101-009912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_piControl_r1i1p1f1_gn_010001-019912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_piControl_r1i1p1f1_gn_020001-029912.nc')
control=xr.concat([control1,control2,control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_abrupt-4xCO2_r1i1p1f1_gn_000101-015012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_abrupt-4xCO2_r1i1p1f1_gn_015101-019912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_abrupt-4xCO2_r1i1p1f1_gn_020001-024912.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2_abrupt-4xCO2_r1i1p1f1_gn_025001-029912.nc')
forcing=xr.concat([forcing1,forcing2,forcing3,forcing4],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')
                  
control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_piControl_r1i1p1f1_gn_000101-005012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_piControl_r1i1p1f1_gn_005101-010012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_piControl_r1i1p1f1_gn_010101-015012.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-FV2_piControl_r1i1p1f1_gn_015101-020012.nc')
control=xr.concat([control1,control2,control3,control4],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_piControl_r1i1p1f1_gn_000101-009912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_piControl_r1i1p1f1_gn_010001-019912.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_abrupt-4xCO2_r1i1p1f1_gn_000101-004912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_abrupt-4xCO2_r1i1p1f1_gn_005001-009912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM_abrupt-4xCO2_r1i1p1f1_gn_010001-015012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-004912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_005001-009912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_010001-014912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_015001-019912.nc')
control=xr.concat([control1,control2,control3,control4],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CIESM_piControl_r1i1p1f1_gr_000101-005012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CIESM_piControl_r1i1p1f1_gr_005101-010012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CIESM_piControl_r1i1p1f1_gr_010101-015012.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CIESM_piControl_r1i1p1f1_gr_015101-020012.nc')
control=xr.concat([control1,control2,control3,control4],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CMCC-CM2-SR5_piControl_r1i1p1f1_gn_185001-209912.nc')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1_piControl_r1i1p1f2_gr_185001-234912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1_abrupt-4xCO2_r1i1p1f2_gr_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-CM6-1-HR_piControl_r1i1p1f2_gr_185001-214912.nc')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-ESM2-1_piControl_r1i1p1f2_gr_185001-234912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_CNRM-ESM2-1_abrupt-4xCO2_r1i1p1f2_gr_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_010101-012512.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_012601-015012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_015101-017512.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_017601-020012.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_020101-022512.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_022601-025012.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-f3-L_piControl_r1i1p1f1_gr_060001-116012.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-f3-L_abrupt-4xCO2_r1i1p1f1_gr_185001-200912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_020001-020912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_021001-021912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_022001-022912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_023001-023912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_024001-024912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_025001-025912.nc')
control7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_026001-026912.nc')
control8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_027001-027912.nc')
control9=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_029001-029912.nc')
control10=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_030001-030912.nc')
control11=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_031001-031912.nc')
control12=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_032001-032912.nc')
control13=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_033001-033912.nc')
control14=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_034001-034912.nc')
control15=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_035001-035912.nc')
control16=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_036001-036912.nc')
control17=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_037001-037912.nc')
control18=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_038001-038912.nc')
control19=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_FGOALS-g3_piControl_r1i1p1f1_gn_039001-039912.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6, control7, control8, control9, control10, control11,\
                  control12, control13, control14, control15, control16, control17, control18, control19],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_055101-065012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_055101-065012.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_000101-010012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_010101-015012.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_piControl_r1i1p1f1_gr1_040101-050012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_piControl_r1i1p1f1_gr1_040101-050012.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_abrupt-4xCO2_r1i1p1f1_gr1_000101-010012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GFDL-ESM4_abrupt-4xCO2_r1i1p1f1_gr1_010101-015012.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output


control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_piControl_r1i1p1f1_gn_415001-420012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_piControl_r1i1p1f1_gn_420101-425012.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_piControl_r1i1p1f1_gn_425101-430012.nc')
control=xr.concat([control1,control2, control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_abrupt-4xCO2_r1i1p1f1_gn_190101-195012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_abrupt-4xCO2_r1i1p1f1_gn_195101-200012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_abrupt-4xCO2_r1i1p1f3_gn_290001-294912.nc')
forcing4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-G_abrupt-4xCO2_r1i1p1f3_gn_295001-299912.nc')

forcing=xr.concat([forcing1,forcing2, forcing3, forcing4],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_piControl_r1i1p1f1_gn_318001-323012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_piControl_r1i1p1f1_gn_323101-328012.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_piControl_r1i1p1f1_gn_328101-333012.nc')
control=xr.concat([control1,control2, control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_abrupt-4xCO2_r1i1p1f1_gn_185001-190012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_abrupt-4xCO2_r1i1p1f1_gn_190101-195012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-1-H_abrupt-4xCO2_r1i1p1f1_gn_195101-200012.nc')

forcing=xr.concat([forcing1,forcing2, forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output


#control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_piControl_r1i1p1f1_gn_200001-202512.nc')
#control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_piControl_r1i1p1f1_gn_202601-205012.nc')
#control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_piControl_r1i1p1f1_gn_205101-207512.nc')
#control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_piControl_r1i1p1f1_gn_207601-210012.nc')
#control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_piControl_r1i1p1f1_gn_210101-212512.nc')
#control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_GISS-E2-2-G_piControl_r1i1p1f1_gn_212601-215012.nc')
#control=xr.concat([control1,control2,control3,control4,control5,control6],'time')
#output=regrid_anomaly(control)
#mylist=xr.concat([mylist,output],'new_dim')
#del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_225001-234912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_225001-234912.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_abrupt-4xCO2_r1i1p1f3_gn_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-LL_abrupt-4xCO2_r1i1p1f3_gn_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_185001-186912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_187001-188912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_189001-190912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_191001-192912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_193001-194912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_195001-196912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_197001-198912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_201001-202912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_203001-204912.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_piControl_r1i1p1f1_gr1_185001-194912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_piControl_r1i1p1f1_gr1_195001-204912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_piControl_r1i1p1f1_gr1_205001-214912.nc')
control=xr.concat([control1,control2,control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_abrupt-4xCO2_r1i1p1f1_gr1_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM4-8_abrupt-4xCO2_r1i1p1f1_gr1_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_piControl_r1i1p1f1_gr1_199601-209512.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_piControl_r1i1p1f1_gr1_209601-219512.nc')
control=xr.concat([control1,control2],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_abrupt-4xCO2_r1i1p1f1_gr1_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_INM-CM5-0_abrupt-4xCO2_r1i1p1f1_gr1_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_285001-304912.nc')
#control = xr.decode_cf(control)
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_IPSL-CM6A-LR_abrupt-4xCO2_r1i1p1f1_gr_185001-214912.nc',use_cftime=True)
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_KACE-1-0-G_piControl_r1i1p1f1_gr_200001-209912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_KACE-1-0-G_piControl_r1i1p1f1_gr_210001-219912.nc')
control=xr.concat([control1,control2],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_piControl_r1i1p1f1_gn_000101-010012.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_piControl_r1i1p1f1_gn_010101-020012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_piControl_r1i1p1f1_gn_020101-030012.nc')
control=xr.concat([control1,control2,control3],'time')
control = control.rename({'longitude': 'lon', 'latitude': 'lat'})
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_abrupt-4xCO2_r1i1p1f1_gn_000101-010012.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_abrupt-4xCO2_r1i1p1f1_gn_010101-020012.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MCM-UA-1-0_abrupt-4xCO2_r1i1p1f1_gn_020101-030012.nc')
forcing=xr.concat([forcing1,forcing2,forcing3],'time')
forcing = forcing.rename({'longitude': 'lon', 'latitude': 'lat'})
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC-ES2L_piControl_r1i1p1f2_gn_225001-234912.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC-ES2L_abrupt-4xCO2_r1i1p1f2_gn_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_piControl_r1i1p1f1_gn_330001-339912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_piControl_r1i1p1f1_gn_340001-349912.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_piControl_r1i1p1f1_gn_390001-399912.nc')
control=xr.concat([control1, control2,control3],'time')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_abrupt-4xCO2_r1i1p1f1_gn_320001-329912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_abrupt-4xCO2_r1i1p1f1_gn_330001-334912.nc')
forcing3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MIROC6_abrupt-4xCO2_r1i1p1f1_gn_335001-344912.nc')
forcing=xr.concat([forcing1,forcing2, forcing3],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_185001-186912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_187001-188912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_189001-190912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_191001-192912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_193001-194912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_195001-196912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_197001-198912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_199001-200912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_201001-202912.nc')
control7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM-1-2-HAM_piControl_r1i1p1f1_gn_203001-204912.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6, control7],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_185001-186912.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_187001-188912.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_189001-190912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_191001-192912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_193001-194912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_195001-196912.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_197001-198912.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_199001-200912.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_201001-202912.nc')
control7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_203001-204912.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6, control7],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')
del output

control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MRI-ESM2-0_piControl_r1i1p1f1_gn_185001-255012.nc')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_MRI-ESM2-0_abrupt-4xCO2_r10i1p1f1_gn_185001-200012.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_piControl_r1i1p1f1_gn_050001-059912.nc')
control2 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_piControl_r1i1p1f1_gn_060001-069912.nc')
control3 = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_piControl_r1i1p1f1_gn_070001-079912.nc')
control=xr.concat([control1, control2,control3],'time')
forcing=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NESM3_abrupt-4xCO2_r1i1p1f1_gn_185001-199912.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NorCPM1_piControl_r1i1p1f1_gn_000101-010012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NorCPM1_piControl_r1i1p1f1_gn_010101-020012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_NorCPM1_piControl_r1i1p1f1_gn_020101-030012.nc')
control=xr.concat([control1,control2,control3],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_000101-001012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_001101-002012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_002101-003012.nc')
control4=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_003101-004012.nc')
control5=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_004101-005012.nc')
control6=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_005101-006012.nc')
control7=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_007101-008012.nc')
control8=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_008101-009012.nc')
control9=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_009101-010012.nc')
control10=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_010101-011012.nc')
control11=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_011101-012012.nc')
control12=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_012101-013012.nc')
control13=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_013101-014012.nc')
control14=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_014101-015012.nc')
control15=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_015101-016012.nc')
control16=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_016101-017012.nc')
control17=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_017101-018012.nc')
control18=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_018101-019012.nc')
control19=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_SAM0-UNICON_piControl_r1i1p1f1_gn_019101-020012.nc')
control=xr.concat([control1,control2,control3,control4,control5,control6, control7, control8, control9, control10, control11,\
                  control12, control13, control14, control15, control16, control17, control18, control19],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')
del output

control1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_TaiESM1_piControl_r1i1p1f1_gn_020101-030012.nc')
control2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_TaiESM1_piControl_r1i1p1f1_gn_030101-040012.nc')
control3=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_TaiESM1_piControl_r1i1p1f1_gn_040101-050012.nc')
control=xr.concat([control1,control2,control3],'time')
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output],'new_dim')


control = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_UKESM1-0-LL_piControl_r1i1p1f2_gn_255001-264912.nc')
forcing1=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_UKESM1-0-LL_abrupt-4xCO2_r1i1p1f2_gn_185001-194912.nc')
forcing2=xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/TS/ts_Amon_UKESM1-0-LL_abrupt-4xCO2_r1i1p1f2_gn_195001-199912.nc')
forcing=xr.concat([forcing1,forcing2],'time')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(control)
mylist=xr.concat([mylist,output], 'new_dim')
del output

xr.concat([mylist], pd.Index(list(model_names), name='new_dim'))


mylist.to_netcdf('/Volumes/Armor_CMIP6/control_timemean_ts_1deg.nc')