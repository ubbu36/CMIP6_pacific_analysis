#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:33:20 2020

@author: ullaheede
"""

# module import
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd

def regrid_anomaly(forcing,a):
#control 

  #  uas_control= control['uas']
 #   uas_control= control['U']

#4xCO2

    uas_4xCO2=forcing['vas']
   # uas_4xCO2=forcing['U']
   # control_timemean=uas_control.mean("time")
    uas_4xCO2_anom=uas_4xCO2#-control_timemean
  

   #uas_4xCO2_anom_an=uas_4xCO2_anom
    uas_4xCO2_anom_an=uas_4xCO2_anom.groupby('time.year').mean('time')

    ds_out = xr.Dataset({'lat': (['lat'], np.arange(-90, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 359, 1)),
                    }
                   )

    regridder = xe.Regridder(uas_4xCO2_anom_an, ds_out, 'bilinear')

    uas_regrid = regridder(uas_4xCO2_anom_an)
  
    
    return uas_regrid
#%%
def regrid_anomaly1(forcing,a):
#control 

  #  uas_control= control['uas']
 #   uas_control= control['U']

#4xCO2
    forcing=forcing.rename(va='vas')
    uas_4xCO2=forcing['vas']

   # uas_4xCO2=forcing['U']
   # control_timemean=uas_control.mean("time")
    uas_4xCO2_anom=uas_4xCO2#-control_timemean
  

   #uas_4xCO2_anom_an=uas_4xCO2_anom
    uas_4xCO2_anom_an=uas_4xCO2_anom.groupby('time.year').mean('time')

    ds_out = xr.Dataset({'lat': (['lat'], np.arange(-90, 90, 1.0)),
                     'lon': (['lon'], np.arange(0, 359, 1)),
                    }
                   )

    regridder = xe.Regridder(uas_4xCO2_anom_an, ds_out, 'bilinear')

    uas_regrid = regridder(uas_4xCO2_anom_an)
  
    
    return uas_regrid 


#%%
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r2i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

ens_number=['r1','r2']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_ACCESS-ESM1-5.nc')

#%%

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_BCC-CSM2-MR_historical_r1i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_BCC-CSM2-MR_historical_r2i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_BCC-CSM2-MR_historical_r3i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

ens_number=['r1','r2','r3']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_BCC-CSM2-MR.nc')

#%%

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_CanESM5_historical_r7i1p2f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_CanESM5_historical_r10i1p2f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_CanESM5_historical_r35i1p2f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

ens_number=['r1','r2','r3']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_CanESM5.nc')

#%%

forcing=xr.open_mfdataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/va_Amon_CESM2_historical_r11i1p1f1_gn_*.nc', concat_dim="time",
                  data_vars='minimal', coords='minimal', compat='override')
forcing=forcing.isel(plev=0)

a=int(forcing.sizes['time']/12)
output=regrid_anomaly1(forcing,a)
mylist=output


# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

mylist=xr.concat([output], 'ens_member')
ens_number=['r1']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_CESM2.nc')

#%%

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_CNRM-CM6-1_historical_r1i1p1f2_gr_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_CNRM-CM6-1_historical_r2i1p1f2_gr_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_CNRM-CM6-1_historical_r6i1p1f2_gr_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

ens_number=['r1','r2','r3']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_CNRM-CM6-1.nc')

#%%

forcing=xr.open_mfdataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/va_Amon_FGOALS-g3_historical_r1i1p1f1_gn_*.nc', concat_dim="time",
                  data_vars='minimal', coords='minimal', compat='override')
forcing=forcing.isel(plev=0)

a=int(forcing.sizes['time']/12)
output=regrid_anomaly1(forcing,a)
mylist=output


# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

mylist=xr.concat([output], 'ens_member')
ens_number=['r1']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_FGOALS.nc')
#%%
forcing = xr.open_mfdataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_GFDL-ESM4_historical_r1i1p1f1_gr1_*.nc', concat_dim="time",
                  data_vars='minimal', coords='minimal', compat='override')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output



mylist=xr.concat([output], 'ens_member')
ens_number=['r1']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_GFDL-ESM4.nc')
#%%

forcing = xr.open_mfdataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_GISS-E2-1-G_historical_r6i1p1f1_gn_*.nc', concat_dim="time",
                  data_vars='minimal', coords='minimal', compat='override')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r2i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

forcing = xr.open_mfdataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_GISS-E2-1-G_historical_r8i1p3f1_gn_*.nc', concat_dim="time",
                  data_vars='minimal', coords='minimal', compat='override')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

ens_number=['r1','r2']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_GISS-E2-1-G.nc')

#%%
forcing = xr.open_mfdataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_HadGEM3-GC31-LL_historical_r3i1p1f3_gn_*.nc', concat_dim="time",
                  data_vars='minimal', coords='minimal', compat='override')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output



mylist=xr.concat([output], 'ens_member')
ens_number=['r1']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_HadGEMGC31-LL.nc')

#%%
forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_IPSL-CM6A-LR_historical_r24i1p1f1_gr_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_IPSL-CM6A-LR_historical_r25i1p1f1_gr_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_IPSL-CM6A-LR_historical_r27i1p1f1_gr_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

ens_number=['r1','r2','r3']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_IPSL_CM6A-LR.nc')
#%%

forcing = xr.open_mfdataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_MIROC6_historical_r10i1p1f1_gn_*.nc', concat_dim="time",
                  data_vars='minimal', coords='minimal', compat='override')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output



mylist=xr.concat([output], 'ens_member')
ens_number=['r1']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_MIROC6.nc')

#%%

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=output

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_MRI-ESM2-0_historical_r1i2p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/vas_Amon_MRI-ESM2-0_historical_r2i1p1f1_gn_185001-201412.nc')
a=int(forcing.sizes['time']/12)
output=regrid_anomaly(forcing,a)
mylist=xr.concat([mylist,output], 'ens_member')

# forcing = xr.open_dataset('/Volumes/Armor_CMIP6/CMIP6_project/WIND_historical/uas_Amon_ACCESS-ESM1-5_historical_r10i1p1f1_gn_185001-201412.nc')
# a=int(forcing.sizes['time']/12)
# output=regrid_anomaly(forcing,a)
# mylist=xr.concat([mylist,output], 'ens_member')

ens_number=['r1','r2','r3']
mylist=mylist.assign_coords(ens_member=ens_number)


mylist.to_netcdf('/Volumes/Armor_CMIP6/vas_historical_MRI-ESM2-0.nc')