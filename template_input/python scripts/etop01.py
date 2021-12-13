#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np 
import matplotlib.pyplot as plt 
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter

copy_of_ETOPO1 = '/Users/benfernando/Desktop/ETOPO1_Ice_g_gdal.grd' #link to your copy of ETOPO1

#### READ, WRITE, PLOT NETCDF FILES


def write_data_to_netcdf_file(average= True, do_filter = True, sigma = 1, write_file = 'etop01.nc'):
    
    latitude_points, longitude_points, elevation_points = load_raw_etop01()
        
    if average == True: 
       latitude_points, longitude_points, elevation_points = average_data(latitude_points, longitude_points, elevation_points)
       
    if do_filter == True: 
        latitude_points, longitude_points, elevation_points = filter_data(latitude_points, longitude_points, elevation_points, sigma)
    
    try:    
        ncfile = Dataset(write_file,'w',format='NETCDF4_CLASSIC')
    except RuntimeError:
        Dataset(write_file).close()
        ncfile = Dataset(write_file,'w',format='NETCDF4_CLASSIC')
    
    ncfile.createDimension('latitude',180) #creates a dimension within the netcdf file with length 180*spd
    ncfile.createDimension('longitude',360)
    ncfile.createDimension('elevation',0) #means that this dimension can grow as needed.

    latitude = ncfile.createVariable('latitude','d',('latitude',)) #creates variables of type 'float' within the netcdf file, with lengths and headers
    longitude = ncfile.createVariable('longitude','d',('longitude',))
    elevation = ncfile.createVariable('elevation','d',('latitude','longitude'))

    longitude[:] = longitude_points #actually fill data in 
    latitude[:] = latitude_points
    elevation[:] = elevation_points

    ncfile.close()
    
    return longitude_points, latitude_points, elevation_points 
    
def read_data_from_netcdf_file(filename='etop01.nc', plot=False):
    
    model = Dataset(filename)
    
    latitude = model.variables['latitude']
    longitude = model.variables['longitude']
    elevation = model.variables['elevation']

    if plot==True:
        plt.figure()        
        p,t = np.meshgrid(longitude, latitude) #create a meshgrid of points to plot
        plt.pcolor(p,t,elevation) 

    return latitude, longitude, elevation


#### RESAMPLE, AVERAGE, SMOOTH ETOP01 data

def average_data(latitude, longitude, elevation): #This function averages each 1 x 1 degree cell 
        
    averaged_grid = np.empty((180,360),'d')
    
    for i in range(0,180):
        for j in range(0,360):
            averaged_grid[i,j]=np.mean(elevation[i*60:(i+1)*60, j*60:(j+1)*60])
            
    latitude_points = latitude[::60]
    longitude_points = longitude[::60]  

    return latitude_points, longitude_points, averaged_grid

def filter_data(latitude, longitude, elevation, sigma):
            
    filtered_data = gaussian_filter(elevation,sigma) #I think this should be ok in terms of units 
    
    return latitude, longitude, filtered_data

####

def load_raw_etop01():
    
    file = Dataset(copy_of_ETOPO1)
    topography_data = file.variables['z']
    
    #create strings for longitude and latitude. Do NOT include wrap-around points (180*60 rather than 180*60 +1)
    latitude = np.linspace(-90,90.001,num=180*60)
    longitude = np.linspace(-180,180.001,num=360*60)

    topographic_grid = np.reshape(topography_data, (10801,21601))
    
    topographic_grid = np.flip(topographic_grid ,0) #reverse latitudes 
    topographic_grid  = np.delete(topographic_grid ,-1,axis=0) #delete wrap-around latitude point
    topographic_grid  = np.delete(topographic_grid ,-1,axis=1) #delete wrap-around longitude point
    #note that the deletion of points is possibly not quite legitimate, as the data is grid-centred. Maybe this causes an issue at poles?

    file.close()
            
    return latitude, longitude, topographic_grid

