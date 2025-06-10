#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np 
import matplotlib.pyplot as plt 
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter, zoom


moho_raw_datafile = '/Users/benfernando/Desktop/Leeds/data/depthtomoho.xyz'
nc_moho_file = '/Users/benfernando/Desktop/Leeds/data/moho.nc'


def load_moho_data(plot=False):
    data = np.hsplit(np.loadtxt(moho_raw_datafile),3)
    
    latitude = np.linspace(-89.999,89.999,num=180) #don't flip this coordinate as else not ascendingly sorted
    longitude = np.linspace(-179.999,179.999,num=360)
    
    depth = np.flip(np.reshape(data[2],(180,360)),0)

    if plot ==True:
        plot_moho(latitude, longitude, depth)

    return latitude, longitude, depth

def write_data_to_netcdf_file(smooth_overall=True, smooth_poles=False,resample=True, sigma=1, resample_factor =10, write_file = 'moho_rf.nc'):
    latitude_points, longitude_points, depth_points = load_moho_data()
    depth_points = depth_points*1000 #convert to metres here 
    
    ###I THOUGHT SMOOTHING AT POLES MIGHT BE NEEDED; ROUTINE BELOW FOR LINEAR INTERPOLATION BETWEEN 80/90 N/S - DOESN'T FIX ISSUES BUT HERE 
    ### IN CASE SOMEONE ELSE WANTS TO TRY IT. CAN ALSO SET ALL POINTS AT 90 DEGREES TO BE THE SAME: THIS IS AN ARTEFACT OF GOING FROM GRID-CENTRED
    ###TO CELL CENTRED (POINTS AT 89.5 DEGREES MOVE TO 90 DEGREES).
    
    #set points at poles to smooth? 
    mean_north_pole_moho_depth = np.mean(depth_points[-1,:]) #origin is in the South-West corner as coordinates must be ascendingly sorted 
    mean_south_pole_moho_depth = np.mean(depth_points[0,:])
    mean_80N_moho_depth = np.mean(depth_points[-9,:])
    mean_80S_moho_depth = np.mean(depth_points[9,:])
    
    grad_in_depth_N=(mean_north_pole_moho_depth - mean_80N_moho_depth)/10.0 #not sure what way round signs should go but it's consistent!
    grad_in_depth_S=(mean_south_pole_moho_depth - mean_80S_moho_depth)/-10.0
               
    if smooth_poles == True:  
        
        for i in range(0,360):
            depth_points[0,i] = mean_north_pole_moho_depth
            depth_points[-1,i] = mean_south_pole_moho_depth
            
            for j in range(170,179): #area south of north pole
                depth_points[j,i] = mean_north_pole_moho_depth - grad_in_depth_N*(180-j)
            
            for k in range(0,9): #area north of south pole
                depth_points[k,i] = mean_south_pole_moho_depth - grad_in_depth_S*(k)

    if smooth_overall ==True: 
    #THIS FUNCTION SMOOTHS THE WHOLE GRID: BEAR IN MIND THAT IF YOU USE THIS AFTER SMOOTHING THE POLES, IT WILL MAKE ALL THE POINTS AT 
    #90 N/S HAVE SLIGHTLY DIFFERENT VALUES DUE TO THE TAILS OF THE GAUSSIAN
        
        depth_points = gaussian_filter(depth_points,sigma)
        
    
    if resample==True:
        #JUST IN CASE RESAMPLING THE MODEL TO HIGHER RESOLUTION REMOVED THE INSTABILITY (DIDN'T SEEM TO?)
        depth_points = zoom(depth_points, resample_factor)
        latitude_points = np.linspace(-90,90,num=180*resample_factor)
        longitude_points = np.linspace(-180, 180,num=360*resample_factor)
        

    try:    
        ncfile = Dataset(write_file,'w',format='NETCDF4_CLASSIC')
    except PermissionError:
        Dataset(write_file).close()
        ncfile = Dataset(write_file,'w',format='NETCDF4_CLASSIC')
    
    ncfile.createDimension('latitude',180*resample_factor) #creates a dimension within the netcdf file with length 180*spd
    ncfile.createDimension('longitude',360*resample_factor)
    ncfile.createDimension('relative_moho_radius',0) #means that this dimension can grow as needed.

    latitude = ncfile.createVariable('latitude','d',('latitude',)) #creates variables of type 'float' within the netcdf file, with lengths and headers
    longitude = ncfile.createVariable('longitude','d',('longitude',))
    relative_moho_radius = ncfile.createVariable('relative_moho_radius','d',('latitude','longitude'))

    longitude[:] = longitude_points #actually fill data in 
    latitude[:] = latitude_points
    
    radius_points = 6371 *1e3 + depth_points #depth points are negative - UPDATED TO CORRECT UNITS 
    
    relative_moho_radius[:] = radius_points - 6346600
    
    ncfile.close()

def read_data_from_netcdf_file(filename='moho.nc', plot=False):
    
    model = Dataset(filename)
    
    latitude = model.variables['latitude']
    longitude = model.variables['longitude']
    rmr = model.variables['relative_moho_radius'] 

    if plot ==True:
        plot_moho(latitude, longitude, rmr)

    return model, latitude, longitude, rmr

def plot_moho(latitude,longitude, rmr):
        p,t = np.meshgrid(longitude, latitude)
        plt.figure()
        plt.pcolormesh(p,t,rmr)
    
