import os
from scipy.io.netcdf import *
import datetime as dt 
import time as timer
import numpy as np
import pandas as pd
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        
        try:
            #print ("type:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
            	pass
                # print ('\t\t%s:' % ncattr,\
                      # repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
        	pass
            # ÷\print ("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        # print ("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
        	pass
            #print ('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        #print ("NetCDF dimension information:")
        for dim in nc_dims:
        	pass
            #print ("\tName:", dim )
            #print ("\t\tsize:", len(nc_fid.dimensions[dim]))
            #print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:

        #print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                pass
                #print ('\tName:', var)
                #print ("\t\tdimensions:", nc_fid.variables[var].dimensions)
                #print ("\t\tsize:", nc_fid.variables[var].size)
                #print_ncattr(var)
    return nc_dims, nc_vars

#IASI : time is 2, longitude is 1, latitude is 0, ch4 is 6

# EXTRACTING FROM METHANE DATASET
img_list=[os.path.join('tanso_fts', x) for x in os.listdir('tanso_fts')]
img_list.sort()


ch4_list=[]
lat_list=[]
lon_list=[]
time_list=[]
day_list=[]
lst_day_list = []
df =pd.DataFrame()

fig, ax = plt.subplots(figsize = (15,10))
plt.xlim(-180, 180)
plt.ylim(-90, 90)
fig.subplots_adjust(left=0., right=1., bottom=0.2, top=1)

for image in range(len(img_list)):
    nc_f = img_list[image]
    nc_fid = Dataset(nc_f, 'r')
    nc_dims, nc_vars = ncdump(nc_fid)

    # EXTRACT DATA 
    xch4= nc_fid.variables[nc_vars[7]][:]
    time= nc_fid.variables[nc_vars[2]][:]
    latitude= nc_fid.variables[nc_vars[4]][:]
    longitude= nc_fid.variables[nc_vars[3]][:]
    
    ch4_list.append(xch4)
    lat_list.append(latitude)
    lon_list.append(longitude)
    day_list.append(time)


    INDEX = image 

    #LOCAL TIME
    #time_list: ALL THE TIMESTEPS FROM 1 DAY
    time_list=([dt.datetime.fromtimestamp(x) for x in day_list[INDEX]])


    time_list_pd = pd.to_datetime(time_list, infer_datetime_format = True)
    df=pd.DataFrame({"timesteps": time_list_pd, "ch4": ch4_list[INDEX],\
     "longitude":lon_list[INDEX], "latitude":lat_list[INDEX]})
    #df = df.append(adder)

#Made dataframe with all the timesteps, ch4, longitude and latitude data for as many days as in the range

    df = df.sort_values(by = "timesteps", ascending =True)




#print(time_first)
#df =(df[((df.timesteps.dt.hour >= 0 ) & (df.timesteps.dt.hour <= 7)) | ((df.timesteps.dt.hour >= 19 ) & (df.timesteps.dt.hour <= 23))])
#df = df[((df.timesteps.dt.day == 5))]
#df = df[((df.ch4 == 0))]
#df = (df[((df.latitude <= 16) & (df.latitude >= 15)) & ((df.longitude <= 15) & (df.longitude >= 14))])
#print(df)
#first =(df.iloc[0,0])
#time_first = ("%s:%s" % (first.hour, first.minute))

#GEOPANDAS BELOW

#Creating Shapely points using longitude and latitude data
    points = df.apply(lambda row: Point(row.longitude, row.latitude), axis =1)
#print(points.head())

#Making a geodataframe & giving it geometry using shapely points from above.
    geo_df = gpd.GeoDataFrame(df, geometry =points)

#Notifying the geodataframe that we are using longitude and latitude as points.
    geo_df.crs = {"init": "epsg:4326"}

    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))


    base = world.plot(ax =ax, color="k")

    plot1 = geo_df.plot(ax = base, column = "ch4", cmap ="OrRd", alpha =0.9, markersize = 9,scheme = "natural_breaks")
#world.boundary.plot()
#geo_df.plot()
#print(geo_df.head())
#print(type(geo_df))
    plt.pause(0.0003)
    #plot1.remove()
    #plt.show(block = True)
'''
#PLOTS BELOW USING BASEMAP
fig, ax= plt.subplots()
plt.xlim(-180, 180)
plt.ylim(-90, 90)
fig.subplots_adjust(left=0., right=1., bottom=0., top=0.9)


m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,\
            llcrnrlon=-180, urcrnrlon=180, resolution='c',lon_0=0)
m.drawcoastlines()
m.drawmapboundary()
parallels = np.arange(-90,90,30)
meridians=np.arange(-180,180,45)
m.drawparallels(parallels)
m.drawmeridians(meridians)

# m IS THUS A BASEMAP OBJECT


#PLOTTING THE MAP ON THE GRAPH

for i in range(len(lon_list)):
    print(str(sum(ch4_list[i])/len(ch4_list[i])) + "\n")

    
    plot1=m.scatter(lon_list[i], lat_list[i],latlon=True,c=np.log(ch4_list[i]),cmap='Reds', alpha=0.5)
    plt.pause(0.0005)
    plot1.remove()

print("dogs")
plt.show()
'''
