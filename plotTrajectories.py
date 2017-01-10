import netCDF4
import os
import glob
import numpy as np
import sys
from IcbFunctions import *
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from shapely.geometry import MultiLineString
from cartopy.feature import NaturalEarthFeature
import gdal
import pylab as m
import matplotlib.colors as mcolors
from matplotlib.colors import colorConverter
import matplotlib as mpl
import matplotlib.cm as cm

nProcs = 800

ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-200, 200, -55, -55],ccrs.PlateCarree())

#drawAntarctica(ax)
pathIcbClass='/Users/imerino/Documents/These/Results/GNM009/icbClassFiles_Dumm/'

trajClassFiles= [pathIcbClass+'/icebergsTrayectoryClass1.nc',pathIcbClass+'icebergsTrayectoryClass2.nc',pathIcbClass+'/icebergsTrayectoryClass5.nc',pathIcbClass+'/icebergsTrayectoryClass7.nc']
#trajClassFiles= [pathIcbClass+'/icebergsTrayectoryClass7.nc']
colors = ['blue','green','red','yellow']
listsct=[]
counter=0
for fileC in trajClassFiles:
    sct=plotTrayectoryClassByColor(fileC,ax,colors[counter],800)
    counter = counter +1
    listsct.append(sct)
    
#plt.scatter([],[],s=10,facecolor='blue',edgecolor='none')
#listsct[1]=ax.scatter([],[],s=10,facecolor='green',edgecolor='none')
#listsct[2]=ax.scatter([],[],s=10,facecolor='red',edgecolor='none')
#listsct[3]=ax.scatter([],[],s=10,facecolor='yellow',edgecolor='none')

ax.legend(listsct,
           ('Class 1', 'Class 2', 'Class 5', 'Class 7'),
           scatterpoints=1,
           loc='upper right',
           ncol=1,
           fontsize=10,
           markerscale=10)

#ncfile = netCDF4.Dataset('../Data/icebergsTrayectoryClass10.nc','a')
#lonVar= np.array(ncfile.variables['lon'])[:,:]
#latVar= np.array(ncfile.variables['lat'])[:,:]
#tsVar= np.array(ncfile.variables['timestep'])[:,:]
#icbNumVar= np.array(ncfile.variables['iceberg_number'])[:,:]
#ncfile.close()
#
#norm = mpl.colors.Normalize(vmin=1500, vmax=3000)
#cmap=cm.autumn
#m = cm.ScalarMappable(norm=norm,cmap=cmap)
#
#nIcb=icbNumVar.shape[0]
#
#for i in np.arange(nIcb):
#    index=np.where(latVar[i,:]<1e30)
#    lat=latVar[i,index]
#    lon=lonVar[i,index]
#    ts=tsVar[i,index]
#    ids=icbNumVar[i]
#    print 'ids',ids
#    #revisa construccion de T
#    print ts
#    print m.to_rgba(ts)
#    process=np.mod(ids[0],256)
#    #if process==6:
#        #ax.scatter(lon,lat,m.to_rgba(ts),cmap=cmap,transform=ccrs.Geodetic(),alpha=1.)
#    ax.scatter(lon,lat,c=ts,cmap=cmap,vmin=15000,vmax=30000,edgecolors='None',transform=ccrs.Geodetic(),alpha=1.)



##drawAntarctica(ax)
#icb =[6,0,0]
#print getIcbClass(icb,path)
#cm = plt.cm.get_cmap('autumn')
#counter=0
#
##createIcbFile('../Data/icebergsTrayectoryClass1.nc')
#ncfile = netCDF4.Dataset('../Data/icebergsTrayectoryClass1.nc','a')
#lonVar= ncfile.variables['lon']
#lonVar2= np.array(ncfile.variables['lon'])[:]
#latVar= ncfile.variables['lat']
#tsVar= ncfile.variables['timestep']
#icbNumVar= ncfile.variables['iceberg_number']
#
#while True:
#    icb=[6+256*counter,0,0]
#
#    lon,lat,ts=extractIcbPath(icb,path)
#    
#    print lonVar2.shape, 'HOLA',lon.size
#    lonVar[counter,:] = lon.reshape(1,lon.size)
#    latVar[counter,:] = lat.reshape(1,lat.size)
#    tsVar[counter,:] = ts.reshape(1,ts.size)
#    icbNumVar[counter,:] = np.array(icb)
#    
#    drawIcebergs(ax,lon,lat)
#    #sc=ax.scatter(lon,lat,transform=ccrs.Geodetic(),edgecolor='black',alpha=0.5,cmap=cm,vmin=1,vmax=60)
#    counter = counter+1
#    
#    if counter==4:
#        break
#
#
ax.coastlines()
ax.gridlines()
plt.show()
