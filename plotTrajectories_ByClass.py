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

colors = ['blue','green','red','yellow']
listsct=[]
counter=0
for fileC in trajClassFiles:
    sct=plotTrayectoryClassByColor(fileC,ax,colors[counter],800)
    counter = counter +1
    listsct.append(sct)


ax.legend(listsct,
           ('Class 1', 'Class 2', 'Class 5', 'Class 7'),
           scatterpoints=1,
           loc='upper right',
           ncol=1,
           fontsize=10,
           markerscale=10)

ax.coastlines()
ax.gridlines()
plt.show()
