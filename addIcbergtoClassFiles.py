import netCDF4
import sys
import os
import glob
import numpy as np
from IcbFunctions import *


#The Class files where icb data will be sorted and stored
classFilespath = '/Users/imerino/Documents/These/Results/GNM017/icbClassFiles/'
#Create 10 class files, to be adapted according to your number of class files
createIcbClassFiles(classFilespath,10)

#Number of Procs being considered (to be adapted to your runs)
nProcs = 800

#Path of the trajectory files 
path='/Users/imerino/Documents/These/Results/GNM017/189960/'

#Path of the trajectory files
path2='/Users/imerino/Documents/These/Results/GNM017/219160/'

#In case Class Files are empty we start to filling up them with data from path, and then from path2
#rangeIcbByClassInitial(classFilespath,path,nProcs)
#rangeIcbByClass(classFilespath,path2,nProcs)

#In case we want to extract full trajectory data from a run in the middle of our simulation 
# we need the list of the last id from path (previous run) to fill up with data from path2
lista=lastIcbergPerProc(path,nProcs) #We get the list of the last icebergs id created by each processor
rangeIcbByClassFromList(lista,classFilespath,path2,nProcs) #We extract the data
