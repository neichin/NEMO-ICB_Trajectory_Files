import netCDF4
import os
import glob
import numpy as np
import sys
import matplotlib.pyplot as plt
import pylab as m
import matplotlib.colors as mcolors
from matplotlib.colors import colorConverter
import matplotlib as mpl
import cartopy.crs as ccrs
import matplotlib.cm as cm

def extractIcbPath(icb,path):
    trajFiles=glob.glob(path+'*.nc')
    lon=[]
    lat=[]
    width=[]
    length=[]
    thickness=[]
    mass=[]
    dt = []
    for j in trajFiles:
        ncfile = netCDF4.Dataset(j,'a')
        ids = np.array(ncfile.variables['iceberg_number'][:,:])
        lonVal = np.array(ncfile.variables['lon'][:])
        latVal = np.array(ncfile.variables['lat'][:])
        #ThickVal = np.array(ncfile.variables['thickness'][:])
        #LenVal = np.array(ncfile.variables['length'][:])
        #WidVal = np.array(ncfile.variables['width'][:])
        #MassVal = np.array(ncfile.variables['mass_scaling'][:])
        ts = np.array(ncfile.variables['timestep'][:])
        arrDetect = ids - icb
        index=np.where(arrDetect[:,0]==0)
        
        lon=np.concatenate((lon,lonVal[index]))
        lat=np.concatenate((lat,latVal[index]))
        #mass=np.concatenate((mass,MassVal[index]))
        #width=np.concatenate((width,WidVal[index]))
        #length=np.concatenate((length,LenVal[index]))
        #thickness=np.concatenate((thickness,ThickVal[index]))
        dt=np.concatenate((dt,ts[index]))
        ncfile.close()
        
    return lon,lat,mass,width,length,thickness,dt

def extractIcbPathFromClassFile(icb,path):
    lon=[]
    lat=[]
    width=[]
    length=[]
    thickness=[]
    mass=[]
    dt = []
    ncfile = netCDF4.Dataset(path,'a')
    ids = np.array(ncfile.variables['iceberg_number'][:,:])
    lonVal = np.array(ncfile.variables['lon'][:])
    latVal = np.array(ncfile.variables['lat'][:])
    #ThickVal = np.array(ncfile.variables['thickness'][:])
    #LenVal = np.array(ncfile.variables['length'][:])
    #WidVal = np.array(ncfile.variables['width'][:])
    #MassVal = np.array(ncfile.variables['mass_scaling'][:])
    ts = np.array(ncfile.variables['timestep'][:])
    arrDetect = ids - icb
    index=np.where(arrDetect[:,0]==0)
    lon=lonVal[index]
    lat=latVal[index]
    #mass=MassVal[index]
    #width=WidVal[index]
    #length=LenVal[index]
    #thickness=ThickVal[index]
    dt=ts[index]
    ncfile.close()
    print 'done'
        
    return lon,lat,dt
    #return lon,lat,mass,width,length,thickness,dt

def IcbExist(icb,path):
    trajFiles=glob.glob(path+'*.nc')
    res2=False
    fileList=[]
    for icbFile in trajFiles:
        res,dummy=IcbInFile(icb,icbFile)
        res2=res2 or res
        if res:
            fileList.append(icbFile)
        
        
    return res2,fileList
        

def IcbInFile(icb,path):
    ncfile = netCDF4.Dataset(path,'a')
    ids = np.array(ncfile.variables['iceberg_number'][:,:])
    arrDetect = abs(ids - icb)
    index=np.where(arrDetect[:,0]==0)
    ncfile.close()
    if len(index[0])>0:
        return True,index[0][0]
    else:
        return False,0
 
def getIcbClass(icb,fileList):
    
    massClass=[8.8e7,4.1e8,3.3e9,1.8e10,3.8e10,7.5e10,1.2e11,2.2e11,3.9e11,7.4e11]
    
    mass=[]
    
    for j in fileList:
        ncfile = netCDF4.Dataset(j,'a')
        ids = np.array(ncfile.variables['iceberg_number'][:,:])
        massVal = np.array(ncfile.variables['mass'][:])
        arrDetect = ids - icb
        index=np.where(arrDetect[:,0]==0)
        temp = massVal[index[0][0]]
        mass.append(temp)
        ncfile.close()
    
    massVal=max(mass)
    temp = abs(massClass - massVal)
    classM = np.argmin(np.array(temp))
        
    return classM+1
    
def createIcbFile(filename):
    ncfile = netCDF4.Dataset(filename,'w',format='NETCDF4') 

    ncfile.createDimension('n')
    ncfile.createDimension('timestep_counter')
    ncfile.createDimension('k',3)
    
    ncfile.createVariable('timestep','i',('n','timestep_counter'))
    ncfile.createVariable('lon','f',('n','timestep_counter'))
    ncfile.createVariable('lat','f',('n','timestep_counter'))
    ncfile.createVariable('length','f',('n','timestep_counter'))
    ncfile.createVariable('width','f',('n','timestep_counter'))
    ncfile.createVariable('thickness','f',('n','timestep_counter'))
    ncfile.createVariable('mass_scaling','f',('n','timestep_counter'))
    ncfile.createVariable('iceberg_number','i',('n','k'))
    ncfile.createVariable('points_number','i',('n'))
    
    ncfile.close()
    print filename , 'created'
    
def createIcbClassFiles(path,nclass):
    for i in np.arange(nclass):
        filename = path+'/icebergsTrayectoryClass'+str(i+1)+'.nc'
        createIcbFile(filename)
        
def drawIcebergs(ax,lon,lat):
    cm = plt.cm.get_cmap('autumn')
    ax.scatter(lon,lat,transform=ccrs.Geodetic(),edgecolor='black',alpha=0.5,cmap=cm,vmin=1,vmax=60)

    
def saveIcbByClass(classFile,classN,path,nProcs):
    trajFiles=glob.glob(path+'*.nc')
    ncfile = netCDF4.Dataset(classFile,'a')
    lonVar= ncfile.variables['lon']
    latVar= ncfile.variables['lat']
    #MassVar= ncfile.variables['mass_scaling']
    #WidVar= ncfile.variables['width']
    #LenVar= ncfile.variables['length']
    #ThickVar= ncfile.variables['thickness']
    tsVar= ncfile.variables['timestep']
    icbNumVar= ncfile.variables['iceberg_number']
    icbCount=0
    
    for icbFile in trajFiles:
        procNum = int(icbFile[len(icbFile)-7:len(icbFile)-3]) + 1
        counter=0
        while True:
            icb = np.array([procNum + nProcs*counter,0,0])
            res,dummyIndex=IcbInFile(icb,icbFile)
            counter=counter+1
            if res:
                icbClass = getIcbClass(icb,path)
                if icbClass==classN:
                    lon,lat,mass,width,length,thickness,ts=extractIcbPath(icb,path)
                    lonVar[icbCount,:] = lon.reshape(1,lon.size)
                    latVar[icbCount,:] = lat.reshape(1,lat.size)
                    #MassVar[icbCount,:] = mass.reshape(1,mass.size)
                    #LenVar[icbCount,:] = length.reshape(1,length.size)
                    #WidVar[icbCount,:] = width.reshape(1,width.size)
                    #ThickVar[icbCount,:] = thickness.reshape(1,thickness.size)
                    tsVar[icbCount,:] = ts.reshape(1,ts.size)
                    icbNumVar[icbCount,:] = np.array(icb)
                    icbCount=icbCount+1
            else:
                break
        
        print procNum
        if procNum > 70:
            print 'break'
            break
            
            
    ncfile.close()


#rangeIcbByClass(classFilespath,path,nProcs)
#    Range individual iceberg data from the raw output files from NEMO-ICB module to a list of 'Class Files'
#       'Class Files' need to be named as icebergsTrayectoryClass+ClassNumber+.nc
#
#INPUTS:
#   classFilespath=Path location of the 'Class Files'    
#   path=Path location of the raw trajectory files (output of NEMO_ICB module
#   nProcs=number of processors
#

def rangeIcbByClass(classFilespath,path,nProcs):
    trajFiles=glob.glob(path+'*.nc') #All the trajectory files
    for icbFile in trajFiles:
        procNum = int(icbFile[len(icbFile)-7:len(icbFile)-3]) + 1 
        print 'lets look for in file ', icbFile 
        counter=0
        while True:
            icb = np.array([procNum + nProcs*counter,0,0]) #Id of icebergs having been launched from the corresponding processor of the current file
            res,fileList=IcbExist(icb,path) #Does iceberg icb exist? and where?
            res2,icbClass,nIcb=icbInClassFiles(classFilespath,icb) #Does icb exist in any of the "Class Files"
            print 'print1',icb,res,res2
            counter=counter+1
            icbCount=0
            if res:
                if not res2:
                    icbClass = getIcbClass(icb,fileList)
                ncfile = netCDF4.Dataset(classFilespath+'/icebergsTrayectoryClass'+str(icbClass)+'.nc','a')
                print 'saving into class file #'+str(icbClass), 'icb num' ,icbCount, 'proc',procNum
                lonVar= ncfile.variables['lon']
                latVar= ncfile.variables['lat']
                #MassVar= ncfile.variables['mass_scaling']
                #WidVar= ncfile.variables['width']
                #LenVar= ncfile.variables['length']
                #ThickVar= ncfile.variables['thickness']
                tsVar= ncfile.variables['timestep']
                icbNumVar= ncfile.variables['iceberg_number']
                icbPointsVar= ncfile.variables['points_number']
                if not res2:
                    nIcb=icbNumVar.shape[0]
                lon,lat,mass,width,length,thickness,ts=extractIcbPath(icb,path) #Get full icb data from trayectory files
                #Filling up Class Files
                if res2:
                    numPoints = icbPointsVar[nIcb]
                    tsVar[nIcb,numPoints:numPoints+len(ts)] = ts.reshape(1,ts.size)
                    lonVar[nIcb,numPoints:numPoints+len(ts)] = lon.reshape(1,lon.size)
                    latVar[nIcb,numPoints:numPoints+len(ts)] = lat.reshape(1,lat.size)
                    #MassVar[nIcb,numPoints:numPoints+len(ts)] = mass.reshape(1,mass.size)
                    #LenVar[nIcb,numPoints:numPoints+len(ts)] = length.reshape(1,length.size)
                    #WidVar[nIcb,numPoints:numPoints+len(ts)] = width.reshape(1,width.size)
                    #ThickVar[nIcb,numPoints:numPoints+len(ts)] = thickness.reshape(1,thickness.size)
                    icbPointsVar[nIcb] = numPoints+lon.size
                else:
                    lonVar[nIcb,:] = lon.reshape(1,lon.size)
                    latVar[nIcb,:] = lat.reshape(1,lat.size)
                    tsVar[nIcb,:] = ts.reshape(1,ts.size)
                    #MassVar[nIcb,:] = mass.reshape(1,mass.size)
                    #LenVar[nIcb,:] = length.reshape(1,length.size)
                    #WidVar[nIcb,:] = width.reshape(1,width.size)
                    #ThickVar[nIcb,:] = thickness.reshape(1,thickness.size)
                    icbPointsVar[nIcb] = lon.size
                    
                icbNumVar[nIcb,:] = np.array(icb)
                ncfile.close()
                icbCount=icbCount+1
            elif not res2:
                print 'no iceberg or no longer floating in the ocean', icb
                break
     
#rangeIcbByClassInitial(classFilespath,path,nProcs)
#    Range individual iceberg data from the raw output files from NEMO-ICB module to a list of 'Class Files'
#       To be used only the first time
#       'Class Files' need to be named as icebergsTrayectoryClass+ClassNumber+.nc
#
#INPUTS:
#   classFilespath=Path location of the 'Class Files'    
#   path=Path location of the raw trajectory files (output of NEMO_ICB module
#   nProcs=number of processors
#

def rangeIcbByClassInitial(classFilespath,path,nProcs):
    trajFiles=glob.glob(path+'*.nc')
    
    for icbFile in trajFiles:
        procNum = int(icbFile[len(icbFile)-7:len(icbFile)-3]) + 1
        if procNum>nProcs/2:
            return
        print 'lets look for in file ', icbFile 
        counter=0
        while True:
            icb = np.array([procNum + nProcs*counter,0,0])
            res,fileList=IcbExist(icb,path)
            res2,icbClass,nIcb=icbInClassFiles(classFilespath,icb)
            res,dummyIndex=IcbInFile(icb,icbFile)
            counter=counter+1
            icbCount=0
            if res:
                icbClass = getIcbClass(icb,fileList)
                ncfile = netCDF4.Dataset(classFilespath+'/icebergsTrayectoryClass'+str(icbClass)+'.nc','a')
                print 'saving into class file #'+str(icbClass), 'icb num' ,icbCount, 'proc',procNum
                lonVar= ncfile.variables['lon']
                latVar= ncfile.variables['lat']
                tsVar= ncfile.variables['timestep']
                #MassVar= ncfile.variables['mass_scaling']
                #WidVar= ncfile.variables['width']
                #LenVar= ncfile.variables['length']
                #ThickVar= ncfile.variables['thickness']
                icbNumVar= ncfile.variables['iceberg_number']
                icbPointsVar= ncfile.variables['points_number']
                nIcb=icbNumVar.shape[0]
                lon,lat,mass,width,length,thickness,ts=extractIcbPath(icb,path)
                lonVar[nIcb,:] = lon.reshape(1,lon.size)
                latVar[nIcb,:] = lat.reshape(1,lat.size)
                #MassVar[nIcb,:] = mass.reshape(1,mass.size)
                #LenVar[nIcb,:] = length.reshape(1,length.size)
                #WidVar[nIcb,:] = width.reshape(1,width.size)
                #ThickVar[nIcb,:] = thickness.reshape(1,thickness.size)
                tsVar[nIcb,:] = ts.reshape(1,ts.size)
                icbPointsVar[nIcb] = lon.size
                icbNumVar[nIcb,:] = np.array(icb)
                ncfile.close()
                icbCount=icbCount+1
            elif not res2:
                print 'no iceberg or no longer floating in the ocean', icb
                break
    
                
#rangeIcbByClassFromList(classFilespath,path,nProcs)
#    Range individual iceberg data from the raw output files from NEMO-ICB module to a list of 'Class Files'
#       To start from the icebergs id specified into lista list
#       'Class Files' need to be named as icebergsTrayectoryClass+ClassNumber+.nc
#
#INPUTS:
#   list1:list of last icberg Id per proc, which is the last of not being considered 
#   classFilespath=Path location of the 'Class Files'    
#   path=Path location of the raw trajectory files (output of NEMO_ICB module
#   nProcs=number of processors
#            
def rangeIcbByClassFromList(lista,classFilespath,path,nProcs):
    trajFiles=glob.glob(path+'*.nc')
    print 'hola',trajFiles
    for icbFile in trajFiles:
        procNum = int(icbFile[len(icbFile)-7:len(icbFile)-3]) + 1
        if procNum>nProcs/2:
            return
        print 'lets look for in file ', icbFile 
        counter=0
        icbNumInit = lista[procNum]
        while True:
            if icbNumInit==0:
                icb = np.array([procNum + nProcs*counter,0,0])
            else:
                icb = np.array([icbNumInit + nProcs*(counter+1),0,0])
            res,fileList=IcbExist(icb,path)
            res2,icbClass,nIcb=icbInClassFiles(classFilespath,icb)
            counter=counter+1
            icbCount=0
            if res:
                icbClass = getIcbClass(icb,fileList)
                ncfile = netCDF4.Dataset(classFilespath+'/icebergsTrayectoryClass'+str(icbClass)+'.nc','a')
                print 'saving into class file #'+str(icbClass), 'icb num' ,icbCount, 'proc',procNum
                lonVar= ncfile.variables['lon']
                latVar= ncfile.variables['lat']
                #MassVar= ncfile.variables['mass_scaling']
                #WidVar= ncfile.variables['width']
                #LenVar= ncfile.variables['length']
                #ThickVar= ncfile.variables['thickness']
                tsVar= ncfile.variables['timestep']
                icbNumVar= ncfile.variables['iceberg_number']
                icbPointsVar= ncfile.variables['points_number']
                nIcb=icbNumVar.shape[0]
                lon,lat,mass,width,length,thickness,ts=extractIcbPath(icb,path)
                lonVar[nIcb,:] = lon.reshape(1,lon.size)
                latVar[nIcb,:] = lat.reshape(1,lat.size)
                #MassVar[nIcb,:] = mass.reshape(1,mass.size)
                #LenVar[nIcb,:] = length.reshape(1,length.size)
                #WidVar[nIcb,:] = width.reshape(1,width.size)
                #ThickVar[nIcb,:] = thickness.reshape(1,thickness.size)
                tsVar[nIcb,:] = ts.reshape(1,ts.size)
                icbPointsVar[nIcb] = lon.size
                icbNumVar[nIcb,:] = np.array(icb)
                ncfile.close()
                icbCount=icbCount+1
            elif not res2:
                print 'no iceberg', icb
                break
                
#rangeIcbByClassFromList(classFilespath,path,nProcs)
#    Range individual iceberg data from the raw output files from NEMO-ICB module to a list of 'Class Files'
#       Only considered icebergs ids between the id sspecified between list1 and list2
#       'Class Files' need to be named as icebergsTrayectoryClass+ClassNumber+.nc
#
#INPUTS:
#   list1:list of last icberg Id per proc, which is the last of not being considered 
#   list2:list of last icberg Id per proc, which is the last of being considered
#   classFilespath=Path location of the 'Class Files'    
#   path=Path location of the raw trajectory files (output of NEMO_ICB module
#   nProcs=number of processors
#            
def rangeIcbByClassBetweenList(list1,list2,classFilespath,path,nProcs):
    trajFiles=glob.glob(path+'*.nc')
    print 'hola',trajFiles
    for icbFile in trajFiles:
        procNum = int(icbFile[len(icbFile)-7:len(icbFile)-3]) + 1
        if procNum>nProcs/2:
            return
        print 'lets look for in file ', icbFile 
        counter=0
        icbNumInit = list1[procNum]
        MaxCounter = list2[procNum]-list1[procNum]
        while counter<MaxCounter:
            if icbNumInit==0:
                icb = np.array([procNum + nProcs*counter,0,0])
            else:
                icb = np.array([icbNumInit + nProcs*(counter+1),0,0])
            res,fileList=IcbExist(icb,path)
            res2,icbClass,nIcb=icbInClassFiles(classFilespath,icb)
            counter=counter+1
            icbCount=0
            if res:
                icbClass = getIcbClass(icb,fileList)
                ncfile = netCDF4.Dataset(classFilespath+'/icebergsTrayectoryClass'+str(icbClass)+'.nc','a')
                print 'saving into class file #'+str(icbClass), 'icb num' ,icbCount, 'proc',procNum
                lonVar= ncfile.variables['lon']
                latVar= ncfile.variables['lat']
                #MassVar= ncfile.variables['mass_scaling']
                #WidVar= ncfile.variables['width']
                #LenVar= ncfile.variables['length']
                #ThickVar= ncfile.variables['thickness']
                tsVar= ncfile.variables['timestep']
                icbNumVar= ncfile.variables['iceberg_number']
                icbPointsVar= ncfile.variables['points_number']
                nIcb=icbNumVar.shape[0]
                lon,lat,mass,width,length,thickness,ts=extractIcbPath(icb,path)
                lonVar[nIcb,:] = lon.reshape(1,lon.size)
                latVar[nIcb,:] = lat.reshape(1,lat.size)
                #MassVar[nIcb,:] = mass.reshape(1,mass.size)
                #LenVar[nIcb,:] = length.reshape(1,length.size)
                #WidVar[nIcb,:] = width.reshape(1,width.size)
                #ThickVar[nIcb,:] = thickness.reshape(1,thickness.size)
                tsVar[nIcb,:] = ts.reshape(1,ts.size)
                icbPointsVar[nIcb] = lon.size
                icbNumVar[nIcb,:] = np.array(icb)
                ncfile.close()
                icbCount=icbCount+1
            elif not res2:
                print 'no iceberg', icb
                break
        

def icbInClassFiles(ClassFilesPath,icb):
    classFiles=glob.glob(ClassFilesPath+'icebergsTrayectoryClass*.nc')
    numClass=1
    res=False
    for f in classFiles:
        numClass = int(f[len(f)-4:len(f)-3])
        if f[len(f)-5] !='s':
            numClass = numClass + int(f[len(f)-5])*10
        
        res,index=IcbInFile(icb,f)
        if res:
            break
        numClass=numClass+1
    
    return res,numClass,index
            
        
def createTrayectoryClassFiles(path,classN):
    pathNewFile= path+'icebergsTrayectoryClass'+str(classN)+'.nc'
    ncfile = netCDF4.Dataset(pathNewFile,'w',format='NETCDF4') 

    ncfile.createDimension('n')
    ncfile.createDimension('timestep_counter')
    ncfile.createDimension('k',3)
    
    ncfile.createVariable('timestep','i',('n','timestep_counter'))
    ncfile.createVariable('lon','f',('n','timestep_counter'))
    ncfile.createVariable('lat','f',('n','timestep_counter'))
    ncfile.createVariable('length','f',('n','timestep_counter'))
    ncfile.createVariable('width','f',('n','timestep_counter'))
    ncfile.createVariable('thickness','f',('n','timestep_counter'))
    ncfile.createVariable('mass_scaling','f',('n','timestep_counter'))
    ncfile.createVariable('iceberg_number','i',('n','k'))
    
    ncfile.close()
    print 'file: ',pathNewFile, 'created'
    return pathNewFile
    
def plotTrayectoryClass(fileName,ax):
    ncfile = netCDF4.Dataset(fileName,'a')
    lonVar= np.array(ncfile.variables['lon'])[:,:]
    latVar= np.array(ncfile.variables['lat'])[:,:]
    tsVar= np.array(ncfile.variables['timestep'])[:,:]
    icbNumVar= np.array(ncfile.variables['iceberg_number'])[:,:]
    ncfile.close()
    
    norm = mpl.colors.Normalize(vmin=1500, vmax=3000)
    cmap=cm.autumn
    m = cm.ScalarMappable(norm=norm,cmap=cmap)
    
    nIcb=icbNumVar.shape[0]
    
    for i in np.arange(nIcb):
        index=np.where(latVar[i,:]<1e30)
        lat=latVar[i,index]
        lon=lonVar[i,index]
        ts=tsVar[i,index]
        ids=icbNumVar[i]
        print 'ids',ids
        #revisa construccion de T
        print ts
        print m.to_rgba(ts)
        process=np.mod(ids[0],256)
        #if process==6:
            #ax.scatter(lon,lat,m.to_rgba(ts),cmap=cmap,transform=ccrs.Geodetic(),alpha=1.)
        ax.scatter(lon,lat,c=ts,cmap=cmap,vmin=15000,vmax=30000,edgecolors='None',transform=ccrs.Geodetic(),alpha=1.)
        
def plotTrayectoryClassByColor(fileName,ax,color,nProcs):
    ncfile = netCDF4.Dataset(fileName,'a')
    lonVar= np.array(ncfile.variables['lon'])[:,:]
    latVar= np.array(ncfile.variables['lat'])[:,:]
    tsVar= np.array(ncfile.variables['timestep'])[:,:]
    icbNumVar= np.array(ncfile.variables['iceberg_number'])[:,:]
    ncfile.close()
    
    nIcb=icbNumVar.shape[0]
    
    for i in np.arange(nIcb):
        index=np.where(latVar[i,:]<1e30)
        lat=latVar[i,index]
        lon=lonVar[i,index]
        ts=tsVar[i,index]
        ids=icbNumVar[i]
        print 'ids',ids
        process=np.mod(ids[0],nProcs)
        #if process==6:
            #ax.scatter(lon,lat,m.to_rgba(ts),cmap=cmap,transform=ccrs.Geodetic(),alpha=1.)
        sct=ax.scatter(lon,lat,facecolor=color,vmin=15000,vmax=30000,edgecolors='None',transform=ccrs.Geodetic(),alpha=1, s=0.6)
    return sct
    
    
def lastIcbergPerProc(path,nProcs):
    files=glob.glob(path+'/*.nc')
    list=np.zeros(nProcs,dtype=int)
    for f in files:
        proc = int(f[len(f)-7:len(f)-3]) + 1
        print f,proc
        ncfile = netCDF4.Dataset(f,'a')
        tsVar= np.array(ncfile.variables['timestep'])[:]
        icbNumVar= np.array(ncfile.variables['iceberg_number'])[:]
        if tsVar.size>0:
            index2=np.where(np.mod(icbNumVar[:,0],nProcs)==proc)
            icbNumTemp=icbNumVar[index2[0],:]
            
            if icbNumTemp.size>0:
                icbNumMax=np.max(icbNumTemp[:,0])
                if np.mod(icbNumMax,nProcs)!=proc:
                    print 'error',np.mod(icbNumMax,nProcs),proc
                print 'nummax',icbNumMax
                list[proc]=icbNumMax
        ncfile.close()
        if proc>nProcs/2:
            break
    
    return list
    