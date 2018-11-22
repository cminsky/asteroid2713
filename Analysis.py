from astropy.io import fits
import numpy as np
import math
from photutils import *
import itertools
import random
import datetime
from scipy import stats
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

class FitsImage:
    """ The FitsImage class loads and holds the data in a FITS image file,
    associated details such as chosen aperture sizes, and includes methods
    for basic photometric analysis and display. """

    def __init__(self,fileName):
        self.fileName = fileName
        hdulist = fits.open(fileName)
        self.data = hdulist[0].data
        self.header = hdulist[0].header
        self.date = hdulist[0].header['DATE-OBS']
        hdulist.close()
        self.refStars = {}

    def setApertures(self,R1,R2,R3):
        self.R1 = R1 # radius of aperture
        self.R2 = R2 # inner radius of annulus
        self.R3 = R3 # outer radius of annulus

    def getData(self):
        return self.data
    def getDate(self):
        return self.date

    def display(self,coords=None,radius=None,aperture=False):
        """ Display the entire image or a square around particular coordinates, 
        optionally with the aperture and annulus size displayed as well. """
        dataSlice = self.data
        
        # trim data to be displayed to square around given coordinates
        if coords:
            x,y = coords[0],coords[1]
            dataSlice = self.data[y-radius:y+radius+1,x-radius:x+radius+1]
        
        # display aperture and annulus circles
        if aperture:
            c1 = plt.Circle((radius,radius),self.R1,color='r',fill=None)
            c2 = plt.Circle((radius,radius),self.R2,color='b',fill=None)
            c3 = plt.Circle((radius,radius),self.R3,color='b',fill=None)
            fig = plt.gcf()
            ax = fig.gca()
            ax.add_artist(c1)
            ax.add_artist(c2)
            ax.add_artist(c3)
            
        plt.imshow(dataSlice,cmap='binary',norm=LogNorm())
        plt.gca().invert_yaxis()
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()

    def getMJD(self):
        """ Modified Julian Date """
        Y,M,D = float(self.date[0:4]),float(self.date[5:7]),float(self.date[8:10])
        h,m,s = float(self.date[11:13]),float(self.date[14:16]),float(self.date[17:])
        # https://en.wikipedia.org/wiki/Julian_day
        JDN = (1461*(Y+4800+(M-14)/12))/4+(367*(M-2-12*((M-14)/12)))/12-(3*((Y+4900+(M-14)/12)/100))/4+D-32075
        JD = JDN + (h-12)/ 24 + m/1440 + s/86400
        MJD = JD -  2400000.5
        return MJD
    
    def getTime(self):
        """ Hours since 00:00 UT """
        h,m,s = float(self.date[11:13]),float(self.date[14:16]),float(self.date[17:])
        hours = h+m/60+s/3600
        return hours

    def getCentroid(self,coords,radius=5):
        """ Coordinates of the centroid of a square within the data."""
        x,y = coords[0],coords[1]
        dataSlice = self.data[y-radius:y+radius+1,x-radius:x+radius+1]
        fromBottomLeft = centroid_2dg(dataSlice)
        newCoords = (x-radius+fromBottomLeft[0],y-radius+fromBottomLeft[1])
        return newCoords

    def getFlux(self,coords,R1=None,R2=None,R3=None):
        """ Calculate the total flux in an aperture minus the background signal in the annulus around it """
        if (R1 or R2 or R3) == None:
            R1 = self.R1
            R2 = self.R2
            R3 = self.R3
            
        # courtesy of MIT 12.410 F18 lecture notes:
        apertures = CircularAperture(coords,r=R1)
        annulus = CircularAnnulus(coords, r_in=R2, r_out=R3)
        apers = [apertures, annulus]
        photTable = aperture_photometry(self.data,apers)
        bkgMean = photTable['aperture_sum_1']/annulus.area()
        bkgSum = bkgMean*apertures.area()
        finalSum = photTable['aperture_sum_0'] - bkgSum
        photTable['residual_aperture_sum'] = finalSum
        
        residuals = [n['residual_aperture_sum'] for n in photTable]
        return residuals

    def curveOfGrowth(self,coords,maxRad=30):
        """ Display a curve of growth to calculate optimal aperture size
        maxRad: int - maximum aperture radius to try
        """
        fluxes = []
        for rad in range(1,maxRad+1):
            fluxes.append(self.getFlux(coords,rad,rad+2,rad+3))
        plt.scatter([rad for rad in range(1,maxRad+1)],fluxes,color='k')
        plt.xlabel('Aperture Radius [pixels]',fontsize=14)
        plt.ylabel('Total Flux in Aperture',fontsize=14)
        plt.ylim(0)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()

class ObservationNight:
    """ Data for a night of observation in the form of a list of FitsImage objects, and
    methods for basic photometric analysis. """
    
    def __init__(self):
        self.files = []
        self.asteroidFluxes = None
        self.asteroidCentroids = None
        self.asteroidRatios = None
        self.magSum = None
        self.asteroidMags = None
        self.refStars = {}

    def loadFiles(self,nameStr,fillNum,startNum,stopNum):
        """ Load FITS files into the list of FitsImage objects.
        nameStr: template for FITS image file names
        fillNum: number of digits in index number
        startNum: first index number
        stopNum: last index number
        returns passes: number of files attempted to load that don't exist
        """
        self.nameStr = nameStr
        self.fillNum = fillNum
        self.startNum = startNum
        self.stopNum = stopNum
        passes = 0
        for n in range(startNum,stopNum+1):
            fileName = (nameStr % str(n).zfill(fillNum))
            try:
                file = FitsImage(fileName)
                self.files.append(file)
            except FileNotFoundError:
                passes += 1
                pass
        return passes

    def popFile(self,index):
        self.files.pop(index)
        
    def setInitialAsteroid(self,x=None,y=None):
        """Set the initial xy coordinates of the asteroid. """
        self.asteroid = (x,y)
        
    def newRefStar(self,name,x,y,mag,centroids=[],fluxes=[]):
        """ Add a star to the reference star dictionary.
        mag: apparent magnitude of star
        """
        self.refStars[name] = {'pos':(x,y),'mag':mag,'centroids':centroids,'fluxes':fluxes}

    def delRefStar(self,name):
        del self.refStars[name]
        
    def setApertures(self,R1):
        """ Set the aperture and annulus sizes for photometry.
        R1 = radius of the aperture
        """
        for file in self.files:
            file.setApertures(R1,R1+2,R1+3)
            
    def getMJDs(self):
        """ Modified Julian Date for each image """
        MJDs = [file.getMJD() for file in self.files]
        return MJDs
    
    def getTimes(self):
        """ Hours since 00:00 UT for each image. """
        times = [file.getTime() for file in self.files]
        return times

    def recursiveCentroids(self,coords,radius=5):
        """ Determine position of an object as it moves throughout the night """
        newPos = coords
        positionList = []
        for file in self.files:
            x,y = newPos[0],newPos[1]
            data = file.getData()
            dataSlice = data[y-radius:y+radius+1,x-radius:x+radius+1]
            delta = centroid_2dg(dataSlice)
            newPos = (newPos[0]-radius+delta[0],newPos[1]-radius+delta[1])
            positionList.append(newPos)
        return positionList
    
    def setAsteroidCentroids(self,radius=5):
        """ Determine position of asteroid """
        asteroidCentroids = self.recursiveCentroids(self.asteroid,radius)
        self.asteroidCentroids = asteroidCentroids
        
    def setStarCentroids(self,radius=5):
        """ Determine position of ref stars """
        if len(self.refStars) == 0:
            self.setRefStars()
        for star in self.refStars:
            initialPos = self.refStars[star]['pos']
            starCentroids = self.recursiveCentroids(initialPos,radius)
            self.refStars[star]['centroids'] = starCentroids

    def recursiveFluxes(self,initialCoords,radius=5,positions=None):
        """ Calculate  flux from an aperture around an object, minus background signal
        from the annulus around it, for each image in a night. """
        if positions == None:
            positions = self.recursiveCentroids(initialCoords,radius)
        objectFluxes = []
        for i in range(len(self.files)):
            file = self.files[i]
            coords = positions[i]
            flux = file.getFlux(coords)
            objectFluxes.append(flux)
        return objectFluxes

    def setAsteroidFluxes(self):
        """ Calculate fluxes of asteroid """
        if self.asteroidCentroids == None:
            self.setAsteroidCentroids()
        asteroidFluxes = self.recursiveFluxes(self.asteroid,positions=self.asteroidCentroids)
        self.asteroidFluxes = asteroidFluxes
        
    def setStarFluxes(self):
        """ Calculate fluxes of ref stars """
        if len(self.refStars) == 0:
            self.setRefStars()
        for star in self.refStars:
            initialPos = self.refStars[star]['pos']
            centroids = self.refStars[star]['centroids']
            if len(centroids) == 0:
                self.setStarCentroids()
                centroids = self.refStars[star]['centroids']
            fluxes = self.recursiveFluxes(initialPos,positions=centroids)
            self.refStars[star]['fluxes'] = fluxes
            
    def setAsteroidRatios(self):
        """Calculate the ratio of flux from the asteroid to combined flux of
        all the reference stars in each image """
        self.setAsteroidFluxes()
        self.setStarFluxes()
        fluxLists = [self.refStars[star]['fluxes'] for star in self.refStars]
        starTotals = [i[0] for i in(np.sum(np.array(fluxLists),axis=0))]
        self.starTotals = starTotals
        ratios = [self.asteroidFluxes[i][0]/starTotals[i] for i in range(len(starTotals))]
        self.asteroidRatios = ratios
        return ratios
    
    def setRefStarMag(self):
        """Calculate the combined magnitude of all the reference stars """
        if len(self.refStars) == 0:
            self.setRefStars()
        magList = [self.refStars[star]['mag'] for star in self.refStars]
        luminosityList = [10**(m/-2.5) for m in magList]
        magSum = -2.5*math.log10(sum(luminosityList))
        self.magSum = magSum
        return magSum
    
    def setAsteroidMags(self):
        """Calculate the magnitude of the asteroid in each image """
        if self.magSum == None:
            self.setRefStarMag()
        ratios = self.setAsteroidRatios()
        magDiffs = [-2.5*math.log10(ratio) for ratio in ratios]
        mags = [mag+self.magSum for mag in magDiffs]
        self.asteroidMags = mags
        return mags
                
    def getFiles(self):
        return self.files
    def getAsteroidMags(self):
        return self.asteroidMags
        
    def graphAsteroidSignal(self):
        """ Time series of the flux from the asteroid throughout the night """
        if self.asteroidFluxes == None:
            self.setAsteroidFluxes()
        plt.scatter(self.getTimes(),self.asteroidFluxes,color='k')
        plt.xlabel('Hours Since 00:00 UT',fontsize=14)
        plt.ylabel('Total Flux in Aperture',fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()
        
    def graphStarFluxes(self):
        """ Time series of the flux from the ref stars throughout the night """
        self.setStarFluxes()
        n = 0
        for star in self.refStars:
            fluxes = self.refStars[star]['fluxes']
            markerStyle = ['o','s','^'][n % 3]
            n += 1
            plt.scatter(self.getTimes(),fluxes,color='k',marker=markerStyle,label=star)
        plt.xlabel('Hours Since 00:00 UT',fontsize=14)
        plt.ylabel('Total Flux in Aperture',fontsize=14)
        plt.legend(fontsize=14,frameon=True,loc='best')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()
        
    def graphStarRatios(self,regress=False):
        """ Time series of the flux ratios of ref stars to each other throughout the night """
        self.setStarFluxes()
        stars = list(self.refStars.keys())
        combos = itertools.combinations(stars,2)
        n = 0
        for pair in combos:
            markerStyle = ['o','s','^'][n % 3]
            fluxes1 = self.refStars[pair[0]]['fluxes']
            fluxes2 = self.refStars[pair[1]]['fluxes']
            ratios = [fluxes1[i][0]/fluxes2[i][0] for i in range(len(fluxes1))]
            meanRatio = np.mean(ratios)
            var = sum([(i-meanRatio)**2 for i in ratios])/(len(ratios)-1)
            stdv = var**(1/2)
            seriesLabel = pair[0]+'/'+pair[1]+', stdv= '+str(stdv)[0:5]
            plt.scatter(self.getTimes(),ratios,label=seriesLabel,color='k',marker=markerStyle)
            if regress:
                plt.plot(self.getTimes(),[meanRatio for i in ratios],color='k',alpha=0.5)
            n += 1
        plt.xlabel('Hours Since 00:00 UT',fontsize=14)
        plt.ylabel('Flux Ratio',fontsize=14)
        plt.legend(fontsize=14,frameon=True,bbox_to_anchor=(0.5,-0.2),loc='upper center')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()
        
    def graphAsteroidMags(self):
        """ Time series of the asteroid magnitude throughout the night """
        if self.asteroidMags == None:
            self.setAsteroidMags()
        plt.scatter(self.getTimes(),self.asteroidMags,color='k')
        plt.gca().invert_yaxis()
        plt.xlabel('Hours Since 00:00 UT',fontsize=14)
        plt.ylabel('Photometric Magnitude',fontsize=14)
        plt.show()