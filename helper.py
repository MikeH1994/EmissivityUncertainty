import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import scipy.optimize
import matplotlib.pyplot as plt
import os

class Model:
    def __init__(self,minW = 7E-6,maxW = 14E-6,minT=273.15,maxT=573.15,responsivityPath = "spectralResponse.txt",gainPath = "tauSignalVsTemp.txt"):
        self.nLookupBins = 1000
        self.c1 = 1.1910E-16
        self.c2 = 0.014388
        self.minW = minW
        self.maxW = maxW
        self.minT = minT
        self.maxT = maxT
        
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.responsivityPath =  os.path.join(self.dir,responsivityPath)
        self.gainPath = os.path.join(self.dir,gainPath)
        
        self.responsivityLUT = None #get the responsivity of a detector at a given wavelength
        self.radianceLUT = None #get the radiance over the detector wavelength range at a given temperature
        self.temperatureLUT = None #get the temperature of a blackbody given a certain radiance
        self.digitalLevelLUT = None #get the digital level of the detector at a given radiance
        self.radianceDLLUT = None #get the radiance that corresponds to a given digital level  
        self.generateLUTs()
        
    def generateLUTs(self):
        #create interpolation function for detector responsivity as a function of wavelength
        responsivityData = np.genfromtxt(self.responsivityPath)
        self.responsivityLUT = interp1d(responsivityData[:,0]*1E-6, responsivityData[:,1],fill_value='extrapolate')
        wavelength = np.arange(self.minW,self.maxW,0.1E-6)
        print("LUT 1/5 generated")
        #create interpolation function for radiance as a function of temperature and vice versa
        step = (self.maxT-self.minT)/self.nLookupBins
        temperature = np.arange(self.minT,self.maxT,step)
        radiance = self.blackbodyRadiance(temperature)
        self.radianceLUT = interp1d(temperature,radiance,fill_value='extrapolate')
        print("LUT 2/5 generated")
        self.temperatureLUT = interp1d(radiance,temperature,fill_value='extrapolate')
        print("LUT 3/5 generated")
        #create interpolation  function for DL as a function of radiance and vice versa
        minSignal,maxSignal = 5000.0,12000.0
        step = (maxSignal-minSignal)/self.nLookupBins
        digitalLevel = np.arange(minSignal,maxSignal,step)
        gain,offset = self.getGainAndOffset()
        radiance = (digitalLevel-offset)/gain #DL = gain*radiance + offset
        self.digitalLevelLUT = interp1d(radiance,digitalLevel,fill_value='extrapolate')
        print("LUT 4/5 generated")
        self.radianceDLLUT = interp1d(digitalLevel,radiance,fill_value='extrapolate') 
        print("LUT 5/5 generated")

    def getGainAndOffset(self):
        data = np.genfromtxt(self.gainPath)
        bbTemp = data[:,0] + 273.15
        radiance = self.blackbodyRadiance(bbTemp)
        digitalLevel = data[:,1]
        par,cov = scipy.optimize.curve_fit(self.linearFit,digitalLevel,radiance,p0=(1.0,1.0))
        gain,offset = par
        return gain,offset
                
    def blackbodyRadiance(self,T):
        #Radiance from a blackbody at a temperature T, scaled by the 
        L = 0
        if type(T) is float:
            return self.integrateFunction(T)
        else:
            arr = np.zeros(T.size)
            for i in range(T.size):
                arr[i] = self.integrateFunction(T[i])
            return arr
            
    def integrateFunction(self,T):
        L = quad(self.spectralRadiance, self.minW, self.maxW, args=T)[0]
        """
        nSections = 2
        step = (self.maxW-self.minW)/nSections
        L = 0

        for i in range(nSections):
            w1 = self.minW + i*step
            w2 = w1 + step
            L+= quad(self.spectralRadiance, w1, w2, args=T)[0]
        """
        return L
        
    def spectralRadiance(self,wavelength,T):
        #Spectral radiance emitted by a blackbody of temperature T
        #at a given wavelength (Scaled by detector response)
        L = self.c1/(wavelength**5)/(np.exp(self.c2/wavelength/T)-1)/1E10

        return L*self.responsivityLUT(wavelength)
    
        
    def linearFit(self,x,m,c):
        return m*x + c
        
    def testResponsivity(self):
        wavelength = np.arange(self.minW,self.maxW,0.1E-6)
        responsivity = self.responsivityLUT(wavelength)
        plt.plot(wavelength,responsivity)
        plt.title("Responsivity vs wavelength")
        plt.show()
    
    def testRadiance(self):
        temperature = np.arange(self.minT,self.maxT,0.1)
        radiance = self.radianceLUT(temperature)
        plt.plot(temperature,radiance)
        plt.title("Radiance vs temperature")
        plt.show()
    
    def testDigitalLevel(self):
        digitalLevel = np.arange(5000,12000)
        radiance = self.radianceDLLUT(digitalLevel)
        plt.plot(radiance,digitalLevel)
        plt.title("Digital level vs radiance") 
        plt.show()   
        
        
    def planckIntegral2(self,T):
        #https://people.physics.tamu.edu/krisciunas/planck.pdf
        #=prefix*integral [x^3/(e^x - 1)]
        #x = hf/kT
        #prefix = pi*(2h/c**2)*(k*T/h)**4
        pass

