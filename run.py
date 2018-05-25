from helper import *
import matplotlib.pyplot as plt
import math

def getUncertainty(model,t,t_ambient,u_t,emissivity,u):
    #get the uncertainty in temperature measurement due to 
    #uncertainty in emissivity when backcorrecting for non-perfect
    #emissivity

    #Assumptions:  (1) background can be represented by a hemispherical
    #                  blackbody
    #              (2) emissivity of the surface is constant over this wavelength
    #                  and temperature range 
    #              (3) imager emissivity is set to 1 
    
    #the backcorrection is equal to the difference between the apparent temperature
    #of the target and the its actual temperature
    
    #the uncertainty is equal to the difference in backcorrections due to the possible
    #range of emissivity values that could be used
  
    corrections = []
    emissivities = [emissivity+u,emissivity-u]
    t_ambients = [t_ambient+u_t,t_ambient-u_t]
    for e in emissivities:
        for t_a in t_ambients:
            #step (1): calculate background radiance from ambient temperature
            #L_background = (1-e)*L_BB(t_ambient)
            L_background = (1-e)*model.radianceLUT(t_a)
            
            #step(2): calculate radiance contribution from surface
            L_signal = e*model.radianceLUT(t)
            
            #step (3) calculate the apparent temperature of the surface
            L_tot = L_signal + L_background
            t_apparent = model.temperatureLUT(L_tot)    
            corrections.append(t_apparent)
    
    return (max(corrections)-min(corrections))/math.sqrt(3)
    
def runEmissivityCalculation(model,t_bkg,u_t,e,u_e,filename):
    print"============================================="
    print"Using file {}".format(filename)
    dir = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(dir,filename)
    data = np.genfromtxt(path)
    epsilon = calcEmissivity(model,data,t_bkg,e) #calc emissivity for each temperature setpoint
    
    
    delta = []
    uncertainties = np.zeros((len(epsilon)))
    emissivities = [e+u_e,e-u_e]
    temperatures = [t_bkg+u_t,t_bkg-u_t]

    for emissivity in emissivities:
        for temperature in temperatures:
            delta.append(calcEmissivity(model,data,temperature,emissivity))
    for i in range(len(delta[0])):
        e1 = max(delta[0][i],delta[1][i],delta[2][i],delta[3][i])
        e2 = min(delta[0][i],delta[1][i],delta[2][i],delta[3][i])
        uncertainties[i] = ((e1-e2)/3.0**0.5)
    for i in range(len(epsilon)):
        print ("HighE temp (degC): {}; emissivity = {}+-{}".format(data[:,0][i],epsilon[i],uncertainties[i]))
    print "Average: {}+-{}: ".format(np.average(epsilon),np.sqrt(np.mean(np.square(uncertainties))))

    print"============================================="

def calcEmissivity(model,data,t_bkg,e):
    highEApparentT = data[:,0] + 273.15
    lowEApparentT = data[:,2] + 273.15
    
    highETotRadiance = model.radianceLUT(highEApparentT)
    L_bkg = model.radianceLUT(t_bkg)
    L_bb = (highETotRadiance-(1-e)*L_bkg)/e
    L_tot = model.radianceLUT(lowEApparentT)
    
    epsilon = (L_tot - L_bkg)/(L_bb - L_bkg)
    #print("{}+-{}".format(np.average(epsilon),np.std(epsilon)))
    return epsilon

def runTemperatureUncertainty(model,emissivity,u,t1 = 10+273.15,t2 = 50+273.15,filename = "out.txt"):
    t_ambient = 22+273.15
    u_ambient = 1
    print"======================================================================="
    
    x = np.arange(t1,t2,0.1)
    y = []
    
    dir = os.path.dirname(os.path.realpath(__file__))
    path = os.path.join(dir,filename)
    f = open(path,'w')
    f.write("#Emissivity = {}+-{}; T_ambient = {}+-{}\n".format(emissivity,u,t_ambient-273.15,u_ambient))
    f.write("#Surface T\t\tT Uncertainty (degC)(k=2)\n")

    for t in x:
        y.append(abs(getUncertainty(model,t,t_ambient,u_ambient,emissivity,u)))
        f.write("{}\t{}\n".format(t,y[-1]))
    f.close()
    plt.plot(x-273.15,y)
    plt.show()
    
if __name__ == "__main__":
    model = Model()
    #model = Model(minW = 2E-6,maxW = 6E-6,responsivityPath = "spectralResponseInfratech.txt")
    runTemperatureUncertainty(model,0.97,0.03,filename = "highET_Uncertainties.txt")
    runTemperatureUncertainty(model,0.8,0.06, filename = "lowET_Uncertainties.txt")
    """
    runEmissivityCalculation(model,22+273.15,1,0.97,0.03,"dataTauNormal.txt")
    runEmissivityCalculation(model,22+273.15,1,0.97,0.03,"dataTauAngle.txt")
    model = Model(minW = 2E-6,maxW = 6E-6,responsivityPath = "spectralResponseInfratech.txt")
    runEmissivityCalculation(model,22+273.15,1,0.97,0.03,"dataInfraNormal.txt")
    runEmissivityCalculation(model,22+273.15,1,0.97,0.03,"dataInfraAngle.txt")
    """
    