from helper import *
import matplotlib.pyplot as plt
import math

def getUncertainty(model,t,t_ambient,emissivity,u):
    #get the uncertainty in temperature measurement due to 
    #uncertainty in emissivity when backcorrecting for non-perfect
    #emissivity

    #Assumptions:  (1) digital level of detector is proportional to flux on it
    #              (2) background can be represented by a hemispherical
    #                  blackbody
    #              (3) emissivity of the surface is constant over this wavelength
    #                  and temperature range 
    #              (4) imager emissivity is set to 1 
    
    #the backcorrection is equal to the difference between the apparent temperature
    #of the target and the its actual temperature
    
    #the uncertainty is equal to the difference in backcorrections due to the possible
    #range of emissivity values that could be used
  
    corrections = []
    emissivities = [emissivity+u,emissivity-u]
    for e in emissivities:
        #step (1): calculate background radiance from ambient temperature
        #L_background = (1-e)*L_BB(t_ambient)
        L_background = (1-e)*model.radianceLUT(t_ambient)
        
        #step(2): calculate radiance contribution from surface
        L_signal = e*model.radianceLUT(t)
        
        #step (3) calculate the apparent temperature of the surface
        L_tot = L_signal + L_background
        t_apparent = model.temperatureLUT(L_tot)    
        corrections.append(t_apparent)
    
    return (corrections[0]-corrections[1])/math.sqrt(3)
    
def getUncertainty2(model,t_app,emissivity,u):
    w = 11E-6
    emissivities = [emissivity+u,emissivity-u]
    corrections = []
    for e in emissivities:
        denom = w*(math.log(e)*np.exp(model.c2/w/t_app)-1)+1
        t = model.c2/denom
        t-=t_app
        corrections.append(t)
    return (corrections[0]-corrections[1])/math.sqrt(3)
  
def run(model,emissivity,u,t1 = 10+273.15,t2 = 50+273.15):


    t_ambient = 22+273.15
    
    x = np.arange(t1,t2,0.1)
    y = []
    
    for t in x:
        y.append(abs(getUncertainty(model,t,t_ambient,emissivity,u)))
        #y.append(abs(getUncertainty2(model,t,emissivity,u)))
    plt.plot(x-273.15,y)
    plt.show()
    
if __name__ == "__main__":
    model = Model()
    run(model,0.97,0.03)
    run(model,0.6,0.1)
