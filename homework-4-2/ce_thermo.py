#Imports all the functions I can't code as well as a few constants
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, eV
import pandas as pd

#Defining my Constants
k_B = k / eV  
#Creating a linespace for temperature in kelvin
T = np.linspace(300, 2000, 1700)

#Setting the directory to our file so I can find the rest of my code, and so I save the graphs to the disired spot
#Dang I must have been really frustrated when microsoft's stupid CoPilot AI told me how to do this
import os
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

#Calculates the partician function for the no split cesium as a function of temperature
def No_split_Ce(T = 500):
    #An array of the relitive energies
    Ei = [0]
    #An array of the degeneracies
    Degen = [14]
    #An empty array will sum together at the end
    Z_unsumed = []
    
    #calculates the effect of each different energy level on Z
    for n in range(len(Ei)):
        Z_unsumed.append(Degen[n]*np.e**(-Ei[n]/(k_B * T)))
    
    #sums up all the energy levels and returns them as Z
    return sum(Z_unsumed)

#Calculates the partician function for the Spin Orbit Coupling cesium as a function of temperature
def SOC_Ce(T = 500):
    #An array of the relitive energies
    Ei = [0,0.28]
    #An array of the degeneracies
    Degen = [6,8]
    #An empty array will sum together at the end
    Z_unsumed = []
    
    #calculates the effect of each different energy level on Z
    for n in range(len(Ei)):
        Z_unsumed.append(Degen[n]*np.e**(-Ei[n]/(k_B * T)))
    
    #sums up all the energy levels and returns them as Z
    return sum(Z_unsumed)

#Calculates the partician function for the Crystal field spliting cesium as a function of temperature
def SOC_CFS_Ce(T = 500):
    #An array of the relitive energies
    Ei = [0,0.12,0.25,0.32,0.46]
    #An array of the degeneracies
    Degen = [4,2,2,4,2]
    #An empty array will sum together at the end
    Z_unsumed = []
    
    #calculates the effect of each different energy level on Z
    for n in range(len(Ei)):
        Z_unsumed.append(Degen[n]*np.e**(-Ei[n]/(k_B * T)))
    
    #sums up all the energy levels and returns them as Z
    return sum(Z_unsumed)

#Defining the energy functions so I only have to define them once

def Internal_Energy(Z,T):
    #Calculates the internal energy using an np gradient
    U = -np.gradient(np.log(Z), 1 / (k_B * T))
    return U

def Free_Energy(Z,T):
    #Calculates the free energy 
    F = (-k_B * T * np.log(Z))
    return F

def Entropy(F, T):
    #Calculates the entropy using an np gradient
    S = -np.gradient(F, T)
    return S

#This function will create all my graphs so I only need to put the graphing stuff once
def Plot_Relation(name, x_name, y_name, y1_name, x, y1, y2_name = 'null', y2 = [], y3_name = 'null', y3 = []):
    
    #plots y1 with respect to x
    plt.plot(x, y1, linestyle = '-', color= 'blue', label= str(y1_name))
    
    #Ensures y2 and y3 are present before graphing
    if y2_name != 'null':
        plt.plot(x, y2, linestyle = '-', color= 'red', label= str(y2_name))
    if y3_name != 'null':
        plt.plot(x, y3, linestyle = '-', color= 'green', label= str(y3_name))

    #lables the X and Y axis
    plt.xlabel(str(x_name))
    plt.ylabel(str(y_name))

    #Sets the legend to the upper left so it's out of the way
    plt.legend().set_loc('upper left')
    #Titles the graph according to what was defined in the function
    plt.title(str(name))

    #Saves the graph and shows it to ensure accuracy, the arguments bbox_inches and pad_inches ensures that both axis titles are legible
    plt.savefig( str(name) + '.png', bbox_inches= 'tight', pad_inches= 0.3)
    plt.show()

#Calculates the partician at every temp
Z_Ce = No_split_Ce(T)
Z_SOC = SOC_Ce(T)
Z_CFS = SOC_CFS_Ce(T)
Plot_Relation('Partition Function of Cesium 3+ with Respect to Temperature','Temperature (K)','Partition Function (Z)'
              ,'No Splting',T,Z_Ce,'SOC spliting', Z_SOC, 'SOC and CFS splitting', Z_CFS)


#Calculates the Internal energy for every partician function
U_Ce = Internal_Energy(Z_Ce,T)
U_SOC = Internal_Energy(Z_SOC,T)
U_CFS = Internal_Energy(Z_CFS,T)
Plot_Relation('Internal Energy of Cesium 3+ with Respect to Temperature','Temperature (K)','Internal Energy (U) (eV)'
              ,'No Splting',T,U_Ce,'SOC spliting', U_SOC, 'SOC and CFS splitting', U_CFS)

#Calculates the Free energy for every partician function
F_Ce = Free_Energy(Z_Ce,T)
F_SOC = Free_Energy(Z_SOC,T)
F_CFS = Free_Energy(Z_CFS,T)
Plot_Relation('Free Energy of Cesium 3+ with Respect to Temperature','Temperature (K)','Free Energy (F) (eV)'
              ,'No Splting',T,F_Ce,'SOC spliting', F_SOC, 'SOC and CFS splitting', F_CFS)

#Calculates the Entropy for every partician function
S_Ce = Entropy(F_Ce,T)
S_SOC = Entropy(F_SOC,T)
S_CFS = Entropy(F_CFS,T)
Plot_Relation('Entropy of Cesium 3+ with Respect to Temperature','Temperature (K)','Entropy (S) (eV/K)'
              ,'No Splting',T,S_Ce,'SOC spliting', S_SOC, 'SOC and CFS splitting', S_CFS)

#Puts all the values into a Pandas Dataframe
data = {'Temperature' : T, 
        'Partition Function No Spliting' : Z_Ce, 'Partition Function SOC Spliting' : Z_SOC, 'Partition Function SOC + CFS Spliting' : Z_CFS,
        'Internal Energy No Spliting' : U_Ce, 'Internal Energy SOC Spliting': U_SOC, 'Internal Energy SOC + CFS Spliting' : U_CFS,
        'Free Energy No Spliting' : U_Ce, 'Free Energy SOC Spliting' : U_SOC, 'Free Energy SOC + CFS Spliting' : U_CFS,
        'Entropy No Spliting' : S_Ce, 'Entropy SOC Spliting': S_SOC, 'Entropy SOC + CFS Spliting' : S_CFS}
df = pd.DataFrame(data)

#another testing print
print(df)

df.to_csv('isothermic-work.csv', index= False)