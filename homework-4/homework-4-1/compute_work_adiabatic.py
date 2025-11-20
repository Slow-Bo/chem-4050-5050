#Imports all the functions I can't code as well as a few constants
import numpy as np
import scipy.constants as con
from scipy import integrate as int
import pandas as pd

#Setting the directory to our file so the CSVs output in the disired spot
#Dang I must have been really frustrated when microsoft's stupid CoPilot AI told me how to do this
import os
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

#Defining my values
R = con.R
n = 1
T = 300
V_init = 0.1
ADBindx = 1.4
Vf = np.linspace(V_init, 3 * V_init, 50)

#Defining a function that calculates pressure using the ideal gas law PV=nRT
def Ideal_Gas_Law(V,n,R,T):
    P = (n*R*T)/V
    return P

#Defining a function that calculates pressure using adiabatic expansion P*V^gamma = constant
def ADB_Expansion(V,Pi,Vi,gamma):
    ADBConstant = Pi * Vi ** gamma
    P = ADBConstant/(V ** gamma)
    return P

#Defining an equation to calculate work adiabaticly
def work_calculator_adb(Vi,Vf,n,R,T,N,ADB):
    """
    Calculates the work done ON an ideal gas assuming a system is isothermic
    parameters:
    Vi: initial Volume
    Vf: final Voulme,
    n = # of mols
    R = Ideal gas constant in J/K*mols
    T = Temperature in Kelvin
    N = number of desired steps
    ADB = The Adabatic index
    outputs: 
    w = Work done on the ideal gas
    """
    #Calculates the inital presure of the system
    P_init = Ideal_Gas_Law(Vi,n,R,T)
    #Calulates the 1D length of the integral
    grid_range = (Vf-Vi)
    #Creates the x and y values and calculate step size
    x = np.linspace(Vi, Vf, N)
    y = ADB_Expansion(x,P_init,Vi,ADB)
    dx = grid_range / N
    #calculates the integrals using the defined values and saves it as work
    w = -int.trapezoid(y,x,dx)
    return w


#A print for testing
print(work_calculator_adb(V_init,3*V_init,n,R,T,1000,ADBindx))

#Calculates the work done at every final volume
work_done = []
#Uses a for loop to ensure theres no wonkyness and it puts one VF in at a time
for v in range(len(Vf)):
    work_done.append(work_calculator_adb(V_init,Vf[v],n,R,T,1000,ADBindx))


#Puts the work done as a function of V final into a Pandas dataframe
data = {'Final Volume (m^3)' : Vf, 'Work Done to Gas (J)' : work_done}
df = pd.DataFrame(data)

#another testing print
print(df)

df.to_csv('adiabatic-work.csv', index= False)