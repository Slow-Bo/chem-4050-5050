#importing all the functions I can't code
from scipy.integrate import trapezoid
import numpy as np
import scipy.constants as con
import pandas as pd

#Setting the directory to our file so I can find the rest of my code, and so I save the graphs to the disired spot
#Dang I must have been really frustrated when microsoft's stupid CoPilot AI told me how to do this
import os
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

#defining my variables
h = con.Planck/con.eV #eV*s
epsilon = 0.0103 #angstroms
sigma = 3.4 #angstroms
Dim = 10 #One side of the box in Angstroms
Tmin = 10 #Kelvin
Tmax = 1000 #Kelvin
k_B = con.k / con.eV #eV/Kelvin
N = 2 #Number of particles
m = 39.95 #Mass of argon (AMU)

#Defines the lenard Jones Function
def lennard_jones(r, epsilon=epsilon, sigma=sigma):
    #I added a small value so that if r is 0 it won't ruin everything
    Atomic_Energy = 4*epsilon*((sigma/(r+1E-30))**12-(sigma/(r+1E-30))**6)
    
    return Atomic_Energy

#Defines the function at the heart of the integral
def partition_core(T,x):
    #Calls the lenard jones function and spits out the result at the core of the partician function
    Result = np.e**(-((1 / (k_B * T))*lennard_jones(np.sqrt(x))))
    return Result

def partician_integral_trapezoid(T, N=20, grid_range=10):
    """
    Computes the partician function of two particles in a box.

    Parameters:
    T temperature of the box
    N (int): Number of grid points along each axis. CANNOT EXCEED 20, USES TOO MUCH MEMORY FOR GRIDSPACE ON MY LAPTOP
    grid_range (float): The size .

    Returns:
    float: The computed overlap integral.
    """
    #creates a linespace for x
    x = np.linspace(0, grid_range, N)

    #claculates the thermal wave value for argon at temp
    therm_wav = np.sqrt((1 / (k_B * T))*h**2/(2*np.pi*m))

    #creates a 6d meshgrid of our x, y, and z values. This works, but I know there is a better way to do this. That dosn't fry my hard-drive.
    X1, Y1, Z1, X2, Y2, Z2 = np.meshgrid(x, x, x, x, x, x, indexing='ij')
    #calculates the radius
    r = np.sqrt((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)

    #calculates the integrand
    integrand = partition_core(T,r)

    dx = grid_range / N
    integrals = []
    # prepares for multi integral integration
    for i in range(6):
        #The first time it uses the outside integrand
        if i == 0:
            integrals.append(trapezoid(integrand,x,dx, axis=0))
        else:
            #every other time it uses the preveous integrand
            old_integrand = integrals[i-1]
            integrals.append(trapezoid(old_integrand,x,dx, axis=0))

    #Calculates the partician function at temperature
    Z = 1/(therm_wav**6 * h**6) * integrals[5]


    return Z

#A Print for Testing
#print(partician_integral_trapezoid(10))

#Creating a linespace to calculate U and Cv over creating 100 points will put my PC to it's maximum limits
T = np.linspace(Tmin,Tmax, 100)

#Collecting the values of Z from the functions above
Z = []
#numbers have to be fed into the partition function one at a time because of how I coded it
for i in T:
    Z.append(partician_integral_trapezoid(i))

#Creating a Pandas Dataframe to compile the results
data = {'Temperature' : T, 
        'Partition Function' : Z}
df = pd.DataFrame(data)

#another testing print
print(df)

#Sending the results to a CV
df.to_csv('Partition_Function_Vs_Temp.csv', index= False)