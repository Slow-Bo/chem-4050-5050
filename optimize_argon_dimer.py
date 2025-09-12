#Importing all the cool functions I couldn't code myself
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import pandas as pd

#A simple function that takes an input of r (and potentially epsilon and sigma) and outputs the atomic energy
def lennard_jones(r, epsilon=0.01, sigma=3.4):
    
    Atomic_Energy = 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    return Atomic_Energy

distance_between_two = opt.minimize(
    #This is the function
    fun=lennard_jones,  
    #This is the initial guess
    x0=0,                    
    #The arguments of the function scipy is not allowed to change
    args=(0.01, 3.4),
    #This is the method, I don't remeber the difference between them so I'm using the first one taught in class
    method="Nelder-Mead",    
    #This is the tolernace because computers have yet to acheve perfection
    tol=1e-6                 
)

#A print used for testing
#print(distance_between_two["x"][0])

#graphs the lennard jones function so that the y value is the function at the given X value.
x = np.linspace(3, 6, 100)
y = lennard_jones(x)

#Marks the calculated minimum
plt.scatter(distance_between_two["x"][0], lennard_jones(distance_between_two["x"][0]), color='red', edgecolor='black', s=50, zorder=5, label="Calculated Minimum")

#creates a plot of the polynomial and the original data.
plt.plot(x, y, linestyle = '-',marker='', color='blue', label='Lennard Jones Function')
plt.xlabel('Distance in Angstroms')
plt.ylabel('Potential Energy')
plt.legend()
plt.title('The relationship of energy and distance')
plt.show()